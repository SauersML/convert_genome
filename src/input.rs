use noodles::vcf::variant::record_buf::samples::sample::Value;
use noodles::vcf::variant::record_buf::{AlternateBases, RecordBuf, Samples, samples::Keys};
use std::io;
use tracing;

use crate::conversion::{ConversionConfig, Ploidy, determine_ploidy, format_genotype};
use crate::dtc::{self, Allele as DtcAllele, Record as DtcRecord};
use crate::reference::ReferenceGenome;
use clap::ValueEnum;

/// Get optional max records limit from environment variable.
/// Set `CONVERT_GENOME_MAX_RECORDS=N` to limit buffered records (for OOM protection).
fn get_max_records_limit() -> Option<usize> {
    std::env::var("CONVERT_GENOME_MAX_RECORDS")
        .ok()
        .and_then(|s| s.parse().ok())
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum InputFormat {
    /// Direct-to-consumer text output (23andMe, etc.)
    Dtc,
    /// Variant Call Format
    Vcf,
    /// Binary Call Format
    Bcf,
    /// Detect format automatically
    Auto,
}

impl InputFormat {
    pub fn detect(path: &std::path::Path) -> Self {
        // First check file extension
        if let Some(filename) = path.file_name().map(|n| n.to_string_lossy().to_lowercase()) {
            if filename.ends_with(".vcf.gz") || filename.ends_with(".vcf") {
                return Self::Vcf;
            } else if filename.ends_with(".bcf") || filename.ends_with(".bcf.gz") {
                return Self::Bcf;
            } else if filename.ends_with(".txt")
                || filename.ends_with(".csv")
                || filename.ends_with(".tsv")
            {
                return Self::Dtc;
            }
        }

        // Check magic bytes for BCF (starts with 'BCF' = 0x42 0x43 0x46)
        if let Ok(file) = std::fs::File::open(path) {
            use std::io::Read;
            let mut reader = std::io::BufReader::new(file);
            let mut magic = [0u8; 3];
            if reader.read_exact(&mut magic).is_ok() {
                if &magic == b"BCF" {
                    return Self::Bcf;
                }
                // VCF files typically start with '#' (0x23) for header
                if magic[0] == b'#' {
                    return Self::Vcf;
                }
            }
        }

        Self::Dtc // Default fallback for text formats
    }
}

/// Trait for a source of genomic variants.
/// Yields `RecordBuf` (generic VCF record wrapper).
pub trait VariantSource {
    fn next_variant(
        &mut self,
        summary: &mut crate::ConversionSummary,
    ) -> Option<io::Result<RecordBuf>>;
}

impl<T: VariantSource + ?Sized> VariantSource for Box<T> {
    fn next_variant(
        &mut self,
        summary: &mut crate::ConversionSummary,
    ) -> Option<io::Result<RecordBuf>> {
        (**self).next_variant(summary)
    }
}

/// Adapter for reading DTC (23andMe) files.
/// Buffers and sorts raw records to ensure cache locality and minimize memory usage.
pub struct DtcSource {
    records: std::vec::IntoIter<DtcRecord>,
    reference: ReferenceGenome,
    config: ConversionConfig,
    header_keys: Keys,
    initial_parse_errors: usize,
    initial_stats_synced: bool,
}

impl DtcSource {
    pub fn new<R: std::io::BufRead>(
        reader: dtc::Reader<R>,
        reference: ReferenceGenome,
        config: ConversionConfig,
    ) -> Self {
        let keys: Keys = vec![String::from("GT")].into_iter().collect();

        // Get optional limit from environment
        let max_records = get_max_records_limit();

        // Read records (with optional limit for OOM protection)
        let mut raw_records = Vec::new();
        let mut parse_errors = 0;
        let mut records_read = 0;

        for res in reader {
            if let Some(limit) = max_records {
                if records_read >= limit {
                    tracing::info!(limit, "reached max_records limit, stopping read");
                    break;
                }
            }
            match res {
                Ok(rec) => {
                    raw_records.push(rec);
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read input line: {}", e);
                }
            }
        }

        // Sort by chromosome index and position to optimize ReferenceGenome cache usage
        raw_records.sort_by_cached_key(|r| {
            let idx = reference.contig_index(&r.chromosome).unwrap_or(usize::MAX);
            (idx, r.position)
        });

        Self {
            records: raw_records.into_iter(),
            reference,
            config,
            header_keys: keys,
            initial_parse_errors: parse_errors,
            initial_stats_synced: false,
        }
    }

    fn convert_record(
        &self,
        record: DtcRecord,
        summary: &mut crate::ConversionSummary,
    ) -> io::Result<Option<RecordBuf>> {
        // Validation and normalization logic moved from conversion.rs

        let canonical_name =
            if let Some(name) = self.reference.resolve_contig_name(&record.chromosome) {
                name.to_string()
            } else {
                // Skipping unknown chromosome
                summary.unknown_chromosomes += 1;
                return Ok(None);
            };

        // Validate position
        let position = record.position;

        // Use reference.base instead of get_context
        let reference_base = match self.reference.base(&canonical_name, position) {
            Ok(base) => base,
            Err(e) => {
                summary.reference_failures += 1;
                tracing::warn!(
                    "reference lookup failed for {}:{}: {}",
                    canonical_name,
                    position,
                    e
                );
                return Ok(None);
            }
        };

        // Parse alleles
        let allele_states = match record.parse_alleles() {
            Ok(states) => states,
            Err(e) => {
                summary.invalid_genotypes += 1;
                tracing::warn!("invalid alleles for {}:{}: {}", canonical_name, position, e);
                return Ok(None);
            }
        };

        // Biological Ploidy
        let ploidy = determine_ploidy(
            &canonical_name,
            position,
            self.config.sex,
            self.config.par_boundaries.as_ref(),
        );

        // Normalize
        let normalized_alleles = match ploidy {
            Ploidy::Zero => {
                summary.missing_genotype_records += 1;
                return Ok(None);
            }
            Ploidy::Haploid => {
                let unique: Vec<_> = allele_states
                    .iter()
                    .filter(|&a| *a != DtcAllele::Missing)
                    .cloned()
                    .collect();

                let mut uniq_vals = unique.clone();
                uniq_vals.dedup();

                if uniq_vals.is_empty() {
                    vec![DtcAllele::Missing]
                } else if uniq_vals.len() == 1 {
                    vec![uniq_vals[0].clone()]
                } else {
                    vec![DtcAllele::Missing]
                }
            }
            Ploidy::Diploid => allele_states,
        };

        // Build Alt Bases
        let mut alt_bases = Vec::new();
        for allele in normalized_alleles.iter() {
            match allele {
                DtcAllele::Base(base) => {
                    let b = base.to_uppercase();
                    let r = reference_base.to_string().to_uppercase();
                    if b != r && !alt_bases.contains(&b) {
                        alt_bases.push(b);
                    }
                }
                DtcAllele::Deletion => {
                    let del = String::from("DEL");
                    if !alt_bases.contains(&del) {
                        alt_bases.push(del);
                    }
                }
                DtcAllele::Insertion => {
                    let ins = String::from("INS");
                    if !alt_bases.contains(&ins) {
                        alt_bases.push(ins);
                    }
                }
                DtcAllele::Missing => {}
            }
        }

        if alt_bases.is_empty() && !self.config.include_reference_sites {
            summary.skipped_reference_sites += 1;
            return Ok(None);
        }

        // Formatted Alts (add < > if needed)
        let formatted_alts: Vec<String> = alt_bases
            .iter()
            .map(|a| {
                if a == "DEL" || a == "INS" {
                    format!("<{}>", a)
                } else {
                    a.clone()
                }
            })
            .collect();

        let alternate_bases = AlternateBases::from(formatted_alts);

        // Genotype string
        let genotype_string = format_genotype(&normalized_alleles, reference_base, &alt_bases)
            .map_err(|e| {
                // Do not increment summary here if we return Err, process_records will count it as parse_error?
                // But wait, user noted double counting.
                // If we return Err, process_records increments parse_errors.
                // If we want to skip it as invalid_genotype, we should return Ok(None) and increment here.
                summary.invalid_genotypes += 1;
                io::Error::new(io::ErrorKind::InvalidData, e.to_string())
            });

        let genotype_string = match genotype_string {
            Ok(s) => s,
            Err(e) => {
                tracing::warn!("genotype formatting failed: {}", e);
                return Ok(None); // Skip, don't return Err to avoid double count
            }
        };

        // Build RecordBuf with IDs
        let ids = record.id.as_ref().map(|id| {
            [id.clone()]
                .into_iter()
                .collect::<noodles::vcf::variant::record_buf::Ids>()
        });

        let samples = Samples::new(
            self.header_keys.clone(),
            vec![vec![Some(Value::String(genotype_string))]],
        );

        // INFO
        let mut info_map: noodles::vcf::variant::record_buf::Info = Default::default();

        let has_del = alt_bases.iter().any(|b| b == "DEL");
        let has_ins = alt_bases.iter().any(|b| b == "INS");

        if has_del || has_ins {
            // Mark as imprecise structural variant
            info_map.insert(
                String::from("IMPRECISE"),
                Some(noodles::vcf::variant::record_buf::info::field::Value::Flag),
            );
            // Set structural variant type
            let svtype = if has_del && has_ins {
                "COMPLEX"
            } else if has_del {
                "DEL"
            } else {
                "INS"
            };
            info_map.insert(
                String::from("SVTYPE"),
                Some(noodles::vcf::variant::record_buf::info::field::Value::String(svtype.into())),
            );
        }

        Ok(Some(
            RecordBuf::builder()
                .set_reference_sequence_name(canonical_name)
                .set_variant_start(noodles::core::Position::new(position as usize).ok_or(
                    io::Error::new(io::ErrorKind::InvalidData, "invalid position"),
                )?)
                .set_ids(ids.unwrap_or_default()) // Use unwrap_or_default if ids can be None
                .set_reference_bases(reference_base.to_string())
                .set_alternate_bases(alternate_bases)
                .set_info(info_map)
                .set_samples(samples)
                .build(),
        ))
    }
}

impl VariantSource for DtcSource {
    fn next_variant(
        &mut self,
        summary: &mut crate::ConversionSummary,
    ) -> Option<io::Result<RecordBuf>> {
        // Sync initial stats once
        if !self.initial_stats_synced {
            summary.parse_errors += self.initial_parse_errors;
            // Note: total_records should typically track *valid* or *attemped* records.
            // If we want total to include malformed lines, add them here:
            // summary.total_records += self.initial_parse_errors as u64;
            // The previous behavior was to count lines read.
            self.initial_stats_synced = true;
        }

        // We iterate over the pre-sorted raw records
        while let Some(record) = self.records.next() {
            summary.total_records += 1;
            match self.convert_record(record, summary) {
                Ok(Some(buf)) => return Some(Ok(buf)),
                Ok(None) => continue, // Filtered
                Err(e) => return Some(Err(e)),
            }
        }
        None
    }
}

// ============================================================================
// VCF Source Adapter
// ============================================================================

use noodles::vcf;
use std::io::BufRead;

/// Adapter for reading VCF files.
/// Buffers and sorts records by contig index and position.
pub struct VcfSource {
    records: std::vec::IntoIter<RecordBuf>,
    initial_parse_errors: usize,
    initial_stats_synced: bool,
}

impl VcfSource {
    pub fn new<R: BufRead>(
        mut reader: vcf::io::Reader<R>,
        reference: &ReferenceGenome,
    ) -> io::Result<Self> {
        let header = reader.read_header()?;

        // Get optional limit from environment
        let max_records = get_max_records_limit();

        let mut raw_records = Vec::new();
        let mut parse_errors = 0;
        let mut records_read = 0;

        for result in reader.record_bufs(&header) {
            if let Some(limit) = max_records {
                if records_read >= limit {
                    tracing::info!(limit, "reached max_records limit, stopping read");
                    break;
                }
            }
            match result {
                Ok(record) => {
                    raw_records.push(record);
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read VCF record: {}", e);
                }
            }
        }

        // Sort by contig index and position
        raw_records.sort_by_cached_key(|r| {
            let idx = reference
                .contig_index(r.reference_sequence_name())
                .unwrap_or(usize::MAX);
            let pos = r.variant_start().map(usize::from).unwrap_or(0);
            (idx, pos)
        });

        Ok(Self {
            records: raw_records.into_iter(),
            initial_parse_errors: parse_errors,
            initial_stats_synced: false,
        })
    }
}

impl VariantSource for VcfSource {
    fn next_variant(
        &mut self,
        summary: &mut crate::ConversionSummary,
    ) -> Option<io::Result<RecordBuf>> {
        if !self.initial_stats_synced {
            summary.parse_errors += self.initial_parse_errors;
            self.initial_stats_synced = true;
        }

        self.records.next().map(|record| {
            summary.total_records += 1;
            Ok(record)
        })
    }
}

// ============================================================================
// BCF Source Adapter
// ============================================================================

use noodles::bcf;

/// Adapter for reading BCF files.
/// Buffers and sorts records by contig index and position.
pub struct BcfSource {
    records: std::vec::IntoIter<RecordBuf>,
    initial_parse_errors: usize,
    initial_stats_synced: bool,
}

impl BcfSource {
    pub fn new<R: std::io::Read>(
        mut reader: bcf::io::Reader<R>,
        reference: &ReferenceGenome,
    ) -> io::Result<Self> {
        let header = reader.read_header()?;

        // Get optional limit from environment
        let max_records = get_max_records_limit();

        let mut raw_records = Vec::new();
        let mut parse_errors = 0;
        let mut records_read = 0;

        for result in reader.record_bufs(&header) {
            if let Some(limit) = max_records {
                if records_read >= limit {
                    tracing::info!(limit, "reached max_records limit, stopping read");
                    break;
                }
            }
            match result {
                Ok(record) => {
                    raw_records.push(record);
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read BCF record: {}", e);
                }
            }
        }

        // Sort by contig index and position
        raw_records.sort_by_cached_key(|r| {
            let idx = reference
                .contig_index(r.reference_sequence_name())
                .unwrap_or(usize::MAX);
            let pos = r.variant_start().map(usize::from).unwrap_or(0);
            (idx, pos)
        });

        Ok(Self {
            records: raw_records.into_iter(),
            initial_parse_errors: parse_errors,
            initial_stats_synced: false,
        })
    }
}

impl VariantSource for BcfSource {
    fn next_variant(
        &mut self,
        summary: &mut crate::ConversionSummary,
    ) -> Option<io::Result<RecordBuf>> {
        if !self.initial_stats_synced {
            summary.parse_errors += self.initial_parse_errors;
            self.initial_stats_synced = true;
        }

        self.records.next().map(|record| {
            summary.total_records += 1;
            Ok(record)
        })
    }
}
