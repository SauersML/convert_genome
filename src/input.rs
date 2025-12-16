use std::io;
use tracing;
use noodles::vcf::variant::record_buf::{RecordBuf, AlternateBases, Samples, samples::Keys};
use noodles::vcf::variant::record_buf::samples::sample::Value;
// Position removed


use crate::dtc::{self, Record as DtcRecord, Allele as DtcAllele};
use crate::reference::ReferenceGenome;
use crate::conversion::{ConversionConfig, Ploidy, determine_ploidy, format_genotype};
use clap::ValueEnum;

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
        if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
            match ext.to_lowercase().as_str() {
                "vcf" | "vcf.gz" => return Self::Vcf,
                "bcf" => return Self::Bcf,
                "txt" | "csv" | "tsv" => return Self::Dtc,
                _ => {}
            }
        }
        // TODO: magic bytes check
        Self::Dtc // Default fallback
    }
}

/// Trait for a source of genomic variants.
/// Yields `RecordBuf` (generic VCF record wrapper).
pub trait VariantSource {
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>>;
}

impl<T: VariantSource + ?Sized> VariantSource for Box<T> {
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>> {
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
    pub fn new<R: std::io::BufRead>(mut reader: dtc::Reader<R>, reference: ReferenceGenome, config: ConversionConfig) -> Self {
        let keys: Keys = vec![String::from("GT")].into_iter().collect();
        
        // Read all records aggressively
        let mut raw_records = Vec::new();
        let mut parse_errors = 0;
        while let Some(res) = reader.next() {
            match res {
                Ok(rec) => raw_records.push(rec),
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read input line: {}", e);
                    // Continue/Skip bad lines
                }
            }
        }
        
        // Sort by chromosome index and position to optimize ReferenceGenome cache usage
        // This is crucial for performance.
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

    fn convert_record(&self, record: DtcRecord, summary: &mut crate::ConversionSummary) -> io::Result<Option<RecordBuf>> {
        // Validation and normalization logic moved from conversion.rs
        
        let canonical_name = if let Some(name) = self.reference.resolve_contig_name(&record.chromosome) {
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
                    canonical_name, position, e
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
        let ploidy = determine_ploidy(&canonical_name, position, self.config.sex, self.config.par_boundaries.as_ref());
        
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
             },
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
                 },
                 DtcAllele::Deletion => {
                     let del = String::from("DEL");
                     if !alt_bases.contains(&del) { alt_bases.push(del); }
                 },
                 DtcAllele::Insertion => {
                     let ins = String::from("INS");
                     if !alt_bases.contains(&ins) { alt_bases.push(ins); }
                 },
                 DtcAllele::Missing => {}
             }
        }
        
        if alt_bases.is_empty() && !self.config.include_reference_sites {
             summary.skipped_reference_sites += 1;
             return Ok(None);
        }
        
        // Formatted Alts (add < > if needed)
        let formatted_alts: Vec<String> = alt_bases.iter().map(|a| {
            if a == "DEL" || a == "INS" { format!("<{}>", a) } else { a.clone() }
        }).collect();
        
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
        let ids = record.id.as_ref().map(|id| [id.clone()].into_iter().collect::<noodles::vcf::variant::record_buf::Ids>());
        
        let samples = Samples::new(
            self.header_keys.clone(),
            vec![vec![Some(Value::String(genotype_string))]]
        );
        
        // INFO
        let mut info_map: noodles::vcf::variant::record_buf::Info = Default::default();
        
        let has_del = alt_bases.iter().any(|b| b == "DEL");
        let has_ins = alt_bases.iter().any(|b| b == "INS");
        
        if has_del || has_ins {
             // IMPRECISE
             info_map.insert(String::from("IMPRECISE"), Some(noodles::vcf::variant::record_buf::info::field::Value::Flag));
             // SVTYPE
             let svtype = if has_del && has_ins { "COMPLEX" } else if has_del { "DEL" } else { "INS" };
             info_map.insert(String::from("SVTYPE"), Some(noodles::vcf::variant::record_buf::info::field::Value::String(svtype.into())));
        }

        Ok(Some(RecordBuf::builder()
            .set_reference_sequence_name(canonical_name)
            .set_variant_start(noodles::core::Position::new(position as usize).ok_or(
                 io::Error::new(io::ErrorKind::InvalidData, "invalid position")
            )?)
            .set_ids(ids.unwrap_or_default()) // Use unwrap_or_default if ids can be None
            .set_reference_bases(reference_base.to_string())
            .set_alternate_bases(alternate_bases)
            .set_info(info_map)
            .set_samples(samples)
            .build()))
    }
}

impl VariantSource for DtcSource {
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>> {
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
    pub fn new<R: BufRead>(mut reader: vcf::io::Reader<R>, reference: &ReferenceGenome) -> io::Result<Self> {
        let header = reader.read_header()?;
        
        let mut raw_records = Vec::new();
        let mut parse_errors = 0;
        
        for result in reader.record_bufs(&header) {
            match result {
                Ok(record) => raw_records.push(record),
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read VCF record: {}", e);
                }
            }
        }
        
        // Sort by contig index and position
        raw_records.sort_by_cached_key(|r| {
            let idx = reference.contig_index(r.reference_sequence_name()).unwrap_or(usize::MAX);
            let pos = r.variant_start().map(|p| usize::from(p)).unwrap_or(0);
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
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>> {
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
    pub fn new<R: std::io::Read>(mut reader: bcf::io::Reader<R>, reference: &ReferenceGenome) -> io::Result<Self> {
        let header = reader.read_header()?;
        
        let mut raw_records = Vec::new();
        let mut parse_errors = 0;
        
        for result in reader.record_bufs(&header) {
            match result {
                Ok(record) => raw_records.push(record),
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read BCF record: {}", e);
                }
            }
        }
        
        // Sort by contig index and position
        raw_records.sort_by_cached_key(|r| {
            let idx = reference.contig_index(r.reference_sequence_name()).unwrap_or(usize::MAX);
            let pos = r.variant_start().map(|p| usize::from(p)).unwrap_or(0);
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
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>> {
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
