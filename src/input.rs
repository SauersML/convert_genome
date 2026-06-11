use noodles::vcf::variant::record_buf::samples::sample::Value;
use noodles::vcf::variant::record_buf::{AlternateBases, RecordBuf, Samples, samples::Keys};
use std::io;
use tracing;

use crate::conversion::{ConversionConfig, Ploidy, determine_ploidy, format_genotype};
use crate::dtc::{Allele as DtcAllele, Record as DtcRecord};
use crate::external_sort::{
    DtcExternalSorter, DtcOrder, RecordExternalSorter, RecordOrder, SortFormat, SortedRecordIter,
};
use crate::reference::ReferenceGenome;
use clap::ValueEnum;

/// Get optional max records limit from environment variable.
/// Set `CONVERT_GENOME_MAX_RECORDS=N` to limit buffered records (for OOM protection).
pub fn get_max_records_limit() -> Option<usize> {
    std::env::var("CONVERT_GENOME_MAX_RECORDS")
        .ok()
        .and_then(|s| s.parse().ok())
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum InputFormat {
    /// Direct-to-consumer text output (23andMe, etc.)
    Dtc,
    /// Illumina GenomeStudio Genotyping ("GSGT") Final Report.
    GenomeStudio,
    /// Variant Call Format
    Vcf,
    /// Binary Call Format
    Bcf,
    /// Detect format automatically
    Auto,
}

impl InputFormat {
    /// Whether this format flows through the DTC conversion pipeline
    /// (record-oriented genotype text that needs build/sex inference and
    /// reference-based REF/ALT resolution). Both 23andMe-style DTC files and
    /// Illumina GenomeStudio Final Reports qualify.
    pub fn is_dtc_like(self) -> bool {
        matches!(self, Self::Dtc | Self::GenomeStudio)
    }

    pub fn detect(path: &std::path::Path) -> Self {
        use std::io::BufRead;

        // Use smart reader to peel off compression layers (GZIP, ZIP, etc.)
        if let Ok(mut reader) = crate::smart_reader::open_input(path) {
            // Peek at the beginning of the stream
            if let Ok(buf) = reader.fill_buf() {
                if buf.len() >= 3 && &buf[..3] == b"BCF" {
                    return Self::Bcf;
                }
                // Illumina GenomeStudio Final Report: starts with a [Header]
                // (or [Data]) section marker, or a "GSGT Version" line.
                if crate::genomestudio::looks_like_genome_studio(buf) {
                    return Self::GenomeStudio;
                }
                // VCF files typically start with '#' (0x23) for header
                if !buf.is_empty() && buf[0] == b'#' {
                    // It could be a comment in DTC, but usually VCF starts with ##fileformat=VCF
                    // We only confidently return VCF if we see the strict header.
                    if buf.starts_with(b"##fileformat=VCF") {
                        return Self::Vcf;
                    }
                    // If it starts with '#' but not the strict header, it might be DTC with comments.
                    // We fall through to extension check or DTC default.
                }
            }
        }

        // Fallback: Check file extension if content detection failed or ambiguous
        // (Though smart_reader failure usually means file error)
        if let Some(filename) = path.file_name().map(|n| n.to_string_lossy().to_lowercase()) {
            if filename.ends_with(".vcf.gz") || filename.ends_with(".vcf") {
                return Self::Vcf;
            } else if filename.ends_with(".bcf") || filename.ends_with(".bcf.gz") {
                return Self::Bcf;
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
    records: crate::external_sort::DtcSortedIter,
    reference: Option<ReferenceGenome>,
    config: ConversionConfig,
    header_keys: Keys,
    /// When true, per-sample [`crate::dtc::ArrayMetrics`] are emitted as the
    /// BAF/LRR/IGC/GTS FORMAT fields (GenomeStudio inputs). The key order here
    /// must match [`DtcSource::header_keys`] and the FORMAT lines in the header.
    emit_metrics: bool,
    initial_parse_errors: usize,
    initial_stats_synced: bool,
    inferred_strand: crate::source_ref::InferredStrand,
}

/// FORMAT field IDs for the GenomeStudio array metrics, in the order they are
/// written into each sample (after `GT`). Shared with the VCF header builder so
/// the declared FORMAT lines and the emitted values never drift apart.
pub const ARRAY_METRIC_KEYS: [&str; 4] = ["BAF", "LRR", "IGC", "GTS"];

impl DtcSource {
    /// Build a source from any iterator of parsed [`DtcRecord`]s.
    ///
    /// Generic over the record iterator so it can consume both the 23andMe-style
    /// [`crate::dtc::Reader`] and the Illumina [`crate::genomestudio::Reader`]; both
    /// yield `Result<DtcRecord, _>` with a `Display`able error.
    pub fn new<I, E>(
        records: I,
        reference: Option<ReferenceGenome>,
        config: ConversionConfig,
        inferred_strand: Option<crate::source_ref::InferredStrand>,
    ) -> io::Result<Self>
    where
        I: IntoIterator<Item = Result<DtcRecord, E>>,
        E: std::fmt::Display,
    {
        // GenomeStudio inputs carry per-sample array metrics; declare them as
        // FORMAT fields after GT. Other DTC inputs emit GT only.
        let emit_metrics = matches!(config.input_format, crate::input::InputFormat::GenomeStudio);
        let keys: Keys = if emit_metrics {
            std::iter::once("GT".to_string())
                .chain(ARRAY_METRIC_KEYS.iter().map(|s| s.to_string()))
                .collect()
        } else {
            vec![String::from("GT")].into_iter().collect()
        };

        // Get optional limit from environment
        let max_records = get_max_records_limit();

        let order = if let Some(ref r) = reference {
            DtcOrder::from_reference(r)
        } else {
            DtcOrder::Natural
        };
        let mut sorter = DtcExternalSorter::new(order);

        // Read records (with optional limit for OOM protection)
        let mut parse_errors = 0;
        let mut records_read = 0;

        for res in records {
            if let Some(limit) = max_records
                && records_read >= limit
            {
                tracing::info!(limit, "reached max_records limit, stopping read");
                break;
            }
            match res {
                Ok(rec) => {
                    sorter.push(rec)?;
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read input line: {}", e);
                }
            }
        }

        let records = sorter.finish()?;

        Ok(Self {
            records,
            reference,
            config,
            header_keys: keys,
            emit_metrics,
            initial_parse_errors: parse_errors,
            initial_stats_synced: false,
            inferred_strand: inferred_strand.unwrap_or(crate::source_ref::InferredStrand::Forward),
        })
    }

    fn convert_record(
        &self,
        record: DtcRecord,
        summary: &mut crate::ConversionSummary,
    ) -> io::Result<Option<RecordBuf>> {
        // Validation and normalization logic moved from conversion.rs

        let canonical_name = if let Some(ref r) = self.reference {
            if let Some(name) = r.resolve_contig_name(&record.chromosome) {
                name.to_string()
            } else {
                // Skipping unknown chromosome
                summary.unknown_chromosomes += 1;
                return Ok(None);
            }
        } else {
            // No reference: pass through chromosome name (normalized slightly)
            record.chromosome.clone()
        };

        // Validate position
        let position = record.position;

        // Use reference.base instead of get_context
        let reference_base = if let Some(ref r) = self.reference {
            match r.base(&canonical_name, position) {
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
            }
        } else {
            // No reference: use placeholder 'N'
            'N'
        };

        // Parse alleles
        let mut allele_states = match record.parse_alleles() {
            Ok(states) => states,
            Err(e) => {
                summary.invalid_genotypes += 1;
                tracing::warn!("invalid alleles for {}:{}: {}", canonical_name, position, e);
                return Ok(None);
            }
        };

        // Apply Strand Lock
        if self.inferred_strand == crate::source_ref::InferredStrand::Reverse {
            for allele in &mut allele_states {
                if let DtcAllele::Base(s) = allele {
                    *s = reverse_complement(s);
                }
            }
        }

        // Biological Ploidy
        let ploidy = determine_ploidy(
            &canonical_name,
            position,
            self.config.sex.unwrap_or(crate::cli::Sex::Unknown),
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

        // Genotype first, then the array metrics (when this is a GenomeStudio
        // source). Order must match `header_keys` / `ARRAY_METRIC_KEYS`.
        let mut sample_values: Vec<Option<Value>> = Vec::with_capacity(1 + ARRAY_METRIC_KEYS.len());
        sample_values.push(Some(Value::String(genotype_string)));
        if self.emit_metrics {
            let m = record.metrics.unwrap_or_default();
            sample_values.push(m.baf.map(Value::Float));
            sample_values.push(m.lrr.map(Value::Float));
            sample_values.push(m.gencall.map(Value::Float));
            sample_values.push(m.gentrain.map(Value::Float));
        }
        let samples = Samples::new(self.header_keys.clone(), vec![sample_values]);

        // Build info section
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
pub struct VcfSource {
    records: SortedRecordIter,
    initial_parse_errors: usize,
    initial_stats_synced: bool,
}

impl VcfSource {
    /// Create a new VcfSource using natural chromosome ordering (no reference required)
    pub fn new_without_reference<R: BufRead>(mut reader: vcf::io::Reader<R>) -> io::Result<Self> {
        let header = reader.read_header()?;
        let max_records = get_max_records_limit();
        let mut parse_errors = 0;
        let mut records_read = 0;
        let mut sorter =
            RecordExternalSorter::new(header.clone(), SortFormat::Vcf, RecordOrder::Natural);

        for result in reader.record_bufs(&header) {
            if let Some(limit) = max_records
                && records_read >= limit
            {
                tracing::info!(limit, "reached max_records limit, stopping read");
                break;
            }
            match result {
                Ok(record) => {
                    sorter.push(record)?;
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read VCF record: {}", e);
                }
            }
        }

        let records = sorter.finish()?;

        Ok(Self {
            records,
            initial_parse_errors: parse_errors,
            initial_stats_synced: false,
        })
    }

    pub fn new<R: BufRead>(
        mut reader: vcf::io::Reader<R>,
        reference: &ReferenceGenome,
    ) -> io::Result<Self> {
        let header = reader.read_header()?;
        let max_records = get_max_records_limit();
        let mut parse_errors = 0;
        let mut records_read = 0;
        let order = RecordOrder::from_reference(reference);
        let mut sorter = RecordExternalSorter::new(header.clone(), SortFormat::Vcf, order);

        for result in reader.record_bufs(&header) {
            if let Some(limit) = max_records
                && records_read >= limit
            {
                tracing::info!(limit, "reached max_records limit, stopping read");
                break;
            }
            match result {
                Ok(record) => {
                    sorter.push(record)?;
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read VCF record: {}", e);
                }
            }
        }

        let records = sorter.finish()?;

        Ok(Self {
            records,
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
            record
        })
    }
}

// ============================================================================
// BCF Source Adapter
// ============================================================================

use noodles::bcf;

/// Adapter for reading BCF files.
pub struct BcfSource {
    records: SortedRecordIter,
    initial_parse_errors: usize,
    initial_stats_synced: bool,
}

impl BcfSource {
    /// Create a new BcfSource using natural chromosome ordering (no reference required)
    pub fn new_without_reference<R: std::io::Read>(
        mut reader: bcf::io::Reader<R>,
    ) -> io::Result<Self> {
        let header = reader.read_header()?;
        let max_records = get_max_records_limit();
        let mut parse_errors = 0;
        let mut records_read = 0;
        let mut sorter =
            RecordExternalSorter::new(header.clone(), SortFormat::Bcf, RecordOrder::Natural);

        for result in reader.record_bufs(&header) {
            if let Some(limit) = max_records
                && records_read >= limit
            {
                tracing::info!(limit, "reached max_records limit, stopping read");
                break;
            }
            match result {
                Ok(record) => {
                    sorter.push(record)?;
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read BCF record: {}", e);
                }
            }
        }

        let records = sorter.finish()?;

        Ok(Self {
            records,
            initial_parse_errors: parse_errors,
            initial_stats_synced: false,
        })
    }

    pub fn new<R: std::io::Read>(
        mut reader: bcf::io::Reader<R>,
        reference: &ReferenceGenome,
    ) -> io::Result<Self> {
        let header = reader.read_header()?;
        let max_records = get_max_records_limit();
        let mut parse_errors = 0;
        let mut records_read = 0;
        let order = RecordOrder::from_reference(reference);
        let mut sorter = RecordExternalSorter::new(header.clone(), SortFormat::Bcf, order);

        for result in reader.record_bufs(&header) {
            if let Some(limit) = max_records
                && records_read >= limit
            {
                tracing::info!(limit, "reached max_records limit, stopping read");
                break;
            }
            match result {
                Ok(record) => {
                    sorter.push(record)?;
                    records_read += 1;
                }
                Err(e) => {
                    parse_errors += 1;
                    tracing::warn!("failed to read BCF record: {}", e);
                }
            }
        }

        let records = sorter.finish()?;

        Ok(Self {
            records,
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
            record
        })
    }
}

// ============================================================================
// Natural Chromosome Ordering (no reference required)
// ============================================================================

/// Returns a sort key for chromosome names using natural ordering.
/// This allows sorting without a reference genome.
///
/// Order: 1-22, X, Y, MT/M, then other contigs alphabetically
pub fn natural_contig_order(name: &str) -> (u32, String) {
    let normalized = name
        .trim_start_matches("chr")
        .trim_start_matches("Chr")
        .trim_start_matches("CHR");

    // Try to parse as a number (autosomes 1-22)
    if let Ok(num) = normalized.parse::<u32>() {
        return (num, String::new());
    }

    // Handle sex chromosomes and mitochondrial
    match normalized.to_uppercase().as_str() {
        "X" => (23, String::new()),
        "Y" => (24, String::new()),
        "M" | "MT" => (25, String::new()),
        other => (1000, other.to_string()), // Other contigs sorted alphabetically after standard ones
    }
}

fn reverse_complement(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            'N' => 'N',
            'n' => 'n',
            _ => c,
        })
        .collect()
}
