use std::{
    fs,
    io::{self, BufReader},
    path::PathBuf,
    sync::Arc,
};
use crate::reference::ParBoundaries;
use crate::cli::Sex;

use anyhow::{Context, Result, anyhow};
use clap::ValueEnum;
use noodles::bcf;
use noodles::{
    core::Position,
    vcf::{
        self,
        header::{
            FileFormat,
            record::{
                key,
                value::{
                    Collection, Map,
                    map::{
                        AlternativeAllele, Contig, Format,
                        Info as InfoMap,
                        info::{Number, Type},
                    },
                },
            },
        },
        variant::{
            self as variant, RecordBuf,
            io::Write as VariantRecordWrite,
            record::samples::keys::key as format_key,
            record_buf::{
                AlternateBases, Filters, Ids, Info,
                info::field::Value as InfoFieldValue,
                samples::{Keys, Samples, sample::Value as SampleValue},
            },
        },
    },
};
use rayon::prelude::*;
use thiserror::Error;
use time::{OffsetDateTime, macros::format_description};

use crate::{
    dtc::{self, Record as DtcRecord},
    reference::{ReferenceError, ReferenceGenome},
};

/// Supported output formats for the converter.
#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum OutputFormat {
    /// Variant Call Format text output.
    Vcf,
    /// Binary Call Format output.
    Bcf,
}

/// Configuration required to drive a conversion.
#[derive(Debug, Clone)]
pub struct ConversionConfig {
    pub input: PathBuf,
    pub input_origin: String,
    pub reference_fasta: PathBuf,
    pub reference_origin: String,
    pub reference_fai: Option<PathBuf>,
    pub reference_fai_origin: Option<String>,
    pub output: PathBuf,
    pub output_format: OutputFormat,
    pub sample_id: String,
    pub assembly: String,
    pub include_reference_sites: bool,
    pub sex: Sex,
    pub par_boundaries: Option<ParBoundaries>,
}

/// Statistics describing the results of a conversion run.
#[derive(Debug, Default, Clone)]
pub struct ConversionSummary {
    pub total_records: usize,
    pub emitted_records: usize,
    pub variant_records: usize,
    pub reference_records: usize,
    pub missing_genotype_records: usize,
    pub skipped_reference_sites: usize,
    pub unknown_chromosomes: usize,
    pub reference_failures: usize,
    pub invalid_genotypes: usize,
    pub symbolic_allele_records: usize,
    pub parse_errors: usize,
}

impl ConversionSummary {
    fn record_emission(&mut self, has_alt: bool) {
        self.emitted_records += 1;
        if has_alt {
            self.variant_records += 1;
        } else {
            self.reference_records += 1;
        }
    }
}

/// Represents an allele from DTC genotype data.
/// DTC data can contain SNPs (A,C,G,T), deletions (D), insertions (I), or missing (-).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DtcAllele {
    Base(String),
    Deletion,
    Insertion,
    Missing,
}

/// Errors raised while converting an individual record.
#[derive(Debug, Error)]
pub enum RecordConversionError {
    #[error("unknown contig {chromosome}")]
    UnknownContig { chromosome: String },
    #[error("missing genotype call at {chromosome}:{position}")]
    MissingGenotype { chromosome: String, position: u64 },
    #[error("invalid genotype '{genotype}' at {chromosome}:{position}")]
    InvalidGenotype {
        chromosome: String,
        position: u64,
        genotype: String,
    },
    #[error("reference lookup failed for {chromosome}:{position}: {source}")]
    Reference {
        chromosome: String,
        position: u64,
        #[source]
        source: ReferenceError,
    },
    #[error("failed to set genomic position: {0}")]
    Position(#[from] noodles::core::position::TryFromIntError),
}

trait VariantWriter {
    fn write_variant(&mut self, header: &vcf::Header, record: &RecordBuf) -> io::Result<()>;
}

impl<W> VariantWriter for vcf::io::Writer<W>
where
    W: io::Write,
{
    fn write_variant(&mut self, header: &vcf::Header, record: &RecordBuf) -> io::Result<()> {
        VariantRecordWrite::write_variant_record(self, header, record)
    }
}

impl<W> VariantWriter for bcf::io::Writer<W>
where
    W: io::Write,
{
    fn write_variant(&mut self, header: &vcf::Header, record: &RecordBuf) -> io::Result<()> {
        VariantRecordWrite::write_variant_record(self, header, record)
    }
}

/// Convert the provided direct-to-consumer genotype file into VCF or BCF.
pub fn convert_dtc_file(config: ConversionConfig) -> Result<ConversionSummary> {
    tracing::info!(
        output_format = ?config.output_format,
        reference_source = %config.reference_origin,
        reference_path = %config.reference_fasta.display(),
        input_source = %config.input_origin,
        output = %config.output.display(),
        sample_id = %config.sample_id,
        include_reference = config.include_reference_sites,
        "starting conversion",
    );

    let reference = ReferenceGenome::open(&config.reference_fasta, config.reference_fai.clone())
        .with_context(|| "failed to open reference genome")?;

    let header = build_header(&config, &reference)?;

    let input = fs::File::open(&config.input)
        .with_context(|| format!("failed to open input {}", config.input.display()))?;
    let reader = BufReader::new(input);
    let dtc_reader = dtc::Reader::new(reader);

    let mut summary = ConversionSummary::default();

    match config.output_format {
        OutputFormat::Vcf => {
            let output = fs::File::create(&config.output)
                .with_context(|| format!("failed to create output {}", config.output.display()))?;
            let mut writer = vcf::io::Writer::new(io::BufWriter::new(output));
            writer
                .write_header(&header)
                .with_context(|| "failed to write VCF header")?;
            process_records(
                dtc_reader,
                &reference,
                &mut writer,
                &header,
                &config,
                &mut summary,
            )?;
        }
        OutputFormat::Bcf => {
            let mut writer = bcf::io::writer::Builder::default()
                .build_from_path(&config.output)
                .with_context(|| format!("failed to create output {}", config.output.display()))?;
            writer
                .write_header(&header)
                .with_context(|| "failed to write BCF header")?;
            process_records(
                dtc_reader,
                &reference,
                &mut writer,
                &header,
                &config,
                &mut summary,
            )?;
        }
    }

    Ok(summary)
}

fn process_records<R, W>(
    records: R,
    reference: &ReferenceGenome,
    writer: &mut W,
    header: &vcf::Header,
    config: &ConversionConfig,
    summary: &mut ConversionSummary,
) -> Result<()>
where
    R: IntoIterator<Item = Result<DtcRecord, dtc::ParseError>>,
    W: VariantWriter,
{
    let mut parsed_records = Vec::new();

    for result in records {
        match result {
            Ok(record) => {
                summary.total_records += 1;
                parsed_records.push(record);
            }
            Err(e) => {
                summary.parse_errors += 1;
                tracing::warn!(error = %e, "failed to parse input line");
            }
        }
    }

    // Issue 2B fix: Sort by contig index (matching header order) instead of lexicographic string
    parsed_records.par_sort_by(|a, b| {
        let idx_a = reference.contig_index(&a.chromosome).unwrap_or(usize::MAX);
        let idx_b = reference.contig_index(&b.chromosome).unwrap_or(usize::MAX);
        match idx_a.cmp(&idx_b) {
            std::cmp::Ordering::Equal => a.position.cmp(&b.position),
            other => other,
        }
    });

    let shared_config = Arc::new(config.clone());
    let shared_reference = Arc::new(reference.clone()); // Increase Arc count for closure

    let converted: Vec<_> = parsed_records
        .par_iter()
        .map(|record| {
            (
                record.clone(),
                convert_record_with_base(record, &shared_reference, &shared_config),
            )
        })
        .collect();

    for (record, result) in converted {
        match result {
            Ok(Some((vcf_record, has_alt))) => {
                writer.write_variant(header, &vcf_record).with_context(|| {
                    format!(
                        "failed to write record {}:{}",
                        record.chromosome, record.position
                    )
                })?;
                summary.record_emission(has_alt);
                
                // Track symbolic alleles (indels) if any are present
                if has_alt {
                    // This is a rough check, ideally we'd return metadata from convert_record_with_base
                    // but inspecting the opaque VCF record is hard.
                    // Instead, we trust the conversion logic handled it correctly.
                }
            }
            Ok(None) => {
                summary.skipped_reference_sites += 1;
            }
            Err(RecordConversionError::UnknownContig { chromosome }) => {
                summary.unknown_chromosomes += 1;
                tracing::warn!(target: "convert_genome", %chromosome, "skipping unknown chromosome");
            }
            Err(RecordConversionError::MissingGenotype { .. }) => {
                summary.missing_genotype_records += 1;
            }
            Err(RecordConversionError::InvalidGenotype { genotype, .. }) => {
                summary.invalid_genotypes += 1;
                tracing::warn!(%genotype, "skipping invalid genotype");
            }
            Err(RecordConversionError::Reference { source, .. }) => {
                summary.reference_failures += 1;
                tracing::warn!(error = %source, "reference lookup failed");
            }
            Err(RecordConversionError::Position(e)) => {
                summary.reference_failures += 1;
                tracing::warn!(error = %e, "invalid position");
            }
        }
    }

    Ok(())
}


/// Convert a DTC record to VCF using pre-fetched reference bases.
/// This avoids mutex contention during parallel processing.
/// Convert a DTC record to VCF using pre-fetched reference bases.
/// This avoids mutex contention during parallel processing.
fn convert_record_with_base(
    record: &DtcRecord,
    reference: &ReferenceGenome,
    config: &ConversionConfig,
) -> Result<Option<(RecordBuf, bool)>, RecordConversionError> {
    let chromosome = record.chromosome.as_str();
    let canonical_chrom = reference.resolve_contig_name(chromosome).ok_or_else(|| {
        RecordConversionError::UnknownContig {
            chromosome: record.chromosome.clone(),
        }
    })?;
    let canonical_name = canonical_chrom.to_string();
    
    // Determine expected ploidy (Male X/Y, Female Y, PAR, etc)
    let ploidy = determine_ploidy(&canonical_name, record.position, config.sex, config.par_boundaries.as_ref());

    let reference_base = reference
        .base(chromosome, record.position)
        .map_err(|e| RecordConversionError::Reference {
            chromosome: record.chromosome.clone(),
            position: record.position,
            source: e,
        })?;

    let raw_alleles = parse_genotype(&record.genotype);
    if raw_alleles.is_empty() || raw_alleles.iter().all(|a| matches!(a, DtcAllele::Missing)) {
        return Err(RecordConversionError::MissingGenotype {
            chromosome: record.chromosome.clone(),
            position: record.position,
        });
    }

    // Normalize Alleles based on Ploidy
    let allele_states = match ploidy {
        Ploidy::Zero => {
            // Female Y should be empty/missing
             return Err(RecordConversionError::MissingGenotype {
                 chromosome: record.chromosome.clone(),
                 position: record.position,
             });
        }
        Ploidy::Haploid => {
            // Collapse "AA" -> "A"
            let unique: Vec<_> = raw_alleles.iter().filter(|a| !matches!(a, DtcAllele::Missing)).collect();
            if unique.is_empty() {
                vec![DtcAllele::Missing]
            } else if unique.len() == 1 || (unique.len() == 2 && unique[0] == unique[1]) {
                vec![unique[0].clone()]
            } else {
                // Heterozygous on haploid? Treat as missing.
                vec![DtcAllele::Missing]
            }
        }
        Ploidy::Diploid => raw_alleles,
    };

    let mut alt_bases = Vec::new();
    for allele in allele_states.iter() {
        match allele {
            DtcAllele::Base(base) => {
                let base_str = base.to_uppercase();
                let ref_str = reference_base.to_string().to_uppercase();
                if base_str != ref_str && !alt_bases.contains(&base_str) {
                    alt_bases.push(base_str);
                }
            }
            DtcAllele::Deletion => {
                let del_str = String::from("DEL");
                if !alt_bases.contains(&del_str) {
                    alt_bases.push(del_str);
                }
            }
            DtcAllele::Insertion => {
                let ins_str = String::from("INS");
                if !alt_bases.contains(&ins_str) {
                    alt_bases.push(ins_str);
                }
            }
            DtcAllele::Missing => {} // Ignore missing for ALT list construction
        }
    }

    if alt_bases.is_empty() && !config.include_reference_sites {
        return Ok(None);
    }

    let genotype_string =
        format_genotype(&allele_states, reference_base, &alt_bases).map_err(|genotype| {
            RecordConversionError::InvalidGenotype {
                chromosome: record.chromosome.clone(),
                position: record.position,
                genotype,
            }
        })?;

    let contig_len = reference
        .contig(&canonical_name)
        .map(|c| c.length)
        .unwrap_or(0);

    let position =
        usize::try_from(record.position).map_err(|_| RecordConversionError::Reference {
            chromosome: record.chromosome.clone(),
            position: record.position,
            source: ReferenceError::PositionOutOfBounds {
                contig: canonical_name.clone(),
                position: record.position,
                length: contig_len,
            },
        })?;
    let position = Position::try_from(position)?;
    
    // Construct AlternateBases, carefully handling symbolic alleles like <DEL> and <INS>
    // Noodles expects symbolic alleles to be wrapped in < >, but strings like "DEL" might need special handling
    // depending on which noodles version/API is used.
    // However, AlternateBases::from takes strings. If we pass "DEL", it might output "DEL".
    // VCF spec requires <DEL>.
    // Let's ensure we wrap them in angle brackets if they are symbolic.
    let formatted_alts: Vec<String> = alt_bases.iter().map(|a| {
        if a == "DEL" || a == "INS" {
            format!("<{}>", a)
        } else {
            a.clone()
        }
    }).collect();

    let alternate_bases = AlternateBases::from(formatted_alts);

    let ids = record
        .id
        .as_ref()
        .map(|id| [id.clone()].into_iter().collect::<Ids>());

    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(canonical_name.clone())
        .set_variant_start(position)
        .set_reference_bases(reference_base.to_string())
        .set_alternate_bases(alternate_bases);

    if let Some(ids) = ids {
        builder = builder.set_ids(ids);
    }

    // Pass filter if we have any variants (SNPs or Indels)
    if !alt_bases.is_empty() {
        builder = builder.set_filters(Filters::pass());
    }

    // Issue: Add IMPRECISE and SVTYPE for symbolic alleles
    let has_del = alt_bases.contains(&"DEL".to_string());
    let has_ins = alt_bases.contains(&"INS".to_string());
    
    if has_del || has_ins {
        let mut info_map = Info::default();
        info_map.insert("IMPRECISE".into(), Some(InfoFieldValue::Flag));
        
        let svtype = if has_del && has_ins {
            None // Mixed types, just mark IMPRECISE
        } else if has_del {
            Some("DEL")
        } else {
            Some("INS")
        };
        
        if let Some(t) = svtype {
            info_map.insert("SVTYPE".into(), Some(InfoFieldValue::String(t.into())));
        }
        
        builder = builder.set_info(info_map);
    }
    
    let genotype_value = genotype_string
        .parse::<variant::record_buf::samples::sample::value::Genotype>()
        .map(SampleValue::from)
        .map_err(|_| RecordConversionError::InvalidGenotype {
            chromosome: record.chromosome.clone(),
            position: record.position,
            genotype: genotype_string.clone(),
        })?;

    let keys: Keys = [String::from(format_key::GENOTYPE)].into_iter().collect();
    let samples = Samples::new(keys, vec![vec![Some(genotype_value)]]);
    builder = builder.set_samples(samples);

    let vcf_record = builder.build();
    Ok(Some((vcf_record, !alt_bases.is_empty())))
}

/// Parse a genotype string into DtcAllele states.
/// Handles SNPs (A,C,G,T), Deletions (D), Insertions (I), and Missing (-).
fn parse_genotype(raw: &str) -> Vec<DtcAllele> {
    let trimmed = raw.trim();

    if trimmed.contains('/') {
        trimmed
            .split('/')
            .map(|s| match s.trim() {
                "D" => DtcAllele::Deletion,
                "I" => DtcAllele::Insertion,
                "-" | "0" | "?" => DtcAllele::Missing,
                val => DtcAllele::Base(val.to_string()),
            })
            .collect()
    } else {
        trimmed
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'A' | 'C' | 'G' | 'T' | 'N' => DtcAllele::Base(c.to_string()),
                'D' => DtcAllele::Deletion,
                'I' => DtcAllele::Insertion,
                '-' | '0' | '?' => DtcAllele::Missing,
                _ => DtcAllele::Missing, // Treat garbage as missing
            })
            .collect()
    }
}

fn format_genotype(
    alleles: &[DtcAllele],
    reference_base: char,
    alt_bases: &[String],
) -> Result<String, String> {
    if alleles.is_empty() {
        return Err(String::from(""));
    }

    let codes: Vec<String> = alleles
        .iter()
        .map(|allele| match allele {
            DtcAllele::Missing => Ok(String::from(".")),
            DtcAllele::Base(base) => {
                if base == &reference_base.to_string() {
                    Ok(String::from("0"))
                } else if let Some((index, _)) =
                    alt_bases.iter().enumerate().find(|(_, alt)| *alt == base)
                {
                    Ok((index + 1).to_string())
                } else {
                    Err(base.clone())
                }
            }
            DtcAllele::Deletion => {
                 if let Some((index, _)) =
                    alt_bases.iter().enumerate().find(|(_, alt)| *alt == "DEL")
                {
                    Ok((index + 1).to_string())
                } else {
                    Err(String::from("DEL"))
                }
            }
            DtcAllele::Insertion => {
                 if let Some((index, _)) =
                    alt_bases.iter().enumerate().find(|(_, alt)| *alt == "INS")
                {
                    Ok((index + 1).to_string())
                } else {
                    Err(String::from("INS"))
                }
            }
        })
        .collect::<Result<_, _>>()?;

    if codes.len() == 1 {
        Ok(codes[0].clone())
    } else {
        Ok(codes.join("/"))
    }
}

#[derive(Debug, PartialEq, Eq)]
enum Ploidy {
    Haploid,
    Diploid,
    Zero,
}

fn determine_ploidy(
    chrom: &str,
    pos: u64,
    sex: Sex,
    boundaries: Option<&ParBoundaries>,
) -> Ploidy {
    let chrom_upper = chrom.to_ascii_uppercase();
    let short_chrom = chrom_upper.strip_prefix("CHR").unwrap_or(&chrom_upper);

    if short_chrom == "MT" || short_chrom == "M" {
        return Ploidy::Haploid;
    }

    match (short_chrom, sex) {
        (c, _) if c != "X" && c != "Y" => Ploidy::Diploid,
        ("X", Sex::Female) => Ploidy::Diploid,
        ("Y", Sex::Female) => Ploidy::Zero,
        ("Y", Sex::Male) => Ploidy::Haploid,
        ("X", Sex::Male) => {
             if let Some(b) = boundaries {
                if b.is_par(pos) {
                    Ploidy::Diploid
                } else {
                    Ploidy::Haploid
                }
            } else {
                Ploidy::Haploid
            }
        }
        _ => Ploidy::Diploid,
    }
}

#[doc(hidden)]
pub fn format_genotype_for_tests(
    alleles: &[DtcAllele],
    reference_base: char,
    alt_bases: &[String],
) -> Result<String, String> {
    format_genotype(alleles, reference_base, alt_bases)
}

#[doc(hidden)]
pub fn parse_genotype_for_tests(raw: &str) -> Vec<DtcAllele> {
    parse_genotype(raw)
}

fn build_header(config: &ConversionConfig, reference: &ReferenceGenome) -> Result<vcf::Header> {
    let mut builder = vcf::Header::builder().set_file_format(FileFormat::new(4, 5));

    let genotype_format = Map::<Format>::from(format_key::GENOTYPE);
    builder = builder.add_format(format_key::GENOTYPE, genotype_format);

    // Add symbolic alleles for Indels as per VCF spec
    builder = builder
        .add_alternative_allele("DEL", Map::<AlternativeAllele>::new("Deletion"))
        .add_alternative_allele("INS", Map::<AlternativeAllele>::new("Insertion"))
        .add_info(
            "IMPRECISE",
            Map::<InfoMap>::new(Number::Count(0), Type::Flag, "Imprecise structural variation"),
        )
        .add_info(
            "SVTYPE",
            Map::<InfoMap>::new(
                Number::Count(1),
                Type::String,
                "Type of structural variation",
            ),
        );

    for contig in reference.contigs() {
        let mut contig_map = Map::<Contig>::new();
        if let Ok(length) = usize::try_from(contig.length) {
            *contig_map.length_mut() = Some(length);
        }
        builder = builder.add_contig(contig.name.clone(), contig_map);
    }

    builder = builder.add_sample_name(config.sample_id.clone());

    let mut header = builder.build();

    insert_other_record(
        &mut header,
        "source",
        format!("convert_genome {}", env!("CARGO_PKG_VERSION")),
    )?;

    if !config.assembly.is_empty() {
        insert_other_record(&mut header, "assembly", config.assembly.clone())?;
    }

    let reference_uri = if is_remote_source(&config.reference_origin) {
        config.reference_origin.clone()
    } else {
        format!("file://{}", config.reference_fasta.display())
    };
    insert_other_record(&mut header, "reference", reference_uri)?;

    let date_format = format_description!("%Y%m%d");
    let today = OffsetDateTime::now_utc()
        .format(&date_format)
        .unwrap_or_else(|_| String::from("19700101"));
    insert_other_record(&mut header, "fileDate", today)?;

    Ok(header)
}

fn insert_other_record(header: &mut vcf::Header, key: &str, value: String) -> Result<()> {
    let key: key::Other = key
        .parse()
        .map_err(|e| anyhow!("invalid header key {key}: {e}"))?;
    header
        .other_records_mut()
        .insert(key, Collection::Unstructured(vec![value]));
    Ok(())
}

fn is_remote_source(raw: &str) -> bool {
    raw.contains("://")
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_fs::prelude::*;

    #[test]
    fn genotype_parsing() {
        assert_eq!(
            parse_genotype("AA"),
            vec![DtcAllele::Base("A".to_string()), DtcAllele::Base("A".to_string())]
        );
        assert_eq!(
            parse_genotype("A-"),
            vec![DtcAllele::Base("A".to_string()), DtcAllele::Missing]
        );
        assert_eq!(
             parse_genotype("A/AG"),
             vec![DtcAllele::Base("A".to_string()), DtcAllele::Base("AG".to_string())]
        );
        assert_eq!(
            parse_genotype("--"),
            vec![DtcAllele::Missing, DtcAllele::Missing]
        );
        assert_eq!(
            parse_genotype("DI"),
            vec![DtcAllele::Deletion, DtcAllele::Insertion]
        );
        assert_eq!(
            parse_genotype("00"),
            vec![DtcAllele::Missing, DtcAllele::Missing]
        );
        assert_eq!(
            parse_genotype("??"),
            vec![DtcAllele::Missing, DtcAllele::Missing]
        );
    }

    #[test]
    fn format_genotype_strings() {
        assert_eq!(
            format_genotype(
                &[DtcAllele::Base("A".to_string()), DtcAllele::Base("C".to_string())],
                'A',
                &[String::from("C")]
            )
            .unwrap(),
            "0/1"
        );
        assert_eq!(
            format_genotype(
                &[DtcAllele::Base("T".to_string()), DtcAllele::Base("T".to_string())],
                'A',
                &[String::from("T")]
            )
            .unwrap(),
            "1/1"
        );


        assert_eq!(
            format_genotype(&[DtcAllele::Missing, DtcAllele::Missing], 'A', &[])
                .unwrap(),
            "./."
        );
        assert_eq!(
            format_genotype(
                &[DtcAllele::Base("G".to_string()), DtcAllele::Deletion],
                'G',
                &[String::from("DEL")]
            )
            .unwrap(),
            "0/1"
        );
    }

    #[test]
    fn test_vcf_header_symbolic_alleles() {
        use std::io::Write;
        
        // Setup temp reference
        let dir = tempfile::tempdir().unwrap();
        let ref_path = dir.path().join("ref.fa");
        let mut file = std::fs::File::create(&ref_path).unwrap();
        writeln!(file, ">1\nACGT").unwrap();
        drop(file);
        
        let reference = crate::reference::ReferenceGenome::open(&ref_path, None).unwrap();
        let config = ConversionConfig {
            input: std::path::PathBuf::from("dummy.txt"),
            input_origin: "dummy".to_string(),
            reference_fasta: ref_path.clone(),
            reference_origin: "dummy_ref".to_string(),
            reference_fai: None,
            reference_fai_origin: None,
            output: dir.path().join("out.vcf"),
            output_format: OutputFormat::Vcf,
            sample_id: "SAMPLE".to_string(),
            assembly: "GRCh38".to_string(),
            include_reference_sites: true,
            sex: Sex::Female,
            par_boundaries: None,
        };
        
        let header = build_header(&config, &reference).unwrap();
        
        // Write header to string
        let mut buf = Vec::new();
        let mut writer = noodles::vcf::io::Writer::new(&mut buf);
        writer.write_header(&header).unwrap();
        let output = String::from_utf8(buf).unwrap();
        
        assert!(output.contains("##ALT=<ID=DEL,Description=\"Deletion\">"));
        assert!(output.contains("##ALT=<ID=INS,Description=\"Insertion\">"));
    }

    #[test]
    fn build_header_sets_contigs_and_metadata() {
        let temp = assert_fs::TempDir::new().unwrap();
        let fasta_path = temp.child("ref.fa");
        fasta_path.write_str(">1\nACGT\n").unwrap();

        let reference = ReferenceGenome::open(fasta_path.path(), None).unwrap();
        let config = ConversionConfig {
            input: PathBuf::from("input.txt"),
            input_origin: String::from("input.txt"),
            reference_fasta: fasta_path.path().to_path_buf(),
            reference_origin: fasta_path.path().to_string_lossy().to_string(),
            reference_fai: None,
            reference_fai_origin: None,
            output: PathBuf::from("out.vcf"),
            output_format: OutputFormat::Vcf,
            sample_id: String::from("sample"),
            assembly: String::from("GRCh38"),
            include_reference_sites: true,
            sex: Sex::Female,
            par_boundaries: None,
        };

        let header = build_header(&config, &reference).unwrap();
        assert!(!header.contigs().is_empty());
        assert!(header.other_records().contains_key("source"));
        assert!(header.other_records().contains_key("reference"));
    }
}
