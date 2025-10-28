use std::{
    fs,
    io::{self, BufReader},
    path::PathBuf,
};

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
                    map::{Contig, Format},
                },
            },
        },
        variant::{
            self as variant, RecordBuf,
            io::Write as VariantRecordWrite,
            record::samples::keys::key as format_key,
            record_buf::{
                AlternateBases, Filters, Ids,
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
    let mut records_vec = Vec::new();

    for result in records {
        match result {
            Ok(record) => {
                summary.total_records += 1;
                records_vec.push(record);
            }
            Err(e) => {
                summary.parse_errors += 1;
                tracing::warn!(error = %e, "failed to parse input line");
            }
        }
    }

    if records_vec.is_empty() {
        return Ok(());
    }

    records_vec.sort_by(|a, b| {
        let idx_a = reference.contig_index(&a.chromosome).unwrap_or(usize::MAX);
        let idx_b = reference.contig_index(&b.chromosome).unwrap_or(usize::MAX);
        idx_a
            .cmp(&idx_b)
            .then_with(|| a.position.cmp(&b.position))
            .then_with(|| a.chromosome.cmp(&b.chromosome))
    });

    let results: Vec<_> = records_vec
        .into_par_iter()
        .map_init(
            || reference.clone(),
            |local_reference, record| {
                let outcome = convert_record(&record, local_reference, config);
                (record, outcome)
            },
        )
        .collect();

    for (record, outcome) in results {
        match outcome {
            Ok(Some((vcf_record, has_alt))) => {
                writer.write_variant(header, &vcf_record).with_context(|| {
                    format!(
                        "failed to write record {}:{}",
                        record.chromosome, record.position
                    )
                })?;
                summary.record_emission(has_alt);
            }
            Ok(None) => {
                summary.skipped_reference_sites += 1;
            }
            Err(RecordConversionError::UnknownContig { chromosome }) => {
                summary.unknown_chromosomes += 1;
                tracing::warn!(
                    target: "convert_genome",
                    %chromosome,
                    "skipping unknown chromosome"
                );
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

fn convert_record(
    record: &DtcRecord,
    reference: &ReferenceGenome,
    config: &ConversionConfig,
) -> Result<Option<(RecordBuf, bool)>, RecordConversionError> {
    let chromosome = record.chromosome.as_str();
    let chromosome_owned = record.chromosome.clone();

    let canonical_chrom = reference.resolve_contig_name(chromosome).ok_or_else(|| {
        RecordConversionError::UnknownContig {
            chromosome: chromosome_owned.clone(),
        }
    })?;
    let canonical_chrom_owned = canonical_chrom.to_string();

    let reference_base = reference
        .base(chromosome, record.position)
        .map_err(|source| RecordConversionError::Reference {
            chromosome: chromosome_owned.clone(),
            position: record.position,
            source,
        })?;

    let allele_states = parse_genotype(&record.genotype);
    if allele_states.is_empty() {
        return Err(RecordConversionError::MissingGenotype {
            chromosome: chromosome_owned.clone(),
            position: record.position,
        });
    }

    let mut alt_bases = Vec::new();
    for allele in allele_states.iter().flatten() {
        let base = allele.to_ascii_uppercase();
        if base != reference_base && !alt_bases.contains(&base) {
            alt_bases.push(base);
        }
    }

    if alt_bases.is_empty() && !config.include_reference_sites {
        return Ok(None);
    }

    let genotype_string =
        format_genotype(&allele_states, reference_base, &alt_bases).map_err(|genotype| {
            RecordConversionError::InvalidGenotype {
                chromosome: chromosome_owned.clone(),
                position: record.position,
                genotype,
            }
        })?;

    let contig_len = reference
        .contig(canonical_chrom)
        .map(|c| c.length)
        .unwrap_or(0);

    let position =
        usize::try_from(record.position).map_err(|_| RecordConversionError::Reference {
            chromosome: chromosome_owned.clone(),
            position: record.position,
            source: ReferenceError::PositionOutOfBounds {
                contig: canonical_chrom_owned.clone(),
                position: record.position,
                length: contig_len,
            },
        })?;
    let position = Position::try_from(position)?;

    let alternate_bases =
        AlternateBases::from(alt_bases.iter().map(|c| c.to_string()).collect::<Vec<_>>());

    let ids = record
        .id
        .as_ref()
        .map(|id| [id.clone()].into_iter().collect::<Ids>());

    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(canonical_chrom_owned.clone())
        .set_variant_start(position)
        .set_reference_bases(reference_base.to_string())
        .set_alternate_bases(alternate_bases);

    if let Some(ids) = ids {
        builder = builder.set_ids(ids);
    }

    if !alt_bases.is_empty() {
        builder = builder.set_filters(Filters::pass());
    }

    let genotype_value = genotype_string
        .parse::<variant::record_buf::samples::sample::value::Genotype>()
        .map(SampleValue::from)
        .map_err(|_| RecordConversionError::InvalidGenotype {
            chromosome: chromosome_owned.clone(),
            position: record.position,
            genotype: genotype_string.clone(),
        })?;

    let keys: Keys = [String::from(format_key::GENOTYPE)].into_iter().collect();
    let samples = Samples::new(keys, vec![vec![Some(genotype_value)]]);
    builder = builder.set_samples(samples);

    let vcf_record = builder.build();
    Ok(Some((vcf_record, !alt_bases.is_empty())))
}

fn parse_genotype(raw: &str) -> Vec<Option<char>> {
    raw.trim()
        .chars()
        .map(|c| match c {
            '-' => None,
            other => Some(other.to_ascii_uppercase()),
        })
        .collect()
}

fn format_genotype(
    alleles: &[Option<char>],
    reference_base: char,
    alt_bases: &[char],
) -> Result<String, String> {
    if alleles.is_empty() {
        return Err(String::from(""));
    }

    let codes: Vec<String> = alleles
        .iter()
        .map(|allele| match allele {
            None => Ok(String::from(".")),
            Some(base) => {
                if *base == reference_base {
                    Ok(String::from("0"))
                } else if let Some((index, _)) =
                    alt_bases.iter().enumerate().find(|(_, alt)| *alt == base)
                {
                    Ok((index + 1).to_string())
                } else {
                    Err(base.to_string())
                }
            }
        })
        .collect::<Result<_, _>>()
        .map_err(|missing| missing)?;

    if codes.len() == 1 {
        Ok(codes[0].clone())
    } else {
        Ok(codes.join("/"))
    }
}

fn build_header(config: &ConversionConfig, reference: &ReferenceGenome) -> Result<vcf::Header> {
    let mut builder = vcf::Header::builder().set_file_format(FileFormat::new(4, 5));

    let genotype_format = Map::<Format>::from(format_key::GENOTYPE);
    builder = builder.add_format(format_key::GENOTYPE, genotype_format);

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
        assert_eq!(parse_genotype("AA"), vec![Some('A'), Some('A')]);
        assert_eq!(parse_genotype("A-"), vec![Some('A'), None]);
        assert_eq!(parse_genotype("--"), vec![None, None]);
    }

    #[test]
    fn format_genotype_strings() {
        assert_eq!(
            format_genotype(&[Some('A'), Some('C')], 'A', &['C']).unwrap(),
            "0/1"
        );
        assert_eq!(
            format_genotype(&[Some('T'), Some('T')], 'A', &['T']).unwrap(),
            "1/1"
        );
        assert_eq!(format_genotype(&[None, None], 'A', &[]).unwrap(), "./.");
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
        };

        let header = build_header(&config, &reference).unwrap();
        assert!(header.contigs().len() >= 1);
        assert!(header.other_records().contains_key("source"));
        assert!(header.other_records().contains_key("reference"));
    }
}
