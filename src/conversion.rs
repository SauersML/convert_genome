use std::{
    fs,
    io::{self, BufReader},
    path::PathBuf,
};
use crate::reference::ParBoundaries;
use crate::cli::Sex;

use anyhow::{Context, Result, anyhow};
use clap::ValueEnum;
use noodles::bcf;
use noodles::{
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
            io::Write as VariantRecordWrite,
            record::samples::keys::key as format_key,
            record_buf::RecordBuf,
        },
    },
};
// rayon removed

use thiserror::Error;
use time::{OffsetDateTime, macros::format_description};

use crate::{
    dtc::{self, Allele as DtcAllele, parse_genotype},
    reference::{ReferenceError, ReferenceGenome},
    plink::PlinkWriter,
    ConversionSummary,
};

/// Supported output formats for the converter.
#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum OutputFormat {
    /// Variant Call Format text output.
    Vcf,
    /// Binary Call Format output.
    Bcf,
    /// PLINK 1.9 binary output (.bed, .bim, .fam).
    Plink,
}

/// Configuration required to drive a conversion.
#[derive(Debug, Clone)]
pub struct ConversionConfig {
    pub input: PathBuf,
    pub input_format: crate::input::InputFormat,
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

// ConversionSummary moved to crate root


// DtcAllele moved to dtc.rs

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

impl VariantWriter for PlinkWriter {
    fn write_variant(&mut self, _header: &vcf::Header, record: &RecordBuf) -> io::Result<()> {
        self.write_variant(record)
    }
}

/// Convert the provided direct-to-consumer genotype file into VCF or BCF.
/// Convert the provided input file into VCF, BCF, or PLINK.
pub fn convert_dtc_file(config: ConversionConfig) -> Result<ConversionSummary> {
    tracing::info!(
        input_format = ?config.input_format,
        output_format = ?config.output_format,
        reference = %config.reference_fasta.display(),
        input = %config.input.display(),
        output = %config.output.display(),
        sample_id = %config.sample_id,
        "starting conversion",
    );

    let reference = ReferenceGenome::open(&config.reference_fasta, config.reference_fai.clone())
        .with_context(|| "failed to open reference genome")?;

    let header = build_header(&config, &reference)?;

    // Instantiate Source Iterator
    let mut source: Box<dyn crate::input::VariantSource> = match config.input_format {
        crate::input::InputFormat::Dtc => {
            let input = fs::File::open(&config.input)
                .with_context(|| format!("failed to open input {}", config.input.display()))?;
            let reader = BufReader::new(input);
            let dtc_reader = dtc::Reader::new(reader);
            let source = crate::input::DtcSource::new(dtc_reader, reference.clone(), config.clone());
            Box::new(source)
        },
        crate::input::InputFormat::Vcf | crate::input::InputFormat::Bcf => {
             // TODO: Implement VcfSource and BcfSource
             return Err(anyhow!("VCF/BCF input not yet implemented"));
        }
        crate::input::InputFormat::Auto => return Err(anyhow!("Auto format detection failed or not resolved")),
    };

    let mut summary = crate::ConversionSummary::default();

    match config.output_format {
        OutputFormat::Vcf => {
            let output = fs::File::create(&config.output)
                .with_context(|| format!("failed to create output {}", config.output.display()))?;
            let mut writer = vcf::io::Writer::new(io::BufWriter::new(output));
            writer
                .write_header(&header)
                .with_context(|| "failed to write VCF header")?;
            process_records(
                source,
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
                source,
                &reference,
                &mut writer,
                &header,
                &config,
                &mut summary,
            )?;
        }
        OutputFormat::Plink => {
            let mut writer = PlinkWriter::new(&config.output)
                .context("failed to create PLINK writer")?;

            writer.write_fam(&config.sample_id, config.sex)
                .context("failed to write FAM file")?;

            process_records(
                source,
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

fn process_records<S, W>(
    mut source: S,
    _reference: &ReferenceGenome,
    writer: &mut W,
    header: &vcf::Header,
    _config: &ConversionConfig,
    summary: &mut crate::ConversionSummary,
) -> Result<()>
where
    S: crate::input::VariantSource,
    W: VariantWriter,
{
    while let Some(result) = source.next_variant(summary) {
        match result {
            Ok(record) => {
                summary.record_emission(!record.alternate_bases().as_ref().is_empty()); 

                writer.write_variant(header, &record)
                    .context("failed to write variant record")?;
            }
            Err(e) => {
                summary.parse_errors += 1;
                // Log but continue (unless fatal?)
                // If it's IO error from source, it might be fatal.
                // But DtcSource returns io::ErrorKind::InvalidData for format errors.
                tracing::warn!(error = %e, "failed to parse/convert input record");
            }
        }
    }

    Ok(())
}

// parse_genotype moved to dtc.rs

pub fn format_genotype(
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
pub enum Ploidy {
    Haploid,
    Diploid,
    Zero,
}

pub fn determine_ploidy(
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
        ("Y", Sex::Male) => {
            if let Some(b) = boundaries {
                if b.is_par(short_chrom, pos) {
                    Ploidy::Diploid
                } else {
                    Ploidy::Haploid
                }
            } else {
                Ploidy::Haploid
            }
        }
        ("X", Sex::Male) => {
             if let Some(b) = boundaries {
                if b.is_par(short_chrom, pos) {
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
            input_format: crate::input::InputFormat::Dtc,
            input_origin: "test_input.txt".into(),
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
            input_format: crate::input::InputFormat::Dtc,
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
