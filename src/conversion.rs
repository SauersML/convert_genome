use crate::cli::Sex;
use crate::reference::ParBoundaries;
use std::{
    fs,
    io::{self, BufReader},
    path::PathBuf,
};

use anyhow::{Context, Result, anyhow};
use clap::ValueEnum;
use noodles::bcf;
use noodles::vcf::{
    self,
    header::{
        FileFormat,
        record::{
            key,
            value::{
                Collection, Map,
                map::{
                    AlternativeAllele, Contig, Format, Info as InfoMap,
                    info::{Number, Type},
                },
            },
        },
    },
    variant::{
        io::Write as VariantRecordWrite, record::samples::keys::key as format_key,
        record_buf::RecordBuf,
    },
};
// rayon removed

use thiserror::Error;
use time::{OffsetDateTime, macros::format_description};

use crate::{
    ConversionSummary,
    dtc::{self, Allele as DtcAllele, parse_genotype},
    plink::PlinkWriter,
    reference::{ReferenceError, ReferenceGenome},
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
    pub standardize: bool,
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
    fn write_variant(&mut self, header: &vcf::Header, record: &RecordBuf) -> io::Result<()> {
        // Header not used for PLINK format but required by trait
        let _ = header;
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
    let source: Box<dyn crate::input::VariantSource> = match config.input_format {
        crate::input::InputFormat::Dtc => {
            let input = fs::File::open(&config.input)
                .with_context(|| format!("failed to open input {}", config.input.display()))?;
            let reader = BufReader::new(input);
            let dtc_reader = dtc::Reader::new(reader);
            let source =
                crate::input::DtcSource::new(dtc_reader, reference.clone(), config.clone());
            Box::new(source)
        }
        crate::input::InputFormat::Vcf => {
            let input = fs::File::open(&config.input)
                .with_context(|| format!("failed to open input {}", config.input.display()))?;
            let reader = BufReader::new(input);
            let vcf_reader = vcf::io::Reader::new(reader);
            let source = crate::input::VcfSource::new(vcf_reader, &reference)
                .with_context(|| "failed to initialize VCF source")?;
            Box::new(source)
        }
        crate::input::InputFormat::Bcf => {
            let input = fs::File::open(&config.input)
                .with_context(|| format!("failed to open input {}", config.input.display()))?;
            let bcf_reader = bcf::io::Reader::new(input);
            let source = crate::input::BcfSource::new(bcf_reader, &reference)
                .with_context(|| "failed to initialize BCF source")?;
            Box::new(source)
        }
        crate::input::InputFormat::Auto => {
            return Err(anyhow!("Auto format detection failed or not resolved"));
        }
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
            let mut writer =
                PlinkWriter::new(&config.output).context("failed to create PLINK writer")?;

            writer
                .write_fam(&config.sample_id, config.sex)
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
    reference: &ReferenceGenome,
    writer: &mut W,
    header: &vcf::Header,
    config: &ConversionConfig,
    summary: &mut crate::ConversionSummary,
) -> Result<()>
where
    S: crate::input::VariantSource,
    W: VariantWriter,
{
    while let Some(result) = source.next_variant(summary) {
        match result {
            Ok(record) => {
                // Apply standardization if requested
                let final_record = if config.standardize {
                    match standardize_record(&record, reference, config) {
                        Ok(Some(standardized)) => standardized,
                        Ok(None) => continue, // Filtered out
                        Err(e) => {
                            tracing::warn!(error = %e, "failed to standardize record, skipping");
                            summary.reference_failures += 1;
                            continue;
                        }
                    }
                } else {
                    record
                };

                summary.record_emission(!final_record.alternate_bases().as_ref().is_empty());

                writer
                    .write_variant(header, &final_record)
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

/// Standardize a VCF record by:
/// 1. Normalizing chromosome name to canonical form
/// 2. Polarizing alleles against reference genome (swap REF/ALT if needed)
/// 3. Generating synthetic ID if missing
///
/// Returns the standardized record, or None if the record should be skipped.
pub fn standardize_record(
    record: &RecordBuf,
    reference: &ReferenceGenome,
    config: &ConversionConfig,
) -> Result<Option<RecordBuf>, RecordConversionError> {
    use noodles::vcf::variant::record_buf::{AlternateBases, Ids};

    let chrom = record.reference_sequence_name();
    let pos = record
        .variant_start()
        .map(|p| usize::from(p) as u64)
        .unwrap_or(0);

    // 1. Normalize chromosome name
    let canonical_name = match reference.resolve_contig_name(chrom) {
        Some(name) => name.to_string(),
        None => {
            tracing::warn!("unknown chromosome: {}", chrom);
            return Ok(None);
        }
    };

    // 2. Get reference base at this position
    let ref_base = match reference.base(&canonical_name, pos) {
        Ok(base) => base.to_ascii_uppercase(),
        Err(e) => {
            tracing::warn!(
                "reference lookup failed at {}:{}: {}",
                canonical_name,
                pos,
                e
            );
            return Err(RecordConversionError::Reference {
                chromosome: canonical_name,
                position: pos,
                source: e,
            });
        }
    };

    let input_ref = record.reference_bases().to_uppercase();
    let ref_base_str = ref_base.to_string();

    // 3. Check if allele polarization is needed
    let (final_ref, final_alts, needs_flip) = if input_ref.len() == 1 && input_ref != ref_base_str {
        // Single-base REF doesn't match reference - need to polarize
        let alt_bases: Vec<String> = record
            .alternate_bases()
            .as_ref()
            .iter()
            .map(|s| s.to_string())
            .collect();

        // Check if reference base is in ALTs
        if let Some(flip_idx) = alt_bases.iter().position(|a| a == &ref_base_str) {
            // Swap: new REF = ref_base, new ALTs = [old_ref] + (old_alts - ref_base)
            let mut new_alts = vec![input_ref.clone()];
            for (i, alt) in alt_bases.iter().enumerate() {
                if i != flip_idx {
                    new_alts.push(alt.clone());
                }
            }
            (ref_base_str.clone(), new_alts, Some(flip_idx + 1)) // +1 because ALT indices are 1-based
        } else {
            // Cannot polarize - REF/ALT don't contain reference base
            tracing::warn!(
                "cannot polarize alleles at {}:{}: REF={} but reference={}, ALTs={:?}",
                canonical_name,
                pos,
                input_ref,
                ref_base_str,
                alt_bases
            );
            (input_ref.clone(), alt_bases, None)
        }
    } else {
        // REF matches or is multi-base (indel) - keep as-is
        let alt_bases: Vec<String> = record
            .alternate_bases()
            .as_ref()
            .iter()
            .map(|s| s.to_string())
            .collect();
        (input_ref.clone(), alt_bases, None)
    };

    // 4. Generate synthetic ID if missing
    let ids = if record.ids().as_ref().is_empty() {
        let alt_str = if final_alts.is_empty() {
            ".".to_string()
        } else {
            final_alts.join(",")
        };
        let synthetic_id = format!("{}:{}:{}:{}", canonical_name, pos, final_ref, alt_str);
        Ids::from_iter(vec![synthetic_id])
    } else {
        record.ids().clone()
    };

    // 5. Apply GT flipping if allele polarization occurred
    let samples = if needs_flip.is_some() {
        tracing::debug!(
            chrom = %canonical_name, pos = pos,
            "Allele polarization applied, flipping GT indices"
        );
        flip_sample_genotypes(record.samples())
    } else {
        record.samples().clone()
    };

    // 6. Check ploidy and enforce if needed
    let ploidy = determine_ploidy(
        &canonical_name,
        pos,
        config.sex,
        config.par_boundaries.as_ref(),
    );
    if ploidy == Ploidy::Zero {
        // Skip this record for this sex
        return Ok(None);
    }
    // Note: Haploid enforcement is already handled by the source for DTC files.
    // For VCF/BCF standardization, we preserve the original ploidy.

    // Build standardized record
    let pos_val = noodles::core::Position::new(pos as usize);
    let pos_val = match pos_val {
        Some(p) => p,
        None => return Ok(None), // Invalid position
    };

    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(canonical_name)
        .set_variant_start(pos_val)
        .set_ids(ids)
        .set_reference_bases(final_ref)
        .set_alternate_bases(AlternateBases::from(final_alts))
        .set_samples(samples);

    // Preserve quality and filters if present
    if let Some(qual) = record.quality_score() {
        builder = builder.set_quality_score(qual);
    }

    Ok(Some(builder.build()))
}

/// Flip genotype indices in all samples (0↔1 swap for biallelic sites)
/// For biallelic sites where ref/alt were swapped, we need to flip 0↔1 in GT
fn flip_sample_genotypes(
    samples: &noodles::vcf::variant::record_buf::Samples,
) -> noodles::vcf::variant::record_buf::Samples {
    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record_buf::Samples;
    use noodles::vcf::variant::record_buf::samples::sample::Value;

    let keys = samples.keys().clone();
    let mut new_values: Vec<Vec<Option<Value>>> = Vec::new();

    for sample in samples.values() {
        // Try to get the GT field and flip it
        if let Some(gt_val) = sample.get(key::GENOTYPE) {
            let mut new_sample_vals: Vec<Option<Value>> = Vec::new();

            // Preserve all fields but flip GT
            for field_val in sample.values() {
                new_sample_vals.push(field_val.clone());
            }

            // Find and flip GT (it should be the first field typically)
            if let Some(Value::String(gt_str)) = gt_val {
                // Replace the GT value with flipped version
                if !new_sample_vals.is_empty() {
                    new_sample_vals[0] = Some(Value::String(flip_gt_string(&gt_str)));
                }
            }

            new_values.push(new_sample_vals);
        } else {
            // No GT field, keep sample as-is
            new_values.push(sample.values().iter().map(|v| v.clone()).collect());
        }
    }

    Samples::new(keys, new_values)
}

/// Flip genotype indices in a GT string (e.g., "0/1" -> "1/0", "0|0" -> "1|1")
fn flip_gt_string(gt: &str) -> String {
    gt.chars()
        .map(|c| match c {
            '0' => '1',
            '1' => '0',
            other => other, // Preserve separators (/, |) and other indices (2, 3, .)
        })
        .collect()
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
            Map::<InfoMap>::new(
                Number::Count(0),
                Type::Flag,
                "Imprecise structural variation",
            ),
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
            vec![
                DtcAllele::Base("A".to_string()),
                DtcAllele::Base("A".to_string())
            ]
        );
        assert_eq!(
            parse_genotype("A-"),
            vec![DtcAllele::Base("A".to_string()), DtcAllele::Missing]
        );
        assert_eq!(
            parse_genotype("A/AG"),
            vec![
                DtcAllele::Base("A".to_string()),
                DtcAllele::Base("AG".to_string())
            ]
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
                &[
                    DtcAllele::Base("A".to_string()),
                    DtcAllele::Base("C".to_string())
                ],
                'A',
                &[String::from("C")]
            )
            .unwrap(),
            "0/1"
        );
        assert_eq!(
            format_genotype(
                &[
                    DtcAllele::Base("T".to_string()),
                    DtcAllele::Base("T".to_string())
                ],
                'A',
                &[String::from("T")]
            )
            .unwrap(),
            "1/1"
        );

        assert_eq!(
            format_genotype(&[DtcAllele::Missing, DtcAllele::Missing], 'A', &[]).unwrap(),
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
            standardize: false,
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
            standardize: false,
        };

        let header = build_header(&config, &reference).unwrap();
        assert!(!header.contigs().is_empty());
        assert!(header.other_records().contains_key("source"));
        assert!(header.other_records().contains_key("reference"));
    }
}
