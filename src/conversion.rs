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
        io::Write as VariantRecordWrite,
        record::samples::keys::key as format_key,
        record_buf::{RecordBuf, Samples, samples::sample::Value},
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
    pub reference_fasta: Option<PathBuf>,
    pub reference_origin: Option<String>,
    pub reference_fai: Option<PathBuf>,
    pub reference_fai_origin: Option<String>,
    pub output: PathBuf,
    pub output_dir: Option<PathBuf>,
    pub output_format: OutputFormat,
    pub sample_id: String,
    pub assembly: String,
    pub include_reference_sites: bool,
    pub sex: Option<Sex>,
    pub par_boundaries: Option<ParBoundaries>,
    pub standardize: bool,
    pub panel: Option<PathBuf>,
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
        self.write_variant(record)
    }
}

/// Pre-scan a DTC file to collect records for inference.
/// This reads the file once and returns records for both build and sex inference.
fn prescan_dtc_records(input: &PathBuf) -> Result<Vec<dtc::Record>> {
    let file = fs::File::open(input)
        .with_context(|| format!("failed to open input for inference: {}", input.display()))?;
    let reader = BufReader::new(file);
    let dtc_reader = dtc::Reader::new(reader);

    // Collect records (using same max_records limit for safety)
    let max_records = crate::input::get_max_records_limit();
    let mut records = Vec::new();
    let mut records_read = 0;

    for res in dtc_reader {
        if let Some(limit) = max_records
            && records_read >= limit
        {
            break;
        }
        if let Ok(rec) = res {
            records.push(rec);
            records_read += 1;
        }
    }

    Ok(records)
}

/// Convert the provided direct-to-consumer genotype file into VCF or BCF.
/// Convert the provided input file into VCF, BCF, or PLINK.
pub fn convert_dtc_file(config: ConversionConfig) -> Result<ConversionSummary> {
    tracing::info!(
        input_format = ?config.input_format,
        output_format = ?config.output_format,
        reference = ?config.reference_fasta,
        input = %config.input.display(),
        output = %config.output.display(),
        sample_id = %config.sample_id,
        panel = ?config.panel,
        standardize = config.standardize,
        "starting conversion",
    );

    // Determine if reference is required based on input format and options
    let requires_reference = matches!(config.input_format, crate::input::InputFormat::Dtc) 
        || config.standardize
        || config.panel.is_some();
    
    let reference: Option<ReferenceGenome> = if let Some(ref fasta_path) = config.reference_fasta {
        Some(ReferenceGenome::open(fasta_path, config.reference_fai.clone())
            .with_context(|| "failed to open reference genome")?)
    } else if requires_reference {
        return Err(anyhow!(
            "Reference FASTA is required for {} input or when using --standardize/--panel",
            match config.input_format {
                crate::input::InputFormat::Dtc => "DTC",
                _ => "this",
            }
        ));
    } else {
        None
    };

    // Track inference results for the report
    let mut sex_inferred = false;
    let mut build_detection: Option<crate::report::BuildDetection> = None;
    let _original_sex_provided = config.sex.is_some();

    // Auto-detect build and sex if not provided (DTC format only for now)
    let mut config = config;
    if matches!(config.input_format, crate::input::InputFormat::Dtc) {
        // Pre-scan DTC file for inference (done once, used for both)
        let prescan_records = prescan_dtc_records(&config.input)?;

        // Detect build if using default GRCh38 - check if it actually matches
        if config.assembly == "GRCh38" {
            tracing::info!("Verifying genome build with check_build...");
            match crate::inference::detect_build_from_dtc(&prescan_records, reference.as_ref().unwrap()) {
                Ok(detected) => {
                    if detected != config.assembly {
                        tracing::warn!(
                            "Detected build {} differs from default {}. Using detected build.",
                            detected,
                            config.assembly
                        );
                        config.assembly = detected.clone();
                    } else {
                        tracing::info!("Build verified: {}", config.assembly);
                    }
                    // Note: We don't have match rates here, would need to modify inference module
                    build_detection = Some(crate::report::BuildDetection {
                        detected_build: detected,
                        hg19_match_rate: 0.0,
                        hg38_match_rate: 0.0,
                    });
                }
                Err(e) => {
                    tracing::warn!("Build detection failed: {}. Using default GRCh38.", e);
                }
            }
        }

        // Infer sex if not provided
        if config.sex.is_none() {
            tracing::info!("Sex not specified, inferring from input data...");
            match crate::inference::infer_sex_from_records(&prescan_records, &config.assembly) {
                Ok(inferred) => {
                    tracing::info!("Inferred sex: {:?}", inferred);
                    config.sex = Some(inferred);
                    sex_inferred = true;
                }
                Err(e) => {
                    tracing::warn!("Sex inference failed: {}. Defaulting to Female.", e);
                    config.sex = Some(Sex::Female);
                    sex_inferred = true;
                }
            }
        }
    }

    // Load panel if provided
    let padded_panel: Option<std::cell::RefCell<crate::panel::PaddedPanel>> =
        if let Some(panel_path) = &config.panel {
            tracing::info!(panel = %panel_path.display(), "loading reference panel");
            let panel_index = crate::panel::PanelIndex::load(panel_path)
                .with_context(|| format!("failed to load panel {}", panel_path.display()))?;
            tracing::info!(sites = panel_index.len(), "panel loaded");
            Some(std::cell::RefCell::new(crate::panel::PaddedPanel::new(
                panel_index,
            )))
        } else {
            None
        };

    let header = build_header(&config, reference.as_ref())?;

    // Instantiate Source Iterator
    let source: Box<dyn crate::input::VariantSource> = match config.input_format {
        crate::input::InputFormat::Dtc => {
            let input = fs::File::open(&config.input)
                .with_context(|| format!("failed to open input {}", config.input.display()))?;
            let reader = BufReader::new(input);
            let dtc_reader = dtc::Reader::new(reader);
            // DTC format always requires reference (checked above)
            let source =
                crate::input::DtcSource::new(dtc_reader, reference.clone().unwrap(), config.clone());
            Box::new(source)
        }
        crate::input::InputFormat::Vcf => {
            let input = fs::File::open(&config.input)
                .with_context(|| format!("failed to open input {}", config.input.display()))?;
            let reader = BufReader::new(input);
            let vcf_reader = vcf::io::Reader::new(reader);
            let source = if let Some(ref ref_genome) = reference {
                crate::input::VcfSource::new(vcf_reader, ref_genome)
                    .with_context(|| "failed to initialize VCF source")?
            } else {
                crate::input::VcfSource::new_without_reference(vcf_reader)
                    .with_context(|| "failed to initialize VCF source")?
            };
            Box::new(source)
        }
        crate::input::InputFormat::Bcf => {
            let input = fs::File::open(&config.input)
                .with_context(|| format!("failed to open input {}", config.input.display()))?;
            let bcf_reader = bcf::io::Reader::new(input);
            let source = if let Some(ref ref_genome) = reference {
                crate::input::BcfSource::new(bcf_reader, ref_genome)
                    .with_context(|| "failed to initialize BCF source")?
            } else {
                crate::input::BcfSource::new_without_reference(bcf_reader)
                    .with_context(|| "failed to initialize BCF source")?
            };
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
                reference.as_ref(),
                &mut writer,
                &header,
                &config,
                &mut summary,
                padded_panel.as_ref(),
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
                reference.as_ref(),
                &mut writer,
                &header,
                &config,
                &mut summary,
                padded_panel.as_ref(),
            )?;
        }
        OutputFormat::Plink => {
            let mut writer =
                PlinkWriter::new(&config.output).context("failed to create PLINK writer")?;

            writer
                .write_fam(&config.sample_id, config.sex.expect("sex should be set"))
                .context("failed to write FAM file")?;

            process_records(
                source,
                reference.as_ref(),
                &mut writer,
                &header,
                &config,
                &mut summary,
                padded_panel.as_ref(),
            )?;
        }
    }

    // Write padded panel if we have one and output_dir is set
    if let (Some(panel), Some(output_dir)) = (&padded_panel, &config.output_dir) {
        let panel = panel.borrow();
        if panel.modified_site_count() > 0 || panel.novel_site_count() > 0 {
            let panel_output = output_dir.join("panel.vcf");
            tracing::info!(
                path = %panel_output.display(),
                modified = panel.modified_site_count(),
                novel = panel.novel_site_count(),
                "writing padded panel"
            );

            if let Some(original_panel_path) = &config.panel {
                crate::panel_writer::write_padded_panel(original_panel_path, &panel, &panel_output)
                    .with_context(|| "failed to write padded panel")?;
            }
        } else {
            tracing::info!("no panel modifications needed, skipping padded panel output");
        }
    }

    // Build and write the run report
    let panel_info = padded_panel.as_ref().map(|p| {
        let panel = p.borrow();
        crate::report::PanelInfo {
            path: config
                .panel
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_default(),
            total_sites: panel.total_site_count(),
            modified_sites: panel.modified_site_count(),
            novel_sites: panel.novel_site_count(),
        }
    });

    let report_builder = crate::report::RunReportBuilder {
        input_path: config.input.display().to_string(),
        input_format: Some(config.input_format),
        input_origin: config.input_origin.clone(),
        output_path: config.output.display().to_string(),
        output_format: Some(config.output_format),
        reference_path: config.reference_fasta.as_ref().map(|p| p.display().to_string()).unwrap_or_default(),
        reference_origin: config.reference_origin.clone().unwrap_or_default(),
        assembly: config.assembly.clone(),
        standardize: config.standardize,
        panel: panel_info,
        sample_id: config.sample_id.clone(),
        sex: config.sex,
        sex_inferred,
        build_detection,
    };
    let report = report_builder.build(&summary);
    if let Err(e) = report.write(&config.output) {
        tracing::warn!("Failed to write run report: {}", e);
    }

    Ok(summary)
}

fn process_records<S, W>(
    mut source: S,
    reference: Option<&ReferenceGenome>,
    writer: &mut W,
    header: &vcf::Header,
    config: &ConversionConfig,
    summary: &mut crate::ConversionSummary,
    panel: Option<&std::cell::RefCell<crate::panel::PaddedPanel>>,
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
                    match standardize_record(&record, reference.expect("reference required for standardize"), config) {
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

                // Apply panel harmonization if panel is provided
                if let Some(panel_cell) = panel {
                    let chrom = final_record.reference_sequence_name().to_string();
                    let pos = final_record.variant_start().map(usize::from).unwrap_or(0) as u64;
                    let ref_base = final_record.reference_bases().to_string();

                    // Collect alleles from the record
                    let record_alts: Vec<String> = final_record
                        .alternate_bases()
                        .as_ref()
                        .iter()
                        .map(|s| s.to_string())
                        .collect();

                    // Build list of all input alleles: [REF, ALT1, ALT2...]
                    let mut all_input_alleles = vec![ref_base.clone()];
                    all_input_alleles.extend(record_alts);

                    // Harmonize alleles against panel (handles strand flips)
                    let mut panel_borrow = panel_cell.borrow_mut();
                    // Note: harmonize_alleles registers alleles with panel for padding
                    // The mapping result could be used to remap GT indices if needed
                    let _ = crate::harmonize::harmonize_alleles(
                        &all_input_alleles,
                        &ref_base,
                        &chrom,
                        pos,
                        &mut panel_borrow,
                    );
                }
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
    let (final_ref, final_alts, needs_remap) = if input_ref.len() == 1 && input_ref != ref_base_str
    {
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
            let mut mapping = std::collections::HashMap::new();

            // Mapping Logic:
            // Old REF (Index 0) -> New Index 1 (it becomes the first ALT)
            mapping.insert(0, 1);

            // Old ALT that matches Reference (Index flip_idx + 1) -> New Index 0 (REF)
            mapping.insert(flip_idx + 1, 0);

            let mut next_new_idx = 2;
            for (i, alt) in alt_bases.iter().enumerate() {
                if i != flip_idx {
                    new_alts.push(alt.clone());
                    // Old ALT index was i + 1. New index is next_new_idx.
                    mapping.insert(i + 1, next_new_idx);
                    next_new_idx += 1;
                }
            }

            (ref_base_str.clone(), new_alts, Some(mapping))
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
        // REF matches or is multi-base (indel)
        // Validation: For Indels, check at least the first base matches reference
        // (VCF spec requires POS to be the position of the first base of REF)
        if input_ref.len() > 1 {
            let first_char = input_ref.chars().next().unwrap();
            let ref_first_char = ref_base_str.chars().next().unwrap(); // ref_base is usually single char from .base() call

            // Wait, reference.base() returns a single char at `pos`.
            // If input_ref="TGT", first base 'T' must match ref_base 'T'.
            if first_char != ref_first_char {
                tracing::warn!(
                    "Indel REF mismatch at {}:{}: user REF={} but reference base={}. Skipping.",
                    canonical_name,
                    pos,
                    input_ref,
                    ref_first_char
                );
                return Ok(None);
            }
        }

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

    // 5. Apply GT remapping if allele polarization occurred
    let samples = if let Some(mapping) = needs_remap {
        tracing::debug!(
            chrom = %canonical_name, pos = pos,
            "Allele polarization applied, remapping GT indices: {:?}", mapping
        );
        remap_sample_genotypes(record.samples(), &mapping)
    } else {
        record.samples().clone()
    };

    // 6. Check ploidy and enforce if needed
    let ploidy = determine_ploidy(
        &canonical_name,
        pos,
        config.sex.unwrap_or(Sex::Female),
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

/// Remap genotype indices in all samples based on mapping
fn remap_sample_genotypes(
    samples: &noodles::vcf::variant::record_buf::Samples,
    mapping: &std::collections::HashMap<usize, usize>,
) -> noodles::vcf::variant::record_buf::Samples {
    let keys = samples.keys();
    let mut new_keys_vec = Vec::new();
    // Filter keys: keep GT and safe fields (GQ, DP, MIN_DP), drop allele-dependent fields (PL, AD, GP)
    for key in keys.as_ref().iter() {
        if key == format_key::GENOTYPE || key == "GQ" || key == "DP" || key == "MIN_DP" {
            new_keys_vec.push(key.clone());
        } else {
            tracing::debug!("Dropping field {} during polarization", key);
        }
    }

    // Create new keys object
    let new_keys: noodles::vcf::variant::record_buf::samples::Keys =
        new_keys_vec.into_iter().collect();

    let mut new_values = Vec::new();

    for sample in samples.values() {
        let mut new_sample_vals = Vec::new();

        // Iterate through valid keys and extract values from original sample
        for key in new_keys.as_ref().iter() {
            let val_opt = sample.get(key).flatten().cloned();

            if key == format_key::GENOTYPE {
                if let Some(Value::String(gt_str)) = &val_opt {
                    new_sample_vals.push(Some(Value::String(remap_gt_string(gt_str, mapping))));
                } else {
                    new_sample_vals.push(val_opt);
                }
            } else {
                new_sample_vals.push(val_opt);
            }
        }
        new_values.push(new_sample_vals);
    }

    Samples::new(new_keys, new_values)
}

/// Remap genotype indices in a GT string (e.g., "0/1" -> "1/0")
/// Handles standard separators (/, |) and ignores '.'
fn remap_gt_string(gt: &str, mapping: &std::collections::HashMap<usize, usize>) -> String {
    // Split by separators, map, join back?
    // Or char-by-char if single digits?
    // GT might have multi-digit indices (e.g. 10). Char mapping is unsafe.
    // Parse fully.

    // Simple parser for VCF GT string:
    // Split on first separator to detect phase?
    // Actually, VCF GT is usually integers separated by / or |.

    // We can use a regex or manual scan.
    // Manual scan is safer and faster.
    let mut result = String::with_capacity(gt.len());
    let mut current_num = String::new();

    for c in gt.chars() {
        if c.is_ascii_digit() {
            current_num.push(c);
        } else {
            // End of number
            if !current_num.is_empty() {
                if let Ok(idx) = current_num.parse::<usize>() {
                    if let Some(new_idx) = mapping.get(&idx) {
                        result.push_str(&new_idx.to_string());
                    } else {
                        // Keep original validity? Or warn?
                        // If mapping incomplete, keep original (risky) or map to .?
                        result.push_str(&current_num);
                    }
                } else {
                    result.push_str(&current_num);
                }
                current_num.clear();
            }
            result.push(c);
        }
    }
    // Final number
    if !current_num.is_empty() {
        if let Ok(idx) = current_num.parse::<usize>() {
            if let Some(new_idx) = mapping.get(&idx) {
                result.push_str(&new_idx.to_string());
            } else {
                result.push_str(&current_num);
            }
        } else {
            result.push_str(&current_num);
        }
    }

    result
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

fn build_header(config: &ConversionConfig, reference: Option<&ReferenceGenome>) -> Result<vcf::Header> {
    let mut builder = vcf::Header::builder().set_file_format(FileFormat::new(4, 5));

    let genotype_format = Map::<Format>::from(format_key::GENOTYPE);
    builder = builder.add_format(format_key::GENOTYPE, genotype_format);

    // Add FORMAT fields that may be preserved during standardization
    use noodles::vcf::header::record::value::map::format::{Number as FmtNumber, Type as FmtType};
    builder = builder
        .add_format(
            "GQ",
            Map::<Format>::new(FmtNumber::Count(1), FmtType::Integer, "Genotype Quality"),
        )
        .add_format(
            "DP",
            Map::<Format>::new(FmtNumber::Count(1), FmtType::Integer, "Read Depth"),
        )
        .add_format(
            "MIN_DP",
            Map::<Format>::new(
                FmtNumber::Count(1),
                FmtType::Integer,
                "Minimum DP observed within the GVCF block",
            ),
        );

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

    // Add contigs from reference if available
    if let Some(ref_genome) = reference {
        for contig in ref_genome.contigs() {
            let mut contig_map = Map::<Contig>::new();
            if let Ok(length) = usize::try_from(contig.length) {
                *contig_map.length_mut() = Some(length);
            }
            builder = builder.add_contig(contig.name.clone(), contig_map);
        }
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

    let reference_uri = if config.reference_origin.as_ref().map(|s| is_remote_source(s)).unwrap_or(false) {
        config.reference_origin.clone().unwrap_or_default()
    } else {
        config.reference_fasta.as_ref()
            .map(|p| format!("file://{}", p.display()))
            .unwrap_or_else(|| "none".to_string())
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
            reference_fasta: Some(ref_path.clone()),
            reference_origin: "dummy_ref".to_string(),
            reference_fai: None,
            reference_fai_origin: None,
            output: dir.path().join("out.vcf"),
            output_dir: None,
            output_format: OutputFormat::Vcf,
            sample_id: "SAMPLE".to_string(),
            assembly: "GRCh38".to_string(),
            include_reference_sites: true,
            sex: Some(Sex::Female),
            par_boundaries: None,
            standardize: false,
            panel: None,
        };

        let header = build_header(&config, Some(&reference)).unwrap();

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
            reference_origin: Some(fasta_path.path().to_string_lossy().to_string()),
            reference_fai: None,
            reference_fai_origin: None,
            output: PathBuf::from("out.vcf"),
            output_dir: None,
            output_format: OutputFormat::Vcf,
            sample_id: String::from("sample"),
            assembly: String::from("GRCh38"),
            include_reference_sites: true,
            sex: Some(Sex::Female),
            par_boundaries: None,
            standardize: false,
            panel: None,
        };

        let header = build_header(&config, Some(&reference)).unwrap();
        assert!(!header.contigs().is_empty());
        assert!(header.other_records().contains_key("source"));
        assert!(header.other_records().contains_key("reference"));
    }
}
