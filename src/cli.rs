use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use tracing_subscriber::{EnvFilter, fmt};

use crate::{ConversionConfig, ConversionSummary, OutputFormat, convert_dtc_file};

#[derive(Debug, Parser)]
#[command(author, version, about = "Convert DTC genotype text files to VCF or BCF", long_about = None)]
struct Cli {
    /// Input DTC genotype file (23andMe, LivingDNA, etc.)
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Reference genome FASTA (GRCh38)
    #[arg(value_name = "REFERENCE")]
    reference: PathBuf,

    /// Output VCF or BCF path
    #[arg(value_name = "OUTPUT")]
    output: PathBuf,

    /// Output file format
    #[arg(long, value_enum, default_value_t = OutputFormat::Vcf)]
    format: OutputFormat,

    /// Optional explicit FASTA index (.fai) path
    #[arg(long, value_name = "FAI")]
    reference_fai: Option<PathBuf>,

    /// Sample identifier to embed in the VCF header
    #[arg(long, value_name = "SAMPLE")]
    sample: Option<String>,

    /// Assembly label to embed in metadata
    #[arg(long, default_value = "GRCh38")]
    assembly: String,

    /// When set, omit reference-only sites from the output
    #[arg(long)]
    variants_only: bool,

    /// Logging verbosity (e.g. error, warn, info, debug)
    #[arg(long, default_value = "info")]
    log_level: String,
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();
    init_logging(&cli.log_level)?;

    let sample_id = cli
        .sample
        .clone()
        .or_else(|| derive_sample_name(&cli.input))
        .unwrap_or_else(|| String::from("sample"));

    let config = ConversionConfig {
        input: cli.input.clone(),
        reference_fasta: cli.reference.clone(),
        reference_fai: cli.reference_fai.clone(),
        output: cli.output.clone(),
        output_format: cli.format,
        sample_id,
        assembly: cli.assembly.clone(),
        include_reference_sites: !cli.variants_only,
    };

    let summary = convert_dtc_file(config)?;
    print_summary(&summary);

    Ok(())
}

fn init_logging(level: &str) -> Result<()> {
    let filter = EnvFilter::try_new(level).unwrap_or_else(|_| EnvFilter::new("info"));
    fmt()
        .with_env_filter(filter)
        .with_target(false)
        .try_init()
        .ok();
    Ok(())
}

fn derive_sample_name(path: &PathBuf) -> Option<String> {
    path.file_stem()
        .map(|s| s.to_string_lossy().replace('.', "_"))
        .filter(|s| !s.is_empty())
}

fn print_summary(summary: &ConversionSummary) {
    println!(
        "Processed {total} records; emitted {emitted} ({variants} variants, {references} reference).",
        total = summary.total_records,
        emitted = summary.emitted_records,
        variants = summary.variant_records,
        references = summary.reference_records,
    );

    if summary.skipped_reference_sites > 0 {
        println!(
            "Skipped {skipped} reference-only sites due to --variants-only.",
            skipped = summary.skipped_reference_sites
        );
    }

    if summary.missing_genotype_records > 0 {
        println!(
            "Encountered {count} sites with missing genotypes.",
            count = summary.missing_genotype_records
        );
    }

    if summary.unknown_chromosomes > 0
        || summary.reference_failures > 0
        || summary.invalid_genotypes > 0
    {
        println!(
            "Warnings: {chrom} unknown chromosomes, {ref_err} reference lookup failures, {invalid} invalid genotypes.",
            chrom = summary.unknown_chromosomes,
            ref_err = summary.reference_failures,
            invalid = summary.invalid_genotypes
        );
    }

    if summary.parse_errors > 0 {
        println!(
            "Ignored {count} malformed input lines.",
            count = summary.parse_errors
        );
    }
}
