use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use tracing_subscriber::{EnvFilter, fmt};
use url::Url;

use crate::{
    ConversionConfig, ConversionSummary, OutputFormat, convert_dtc_file,
    remote::{self, RemoteResource},
};

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

    let mut resources = ResourceManager::new();
    let input_origin = cli.input.to_string_lossy().to_string();
    let reference_origin = cli.reference.to_string_lossy().to_string();
    let reference_fai_origin = cli
        .reference_fai
        .as_ref()
        .map(|path| path.to_string_lossy().to_string());

    let resolved_input = resources.resolve(&cli.input)?;
    let resolved_reference = resources.resolve(&cli.reference)?;
    let resolved_reference_fai = match &cli.reference_fai {
        Some(path) => Some(resources.resolve(path)?),
        None => None,
    };

    let config = ConversionConfig {
        input: resolved_input,
        input_origin,
        reference_fasta: resolved_reference,
        reference_origin,
        reference_fai: resolved_reference_fai,
        reference_fai_origin,
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

fn derive_sample_name(path: &Path) -> Option<String> {
    if let Some(raw) = path.to_str()
        && raw.contains("://")
        && let Ok(url) = Url::parse(raw)
    {
        return url
            .path_segments()
            .and_then(|mut segments| segments.next_back())
            .filter(|segment| !segment.is_empty())
            .map(|segment| segment.replace('.', "_"));
    }

    path.file_stem()
        .map(|s| s.to_string_lossy().replace('.', "_"))
        .filter(|s| !s.is_empty())
}

struct ResourceManager {
    remotes: Vec<RemoteResource>,
}

impl ResourceManager {
    fn new() -> Self {
        Self {
            remotes: Vec::new(),
        }
    }

    fn resolve<P>(&mut self, path: P) -> Result<PathBuf>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref();
        let raw = path.to_string_lossy();
        if let Some(url) = parse_url(&raw) {
            let resource = remote::fetch_remote_resource(&url)
                .with_context(|| format!("failed to fetch {url}"))?;
            let local_path = resource.local_path().to_path_buf();
            self.remotes.push(resource);
            Ok(local_path)
        } else {
            Ok(path.to_path_buf())
        }
    }
}

fn parse_url(raw: &str) -> Option<Url> {
    if raw.contains("://") {
        Url::parse(raw).ok()
    } else {
        None
    }
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

    if summary.symbolic_allele_records > 0 {
        println!(
            "Skipped {count} indel-like genotypes (represented as symbolic alleles).",
            count = summary.symbolic_allele_records
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
