use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use tracing_subscriber::{EnvFilter, fmt};
use url::Url;

use crate::{
    ConversionConfig, ConversionSummary, OutputFormat, convert_dtc_file,
    input::InputFormat,
    remote::{self, RemoteResource},
};

#[derive(Debug, Clone, Copy, Eq, PartialEq, clap::ValueEnum)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

#[derive(Debug, Parser)]
#[command(author, version, about = "Convert DTC genotype text files to VCF or BCF", long_about = None)]
struct Cli {
    /// Input DTC genotype file (23andMe, LivingDNA, etc.)
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    /// Input file format (auto-detected if not specified)
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    input_format: InputFormat,

    /// Reference genome FASTA (GRCh38). If omitted, a known reference is downloaded as needed.
    #[arg(long, value_name = "REFERENCE")]
    reference: Option<PathBuf>,

    /// Output VCF or BCF path (mutually exclusive with --output-dir)
    #[arg(value_name = "OUTPUT", conflicts_with = "output_dir")]
    output: Option<PathBuf>,

    /// Output directory for panel mode (produces panel.vcf + genotypes.vcf)
    #[arg(long, value_name = "DIR", conflicts_with = "output")]
    output_dir: Option<PathBuf>,

    /// Output file format
    #[arg(long, value_enum, default_value_t = OutputFormat::Vcf)]
    format: OutputFormat,

    /// Optional explicit FASTA index (.fai) path
    #[arg(long, value_name = "FAI")]
    reference_fai: Option<PathBuf>,

    /// Reference panel VCF/BCF for allele harmonization (Beagle compatibility)
    #[arg(long, value_name = "FILE")]
    panel: Option<PathBuf>,

    /// Sample identifier to embed in the VCF header
    #[arg(long, value_name = "SAMPLE")]
    sample: Option<String>,

    /// Assembly label to embed in metadata
    #[arg(long, default_value = "GRCh38")]
    assembly: String,

    /// When set, omit reference-only sites from the output
    #[arg(long)]
    variants_only: bool,

    /// Standardize/normalize the input file without format conversion.
    /// Performs: chromosome naming normalization, allele polarization against
    /// reference, sex chromosome ploidy enforcement, and sorting.
    #[arg(long)]
    standardize: bool,

    /// Logging verbosity (e.g. error, warn, info, debug)
    #[arg(long, default_value = "info")]
    log_level: String,

    /// Sex of the sample (auto-detected if not specified)
    #[arg(long, value_enum)]
    sex: Option<Sex>,
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();
    init_logging(&cli.log_level)?;

    // Validate output arguments
    let output = match (&cli.output, &cli.output_dir, &cli.panel) {
        (Some(output), None, _) => output.clone(),
        (None, Some(dir), Some(_)) => {
            // Panel mode - create output directory and use genotypes.vcf as output
            std::fs::create_dir_all(dir)
                .with_context(|| format!("failed to create output directory: {}", dir.display()))?;
            dir.join("genotypes.vcf")
        }
        (None, Some(_), None) => {
            anyhow::bail!("--output-dir requires --panel");
        }
        (None, None, _) => {
            anyhow::bail!("Either --output or --output-dir (with --panel) is required");
        }
        _ => unreachable!(), // conflicts_with handles other cases
    };

    let sample_id = cli
        .sample
        .clone()
        .or_else(|| derive_sample_name(&cli.input))
        .unwrap_or_else(|| String::from("sample"));

    let mut resources = ResourceManager::new();
    let input_origin = cli.input.to_string_lossy().to_string();
    let reference_origin = cli
        .reference
        .as_ref()
        .map(|path| path.to_string_lossy().to_string());
    let reference_fai_origin = cli
        .reference_fai
        .as_ref()
        .map(|path| path.to_string_lossy().to_string());

    let resolved_input = resources.resolve(&cli.input)?;
    if cli.reference.is_none() && cli.reference_fai.is_some() {
        anyhow::bail!("--reference-fai requires --reference");
    }

    let resolved_reference = match &cli.reference {
        Some(path) => Some(resources.resolve(path)?),
        None => None,
    };
    let resolved_reference_fai = match &cli.reference_fai {
        Some(path) => Some(resources.resolve(path)?),
        None => None,
    };

    // Resolve panel if provided
    let resolved_panel = match &cli.panel {
        Some(path) => Some(resources.resolve(path)?),
        None => None,
    };

    let input_format = if matches!(cli.input_format, InputFormat::Auto) {
        InputFormat::detect(&resolved_input)
    } else {
        cli.input_format
    };

    let config = ConversionConfig {
        input: resolved_input,
        input_format,
        input_origin,
        reference_fasta: resolved_reference,
        reference_origin,
        reference_fai: resolved_reference_fai,
        reference_fai_origin,
        output: output.clone(),
        output_dir: cli.output_dir.clone(),
        output_format: cli.format,
        sample_id,
        assembly: cli.assembly.clone(),
        include_reference_sites: !cli.variants_only,
        sex: cli.sex,
        par_boundaries: crate::reference::ParBoundaries::new(&cli.assembly),
        standardize: cli.standardize,
        panel: resolved_panel,
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn parses_output_without_reference_positional() {
        let cli = Cli::parse_from(["convert_genome", "input.txt", "output.vcf"]);
        assert_eq!(cli.input, PathBuf::from("input.txt"));
        assert_eq!(cli.reference, None);
        assert_eq!(cli.output, Some(PathBuf::from("output.vcf")));
    }

    #[test]
    fn parses_explicit_reference_flag() {
        let cli = Cli::parse_from([
            "convert_genome",
            "input.txt",
            "--reference",
            "ref.fa",
            "output.vcf",
        ]);
        assert_eq!(cli.reference, Some(PathBuf::from("ref.fa")));
        assert_eq!(cli.output, Some(PathBuf::from("output.vcf")));
    }
}
