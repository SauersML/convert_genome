use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

use convert_genome::{
    ConversionOptions, OutputFormat, ReferenceGenome, convert, load_dtc_records, write_records,
};

#[derive(Parser, Debug)]
#[command(author, version, about = "Convert direct-to-consumer genotypes to VCF/BCF", long_about = None)]
struct Cli {
    /// Path to the direct-to-consumer genotype text file.
    #[arg(long)]
    dtc: PathBuf,

    /// Reference FASTA (hg38) used for determining REF alleles.
    #[arg(long)]
    reference: PathBuf,

    /// Output file path.
    #[arg(long)]
    output: PathBuf,

    /// Output format: vcf or bcf.
    #[arg(long, default_value = "vcf")]
    format: OutputFormat,

    /// Sample name to emit in the VCF/BCF header.
    #[arg(long, default_value = "sample")]
    sample_name: String,

    /// Skip homozygous reference calls to reduce output size.
    #[arg(long)]
    skip_reference_calls: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let reference = ReferenceGenome::open(&cli.reference)?;
    let mut options = ConversionOptions::default();
    options.sample_name = cli.sample_name.clone();
    options.skip_reference_calls = cli.skip_reference_calls;

    let records = load_dtc_records(&cli.dtc)?;
    let result = convert(records, reference, &options)?;
    write_records(&cli.output, cli.format, &result.header, &result.records)?;

    Ok(())
}
