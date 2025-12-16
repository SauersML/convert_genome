use std::{fs, path::PathBuf};

use anyhow::Result;
use assert_fs::{TempDir, prelude::*};
use convert_genome::reference::ReferenceGenome;
use convert_genome::{ConversionConfig, ConversionSummary, OutputFormat, convert_dtc_file};
use noodles::bcf::io::reader::Builder as BcfReaderBuilder;
use rayon::ThreadPoolBuilder;

fn write_reference(dir: &TempDir) -> Result<PathBuf> {
    let fasta = dir.child("ref.fa");
    fasta.write_str(">chr1\nACGTACGT\n>chr2\nTTTTCCCC\n")?;
    Ok(fasta.path().to_path_buf())
}

fn write_dtc(dir: &TempDir, contents: &str) -> Result<PathBuf> {
    let input = dir.child("input.txt");
    input.write_str(contents)?;
    Ok(input.path().to_path_buf())
}

fn base_config(input: PathBuf, reference: PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_format: convert_genome::input::InputFormat::Dtc,
        input_origin: input.to_string_lossy().to_string(),
        reference_fasta: reference.clone(),
        reference_origin: reference.to_string_lossy().to_string(),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_format: convert_genome::OutputFormat::Vcf,
        sample_id: "SAMPLE".into(),
        assembly: "GRCh38".into(),
        include_reference_sites: true,
        sex: convert_genome::cli::Sex::Female,
        par_boundaries: None,
    }
}

fn run_conversion_with_threads(
    config: ConversionConfig,
    threads: usize,
) -> Result<(ConversionSummary, Vec<u8>)> {
    let output = config.output.clone();
    let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let summary = pool.install(|| convert_dtc_file(config))?;
    let data = fs::read(output)?;
    Ok((summary, data))
}

#[test]
fn converts_to_vcf_and_bcf() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    let input = write_dtc(&temp, "rs1\t1\t2\tAA\nrs2\t1\t3\tAG\nrs3\t2\t4\tTT\n")?;

    let vcf_path = temp.child("out.vcf");
    let config = base_config(
        input.clone(),
        reference.clone(),
        vcf_path.path().to_path_buf(),
    );
    let summary = convert_dtc_file(config.clone())?;

    assert_eq!(summary.emitted_records, 3);
    assert_eq!(summary.variant_records, 2);
    assert_eq!(summary.reference_records, 1);

    let vcf_data = fs::read_to_string(vcf_path.path())?;
    assert!(vcf_data.contains("#CHROM\tPOS"));
    assert!(vcf_data.contains("1\t2"));

    let bcf_path = temp.child("out.bcf");
    let mut bcf_config = config;
    bcf_config.output = bcf_path.path().to_path_buf();
    bcf_config.output_format = OutputFormat::Bcf;
    let bcf_summary = convert_dtc_file(bcf_config)?;
    assert_eq!(bcf_summary.emitted_records, 3);

    let mut reader = BcfReaderBuilder::default().build_from_path(bcf_path.path())?;
    reader.read_header()?;
    let mut records = reader.records();
    assert!(records.next().transpose()?.is_some());

    Ok(())
}

#[test]
fn handles_empty_input() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    let input = write_dtc(&temp, "")?;

    let vcf_path = temp.child("empty.vcf");
    let config = base_config(input, reference, vcf_path.path().to_path_buf());
    let summary = convert_dtc_file(config)?;
    assert_eq!(summary.total_records, 0);
    let vcf_data = fs::read_to_string(vcf_path.path())?;
    assert!(vcf_data.contains("#CHROM"));
    Ok(())
}

#[test]
fn reports_unknown_chromosomes() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    let input = write_dtc(&temp, "rs1\tUn\t1\tAA\n")?;

    let vcf_path = temp.child("unknown.vcf");
    let config = base_config(input, reference, vcf_path.path().to_path_buf());
    let summary = convert_dtc_file(config)?;
    assert_eq!(summary.unknown_chromosomes, 1);
    let data = fs::read_to_string(vcf_path.path())?;
    let variants: Vec<_> = data.lines().filter(|line| !line.starts_with('#')).collect();
    assert!(variants.is_empty());
    Ok(())
}

#[test]
fn parallel_matches_single_thread() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    let input = write_dtc(
        &temp,
        "rs1\t1\t2\tAA\nrs2\t1\t3\tAG\nrs3\t2\t4\tTT\nrs4\t2\t5\tCT\n",
    )?;

    let single_output = temp.child("single.vcf");
    let parallel_output = temp.child("parallel.vcf");

    let single_config = base_config(
        input.clone(),
        reference.clone(),
        single_output.path().to_path_buf(),
    );
    let parallel_config = base_config(input, reference, parallel_output.path().to_path_buf());

    let (single_summary, single_bytes) = run_conversion_with_threads(single_config, 1)?;
    let (parallel_summary, parallel_bytes) = run_conversion_with_threads(parallel_config, 4)?;

    assert_eq!(
        single_summary.emitted_records,
        parallel_summary.emitted_records
    );
    assert_eq!(single_bytes, parallel_bytes);
    Ok(())
}

#[test]
fn reference_cache_populates() -> Result<()> {
    let temp = TempDir::new()?;
    let fasta = write_reference(&temp)?;
    let reference = ReferenceGenome::open(&fasta, None)?;
    assert_eq!(reference.cache_len(), 0);
    assert_eq!(reference.base("1", 1)?, 'A');
    assert!(reference.cache_len() >= 1);
    Ok(())
}
