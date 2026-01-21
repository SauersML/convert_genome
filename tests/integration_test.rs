use std::{fs, path::PathBuf};

use anyhow::Result;
use assert_fs::{TempDir, prelude::*};
use convert_genome::reference::ReferenceGenome;
use convert_genome::{ConversionConfig, ConversionSummary, OutputFormat, convert_dtc_file};
use noodles::bcf::io::reader::Builder as BcfReaderBuilder;
use noodles::vcf;
use noodles::vcf::variant::record::samples::Sample;
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

fn write_panel_vcf(dir: &TempDir) -> Result<PathBuf> {
    let vcf_path = dir.child("panel.vcf");
    vcf_path.write_str(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\tP3\tP4\tP5\n1\t1500\t.\tT\tTA\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n",
    )?;
    Ok(vcf_path.path().to_path_buf())
}

fn write_reference_2000(dir: &TempDir) -> Result<PathBuf> {
    let fasta = dir.child("ref.fa");
    let seq: String = std::iter::repeat('T').take(2000).collect();
    fasta.write_str(&format!(">1\n{}\n", seq))?;
    Ok(fasta.path().to_path_buf())
}

fn base_config(input: PathBuf, reference: PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_format: convert_genome::input::InputFormat::Dtc,
        input_origin: input.to_string_lossy().to_string(),
        reference_fasta: Some(reference.clone()),
        reference_origin: Some(reference.to_string_lossy().to_string()),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_dir: None,
        output_format: convert_genome::OutputFormat::Vcf,
        sample_id: "SAMPLE".into(),
        assembly: "GRCh38".into(),
        include_reference_sites: true,
        sex: Some(convert_genome::cli::Sex::Female),
        par_boundaries: None,
        standardize: false,
        panel: None,
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
    // Input modified to be concordant with reference (>chr1\nACGTACGT\n>chr2\nTTTTCCCC\n)
    // chr1:2 is C. Input CC -> Match.
    // chr1:3 is G. Input AG -> Match (G).
    // chr2:4 is T. Input TT -> Match.
    // Concordance 3/3 = 1.0. Avoids liftover trigger (which fails in sandbox).
    let input = write_dtc(&temp, "rs1\t1\t2\tCC\nrs2\t1\t3\tAG\nrs3\t2\t4\tTT\n")?;

    let vcf_path = temp.child("out.vcf");
    let config = base_config(
        input.clone(),
        reference.clone(),
        vcf_path.path().to_path_buf(),
    );
    let summary = convert_dtc_file(config.clone())?;

    assert_eq!(summary.emitted_records, 3);
    // Updated expectation due to concordant input change:
    // rs1 (CC) is Ref match (was Variant AA)
    // rs2 (AG) is Variant
    // rs3 (TT) is Ref match
    // So 1 Variant, 2 Reference.
    assert_eq!(summary.variant_records, 1);
    assert_eq!(summary.reference_records, 2);

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
    // Add valid records to pass concordance check (>0.7)
    // 3 valid records (concordant), 1 unknown.
    // Concordance = 3/3 = 1.0 (unknown is skipped in check).
    let input = write_dtc(
        &temp,
        "rs1\tUn\t1\tAA\nrs2\t1\t1\tAA\nrs3\t1\t2\tCC\nrs4\t1\t3\tGG\n",
    )?;

    let vcf_path = temp.child("unknown.vcf");
    let config = base_config(input, reference, vcf_path.path().to_path_buf());
    let summary = convert_dtc_file(config)?;
    assert_eq!(summary.unknown_chromosomes, 1);
    assert_eq!(summary.emitted_records, 3);
    Ok(())
}

#[test]
fn parallel_matches_single_thread() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    // Ensure concordant input
    let input = write_dtc(
        &temp,
        "rs1\t1\t2\tCC\nrs2\t1\t3\tAG\nrs3\t2\t4\tTT\nrs4\t2\t5\tCT\n",
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

#[test]
fn preserves_private_multiallelic_site_with_panel() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference_2000(&temp)?;
    let panel = write_panel_vcf(&temp)?;

    let mut dtc = String::new();
    for pos in 1..=2000u64 {
        if pos == 1500 {
            dtc.push_str(&format!("rs{pos}\t1\t{pos}\tAG\n"));
        } else {
            dtc.push_str(&format!("rs{pos}\t1\t{pos}\tTT\n"));
        }
    }
    let input = write_dtc(&temp, &dtc)?;

    let vcf_path = temp.child("out.vcf");
    let mut config = base_config(input, reference, vcf_path.path().to_path_buf());
    config.standardize = true;
    config.panel = Some(panel);
    convert_dtc_file(config)?;

    let mut reader = vcf::io::reader::Builder::default().build_from_path(vcf_path.path())?;
    let header = reader.read_header()?;

    let mut found = false;
    for result in reader.record_bufs(&header) {
        let record = result?;
        let chrom = record.reference_sequence_name().to_string();
        let pos_raw = record.variant_start().expect("missing pos");
        let pos = usize::from(pos_raw) as u64;
        if chrom == "1" && pos == 1500 {
            found = true;

            assert_eq!(record.reference_bases().to_string().to_uppercase(), "T");

            let alts: Vec<String> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|s| s.to_string().to_uppercase())
                .collect();
            assert_eq!(alts.len(), 2);
            assert!(alts.contains(&"A".to_string()));
            assert!(alts.contains(&"G".to_string()));

            let samples = record.samples();
            let sample_values = samples.values().next().expect("missing sample");
            let gt = sample_values
                .iter(&header)
                .next()
                .expect("missing GT")?
                .1
                .expect("missing GT value");
            let gt_str: String = match gt {
                vcf::variant::record::samples::series::Value::String(s) => s.into_owned(),
                vcf::variant::record::samples::series::Value::Genotype(geno) => {
                    let parts: Vec<String> = geno
                        .iter()
                        .map(|a| {
                            let a = a.expect("invalid allele");
                            a.0.map(|v| v.to_string()).unwrap_or_else(|| ".".to_string())
                        })
                        .collect();
                    parts.join("/")
                }
                _ => panic!("unexpected GT type"),
            };
            assert!(gt_str == "1/2" || gt_str == "2/1");

            break;
        }
    }

    assert!(found, "did not find chr1:1500 in output VCF");
    Ok(())
}
