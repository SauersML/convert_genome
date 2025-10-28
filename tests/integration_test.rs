use std::{
    fs,
    io::{self, Read},
    path::PathBuf,
};

use convert_genome::{
    conversion::{ConversionConfig, OutputFormat, convert_dtc_file},
    reference::ReferenceGenome,
};
use noodles::{
    bcf::io::reader::Builder as BcfReaderBuilder,
    vcf::{self, io::reader::Builder as VcfReaderBuilder},
};
use rayon::ThreadPoolBuilder;
use tempfile::tempdir;

fn write_reference(dir: &tempfile::TempDir) -> io::Result<PathBuf> {
    let reference_path = dir.path().join("reference.fa");
    fs::write(&reference_path, ">chr1\nACGTACGT\n>chr2\nTTTTGGGG\n")?;
    Ok(reference_path)
}

fn write_dtc(dir: &tempfile::TempDir, name: &str, contents: &str) -> io::Result<PathBuf> {
    let path = dir.path().join(name);
    fs::write(&path, contents)?;
    Ok(path)
}

fn base_config(
    input: PathBuf,
    reference: &PathBuf,
    output: PathBuf,
    format: OutputFormat,
    include_reference_sites: bool,
) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_origin: input.display().to_string(),
        reference_fasta: reference.clone(),
        reference_origin: reference.display().to_string(),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_format: format,
        sample_id: String::from("sample"),
        assembly: String::from("GRCh38"),
        include_reference_sites,
    }
}

fn read_vcf_records(path: &PathBuf) -> io::Result<Vec<vcf::Record>> {
    let mut reader = VcfReaderBuilder::default().build_from_path(path)?;
    let _header = reader.read_header()?;
    let mut records = Vec::new();
    for result in reader.records() {
        records.push(result.expect("read VCF record"));
    }
    Ok(records)
}

#[test]
fn full_pipeline_produces_vcf_output() {
    let dir = tempdir().unwrap();
    let reference_path = write_reference(&dir).unwrap();
    let input_path = write_dtc(
        &dir,
        "input.txt",
        "#comment\nrs1\t1\t1\tAA\nrs2\t1\t2\tAG\nrs3\t2\t1\tTT\n",
    )
    .unwrap();
    let output_path = dir.path().join("out.vcf");

    let config = base_config(
        input_path,
        &reference_path,
        output_path.clone(),
        OutputFormat::Vcf,
        true,
    );

    let summary = convert_dtc_file(config).expect("convert to VCF");
    assert_eq!(summary.total_records, 3);
    assert_eq!(summary.emitted_records, 3);
    assert_eq!(summary.variant_records, 1);
    assert_eq!(summary.reference_records, 2);

    let records = read_vcf_records(&output_path).unwrap();
    assert_eq!(records.len(), 3);
    assert!(
        records
            .iter()
            .any(|record| record.reference_sequence_name() == "chr1")
    );
}

#[test]
fn full_pipeline_produces_bcf_output() {
    let dir = tempdir().unwrap();
    let reference_path = write_reference(&dir).unwrap();
    let input_path = write_dtc(&dir, "input.txt", "rs1\t1\t1\tAA\nrs2\t1\t2\tAG\n").unwrap();
    let output_path = dir.path().join("out.bcf");

    let config = base_config(
        input_path,
        &reference_path,
        output_path.clone(),
        OutputFormat::Bcf,
        true,
    );

    let summary = convert_dtc_file(config).expect("convert to BCF");
    assert_eq!(summary.total_records, 2);
    assert_eq!(summary.emitted_records, 2);

    let mut reader = BcfReaderBuilder::default()
        .build_from_path(&output_path)
        .expect("open bcf");
    let _header = reader.read_header().expect("read header");
    let mut records = Vec::new();
    for result in reader.records() {
        records.push(result.expect("read record"));
    }
    assert_eq!(records.len(), 2);
}

#[test]
fn empty_input_yields_header_only() {
    let dir = tempdir().unwrap();
    let reference_path = write_reference(&dir).unwrap();
    let input_path = write_dtc(&dir, "input.txt", "#only\n").unwrap();
    let output_path = dir.path().join("empty.vcf");

    let config = base_config(
        input_path,
        &reference_path,
        output_path.clone(),
        OutputFormat::Vcf,
        false,
    );

    let summary = convert_dtc_file(config).expect("convert empty");
    assert_eq!(summary.total_records, 0);
    assert_eq!(summary.emitted_records, 0);

    let mut contents = String::new();
    fs::File::open(&output_path)
        .unwrap()
        .read_to_string(&mut contents)
        .unwrap();
    assert!(contents.contains("#CHROM"));
}

#[test]
fn malformed_input_increments_parse_errors() {
    let dir = tempdir().unwrap();
    let reference_path = write_reference(&dir).unwrap();
    let input_path = write_dtc(
        &dir,
        "input.txt",
        "rs1\t1\t1\tAA\ninvalid line\nrs2\t1\t2\tAG\n",
    )
    .unwrap();
    let output_path = dir.path().join("malformed.vcf");

    let config = base_config(
        input_path,
        &reference_path,
        output_path,
        OutputFormat::Vcf,
        true,
    );

    let summary = convert_dtc_file(config).expect("convert malformed");
    assert_eq!(summary.total_records, 2);
    assert_eq!(summary.parse_errors, 1);
    assert!(summary.emitted_records >= 1);
}

#[test]
fn unknown_chromosomes_are_reported() {
    let dir = tempdir().unwrap();
    let reference_path = write_reference(&dir).unwrap();
    let input_path = write_dtc(&dir, "input.txt", "rs1\t1\t1\tAA\nrsx\tUnk\t5\tGG\n").unwrap();
    let output_path = dir.path().join("unknown.vcf");

    let config = base_config(
        input_path,
        &reference_path,
        output_path,
        OutputFormat::Vcf,
        true,
    );

    let summary = convert_dtc_file(config).expect("convert unknown");
    assert_eq!(summary.total_records, 2);
    assert_eq!(summary.unknown_chromosomes, 1);
}

#[test]
fn parallel_and_single_thread_outputs_match() {
    let dir = tempdir().unwrap();
    let reference_path = write_reference(&dir).unwrap();
    let mut dtc_content = String::new();
    for i in 1..=50 {
        dtc_content.push_str(&format!("rs{0}\t1\t{0}\tAG\n", i));
    }
    let input_path = write_dtc(&dir, "input.txt", &dtc_content).unwrap();

    let sequential_output = dir.path().join("seq.vcf");
    let parallel_output = dir.path().join("par.vcf");

    let seq_config = base_config(
        input_path.clone(),
        &reference_path,
        sequential_output.clone(),
        OutputFormat::Vcf,
        true,
    );
    let par_config = base_config(
        input_path,
        &reference_path,
        parallel_output.clone(),
        OutputFormat::Vcf,
        true,
    );

    let seq_pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let seq_summary =
        seq_pool.install(|| convert_dtc_file(seq_config).expect("sequential conversion"));
    assert_eq!(seq_summary.total_records, 50);

    let par_pool = ThreadPoolBuilder::new().num_threads(4).build().unwrap();
    let par_summary =
        par_pool.install(|| convert_dtc_file(par_config).expect("parallel conversion"));
    assert_eq!(par_summary.total_records, 50);
    assert_eq!(seq_summary.emitted_records, par_summary.emitted_records);

    let seq_contents = fs::read_to_string(&sequential_output).unwrap();
    let par_contents = fs::read_to_string(&parallel_output).unwrap();
    assert_eq!(seq_contents, par_contents);
}

#[test]
fn reference_cache_is_shared_across_clones() {
    let dir = tempdir().unwrap();
    let reference_path = write_reference(&dir).unwrap();
    let reference = ReferenceGenome::open(&reference_path, None).unwrap();
    assert_eq!(reference.cache_usage(), 0);

    assert_eq!(reference.base("1", 2).unwrap(), 'C');
    assert!(reference.cache_usage() > 0);

    let clone = reference.clone();
    assert_eq!(clone.base("1", 2).unwrap(), 'C');
    assert_eq!(reference.cache_usage(), clone.cache_usage());
}
