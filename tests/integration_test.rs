use std::{
    fs::{self, File},
    io::{BufReader, Write},
    path::PathBuf,
};

use convert_genome::{
    ConversionConfig, OutputFormat, convert_dtc_file, reference::ReferenceGenome,
};
use noodles::{
    bcf::io::reader::Builder as BcfReaderBuilder,
    vcf::{io::Reader as VcfReader, record::Record},
};
use tempfile::tempdir;

fn write_reference(dir: &PathBuf) -> PathBuf {
    let fasta_path = dir.join("reference.fa");
    let mut file = File::create(&fasta_path).expect("create reference");
    writeln!(file, ">1").unwrap();
    writeln!(file, "ACGTACGTACGT").unwrap();
    writeln!(file, ">2").unwrap();
    writeln!(file, "GGGGCCTTAA").unwrap();
    fasta_path
}

fn write_dtc(dir: &PathBuf, contents: &str) -> PathBuf {
    let dtc_path = dir.join("input.txt");
    fs::write(&dtc_path, contents).expect("write dtc");
    dtc_path
}

fn base_config(input: PathBuf, reference: PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_origin: input.to_string_lossy().into_owned(),
        reference_fasta: reference.clone(),
        reference_origin: reference.to_string_lossy().into_owned(),
        reference_fai: None,
        reference_fai_origin: None,
        output: output.clone(),
        output_format: OutputFormat::Vcf,
        sample_id: String::from("sample"),
        assembly: String::from("GRCh38"),
        include_reference_sites: true,
    }
}

#[test]
fn converts_dtc_to_vcf() {
    let temp = tempdir().unwrap();
    let root = temp.path().to_path_buf();
    let fasta = write_reference(&root);
    let dtc = write_dtc(
        &root,
        "#source\nrs1\t1\t2\tAC\nrs2\t1\t5\tTT\nrs3\t2\t3\tGG\n",
    );
    let output = root.join("output.vcf");
    let config = base_config(dtc, fasta.clone(), output.clone());

    let summary = convert_dtc_file(config).expect("conversion succeeds");
    assert_eq!(summary.total_records, 3);
    assert_eq!(
        summary.variant_records + summary.reference_records,
        summary.emitted_records
    );

    let mut reader = VcfReader::new(BufReader::new(File::open(&output).unwrap()));
    let _header = reader.read_header().expect("header");
    let records: Vec<Record> = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("records");
    assert_eq!(records.len(), 3);
    assert_eq!(records[0].reference_sequence_name(), "1");
    let pos = records[0]
        .variant_start()
        .and_then(|p| p.ok())
        .expect("position");
    assert_eq!(usize::from(pos), 2);
}

#[test]
fn converts_dtc_to_bcf() {
    let temp = tempdir().unwrap();
    let root = temp.path().to_path_buf();
    let fasta = write_reference(&root);
    let dtc = write_dtc(&root, "rs1\t1\t2\tAC\n");
    let output = root.join("output.bcf");
    let mut config = base_config(dtc, fasta, output.clone());
    config.output = output.clone();
    config.output_format = OutputFormat::Bcf;

    let summary = convert_dtc_file(config).expect("conversion succeeds");
    assert_eq!(summary.total_records, 1);
    assert_eq!(summary.emitted_records, 1);

    let mut reader = BcfReaderBuilder::default()
        .build_from_path(&output)
        .expect("open bcf");
    let _header = reader.read_header().expect("header");
    let records: Vec<_> = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("records");
    assert_eq!(records.len(), 1);
}

#[test]
fn handles_empty_input() {
    let temp = tempdir().unwrap();
    let root = temp.path().to_path_buf();
    let fasta = write_reference(&root);
    let dtc = write_dtc(&root, "");
    let output = root.join("output.vcf");
    let config = base_config(dtc, fasta, output);

    let summary = convert_dtc_file(config).expect("conversion succeeds");
    assert_eq!(summary.total_records, 0);
    assert_eq!(summary.emitted_records, 0);
    assert_eq!(summary.parse_errors, 0);
}

#[test]
fn counts_malformed_input_and_unknown_chromosomes() {
    let temp = tempdir().unwrap();
    let root = temp.path().to_path_buf();
    let fasta = write_reference(&root);
    let malformed = write_dtc(&root, "bad line\n");
    let output = root.join("malformed.vcf");
    let config = base_config(malformed.clone(), fasta.clone(), output.clone());

    let summary = convert_dtc_file(config).expect("conversion succeeds");
    assert_eq!(summary.total_records, 0);
    assert_eq!(summary.parse_errors, 1);

    let unknown = write_dtc(&root, "rsX\tchrUn\t1\tAA\n");
    let output_unknown = root.join("unknown.vcf");
    let config_unknown = base_config(unknown, fasta, output_unknown);
    let summary_unknown = convert_dtc_file(config_unknown).expect("conversion succeeds");
    assert_eq!(summary_unknown.unknown_chromosomes, 1);
    assert_eq!(summary_unknown.emitted_records, 0);
}

#[test]
fn parallel_matches_single_thread_output() {
    let temp = tempdir().unwrap();
    let root = temp.path().to_path_buf();
    let fasta = write_reference(&root);
    let dtc = write_dtc(
        &root,
        "rs1\t1\t2\tAC\nrs2\t1\t5\tTT\nrs3\t2\t3\tGG\nrs4\t2\t5\tAA\n",
    );

    let sequential_output = root.join("sequential.vcf");
    let sequential = base_config(dtc.clone(), fasta.clone(), sequential_output.clone());
    unsafe {
        std::env::set_var("RAYON_NUM_THREADS", "1");
    }
    let seq_summary = convert_dtc_file(sequential).expect("sequential run succeeds");
    assert_eq!(seq_summary.total_records, 4);

    let parallel_output = root.join("parallel.vcf");
    let parallel = base_config(dtc, fasta, parallel_output.clone());
    unsafe {
        std::env::set_var("RAYON_NUM_THREADS", "4");
    }
    let par_summary = convert_dtc_file(parallel).expect("parallel run succeeds");
    assert_eq!(par_summary.total_records, 4);

    let seq_bytes = fs::read(sequential_output).unwrap();
    let par_bytes = fs::read(parallel_output).unwrap();
    assert_eq!(seq_bytes, par_bytes);

    unsafe {
        std::env::remove_var("RAYON_NUM_THREADS");
    }
}

#[test]
fn cached_reference_is_thread_safe() {
    let temp = tempdir().unwrap();
    let root = temp.path().to_path_buf();
    let fasta = write_reference(&root);
    let reference = ReferenceGenome::open(&fasta, None).expect("open reference");

    let mut handles = Vec::new();
    for _ in 0..8 {
        let clone = reference.clone();
        handles.push(std::thread::spawn(move || {
            for _ in 0..64 {
                assert_eq!(clone.base("1", 2).unwrap(), 'C');
                assert_eq!(clone.base("2", 4).unwrap(), 'G');
            }
        }));
    }

    for handle in handles {
        handle.join().expect("thread");
    }
}
