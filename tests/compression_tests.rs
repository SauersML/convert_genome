#[test]
fn test_dtc_with_comments() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    // DTC file starting with comments
    let dtc_content = "# This is a comment\n# Another comment\nrs1\t1\t1\tAC\n";

    let input = temp.child("test.txt");
    input.write_str(dtc_content).unwrap();

    let output = temp.child("out.vcf");
    let config = create_config(
        input.path().to_path_buf(),
        ref_path,
        output.path().to_path_buf(),
    );

    // Should detect as DTC (because it fails VCF check) and convert successfully
    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

use assert_fs::prelude::*;
use convert_genome::{ConversionConfig, OutputFormat, convert_dtc_file};
use std::io::{Cursor, Write};
use std::path::PathBuf;
use zip::write::FileOptions;

fn create_nested_file(
    dir: &assert_fs::TempDir,
    content: &str,
    layers: &[&str],
    filename: &str,
) -> PathBuf {
    let mut current_data = content.as_bytes().to_vec();
    let mut current_name = filename.to_string();

    for layer in layers.iter().rev() {
        match *layer {
            "gz" => {
                let mut encoder =
                    flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
                encoder.write_all(&current_data).unwrap();
                current_data = encoder.finish().unwrap();
                current_name = format!("{}.gz", current_name);
            }
            "zip" => {
                let mut buf = Vec::new();
                {
                    let mut zip = zip::ZipWriter::new(Cursor::new(&mut buf));
                    zip.start_file::<&str, ()>(&current_name, FileOptions::default())
                        .unwrap();
                    zip.write_all(&current_data).unwrap();
                    zip.finish().unwrap();
                }
                current_data = buf;
                current_name = format!("{}.zip", current_name);
            }
            _ => panic!("Unknown layer type"),
        }
    }

    let path = dir.child(&current_name);
    path.write_binary(&current_data).unwrap();
    path.path().to_path_buf()
}

fn create_reference(dir: &assert_fs::TempDir) -> PathBuf {
    let fasta = dir.child("ref.fa");
    fasta.write_str(">1\nACGTACGTACGT\n").unwrap();
    fasta.path().to_path_buf()
}

fn create_config(input: PathBuf, reference: PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input,
        input_format: convert_genome::input::InputFormat::Auto, // Test auto-detection
        input_origin: "test".into(),
        reference_fasta: Some(reference),
        reference_origin: Some("ref".into()),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_dir: None,
        output_format: OutputFormat::Vcf,
        sample_id: "SAMPLE".into(),
        assembly: "GRCh38".into(),
        include_reference_sites: true,
        sex: Some(convert_genome::cli::Sex::Female),
        par_boundaries: None,
        standardize: false,
        panel: None,
    }
}

#[test]
fn test_plain_vcf() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
    let input = temp.child("test.vcf");
    input.write_str(vcf_content).unwrap();

    let output = temp.child("out.vcf");
    let config = create_config(
        input.path().to_path_buf(),
        ref_path,
        output.path().to_path_buf(),
    );

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

#[test]
fn test_vcf_gz() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
    let input = create_nested_file(&temp, vcf_content, &["gz"], "test.vcf");

    let output = temp.child("out.vcf");
    let config = create_config(input, ref_path, output.path().to_path_buf());

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

#[test]
fn test_vcf_zip() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
    let input = create_nested_file(&temp, vcf_content, &["zip"], "test.vcf");

    let output = temp.child("out.vcf");
    let config = create_config(input, ref_path, output.path().to_path_buf());

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

#[test]
fn test_vcf_gz_zip() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
    // test.vcf -> test.vcf.gz -> test.vcf.gz.zip
    let input = create_nested_file(&temp, vcf_content, &["gz", "zip"], "test.vcf");

    let output = temp.child("out.vcf");
    let config = create_config(input, ref_path, output.path().to_path_buf());

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

#[test]
fn test_lying_extension_gz_as_txt() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";

    // Create GZ content
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(vcf_content.as_bytes()).unwrap();
    let gz_data = encoder.finish().unwrap();

    // Write as .txt
    let input = temp.child("test.txt");
    input.write_binary(&gz_data).unwrap();

    let output = temp.child("out.vcf");
    let config = create_config(
        input.path().to_path_buf(),
        ref_path,
        output.path().to_path_buf(),
    );

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

#[test]
fn test_lying_extension_zip_as_gz() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";

    // Create ZIP content
    let mut buf = Vec::new();
    {
        let mut zip = zip::ZipWriter::new(Cursor::new(&mut buf));
        zip.start_file::<&str, ()>("test.vcf", FileOptions::default())
            .unwrap();
        zip.write_all(vcf_content.as_bytes()).unwrap();
        zip.finish().unwrap();
    }

    // Write as .gz
    let input = temp.child("test.gz");
    input.write_binary(&buf).unwrap();

    let output = temp.child("out.vcf");
    let config = create_config(
        input.path().to_path_buf(),
        ref_path,
        output.path().to_path_buf(),
    );

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

#[test]
fn test_dtc_nested_compression() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let dtc_content = "rs1\t1\t1\tAC\n";

    // Create dtc.txt.gz.zip
    let input = create_nested_file(&temp, dtc_content, &["gz", "zip"], "test.txt");

    let output = temp.child("out.vcf");
    let config = create_config(input, ref_path, output.path().to_path_buf());

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}

#[test]
fn test_double_gz() {
    let temp = assert_fs::TempDir::new().unwrap();
    let ref_path = create_reference(&temp);
    let vcf_content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n";
    // test.vcf.gz.gz
    let input = create_nested_file(&temp, vcf_content, &["gz", "gz"], "test.vcf");

    let output = temp.child("out.vcf");
    let config = create_config(input, ref_path, output.path().to_path_buf());

    let summary = convert_dtc_file(config).unwrap();
    assert_eq!(summary.variant_records, 1);
}
