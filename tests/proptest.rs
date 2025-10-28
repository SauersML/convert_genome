use std::{fs, path::PathBuf};

use convert_genome::{
    conversion::{ConversionConfig, OutputFormat, convert_dtc_file},
    dtc,
    reference::ReferenceGenome,
};
use proptest::prelude::*;
use tempfile::tempdir;

fn setup_reference(sequence: &str) -> (tempfile::TempDir, PathBuf) {
    let dir = tempdir().unwrap();
    let path = dir.path().join("ref.fa");
    fs::write(&path, format!(">chr1\n{}\n", sequence)).unwrap();
    (dir, path)
}

fn base_config(input: PathBuf, reference: &PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_origin: input.display().to_string(),
        reference_fasta: reference.clone(),
        reference_origin: reference.display().to_string(),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_format: OutputFormat::Vcf,
        sample_id: String::from("sample"),
        assembly: String::from(""),
        include_reference_sites: true,
    }
}

proptest! {
    #[test]
    fn dtc_reader_handles_arbitrary_input(raw in prop::collection::vec(any::<u8>(), 0..256)) {
        let mut reader = dtc::Reader::new(&raw[..]);
        while let Some(result) = reader.next() {
            let _ = result;
        }
    }
}

proptest! {
    #[test]
    fn genotype_combinations_convert(genotype in prop::collection::vec(prop::sample::select(vec!['A', 'C', 'G', 'T', '-']), 0..=2)) {
        let genotype_string: String = genotype.iter().collect();
        let (dir, reference_path) = setup_reference("AAAA");
        let input_path = dir.path().join("input.txt");
        fs::write(&input_path, format!("rs1\t1\t1\t{}\n", genotype_string)).unwrap();
        let output_path = dir.path().join("out.vcf");
        let config = base_config(input_path, &reference_path, output_path.clone());

        let summary = convert_dtc_file(config).expect("conversion");

        if genotype_string.trim().is_empty() {
            assert_eq!(summary.total_records, 0);
            assert_eq!(summary.parse_errors, 1);
        } else {
            assert_eq!(summary.total_records, 1);
            let mut expected_alt = Vec::new();
            for allele in genotype_string.chars() {
                if allele == '-' {
                    continue;
                }
                let upper = allele.to_ascii_uppercase();
                if upper != 'A' && !expected_alt.contains(&upper) {
                    expected_alt.push(upper);
                }
            }

            assert_eq!(summary.missing_genotype_records, 0);
            assert_eq!(summary.emitted_records, 1);
            assert_eq!(summary.variant_records, usize::from(!expected_alt.is_empty()));
            assert_eq!(summary.parse_errors, 0);
        }

        let _ = fs::remove_file(output_path);
    }
}

proptest! {
    #[test]
    fn reference_lookup_matches_sequence(pos in 1u64..=32) {
        let sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let (_dir, reference_path) = setup_reference(sequence);
        let reference = ReferenceGenome::open(&reference_path, None).unwrap();
        let expected = sequence.as_bytes()[usize::try_from(pos - 1).unwrap()] as char;
        let observed = reference.base("1", pos).unwrap();
        assert_eq!(observed, expected.to_ascii_uppercase());
        assert!(reference.cache_usage() >= 1);
    }
}

proptest! {
    #[test]
    fn concurrent_reference_access_is_consistent(positions in prop::collection::vec(1u64..=16, 1..16)) {
        let sequence = "ACGTACGTACGTACGT";
        let (_dir, reference_path) = setup_reference(sequence);
        let reference = ReferenceGenome::open(&reference_path, None).unwrap();

        std::thread::scope(|scope| {
            for pos in positions.iter().copied() {
                let clone = reference.clone();
                let expected = sequence.as_bytes()[usize::try_from(pos - 1).unwrap()] as char;
                scope.spawn(move || {
                    let base = clone.base("1", pos).unwrap();
                    assert_eq!(base, expected.to_ascii_uppercase());
                });
            }
        });

        assert!(reference.cache_usage() > 0);
    }
}
