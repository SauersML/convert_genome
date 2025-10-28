use std::{fs::File, io::Write, sync::Arc};

use convert_genome::{
    ConversionConfig, OutputFormat, convert_dtc_file, dtc, reference::ReferenceGenome,
};
use proptest::prelude::*;
use tempfile::tempdir;

fn create_reference() -> (tempfile::TempDir, std::path::PathBuf) {
    let temp = tempdir().unwrap();
    let fasta_path = temp.path().join("ref.fa");
    let mut file = File::create(&fasta_path).unwrap();
    writeln!(file, ">1").unwrap();
    writeln!(file, "ACGTACGT").unwrap();
    drop(file);
    (temp, fasta_path)
}

fn run_conversion(fasta: &std::path::Path, genotype: &str) -> convert_genome::ConversionSummary {
    let dtc_path = fasta.with_file_name("input.txt");
    std::fs::write(&dtc_path, format!("rs1\t1\t2\t{}\n", genotype)).unwrap();
    let output = fasta.with_file_name("out.vcf");
    let config = ConversionConfig {
        input: dtc_path.clone(),
        input_origin: dtc_path.to_string_lossy().into_owned(),
        reference_fasta: fasta.to_path_buf(),
        reference_origin: fasta.to_string_lossy().into_owned(),
        reference_fai: None,
        reference_fai_origin: None,
        output: output.clone(),
        output_format: OutputFormat::Vcf,
        sample_id: String::from("sample"),
        assembly: String::from("test"),
        include_reference_sites: true,
    };
    convert_dtc_file(config).unwrap()
}

proptest! {
    #[test]
    fn dtc_reader_handles_random_input(lines in proptest::collection::vec(
        proptest::string::string_regex("[A-Za-z0-9\\t -]{0,32}").unwrap(),
        0..32
    )) {
        let joined = if lines.is_empty() {
            String::new()
        } else {
            lines.join("\n")
        };
        let cursor = std::io::Cursor::new(joined);
        let mut reader = dtc::Reader::new(cursor);
        while let Some(result) = reader.next() {
            let _ = result;
        }
    }
}

proptest! {
    #[test]
    fn conversion_handles_varied_genotypes(genotype in proptest::collection::vec(
        prop_oneof![
            Just('A'), Just('C'), Just('G'), Just('T'), Just('a'), Just('c'), Just('g'), Just('t'),
            Just('-'), Just('N'), Just('0')
        ], 1..=4
    )) {
        let (temp, fasta) = create_reference();
        let genotype_string: String = genotype.into_iter().collect();
        let summary = run_conversion(&fasta, &genotype_string);

        let valid = genotype_string.chars().all(|c| matches!(c, 'A'|'C'|'G'|'T'|'a'|'c'|'g'|'t'|'-'));
        if valid {
            prop_assert_eq!(summary.invalid_genotypes, 0);
        } else {
            prop_assert!(summary.invalid_genotypes <= 1);
        }

        drop(temp);
    }
}

proptest! {
    #[test]
    fn reference_lookup_matches_sequence(position in 0u64..12) {
        let (temp, fasta) = create_reference();
        let reference = ReferenceGenome::open(&fasta, None).unwrap();
        let sequence: Vec<char> = "ACGTACGT".chars().collect();

        if position == 0 || position as usize > sequence.len() {
            let result = reference.base("1", position);
            prop_assert!(result.is_err());
        } else {
            let result = reference.base("1", position).unwrap();
            prop_assert_eq!(result, sequence[(position as usize) - 1]);
        }

        drop(temp);
    }
}

proptest! {
    #[test]
    fn reference_is_safe_under_parallel_access(positions in proptest::collection::vec(1u64..=8, 8..32)) {
        let (temp, fasta) = create_reference();
        let reference = ReferenceGenome::open(&fasta, None).unwrap();
        let sequence: Arc<Vec<char>> = Arc::new("ACGTACGT".chars().collect());

        std::thread::scope(|scope| {
            for chunk in positions.chunks(4) {
                let clone = reference.clone();
                let seq = Arc::clone(&sequence);
                let chunk = chunk.to_vec();
                scope.spawn(move || {
                    for &pos in &chunk {
                        let base = clone.base("1", pos).unwrap();
                        assert_eq!(base, seq[(pos as usize) - 1]);
                    }
                });
            }
        });

        drop(temp);
    }
}
