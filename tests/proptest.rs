use std::{
    io::Cursor,
    sync::Arc,
    time::{SystemTime, UNIX_EPOCH},
};

use convert_genome::{conversion::format_genotype_for_tests, dtc, reference::ReferenceGenome};
use proptest::prelude::*;
use rayon::{ThreadPoolBuilder, prelude::*};

fn build_reference(sequence: &str) -> ReferenceGenome {
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("convert_genome_prop_{unique}"));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("ref.fa");
    std::fs::write(&path, format!(">chr1\n{}\n", sequence)).unwrap();
    ReferenceGenome::open(&path, None).unwrap()
}

proptest! {
    #[test]
    fn reader_handles_arbitrary_input(data in proptest::collection::vec(any::<u8>(), 0..1024)) {
        let cursor = Cursor::new(data);
        let mut reader = dtc::Reader::new(cursor);
        while let Some(record) = reader.next() {
            let _ = record;
        }
    }
}

proptest! {
    #[test]
    fn genotype_parsing_round_trips(
        reference_base in prop::sample::select(vec!['A', 'C', 'G', 'T']),
        alleles in proptest::collection::vec(prop::sample::select(vec!['A', 'C', 'G', 'T', '-']), 1..3),
    ) {
        let parsed: Vec<Option<char>> = alleles
            .iter()
            .map(|&c| if c == '-' { None } else { Some(c) })
            .collect();

        let mut alts = Vec::new();
        for allele in parsed.iter().flatten() {
            if *allele != reference_base && !alts.contains(allele) {
                alts.push(*allele);
            }
        }

        let formatted = format_genotype_for_tests(&parsed, reference_base, &alts)
            .expect("valid genotype formatting");
        prop_assert!(!formatted.is_empty());
        prop_assert!(formatted
            .chars()
            .all(|c: char| c.is_ascii_digit() || matches!(c, '.' | '/' )));
    }
}

proptest! {
    #[test]
    fn invalid_alleles_are_rejected(
        reference_base in prop::sample::select(vec!['A', 'C', 'G', 'T']),
        invalid in prop::sample::select(vec!['X', 'Y', 'Z']),
    ) {
        let alleles = vec![Some(invalid)];
        let result = format_genotype_for_tests(&alleles, reference_base, &[]);
        prop_assert!(result.is_err());
    }
}

proptest! {
    #[test]
    fn reference_lookup_returns_bases(position in 1u64..=16) {
        let reference = build_reference("ACGTACGTACGTACGT");
        let base = reference.base("1", position).unwrap();
        prop_assert!(matches!(base, 'A' | 'C' | 'G' | 'T'));
    }
}

proptest! {
    #[test]
    fn concurrent_reference_access(positions in proptest::collection::vec(1u64..=16, 2..16)) {
        let reference = Arc::new(build_reference("ACGTACGTACGTACGT"));
        let pool = ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let results = pool.install(|| {
            positions
                .par_iter()
                .map(|pos| reference.base("1", *pos).unwrap())
                .collect::<Vec<_>>()
        });

        for base in results {
            prop_assert!(matches!(base, 'A' | 'C' | 'G' | 'T'));
        }
    }
}
