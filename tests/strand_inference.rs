use convert_genome::dtc::Record as DtcRecord;
use convert_genome::reference::ReferenceGenome;
use convert_genome::source_ref::{InferredStrand, infer_strand_lock};
use std::io::Write;

fn create_source_reference(dir: &tempfile::TempDir) -> ReferenceGenome {
    let fasta_path = dir.path().join("source.fa");
    let mut file = std::fs::File::create(&fasta_path).unwrap();
    // chr1: A at 100, G at 200
    writeln!(file, ">chr1").unwrap();
    // 0123456789...
    // Create a sequence long enough.
    // We'll put A at 100 and G at 200.
    // Let's just make a dummy sequence where we control specific positions.
    // Indexing is 0-based in file, 1-based in DtcRecord? ReferenceGenome::base uses 0-based if I remember correctly or converts?
    // Let's check ReferenceGenome::base implementation again. conversion.rs passes rec.position directly.
    // DtcRecord uses 1-based coordinates usually?
    // Wait, let's verify ReferenceGenome::base checks.
    // It says: "position {position} is outside".
    // "let pos = usize::try_from(position)..."
    // "let start = Position::try_from(pos)?;"
    // noodles Position is 1-based.
    // So ReferenceGenome::base expects 1-based coordinates.

    // We need 200 bases.
    let mut seq = vec!['N'; 300];
    seq[99] = 'A'; // Position 100 (1-based)
    seq[199] = 'G'; // Position 200

    // Write in chunks of 80 chars
    let seq_str: String = seq.into_iter().collect();
    writeln!(file, "{}", seq_str).unwrap();

    ReferenceGenome::open(&fasta_path, None).unwrap()
}

#[test]
fn test_infer_forward_strand() {
    let dir = tempfile::tempdir().unwrap();
    let reference = create_source_reference(&dir);

    let records = vec![
        DtcRecord {
            chromosome: "1".to_string(),
            position: 100,
            genotype: "AA".to_string(), // Matches Ref A
            id: None,
        },
        DtcRecord {
            chromosome: "1".to_string(),
            position: 200,
            genotype: "GG".to_string(), // Matches Ref G
            id: None,
        },
    ];

    let result = infer_strand_lock(&records, &reference).unwrap();
    assert_eq!(result, InferredStrand::Forward);
}

#[test]
fn test_infer_reverse_strand() {
    let dir = tempfile::tempdir().unwrap();
    let reference = create_source_reference(&dir);

    // Ref is A at 100, G at 200.
    // Reverse Complement of A is T.
    // Reverse Complement of G is C.

    // We need > 50 records to trigger inference
    let mut records = Vec::new();
    for _ in 0..30 {
        records.push(DtcRecord {
            chromosome: "1".to_string(),
            position: 100,
            genotype: "TT".to_string(), // Matches RC(A)
            id: None,
        });
        records.push(DtcRecord {
            chromosome: "1".to_string(),
            position: 200,
            genotype: "CC".to_string(), // Matches RC(G)
            id: None,
        });
    }

    let result = infer_strand_lock(&records, &reference).unwrap();
    assert_eq!(result, InferredStrand::Reverse);
}

#[test]
fn test_infer_ambiguous_skips() {
    let dir = tempfile::tempdir().unwrap();
    let reference = create_source_reference(&dir);

    // A/T SNPs are ambiguous (A -> T flip looks like T)
    // The code should skip these.

    let records = vec![DtcRecord {
        chromosome: "1".to_string(),
        position: 100,
        genotype: "AT".to_string(), // Transversion check should return false
        id: None,
    }];

    // If all skipped, it defaults to Forward with a warning (tested < 50)
    let result = infer_strand_lock(&records, &reference).unwrap();
    assert_eq!(result, InferredStrand::Forward);
}
