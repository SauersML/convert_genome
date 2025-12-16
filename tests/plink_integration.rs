use convert_genome::cli::Sex;
use convert_genome::plink::PlinkWriter;
use noodles::core::Position;
use noodles::vcf::variant::record_buf::{
    AlternateBases, Ids, RecordBuf, Samples, samples::Keys, samples::sample::Value,
};
use std::fs;
use std::path::Path;
use tempfile::tempdir;

#[test]
fn test_plink_writer_spec_example() {
    let dir = tempdir().unwrap();
    let prefix = dir.path().join("test");

    let mut writer = PlinkWriter::new(&prefix).expect("failed to create writer");

    // Keys for Samples
    let keys: Keys = vec![String::from("GT")].into_iter().collect();

    // SNP 1: 1 1 0 0 1 0 (Genotypes: G G, A A, 0 0, A A, A A, A A)
    // Ref: G, Alt: A.
    // Spec: 00(HomRef), 11(HomAlt), 01(Missing), 11, 11, 11.
    // Byte: DC (11011100).
    // Byte 2: 0x0F (00001111). S5(11), S6(11).
    // Bits: 11 11.
    // 0000 1111 -> 0F. Correct.

    let samples = vec![
        vec![Some(Value::String(String::from("0/0")))],
        vec![Some(Value::String(String::from("1/1")))],
        vec![Some(Value::String(String::from("./.")))],
        vec![Some(Value::String(String::from("1/1")))],
        vec![Some(Value::String(String::from("1/1")))],
        vec![Some(Value::String(String::from("1/1")))],
    ];
    let samples_buf = Samples::new(keys.clone(), samples);

    let ids: Ids = [String::from("snp1")].into_iter().collect();

    let record1 = RecordBuf::builder()
        .set_reference_sequence_name("1")
        .set_ids(ids)
        .set_variant_start(Position::try_from(1).unwrap())
        .set_reference_bases("G")
        .set_alternate_bases(AlternateBases::from(vec![String::from("A")]))
        .set_samples(samples_buf)
        .build();

    // Snp 2: genotypes (HomAlt, Missing, Het, HomAlt, HomAlt, HomAlt)
    // Expected bytes: E7 0F

    let samples2 = vec![
        vec![Some(Value::String(String::from("1/1")))],
        vec![Some(Value::String(String::from("./.")))],
        vec![Some(Value::String(String::from("0/1")))],
        vec![Some(Value::String(String::from("1/1")))],
        vec![Some(Value::String(String::from("1/1")))],
        vec![Some(Value::String(String::from("1/1")))],
    ];
    let samples_buf2 = Samples::new(keys.clone(), samples2);

    let ids2: Ids = [String::from("snp2")].into_iter().collect();

    let record2 = RecordBuf::builder()
        .set_reference_sequence_name("1")
        .set_ids(ids2)
        .set_variant_start(Position::try_from(2).unwrap())
        .set_reference_bases("1")
        .set_alternate_bases(AlternateBases::from(vec![String::from("2")]))
        .set_samples(samples_buf2)
        .build();

    // SNP 3: 0 0 1 2 1 2
    // Bytes: 6B (01 10 10 11 -> S4 S3 S2 S1) (Missing Het Het HomAlt)

    let samples3 = vec![
        vec![Some(Value::String(String::from("1/1")))],
        vec![Some(Value::String(String::from("0/1")))],
        vec![Some(Value::String(String::from("0/1")))],
        vec![Some(Value::String(String::from("./.")))],
        vec![Some(Value::String(String::from("./.")))],
        vec![Some(Value::String(String::from("0/0")))],
    ];
    let samples_buf3 = Samples::new(keys.clone(), samples3);

    let ids3: Ids = [String::from("snp3")].into_iter().collect();

    let record3 = RecordBuf::builder()
        .set_reference_sequence_name("1")
        .set_ids(ids3)
        .set_variant_start(Position::try_from(3).unwrap())
        .set_reference_bases("A")
        .set_alternate_bases(AlternateBases::from(vec![String::from("C")]))
        .set_samples(samples_buf3)
        .build();

    writer.write_variant(&record1).unwrap();
    writer.write_variant(&record2).unwrap();
    writer.write_variant(&record3).unwrap();

    drop(writer);

    let bed_content = fs::read(prefix.with_extension("bed")).unwrap();
    let expected = vec![
        0x6C, 0x1B, 0x01, // Magic
        0xDC, 0x0F, // SNP1
        0xE7, 0x0F, // SNP2
        0x6B, 0x01, // SNP3
    ];
    assert_eq!(bed_content, expected);
}
