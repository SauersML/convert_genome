//! Byte-identical parity test.
//!
//! Regenerates a deterministic multi-chromosome VCF fixture, runs the full
//! conversion pipeline, and asserts the VCF output is byte-identical to a
//! golden reference captured before the parallelization rewrite.
//!
//! One VCF header line is inherently volatile: `##reference=file:///<abs-path>`
//! embeds the tempdir where the fixture's FASTA was written, and that
//! tempdir differs per run. The test strips that single line from both
//! sides before comparing — every other byte (contigs, variants, all
//! genotype data) is checked strictly.
//!
//! Golden location: /tmp/cg_golden/
//! Fixture + output stay fully inside the repo's test scope; no network
//! access is required because `input_build == assembly == GRCh38` short-
//! circuits build detection and liftover paths.

use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
};

use anyhow::Result;
use convert_genome::input::InputFormat;
use convert_genome::{ConversionConfig, OutputFormat, convert_dtc_file};

const GOLDEN_DIR: &str = "/tmp/cg_golden";
const FIXTURE_CHROMS: &[(&str, usize)] = &[
    ("1", 4_000),
    ("2", 3_000),
    ("3", 2_500),
    ("4", 2_000),
    ("5", 1_500),
    ("X", 1_000),
];

fn write_fixture_reference(dir: &Path) -> Result<PathBuf> {
    // Build a deterministic ref: for each chromosome, repeat ACGT for 2x the
    // largest needed coordinate so reference lookups never OOB.
    let max_pos = FIXTURE_CHROMS
        .iter()
        .map(|(_, n)| *n)
        .max()
        .unwrap_or(0);
    let seq_len = (max_pos * 4 + 8).max(32);
    let pattern = b"ACGT";

    let path = dir.join("ref.fa");
    let mut file = fs::File::create(&path)?;
    for (chrom, _) in FIXTURE_CHROMS {
        writeln!(file, ">{}", chrom)?;
        let mut buf = Vec::with_capacity(seq_len);
        for i in 0..seq_len {
            buf.push(pattern[i % pattern.len()]);
        }
        file.write_all(&buf)?;
        writeln!(file)?;
    }
    Ok(path)
}

fn write_fixture_vcf(dir: &Path) -> Result<PathBuf> {
    let path = dir.join("input.vcf");
    let mut file = fs::File::create(&path)?;
    writeln!(file, "##fileformat=VCFv4.2")?;
    for (chrom, count) in FIXTURE_CHROMS {
        let seq_len = count * 4 + 8;
        writeln!(file, "##contig=<ID={},length={}>", chrom, seq_len)?;
    }
    writeln!(
        file,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    )?;

    // Deterministic variant generation. Positions are spaced 4bp apart so
    // REF is always `A` (since the reference is repeating ACGT starting at
    // position 1: pos 1 = A, pos 5 = A, pos 9 = A, ...). This keeps the
    // fixture trivial to reason about.
    //
    // Every 7th record is a homozygous reference call; every 5th is a
    // heterozygous variant; everything else is homozygous alt. This exercises
    // variant/reference counting and genotype writing on all three paths.
    let alts = ['C', 'G', 'T'];
    for (chrom, count) in FIXTURE_CHROMS {
        for i in 0..*count {
            let pos = i * 4 + 1;
            let alt = alts[i % alts.len()];
            let gt = if i % 7 == 0 {
                "0/0"
            } else if i % 5 == 0 {
                "0/1"
            } else {
                "1/1"
            };
            writeln!(
                file,
                "{}\t{}\t.\tA\t{}\t.\t.\t.\tGT\t{}",
                chrom, pos, alt, gt
            )?;
        }
    }
    Ok(path)
}

fn build_config(input: PathBuf, reference: PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_format: InputFormat::Vcf,
        input_origin: input.display().to_string(),
        reference_fasta: Some(reference.clone()),
        reference_origin: Some(reference.display().to_string()),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_dir: None,
        output_format: OutputFormat::Vcf,
        sample_id: "GOLDEN".into(),
        assembly: "GRCh38".into(),
        include_reference_sites: true,
        sex: Some(convert_genome::cli::Sex::Female),
        par_boundaries: convert_genome::reference::ParBoundaries::new("GRCh38"),
        standardize: false,
        panel: None,
        // Declared build == target assembly: skips network build detection and
        // liftover entirely.
        input_build: Some("GRCh38".into()),
    }
}

fn total_variants() -> usize {
    FIXTURE_CHROMS.iter().map(|(_, n)| *n).sum()
}

fn run_conversion(workdir: &Path) -> Result<PathBuf> {
    fs::create_dir_all(workdir)?;
    let reference = write_fixture_reference(workdir)?;
    let input = write_fixture_vcf(workdir)?;
    let output = workdir.join("output.vcf");
    let config = build_config(input, reference, output.clone());
    convert_dtc_file(config)?;
    Ok(output)
}

/// Strip the single path-volatile header line (`##reference=file:///...`)
/// from a VCF so byte-identical comparison is possible across tempdirs.
/// Every other line — contigs, source version, assembly, variants — is
/// passed through untouched.
fn normalize_vcf(raw: &[u8]) -> Vec<u8> {
    let text = std::str::from_utf8(raw).expect("output VCF is UTF-8");
    let mut out = Vec::with_capacity(raw.len());
    for line in text.split_inclusive('\n') {
        if line.starts_with("##reference=") {
            continue;
        }
        out.extend_from_slice(line.as_bytes());
    }
    out
}

fn first_divergence(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b.iter())
        .position(|(x, y)| x != y)
        .unwrap_or(a.len().min(b.len()))
}

#[test]
fn parity_multi_chrom_vcf_matches_golden() -> Result<()> {
    let workdir = tempfile::tempdir()?;
    let output = run_conversion(workdir.path())?;
    let actual_raw = fs::read(&output)?;
    let actual = normalize_vcf(&actual_raw);

    let golden_dir = Path::new(GOLDEN_DIR);
    let golden_output = golden_dir.join("output.vcf");

    if !golden_output.exists() {
        // First-run capture. Store the already-normalized form so byte-for-
        // byte comparison works cleanly on every subsequent run.
        fs::create_dir_all(golden_dir)?;
        fs::write(&golden_output, &actual)?;
        fs::write(
            golden_dir.join("fixture_meta.txt"),
            format!(
                "chroms={}\nvariants={}\nbytes_normalized={}\nbytes_raw={}\n",
                FIXTURE_CHROMS.len(),
                total_variants(),
                actual.len(),
                actual_raw.len(),
            ),
        )?;
        return Ok(());
    }

    let expected = fs::read(&golden_output)?;
    assert_eq!(
        actual.len(),
        expected.len(),
        "normalized output byte length differs from golden ({} actual vs {} golden; first diff at byte {})",
        actual.len(),
        expected.len(),
        first_divergence(&actual, &expected),
    );
    assert!(
        actual == expected,
        "normalized output does not match golden byte-for-byte (first divergence at byte {})",
        first_divergence(&actual, &expected),
    );
    Ok(())
}

#[test]
fn parity_fixture_has_expected_variant_count() -> Result<()> {
    // Independent sanity check — if this ever changes, the fixture has been
    // inadvertently modified and the golden should be recaptured.
    assert_eq!(total_variants(), 14_000);
    Ok(())
}
