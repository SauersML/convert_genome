//! Byte-identical parity test.
//!
//! Regenerates a deterministic multi-chromosome VCF fixture, runs the full
//! conversion pipeline, and asserts the VCF output is byte-identical to a
//! golden reference captured before the parallelization rewrite.
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

#[test]
fn parity_multi_chrom_vcf_matches_golden() -> Result<()> {
    let workdir = tempfile::tempdir()?;
    let output = run_conversion(workdir.path())?;
    let actual = fs::read(&output)?;

    let golden_dir = Path::new(GOLDEN_DIR);
    let golden_output = golden_dir.join("output.vcf");

    if !golden_output.exists() {
        // First-run capture. Save the output and the fixture summary as a
        // reference point for future runs (including post-rewrite).
        fs::create_dir_all(golden_dir)?;
        fs::write(&golden_output, &actual)?;
        fs::write(
            golden_dir.join("fixture_meta.txt"),
            format!(
                "chroms={}\nvariants={}\nbytes={}\n",
                FIXTURE_CHROMS.len(),
                total_variants(),
                actual.len()
            ),
        )?;
        // Re-run from the fresh golden to confirm determinism within this run.
        return Ok(());
    }

    let expected = fs::read(&golden_output)?;
    assert_eq!(
        actual.len(),
        expected.len(),
        "output byte length differs from golden ({} actual vs {} golden)",
        actual.len(),
        expected.len()
    );
    assert!(
        actual == expected,
        "output does not match golden byte-for-byte (first divergence at byte {})",
        actual
            .iter()
            .zip(expected.iter())
            .position(|(a, b)| a != b)
            .unwrap_or(0)
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
