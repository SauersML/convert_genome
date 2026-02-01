use std::path::PathBuf;

use anyhow::{Context, Result};
use assert_fs::{TempDir, prelude::*};
use convert_genome::conversion::ConversionConfig;
use convert_genome::liftover::ChainRegistry;
use convert_genome::reference::ReferenceGenome;
use convert_genome::{OutputFormat, convert_dtc_file};
use serial_test::serial;

fn ensure_reference_paths() -> Result<(PathBuf, PathBuf)> {
    convert_genome::source_ref::load_source_reference("GRCh37")?;
    convert_genome::source_ref::load_source_reference("GRCh38")?;
    let refs_dir = convert_genome::source_ref::convert_genome_refs_dir()?;
    Ok((refs_dir.join("hg19.fa"), refs_dir.join("hg38.fa")))
}

fn other_bases(base: char) -> Vec<char> {
    ['A', 'C', 'G', 'T']
        .into_iter()
        .filter(|b| *b != base)
        .collect()
}

fn find_positions_with_base_diff(
    hg19: &ReferenceGenome,
    hg38: &ReferenceGenome,
    chain: Option<&convert_genome::liftover::ChainMap>,
    chrom: &str,
    target: usize,
) -> Result<Vec<u64>> {
    let mut positions = Vec::with_capacity(target);
    for pos in 1u64..=1_000_000u64 {
        if positions.len() >= target {
            break;
        }

        if let Some(chain) = chain {
            if chain.lift(chrom, pos - 1).is_err() {
                continue;
            }
        }

        let base_19 = hg19.base(chrom, pos)?.to_ascii_uppercase();
        let base_38 = hg38.base(chrom, pos)?.to_ascii_uppercase();
        if !matches!(base_19, 'A' | 'C' | 'G' | 'T') {
            continue;
        }
        if !matches!(base_38, 'A' | 'C' | 'G' | 'T') {
            continue;
        }
        if base_19 != base_38 {
            positions.push(pos);
        }
    }

    if positions.len() < target {
        anyhow::bail!(
            "Unable to find enough divergent positions on {} (found {})",
            chrom,
            positions.len()
        );
    }

    Ok(positions)
}

fn write_vcf(path: &std::path::Path, chrom: &str, positions: &[u64], refs: &[char]) -> Result<()> {
    let mut contents = String::new();
    contents.push_str("##fileformat=VCFv4.2\n");
    contents.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n");

    for (pos, ref_base) in positions.iter().zip(refs.iter()) {
        let alts: String = other_bases(*ref_base)
            .into_iter()
            .map(|b| b.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let id = format!("rs{}", pos);
        contents.push_str(&format!(
            "{chrom}\t{pos}\t{id}\t{ref_base}\t{alts}\t.\t.\t.\tGT\t0/1\n",
            chrom = chrom,
            pos = pos,
            id = id,
            ref_base = ref_base,
            alts = alts
        ));
    }

    std::fs::write(path, contents)?;
    Ok(())
}

fn base_config(input: PathBuf, reference: PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_format: convert_genome::input::InputFormat::Vcf,
        input_origin: input.to_string_lossy().to_string(),
        reference_fasta: Some(reference.clone()),
        reference_origin: Some(reference.to_string_lossy().to_string()),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_dir: None,
        output_format: OutputFormat::Vcf,
        sample_id: "SAMPLE".into(),
        assembly: "GRCh38".into(),
        include_reference_sites: true,
        sex: None,
        par_boundaries: convert_genome::reference::ParBoundaries::new("GRCh38"),
        standardize: false,
        panel: None,
    }
}

fn find_record_by_id(
    data: &str,
    id: &str,
) -> Option<(String, u64)> {
    for line in data.lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = fields.next()?.to_string();
        let pos = fields.next()?.parse::<u64>().ok()?;
        let record_id = fields.next()?;
        if record_id == id {
            return Some((chrom, pos));
        }
    }
    None
}

#[test]
#[serial]
fn vcf_build_detection_triggers_liftover() -> Result<()> {
    let temp = TempDir::new()?;
    let (hg19_path, hg38_path) = ensure_reference_paths()?;
    let hg19 = ReferenceGenome::open(&hg19_path, None)?;
    let hg38 = ReferenceGenome::open(&hg38_path, None)?;

    let chain = ChainRegistry::new()?.get_chain("GRCh37", "GRCh38")?;
    let chrom = "chr1";
    let positions = find_positions_with_base_diff(&hg19, &hg38, Some(&chain), chrom, 20)?;
    let refs: Vec<char> = positions
        .iter()
        .map(|pos| hg19.base(chrom, *pos).unwrap().to_ascii_uppercase())
        .collect();

    let input_vcf = temp.child("input_hg19.vcf");
    write_vcf(input_vcf.path(), chrom, &positions, &refs)?;

    let output_vcf = temp.child("output_hg38.vcf");
    let config = base_config(
        input_vcf.path().to_path_buf(),
        hg38_path.clone(),
        output_vcf.path().to_path_buf(),
    );
    let summary = convert_dtc_file(config)?;
    assert_eq!(summary.emitted_records, positions.len());

    let first_pos = positions[0];
    let (expected_chrom, expected_pos_0, _) = chain
        .lift(chrom, first_pos - 1)
        .map_err(|e| anyhow::anyhow!("failed to lift expected position: {:?}", e))?;
    let expected_pos = expected_pos_0 + 1;
    let expected_id = format!("rs{}", first_pos);

    let output_data = std::fs::read_to_string(output_vcf.path())?;
    let (actual_chrom, actual_pos) =
        find_record_by_id(&output_data, &expected_id).context("missing lifted record")?;

    assert_eq!(actual_chrom, expected_chrom);
    assert_eq!(actual_pos, expected_pos);

    Ok(())
}

#[test]
#[serial]
fn vcf_build_detection_skips_liftover_when_matching() -> Result<()> {
    let temp = TempDir::new()?;
    let (hg19_path, hg38_path) = ensure_reference_paths()?;
    let hg19 = ReferenceGenome::open(&hg19_path, None)?;
    let hg38 = ReferenceGenome::open(&hg38_path, None)?;

    let chrom = "chr1";
    let positions = find_positions_with_base_diff(&hg19, &hg38, None, chrom, 20)?;
    let refs: Vec<char> = positions
        .iter()
        .map(|pos| hg38.base(chrom, *pos).unwrap().to_ascii_uppercase())
        .collect();

    let input_vcf = temp.child("input_hg38.vcf");
    write_vcf(input_vcf.path(), chrom, &positions, &refs)?;

    let output_vcf = temp.child("output_hg38.vcf");
    let config = base_config(
        input_vcf.path().to_path_buf(),
        hg38_path.clone(),
        output_vcf.path().to_path_buf(),
    );
    let summary = convert_dtc_file(config)?;
    assert_eq!(summary.emitted_records, positions.len());

    let first_pos = positions[0];
    let expected_id = format!("rs{}", first_pos);
    let output_data = std::fs::read_to_string(output_vcf.path())?;
    let (actual_chrom, actual_pos) =
        find_record_by_id(&output_data, &expected_id).context("missing record")?;

    assert_eq!(actual_chrom, chrom);
    assert_eq!(actual_pos, first_pos);

    Ok(())
}
