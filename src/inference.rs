//! Automatic inference of genome build and sample sex.
//!
//! Uses external libraries:
//! - `check_build` for genome build detection (hg19 vs hg38)
//! - `infer_sex` for biological sex inference from variant data

use std::path::Path;

use anyhow::{Context, Result};
use noodles::{bcf, vcf};

use crate::cli::Sex;
use crate::dtc::{self, Allele, Record as DtcRecord};
use crate::input::InputFormat;
use crate::smart_reader;

#[derive(Debug, Clone)]
pub struct BuildDetectionResult {
    pub detected_build: String,
    pub hg19_match_rate: f64,
    pub hg38_match_rate: f64,
}

/// Infer sex from DTC records.
///
/// This function processes DTC records and uses the infer_sex library
/// to determine biological sex based on X/Y chromosome variant patterns.
pub fn infer_sex_from_records(records: &[DtcRecord], build: &str) -> Result<Sex> {
    use infer_sex::{
        Chromosome, DecisionThresholds, GenomeBuild, InferenceConfig, InferredSex,
        PlatformDefinition, SexInferenceAccumulator, VariantInfo,
    };

    let genome_build = if build.contains("37") || build.to_lowercase().contains("hg19") {
        GenomeBuild::Build37
    } else {
        GenomeBuild::Build38
    };

    // Count attempted loci by chromosome type
    let mut n_autosomes = 0u64;
    let mut n_y_nonpar = 0u64;

    let constants = genome_build.algorithm_constants();

    for rec in records {
        let chrom = classify_chromosome(&rec.chromosome);
        match chrom {
            Chromosome::Autosome => n_autosomes += 1,
            Chromosome::Y => {
                if constants.is_in_y_non_par(rec.position) {
                    n_y_nonpar += 1;
                }
            }
            Chromosome::X => {}
        }
    }

    // Need at least some autosomal variants
    if n_autosomes == 0 {
        anyhow::bail!("No autosomal variants found - cannot infer sex");
    }

    let config = InferenceConfig {
        build: genome_build,
        platform: PlatformDefinition {
            n_attempted_autosomes: n_autosomes,
            n_attempted_y_nonpar: n_y_nonpar.max(1), // Avoid division by zero
        },
        thresholds: Some(DecisionThresholds::default()),
    };

    let mut acc = SexInferenceAccumulator::new(config);

    for rec in records {
        if rec.is_missing() {
            continue;
        }

        let chrom = classify_chromosome(&rec.chromosome);
        let is_het = is_heterozygous(&rec.genotype);

        acc.process_variant(&VariantInfo {
            chrom,
            pos: rec.position,
            is_heterozygous: is_het,
        });
    }

    let result = acc
        .finish()
        .map_err(|e| anyhow::anyhow!("Sex inference failed: {:?}", e))?;

    tracing::info!(
        "Sex inference: {:?} (Y density: {:?}, X/auto ratio: {:?})",
        result.final_call,
        result.report.y_genome_density,
        result.report.x_autosome_het_ratio
    );

    match result.final_call {
        InferredSex::Male => Ok(Sex::Male),
        InferredSex::Female => Ok(Sex::Female),
        InferredSex::Indeterminate => Ok(Sex::Unknown),
    }
}

fn classify_build(detection: check_build::BuildResult) -> BuildDetectionResult {
    let detected_build = match detection.better_match() {
        Some(check_build::Reference::Hg19) => "GRCh37".to_string(),
        Some(check_build::Reference::Hg38) | None => "GRCh38".to_string(),
    };

    BuildDetectionResult {
        detected_build,
        hg19_match_rate: detection.hg19_match_rate,
        hg38_match_rate: detection.hg38_match_rate,
    }
}

fn detect_build_from_variants(
    variants: &[check_build::Variant],
) -> Result<Option<BuildDetectionResult>> {
    if variants.is_empty() {
        return Ok(None);
    }

    // Ensure we have cached references (downloads if needed, uses cache otherwise)
    let refs_dir = crate::source_ref::convert_genome_refs_dir()?;
    let hg19_path = refs_dir.join("hg19.fa");
    let hg38_path = refs_dir.join("hg38.fa");

    // Pre-cache both references if they don't exist
    // This downloads once and reuses forever
    if !hg19_path.exists() || !hg38_path.exists() {
        tracing::info!("Caching reference genomes for build detection");
        crate::source_ref::load_source_reference("GRCh37")?;
        crate::source_ref::load_source_reference("GRCh38")?;
    }

    // Call check_build with our cached paths - no downloads, no MD5 checks
    let result = check_build::detect_build_from_positions_with_refs(
        variants,
        hg19_path.to_string_lossy().as_ref(),
        hg38_path.to_string_lossy().as_ref(),
    )?;

    tracing::info!(
        "Build detection: hg19={:.1}%, hg38={:.1}%",
        result.hg19_match_rate,
        result.hg38_match_rate
    );

    Ok(Some(classify_build(result)))
}

/// Detect genome build from a VCF file using check_build.
///
/// Returns "GRCh37" or "GRCh38" based on reference allele matching.
pub fn detect_build_from_vcf(vcf_path: &Path) -> Result<Option<BuildDetectionResult>> {
    let reader = smart_reader::open_input(vcf_path)
        .with_context(|| format!("failed to open input {}", vcf_path.display()))?;
    let mut vcf_reader = vcf::io::Reader::new(reader);
    let header = vcf_reader.read_header()?;

    let mut variants = Vec::with_capacity(1000);
    let mut records_seen = 0usize;

    for result in vcf_reader.record_bufs(&header) {
        let record = match result {
            Ok(record) => record,
            Err(e) => {
                tracing::warn!("failed to read VCF record: {}", e);
                continue;
            }
        };

        records_seen += 1;
        if variants.len() >= 1000 {
            break;
        }

        let chrom = record.reference_sequence_name().to_uppercase();
        if matches!(
            chrom.as_str(),
            "X" | "Y" | "CHRX" | "CHRY" | "MT" | "CHRM"
        ) {
            continue;
        }

        let pos = match record.variant_start().map(|p| usize::from(p) as u64) {
            Some(pos) if pos > 0 => pos,
            _ => continue,
        };

        let ref_base = record.reference_bases().to_ascii_uppercase();
        if ref_base.len() != 1 {
            continue;
        }
        let base = ref_base.chars().next().unwrap();
        if !matches!(base, 'A' | 'C' | 'G' | 'T') {
            continue;
        }

        variants.push(check_build::Variant {
            chrom: record.reference_sequence_name().to_string(),
            pos,
            ref_base,
        });
    }

    if variants.len() < 100 {
        tracing::warn!(
            variants = variants.len(),
            records = records_seen,
            "Low variant count for build detection"
        );
    }

    detect_build_from_variants(&variants)
}

pub fn detect_build_from_bcf(bcf_path: &Path) -> Result<Option<BuildDetectionResult>> {
    let reader = smart_reader::open_input(bcf_path)
        .with_context(|| format!("failed to open input {}", bcf_path.display()))?;
    let mut bcf_reader = bcf::io::Reader::new(reader);
    let header = bcf_reader.read_header()?;

    let mut variants = Vec::with_capacity(1000);
    let mut records_seen = 0usize;

    for result in bcf_reader.record_bufs(&header) {
        let record = match result {
            Ok(record) => record,
            Err(e) => {
                tracing::warn!("failed to read BCF record: {}", e);
                continue;
            }
        };

        records_seen += 1;
        if variants.len() >= 1000 {
            break;
        }

        let chrom = record.reference_sequence_name().to_uppercase();
        if matches!(
            chrom.as_str(),
            "X" | "Y" | "CHRX" | "CHRY" | "MT" | "CHRM"
        ) {
            continue;
        }

        let pos = match record.variant_start().map(|p| usize::from(p) as u64) {
            Some(pos) if pos > 0 => pos,
            _ => continue,
        };

        let ref_base = record.reference_bases().to_ascii_uppercase();
        if ref_base.len() != 1 {
            continue;
        }
        let base = ref_base.chars().next().unwrap();
        if !matches!(base, 'A' | 'C' | 'G' | 'T') {
            continue;
        }

        variants.push(check_build::Variant {
            chrom: record.reference_sequence_name().to_string(),
            pos,
            ref_base,
        });
    }

    if variants.len() < 100 {
        tracing::warn!(
            variants = variants.len(),
            records = records_seen,
            "Low variant count for build detection"
        );
    }

    detect_build_from_variants(&variants)
}

pub fn detect_build_from_variant_file(
    path: &Path,
    format: InputFormat,
) -> Result<Option<BuildDetectionResult>> {
    match format {
        InputFormat::Vcf => detect_build_from_vcf(path),
        InputFormat::Bcf => detect_build_from_bcf(path),
        _ => anyhow::bail!("Build detection requires VCF or BCF input"),
    }
}

/// Detect build from DTC records using check_build.
///
/// Builds Variant list from records (positions only), then checks against hg19/hg38.
pub fn detect_build_from_dtc(records: &[DtcRecord]) -> Result<BuildDetectionResult> {
    // Build Variant list for check_build (max 1000 for speed)
    let mut variants = Vec::with_capacity(1000);

    for rec in records.iter().take(5000) {
        if variants.len() >= 1000 {
            break;
        }

        // Skip sex chromosomes for build detection (less reliable)
        let chrom_upper = rec.chromosome.to_uppercase();
        if chrom_upper == "X"
            || chrom_upper == "Y"
            || chrom_upper == "CHRX"
            || chrom_upper == "CHRY"
            || chrom_upper == "MT"
            || chrom_upper == "CHRM"
        {
            continue;
        }

        // Parse genotype to find homozygous calls to use as candidate ref bases
        let alleles = dtc::parse_genotype(&rec.genotype);

        // We only use homozygous calls (e.g. AA, TT) because statistically
        // the allele is likely to match the reference.
        // We skip heterozygous calls because we don't know which is ref.
        if alleles.len() == 2 {
            if let (Allele::Base(b1), Allele::Base(b2)) = (&alleles[0], &alleles[1]) {
                if b1 == b2 && b1.len() == 1 {
                    variants.push(check_build::Variant {
                        chrom: rec.chromosome.clone(),
                        pos: rec.position,
                        ref_base: b1.clone(),
                    });
                }
            }
        }
    }

    if variants.len() < 100 {
        tracing::warn!("Only {} variants for build detection", variants.len());
    }

    detect_build_from_variants(&variants)?.ok_or_else(|| {
        anyhow::anyhow!("No variants available for build detection")
    })
}

/// Classify chromosome string into Chromosome enum for infer_sex.
fn classify_chromosome(chrom: &str) -> infer_sex::Chromosome {
    let normalized = chrom
        .trim()
        .to_uppercase()
        .trim_start_matches("CHR")
        .to_string();

    match normalized.as_str() {
        "X" | "23" | "25" | "PAR" | "PSEUDOAUTOSOMAL" => infer_sex::Chromosome::X,
        "Y" | "24" => infer_sex::Chromosome::Y,
        "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9" | "10" | "11" | "12" | "13" | "14"
        | "15" | "16" | "17" | "18" | "19" | "20" | "21" | "22" => infer_sex::Chromosome::Autosome,
        "M" | "MT" | "MITO" | "MITOCHONDRIAL" | "26" => infer_sex::Chromosome::Autosome,
        _ => infer_sex::Chromosome::Autosome, // Treat unknown as autosome (conservative)
    }
}

/// Check if a genotype string represents a heterozygous call.
fn is_heterozygous(genotype: &str) -> bool {
    let alleles = dtc::parse_genotype(genotype);

    if alleles.len() != 2 {
        return false;
    }

    // Heterozygous if the two alleles are different (and neither is missing)
    match (&alleles[0], &alleles[1]) {
        (Allele::Base(a), Allele::Base(b)) => !a.eq_ignore_ascii_case(b),
        (Allele::Missing, _) | (_, Allele::Missing) => false,
        (a, b) => a != b,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_chromosome() {
        assert_eq!(classify_chromosome("1"), infer_sex::Chromosome::Autosome);
        assert_eq!(
            classify_chromosome("chr22"),
            infer_sex::Chromosome::Autosome
        );
        assert_eq!(classify_chromosome("X"), infer_sex::Chromosome::X);
        assert_eq!(classify_chromosome("chrY"), infer_sex::Chromosome::Y);
        assert_eq!(classify_chromosome("23"), infer_sex::Chromosome::X);
        assert_eq!(classify_chromosome("24"), infer_sex::Chromosome::Y);
        assert_eq!(classify_chromosome("25"), infer_sex::Chromosome::X);
        assert_eq!(classify_chromosome("26"), infer_sex::Chromosome::Autosome);
        assert_eq!(classify_chromosome("chrM"), infer_sex::Chromosome::Autosome);
    }

    #[test]
    fn test_is_heterozygous() {
        assert!(is_heterozygous("AG"));
        assert!(is_heterozygous("A/G"));
        assert!(!is_heterozygous("AA"));
        assert!(!is_heterozygous("--"));
        assert!(!is_heterozygous("A"));
    }

    #[test]
    fn test_infer_sex_indeterminate_without_sex_chromosomes() {
        let records = vec![
            DtcRecord {
                id: Some("rs1".to_string()),
                chromosome: "1".to_string(),
                position: 100,
                genotype: "AA".to_string(),
            },
            DtcRecord {
                id: Some("rs2".to_string()),
                chromosome: "2".to_string(),
                position: 200,
                genotype: "GG".to_string(),
            },
            DtcRecord {
                id: Some("rs3".to_string()),
                chromosome: "3".to_string(),
                position: 300,
                genotype: "TT".to_string(),
            },
        ];

        let sex = infer_sex_from_records(&records, "GRCh38").unwrap();
        assert_eq!(sex, Sex::Unknown);
    }
}
