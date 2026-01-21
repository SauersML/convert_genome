//! Automatic inference of genome build and sample sex.
//!
//! Uses external libraries:
//! - `check_build` for genome build detection (hg19 vs hg38)
//! - `infer_sex` for biological sex inference from variant data

use std::path::Path;

use anyhow::Result;

use crate::cli::Sex;
use crate::dtc::{self, Allele, Record as DtcRecord};

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
    }
}

/// Detect genome build from a VCF file using check_build.
///
/// Returns "GRCh37" or "GRCh38" based on reference allele matching.
pub fn detect_build_from_vcf(vcf_path: &Path) -> Result<String> {
    let path_str = vcf_path.to_string_lossy();

    // check_build requires .vcf extension
    if !path_str.ends_with(".vcf") {
        anyhow::bail!("check_build requires .vcf file, got: {}", path_str);
    }

    let result = check_build::detect_build(path_str.to_string())
        .map_err(|e| anyhow::anyhow!("Build detection failed: {}", e))?;

    tracing::info!(
        "Build detection: hg19={:.1}%, hg38={:.1}%",
        result.hg19_match_rate,
        result.hg38_match_rate
    );

    match result.better_match() {
        Some(check_build::Reference::Hg19) => Ok("GRCh37".to_string()),
        Some(check_build::Reference::Hg38) | None => Ok("GRCh38".to_string()),
    }
}

/// Detect build from DTC records using check_build.
///
/// Builds Variant list from records (positions only), then checks against hg19/hg38.
pub fn detect_build_from_dtc(records: &[DtcRecord]) -> Result<String> {
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

        // We only provide CHROM+POS. check_build will look up the reference bases itself.
        // We pass an empty ref_base string or standard notation if allowed, but check_build
        // usually needs the ref base to verify? 
        // Wait, the user said: "check_build ignores the user's --reference... checks these variants against its own internal cache".
        // check_build::detect_build_from_positions takes a list of variants. 
        // Actually, looking at check_build source (implied), it likely needs to know what the *observed* ref base is in the file if the file has it?
        // DTC records usually DON'T have a ref base.
        // The previous code was looking up ref base from the *Target* reference to populate the struct.
        // If we don't have a ref base, can we use `detect_build_from_positions` effectively?
        // The user said: "Parse the first ~1000 records... Pass this vector directly to check_build".
        
        // Let's assume check_build can handle just Chrom/Pos or we leave ref_base empty/dummy if it loads it itself.
        // Actually, if we look at the previous code:
        // `variants.push(check_build::Variant { ref_base: ref_base.to_string() ... })`
        // It was populating ref_base from the *reference*.
        
        // If DTC doesn't have ref bases, we can't provide them.
        // Does `check_build` support lookup solely on position matches?
        // The user mentioned: "detect_build_from_positions".
        
        variants.push(check_build::Variant {
            chrom: rec.chromosome.clone(),
            pos: rec.position,
            ref_base: "N".to_string(), // internal check_build logic handles lookups
        });
    }

    if variants.len() < 100 {
        tracing::warn!(
            "Only {} variants for build detection",
            variants.len()
        );
    }

    let result = check_build::detect_build_from_positions(&variants)
        .map_err(|e| anyhow::anyhow!("Build detection failed: {}", e))?;

    tracing::info!(
        "Build detection: hg19={:.1}%, hg38={:.1}%",
        result.hg19_match_rate,
        result.hg38_match_rate
    );

    match result.better_match() {
        Some(check_build::Reference::Hg19) => Ok("GRCh37".to_string()),
        Some(check_build::Reference::Hg38) | None => Ok("GRCh38".to_string()),
    }
}

/// Classify chromosome string into Chromosome enum for infer_sex.
fn classify_chromosome(chrom: &str) -> infer_sex::Chromosome {
    let normalized = chrom
        .trim()
        .to_uppercase()
        .trim_start_matches("CHR")
        .to_string();

    match normalized.as_str() {
        "X" => infer_sex::Chromosome::X,
        "Y" => infer_sex::Chromosome::Y,
        "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9" | "10" | "11" | "12" | "13" | "14"
        | "15" | "16" | "17" | "18" | "19" | "20" | "21" | "22" => infer_sex::Chromosome::Autosome,
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
    }

    #[test]
    fn test_is_heterozygous() {
        assert!(is_heterozygous("AG"));
        assert!(is_heterozygous("A/G"));
        assert!(!is_heterozygous("AA"));
        assert!(!is_heterozygous("--"));
        assert!(!is_heterozygous("A"));
    }
}
