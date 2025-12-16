//! Panel harmonization logic for Beagle compatibility.
//!
//! This module contains functions to harmonize user genotype data against
//! a reference panel, ensuring REF/ALT encoding matches for Beagle imputation.

use crate::panel::{PaddedPanel, PanelSite};

/// Get the complement of a DNA base.
pub fn complement(base: char) -> char {
    match base.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => base,
    }
}

/// Get the complement of a DNA sequence.
pub fn complement_seq(seq: &str) -> String {
    seq.chars().map(complement).collect()
}

/// Check if a base pair is ambiguous (A/T or C/G).
/// These cannot be strand-flipped reliably.
pub fn is_ambiguous_snp(allele1: &str, allele2: &str) -> bool {
    if allele1.len() != 1 || allele2.len() != 1 {
        return false;
    }
    let a = allele1.chars().next().unwrap().to_ascii_uppercase();
    let b = allele2.chars().next().unwrap().to_ascii_uppercase();
    matches!((a, b), ('A', 'T') | ('T', 'A') | ('C', 'G') | ('G', 'C'))
}

/// Result of harmonizing a user genotype against a panel.
#[derive(Debug, Clone)]
pub enum HarmonizationResult {
    /// Successfully mapped to panel alleles.
    /// Contains (allele1_index, allele2_index) in panel encoding.
    Success { gt_indices: (usize, usize) },

    /// User alleles don't match reference genome and can't be strand-flipped.
    InvalidAlleles {
        user_bases: Vec<String>,
        ref_base: String,
    },

    /// Site not processable (e.g., indel mismatch).
    Skip { reason: String },
}

/// Harmonize user bases against a reference panel site.
///
/// # Arguments
/// * `user_bases` - The actual bases the user has (e.g., ["A", "G"])
/// * `ref_base` - The reference genome base at this position
/// * `chrom` - Chromosome name
/// * `pos` - Position
/// * `panel` - Mutable reference to the padded panel for tracking modifications
///
/// # Returns
/// * `HarmonizationResult` indicating success or failure with allele indices
pub fn harmonize_genotype(
    user_bases: &[String],
    ref_base: &str,
    chrom: &str,
    pos: u64,
    panel: &mut PaddedPanel,
) -> HarmonizationResult {
    // Validate user bases against reference
    let mut validated_bases = Vec::with_capacity(user_bases.len());
    let mut needs_flip = false;

    for base in user_bases {
        if base.eq_ignore_ascii_case(ref_base) {
            // Matches reference - good
            validated_bases.push(base.to_uppercase());
        } else if base.len() == 1 && ref_base.len() == 1 {
            // Check if complement matches reference
            let comp = complement_seq(base);
            if comp.eq_ignore_ascii_case(ref_base) {
                needs_flip = true;
                validated_bases.push(base.to_uppercase());
            } else {
                // User has a variant allele (or complement of variant)
                validated_bases.push(base.to_uppercase());
            }
        } else {
            // Indel - keep as-is
            validated_bases.push(base.to_uppercase());
        }
    }

    // Apply strand flip if detected (complement all bases)
    if needs_flip {
        validated_bases = validated_bases.iter().map(|b| complement_seq(b)).collect();
    }

    // Now map each base to panel allele index
    let mut gt_indices = Vec::with_capacity(validated_bases.len());

    for base in &validated_bases {
        let idx = panel.get_or_add_allele_index(chrom, pos, base, ref_base);
        gt_indices.push(idx);
    }

    if gt_indices.len() >= 2 {
        HarmonizationResult::Success {
            gt_indices: (gt_indices[0], gt_indices[1]),
        }
    } else if gt_indices.len() == 1 {
        // Haploid
        HarmonizationResult::Success {
            gt_indices: (gt_indices[0], gt_indices[0]),
        }
    } else {
        HarmonizationResult::Skip {
            reason: "No bases to process".to_string(),
        }
    }
}

/// Get a panel site's alleles merged with any added ALTs.
pub fn get_merged_alts(site: &PanelSite, added: Option<&Vec<String>>) -> Vec<String> {
    let mut alts = site.alt_alleles.clone();
    if let Some(added_alts) = added {
        alts.extend(added_alts.iter().cloned());
    }
    alts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(complement('A'), 'T');
        assert_eq!(complement('T'), 'A');
        assert_eq!(complement('C'), 'G');
        assert_eq!(complement('G'), 'C');
        assert_eq!(complement('a'), 'T');
    }

    #[test]
    fn test_complement_seq() {
        assert_eq!(complement_seq("ATCG"), "TAGC");
        assert_eq!(complement_seq("atcg"), "TAGC");
    }

    #[test]
    fn test_is_ambiguous_snp() {
        assert!(is_ambiguous_snp("A", "T"));
        assert!(is_ambiguous_snp("T", "A"));
        assert!(is_ambiguous_snp("C", "G"));
        assert!(is_ambiguous_snp("G", "C"));
        assert!(!is_ambiguous_snp("A", "G"));
        assert!(!is_ambiguous_snp("A", "C"));
        assert!(!is_ambiguous_snp("AT", "TA")); // Not SNPs
    }
}
