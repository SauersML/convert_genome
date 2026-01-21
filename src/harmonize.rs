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
    /// Diploid genotype mapped to panel alleles.
    Diploid(usize, usize),

    /// Haploid genotype mapped to panel alleles.
    Haploid(usize),

    /// User alleles don't match reference genome and can't be strand-flipped.
    InvalidAlleles {
        user_bases: Vec<String>,
        ref_base: String,
    },

    /// Site not processable (e.g., indel mismatch).
    Skip { reason: String },
}

/// Harmonize a list of alleles against the panel.
/// Returns the panel allele indices for each input allele.
/// Handles strand flipping if necessary.
pub fn harmonize_alleles(
    alleles: &[String],
    ref_base: &str,
    chrom: &str,
    pos: u64,
    panel: &mut PaddedPanel,
) -> Result<Vec<usize>, String> {
    // Validate user bases against reference
    let mut validated_bases = Vec::with_capacity(alleles.len());
    let mut has_ref_match = false;
    let mut has_ref_comp_match = false;

    // Check strict matching first
    for base in alleles {
        if base.eq_ignore_ascii_case(ref_base) {
            has_ref_match = true;
        } else if base.len() == 1 && ref_base.len() == 1 {
            let comp = complement_seq(base);
            if comp.eq_ignore_ascii_case(ref_base) {
                has_ref_comp_match = true;
            }
        }
    }

    // Apply flip decision
    // Priority: No Flip (Direct Match) > Flip (Comp Match)
    let needs_flip = !has_ref_match && has_ref_comp_match;

    let original_site = panel.get_original(chrom, pos);
    let panel_ref = original_site.map(|s| s.ref_allele.as_str());

    for (i, base) in alleles.iter().enumerate() {
        if needs_flip {
            validated_bases.push(complement_seq(base));
            continue;
        }

        let base_upper = base.to_uppercase();
        if i == 0 {
            validated_bases.push(base_upper);
            continue;
        }

        let mut rescued = None;
        if let Some(site) = original_site {
            let direct_matches_alt = site
                .alt_alleles
                .iter()
                .any(|a| a.eq_ignore_ascii_case(&base_upper));
            if direct_matches_alt {
                rescued = Some(base_upper.clone());
            } else {
                let comp = complement_seq(&base_upper);
                let comp_matches_alt = site
                    .alt_alleles
                    .iter()
                    .any(|a| a.eq_ignore_ascii_case(&comp));
                let comp_matches_panel_ref = panel_ref
                    .map(|r| r.eq_ignore_ascii_case(&comp))
                    .unwrap_or(false);
                if comp_matches_alt && !comp_matches_panel_ref {
                    rescued = Some(comp);
                }
            }
        }

        validated_bases.push(rescued.unwrap_or(base_upper));
    }

    // Map to indices
    let mut indices = Vec::with_capacity(validated_bases.len());
    for base in validated_bases {
        let idx = panel.get_or_add_allele_index(chrom, pos, &base, ref_base);
        indices.push(idx);
    }

    Ok(indices)
}

/// Harmonize user GENOTYPE (called alleles) against a reference panel site.
pub fn harmonize_genotype(
    user_bases: &[String],
    ref_base: &str,
    chrom: &str,
    pos: u64,
    panel: &mut PaddedPanel,
) -> HarmonizationResult {
    match harmonize_alleles(user_bases, ref_base, chrom, pos, panel) {
        Ok(indices) => match indices.len() {
            1 => HarmonizationResult::Haploid(indices[0]),
            2 => HarmonizationResult::Diploid(indices[0], indices[1]),
            _ => HarmonizationResult::InvalidAlleles {
                user_bases: user_bases.to_vec(),
                ref_base: ref_base.to_string(),
            },
        },
        Err(e) => HarmonizationResult::Skip { reason: e },
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

    #[test]
    fn test_ambiguous_flip_prevention() {
        let chrom = "1";
        let pos = 100;
        let original_index = crate::panel::PanelIndex::default();
        let mut panel = PaddedPanel::new(original_index);

        // Ref=A. Input=["A", "T"].
        // "A" matches Ref. "T" matches Comp(Ref).
        // Ambiguous. Priority: Fwd match ("A") -> No Flip.
        let res = harmonize_genotype(
            &["A".to_string(), "T".to_string()],
            "A",
            chrom,
            pos,
            &mut panel,
        );
        match res {
            HarmonizationResult::Diploid(idx1, idx2) => {
                // Input [A, T] with A as reference should yield indices (0, 1)
                assert_eq!(idx1, 0, "First allele should be Ref (A)");
                assert_eq!(idx2, 1, "Second allele should be Alt (T)");
            }
            _ => panic!("Expected Diploid result"),
        }
    }

    #[test]
    fn test_flip_when_no_ref_match() {
        let chrom = "1";
        let pos = 100;
        let original_index = crate::panel::PanelIndex::default();
        let mut panel = PaddedPanel::new(original_index);

        // Ref=A. Input=["T", "C"].
        // "T" matches Comp(Ref)=A. "C" matches Comp(G).
        // No direct Ref match ("A").
        // Has Ref Comp match ("T" -> "A").
        // Decision: FLIP.
        // Flipped: T->A, C->G.
        // Result: [A, G].
        // A is Ref (0). G is novel (1).
        // Expect (0, 1).
        let res = harmonize_genotype(
            &["T".to_string(), "C".to_string()],
            "A",
            chrom,
            pos,
            &mut panel,
        );
        match res {
            HarmonizationResult::Diploid(idx1, idx2) => {
                assert_eq!(idx1, 0, "First allele should be flipped to Ref (A)");
                assert_eq!(idx2, 1, "Second allele should be flipped to Alt (G)");
            }
            _ => panic!("Expected Diploid result"),
        }
    }
}
