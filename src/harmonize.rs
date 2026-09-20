//! Site-level harmonization of a sample's SNV genotype against a reference panel.
//!
//! Strand is a property of a SITE, never of one allele. The alleles a sample
//! carries at an SNV are either all on the panel's strand, or all on the other
//! one; the only honest question is which, and for an A/T or C/G site the
//! alleles alone cannot answer it.
//!
//! The earlier per-allele "rescue" complemented any single ALT that missed the
//! panel but whose complement hit it. At a split multi-allelic panel site
//! (A>C and A>G written as two records) whose A>G record was the one kept, a
//! true A/C call became A/G. It also turned genuine novel alleles into panel
//! alleles at biallelic sites. It is gone: [`resolve_snv`] decides per site,
//! with the counts a Michigan-style QC report uses (match, strand flip,
//! mismatch, palindromic), and resolves palindromic homozygotes from the
//! panel allele frequency only when the input's strand is uncertain.

/// Smallest strand-flip prior at which an unambiguous flip is resolved by
/// complementing. Below it (a reference-oriented input: a sequencing VCF, a
/// Plus-strand array export, a 23andMe-style table) a genotype that only fits
/// the panel after complementing is more likely a genotyping error or a novel
/// allele than a strand error, and is left missing for imputation to fill.
pub const FLIP_RESOLVE_MIN_PRIOR: f64 = 0.01;

/// Minor-allele frequency at and above which an A/T or C/G homozygote on a
/// strand-uncertain input is dropped outright. Near 0.5 the two strands give
/// almost the same genotype likelihood, and a global panel frequency is not
/// precise enough across ancestries to tell them apart.
pub const PALINDROME_AMBIGUITY_MAF: f64 = 0.40;

/// A palindromic homozygote is kept (or complemented) only when the posterior
/// probability that it sits on the other strand is at most this (or at least
/// one minus this). Anything in between is dropped.
pub const PALINDROME_POSTERIOR_TOLERANCE: f64 = 0.01;

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

/// What to do with the sample's called alleles at one SNV site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Resolution {
    /// The alleles are on the panel's strand; use them as called.
    Keep,
    /// The alleles are on the other strand; use their complements.
    Complement,
    /// The genotype cannot be placed on the panel's alleles; emit it as missing.
    Missing,
}

/// Why a site was resolved the way it was. One counter per class goes into the
/// run report's `panel_qc` block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteClass {
    /// Called alleles are panel alleles (non-palindromic, one panel ALT).
    Match,
    /// Only the complements are panel alleles, and the input's strand is
    /// uncertain enough that a strand flip is the likely cause: complemented.
    FlipResolved,
    /// Only the complements are panel alleles, on a reference-oriented input:
    /// left missing.
    FlipUnresolved,
    /// Neither the alleles nor their complements are panel alleles.
    AlleleMismatch,
    /// A/T or C/G heterozygote: strand-invariant, kept.
    PalindromeHet,
    /// A/T or C/G homozygote kept as called.
    PalindromeKept,
    /// A/T or C/G homozygote complemented (swapped between the homozygotes).
    PalindromeFlipped,
    /// A/T or C/G homozygote in the ambiguity band: dropped.
    PalindromeAmbiguous,
    /// A/T or C/G homozygote whose strand posterior is neither near 0 nor
    /// near 1, or whose panel frequency is unknown: dropped.
    PalindromeUnresolved,
    /// Several panel ALTs at the position and the call fits them directly.
    MultiAllelicMatch,
    /// Several panel ALTs at the position and the call does not fit them
    /// directly. Never strand-rescued.
    MultiAllelicMismatch,
}

/// Decide one SNV site.
///
/// * `called` — the sample's called alleles in genotype order (a homozygote
///   repeats its allele); missing alleles are excluded by the caller.
/// * `panel_ref` / `panel_alts` — the panel's REF and the union of its SNV
///   ALTs at this position (split records merged).
/// * `alt_frequency` — the panel frequency of the single ALT of a biallelic
///   site, when the panel carries AC/AN.
/// * `strand_prior` — the probability that any one site of this input is
///   reported on the opposite strand. Zero for inputs whose strand is fixed by
///   construction (a sequencing VCF). Measured by the caller from the input's
///   non-palindromic heterozygotes otherwise.
pub fn resolve_snv(
    called: &[String],
    panel_ref: &str,
    panel_alts: &[String],
    alt_frequency: Option<f64>,
    strand_prior: f64,
) -> (Resolution, SiteClass) {
    let upper = |s: &str| s.to_ascii_uppercase();
    let panel_ref = upper(panel_ref);
    let panel_alts: Vec<String> = panel_alts.iter().map(|a| upper(a)).collect();
    let in_panel = |allele: &str| allele == panel_ref || panel_alts.iter().any(|a| a == allele);

    let called: Vec<String> = called.iter().map(|a| upper(a)).collect();
    let direct = called.iter().all(|a| in_panel(a));

    if panel_alts.len() > 1 {
        return if direct {
            (Resolution::Keep, SiteClass::MultiAllelicMatch)
        } else {
            (Resolution::Missing, SiteClass::MultiAllelicMismatch)
        };
    }

    let palindromic = panel_alts.len() == 1 && complement_seq(&panel_ref) == panel_alts[0];
    if !palindromic {
        if direct {
            return (Resolution::Keep, SiteClass::Match);
        }
        let flipped = called.iter().all(|a| in_panel(&complement_seq(a)));
        if !flipped {
            return (Resolution::Missing, SiteClass::AlleleMismatch);
        }
        return if strand_prior >= FLIP_RESOLVE_MIN_PRIOR {
            (Resolution::Complement, SiteClass::FlipResolved)
        } else {
            (Resolution::Missing, SiteClass::FlipUnresolved)
        };
    }

    if !direct {
        return (Resolution::Missing, SiteClass::AlleleMismatch);
    }
    let heterozygous = called.windows(2).any(|w| w[0] != w[1]);
    if heterozygous {
        return (Resolution::Keep, SiteClass::PalindromeHet);
    }
    if strand_prior <= 0.0 {
        return (Resolution::Keep, SiteClass::PalindromeKept);
    }
    let Some(p) = alt_frequency.filter(|p| p.is_finite() && (0.0..=1.0).contains(p)) else {
        return (Resolution::Missing, SiteClass::PalindromeUnresolved);
    };
    if p.min(1.0 - p) >= PALINDROME_AMBIGUITY_MAF {
        return (Resolution::Missing, SiteClass::PalindromeAmbiguous);
    }
    let observed_is_alt = called.first().is_some_and(|a| *a == panel_alts[0]);
    let (hom_alt, hom_ref) = (p * p, (1.0 - p) * (1.0 - p));
    let (as_called, as_flipped) = if observed_is_alt {
        (hom_alt, hom_ref)
    } else {
        (hom_ref, hom_alt)
    };
    let prior = strand_prior.min(1.0);
    let numerator = prior * as_flipped;
    let denominator = numerator + (1.0 - prior) * as_called;
    if denominator <= 0.0 {
        return (Resolution::Missing, SiteClass::PalindromeUnresolved);
    }
    let posterior_flipped = numerator / denominator;
    if posterior_flipped <= PALINDROME_POSTERIOR_TOLERANCE {
        (Resolution::Keep, SiteClass::PalindromeKept)
    } else if posterior_flipped >= 1.0 - PALINDROME_POSTERIOR_TOLERANCE {
        (Resolution::Complement, SiteClass::PalindromeFlipped)
    } else {
        (Resolution::Missing, SiteClass::PalindromeUnresolved)
    }
}

/// Get a panel site's alleles merged with any added ALTs.
pub fn get_merged_alts(site: &crate::panel::PanelSite, added: Option<&Vec<String>>) -> Vec<String> {
    let mut alts = site.alt_alleles.clone();
    if let Some(added_alts) = added {
        alts.extend(added_alts.iter().cloned());
    }
    alts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(v: &[&str]) -> Vec<String> {
        v.iter().map(|a| a.to_string()).collect()
    }

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
    fn split_site_a_c_call_is_never_rescued_onto_g() {
        // Panel A>C and A>G (split), sample A/C: fits directly, kept as A/C.
        assert_eq!(
            resolve_snv(&s(&["A", "C"]), "A", &s(&["C", "G"]), None, 0.08),
            (Resolution::Keep, SiteClass::MultiAllelicMatch)
        );
        // Panel A>G only (the other split record lost), sample A/C: the old
        // per-allele rescue complemented C into G. Now it is a mismatch.
        assert_eq!(
            resolve_snv(&s(&["A", "C"]), "A", &s(&["G"]), None, 0.08),
            (Resolution::Missing, SiteClass::AlleleMismatch)
        );
        // Several panel ALTs and a call that fits only after complementing:
        // never strand-rescued.
        assert_eq!(
            resolve_snv(&s(&["T", "T"]), "A", &s(&["C", "G"]), None, 0.5),
            (Resolution::Missing, SiteClass::MultiAllelicMismatch)
        );
    }

    #[test]
    fn unambiguous_flip_is_resolved_only_on_a_strand_uncertain_input() {
        // Panel A>G; a reverse-strand het reads T/C and a reverse hom-alt CC.
        assert_eq!(
            resolve_snv(&s(&["T", "C"]), "A", &s(&["G"]), None, 0.08),
            (Resolution::Complement, SiteClass::FlipResolved)
        );
        assert_eq!(
            resolve_snv(&s(&["C", "C"]), "A", &s(&["G"]), None, 0.08),
            (Resolution::Complement, SiteClass::FlipResolved)
        );
        assert_eq!(
            resolve_snv(&s(&["C", "C"]), "A", &s(&["G"]), None, 0.0),
            (Resolution::Missing, SiteClass::FlipUnresolved)
        );
        assert_eq!(
            resolve_snv(&s(&["A", "G"]), "A", &s(&["G"]), None, 0.08),
            (Resolution::Keep, SiteClass::Match)
        );
    }

    #[test]
    fn palindromes_follow_the_prior_and_the_ambiguity_band() {
        // Het A/T is strand-invariant.
        assert_eq!(
            resolve_snv(&s(&["A", "T"]), "A", &s(&["T"]), Some(0.45), 0.5),
            (Resolution::Keep, SiteClass::PalindromeHet)
        );
        // Strand fixed by construction: homozygotes kept whatever the frequency.
        assert_eq!(
            resolve_snv(&s(&["T", "T"]), "A", &s(&["T"]), Some(0.01), 0.0),
            (Resolution::Keep, SiteClass::PalindromeKept)
        );
        // Strand uncertain, MAF in the band: dropped.
        assert_eq!(
            resolve_snv(&s(&["A", "A"]), "A", &s(&["T"]), Some(0.45), 0.08),
            (Resolution::Missing, SiteClass::PalindromeAmbiguous)
        );
        // Common homozygote on a Forward-like input (8% reversed): kept.
        assert_eq!(
            resolve_snv(&s(&["A", "A"]), "A", &s(&["T"]), Some(0.05), 0.08),
            (Resolution::Keep, SiteClass::PalindromeKept)
        );
        // Rare homozygote there: as likely a flipped common one; dropped.
        assert_eq!(
            resolve_snv(&s(&["T", "T"]), "A", &s(&["T"]), Some(0.05), 0.08),
            (Resolution::Missing, SiteClass::PalindromeUnresolved)
        );
        // Mostly reversed input (TOP/BOT-like): the rare homozygote is a flip.
        assert_eq!(
            resolve_snv(&s(&["T", "T"]), "A", &s(&["T"]), Some(0.01), 0.5),
            (Resolution::Complement, SiteClass::PalindromeFlipped)
        );
        // No panel frequency on a strand-uncertain input: dropped.
        assert_eq!(
            resolve_snv(&s(&["A", "A"]), "A", &s(&["T"]), None, 0.08),
            (Resolution::Missing, SiteClass::PalindromeUnresolved)
        );
        // An allele outside the pair is a mismatch, not a palindrome question.
        assert_eq!(
            resolve_snv(&s(&["A", "C"]), "A", &s(&["T"]), Some(0.2), 0.0),
            (Resolution::Missing, SiteClass::AlleleMismatch)
        );
    }
}
