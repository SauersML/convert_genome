//! Reference panel loading and harmonization support.
//!
//! This module provides functionality to load a reference panel (VCF/BCF),
//! index its sites for fast lookup, and track modifications (padding) needed
//! to represent novel alleles from user data.

use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;

use anyhow::{Context, Result};
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases as _;

/// A site definition from the reference panel.
#[derive(Debug, Clone)]
pub struct PanelSite {
    pub chrom: String,
    pub pos: u64,
    pub id: Option<String>,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    /// Panel ALT allele counts (INFO/AC, one per ALT), when the panel carries them.
    pub alt_counts: Vec<Option<u32>>,
    /// Panel allele number (INFO/AN), when the panel carries it.
    pub allele_number: Option<u32>,
}

impl PanelSite {
    /// A site with no allele-frequency information.
    pub fn new(chrom: &str, pos: u64, ref_allele: &str, alt_alleles: &[&str]) -> Self {
        Self {
            chrom: chrom.to_string(),
            pos,
            id: None,
            ref_allele: ref_allele.to_string(),
            alt_alleles: alt_alleles.iter().map(|a| a.to_string()).collect(),
            alt_counts: Vec::new(),
            allele_number: None,
        }
    }

    /// Frequency of ALT `alt_index` (0-based into `alt_alleles`) in the panel,
    /// from INFO/AC and INFO/AN. `None` when the panel does not carry them.
    pub fn alt_frequency(&self, alt_index: usize) -> Option<f64> {
        let an = self.allele_number.filter(|&an| an > 0)?;
        let ac = self.alt_counts.get(alt_index).copied().flatten()?;
        Some(f64::from(ac) / f64::from(an))
    }

    /// True when every allele at the site is a single base (an SNV record).
    pub fn is_snv(&self) -> bool {
        self.ref_allele.len() == 1 && self.alt_alleles.iter().all(|a| a.len() == 1)
    }

    /// Find the index of an allele (0 = REF, 1+ = ALT).
    /// Returns None if allele is not present.
    pub fn allele_index(&self, allele: &str) -> Option<usize> {
        if allele.eq_ignore_ascii_case(&self.ref_allele) {
            return Some(0);
        }
        for (i, alt) in self.alt_alleles.iter().enumerate() {
            if allele.eq_ignore_ascii_case(alt) {
                return Some(i + 1);
            }
        }
        None
    }

    /// Check if this site contains the given allele.
    pub fn has_allele(&self, allele: &str) -> bool {
        self.allele_index(allele).is_some()
    }
}

/// One panel record, stored compactly: a genome-wide panel holds tens of
/// millions of these, and a pipeline process imputes several chromosomes at once.
#[derive(Debug, Clone)]
struct CompactRecord {
    pos: u64,
    /// INFO/AN, or `u32::MAX` when absent.
    allele_number: u32,
    /// INFO/AC of the first ALT, or `u32::MAX` when absent.
    first_count: u32,
    /// `REF` and the ALTs, tab then comma separated (`A\tG`, `AT\tA,ATT`).
    alleles: Box<str>,
    /// INFO/AC of the ALTs after the first (`u32::MAX` when absent).
    other_counts: Option<Box<[u32]>>,
}

const MISSING_COUNT: u32 = u32::MAX;

impl CompactRecord {
    fn from_site(site: &PanelSite) -> Self {
        let mut alleles = site.ref_allele.clone();
        alleles.push('\t');
        alleles.push_str(&site.alt_alleles.join(","));
        let count = |i: usize| {
            site.alt_counts
                .get(i)
                .copied()
                .flatten()
                .unwrap_or(MISSING_COUNT)
        };
        let other_counts =
            (site.alt_alleles.len() > 1).then(|| (1..site.alt_alleles.len()).map(count).collect());
        Self {
            pos: site.pos,
            allele_number: site.allele_number.unwrap_or(MISSING_COUNT),
            first_count: count(0),
            alleles: alleles.into_boxed_str(),
            other_counts,
        }
    }

    fn to_site(&self, chrom: &str) -> PanelSite {
        let (reference, alts) = self.alleles.split_once('\t').unwrap_or((&self.alleles, ""));
        let alt_alleles: Vec<String> = if alts.is_empty() {
            Vec::new()
        } else {
            alts.split(',').map(str::to_string).collect()
        };
        let known = |c: u32| (c != MISSING_COUNT).then_some(c);
        let mut alt_counts = Vec::with_capacity(alt_alleles.len());
        if !alt_alleles.is_empty() {
            alt_counts.push(known(self.first_count));
        }
        if let Some(rest) = &self.other_counts {
            alt_counts.extend(rest.iter().map(|&c| known(c)));
        }
        PanelSite {
            chrom: chrom.to_string(),
            pos: self.pos,
            id: None,
            ref_allele: reference.to_string(),
            alt_alleles,
            alt_counts,
            allele_number: known(self.allele_number),
        }
    }
}

/// Index of panel records for lookup by position.
///
/// A position can hold several records: panels built with `bcftools norm -m-`
/// split a multi-allelic site into one biallelic record per ALT (A>C and A>G at
/// the same POS), and an SNV and an indel can share a POS. Every record is kept,
/// in file order. Keying on position alone and letting the last record win lost
/// the other ALTs, and harmonization then "rescued" a true A/C call onto the kept
/// A>G record.
///
/// Records are held per chromosome, sorted by position (stably, so records at
/// one position keep file order), and found by binary search.
#[derive(Default)]
pub struct PanelIndex {
    /// Chromosome name (panel naming) -> its records, sorted by position.
    records: HashMap<String, Vec<CompactRecord>>,
    /// Chromosome order from the panel header for sorting output
    chrom_order: Vec<String>,
    record_count: usize,
}

impl PanelIndex {
    /// Load a panel index from a VCF file.
    /// Only loads site definitions (chrom, pos, ref, alt), not genotypes.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");

        if ext == "bcf" {
            Self::load_bcf(path)
        } else {
            Self::load_vcf(path)
        }
    }

    fn load_vcf<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let reader = crate::smart_reader::open_input(path)
            .with_context(|| format!("failed to open panel {}", path.display()))?;

        let mut vcf_reader = vcf::io::Reader::new(reader);
        let header = vcf_reader.read_header()?;

        // Extract chromosome order from header contigs
        let chrom_order: Vec<String> = header.contigs().keys().map(|k| k.to_string()).collect();

        let mut index = Self {
            chrom_order,
            ..Self::default()
        };
        let mut record_count = 0u64;

        let mut record = vcf::Record::default();
        loop {
            if vcf_reader.read_record(&mut record)? == 0 {
                break;
            }
            record_count += 1;

            let chrom = record.reference_sequence_name().to_string();
            let pos = record
                .variant_start()
                .transpose()?
                .map(|p| usize::from(p) as u64)
                .unwrap_or(0);

            let ref_allele = record.reference_bases().to_string();
            let alt_alleles: Vec<String> = record
                .alternate_bases()
                .iter()
                .map(|a| a.map(|s| s.to_string()))
                .collect::<std::io::Result<Vec<_>>>()?;

            let (alt_counts, allele_number) =
                allele_counts(&record.info(), &header, alt_alleles.len());

            index.insert(PanelSite {
                chrom,
                pos,
                id: None,
                ref_allele,
                alt_alleles,
                alt_counts,
                allele_number,
            });

            if record_count.is_multiple_of(1_000_000) {
                tracing::info!("Loaded {} panel sites...", record_count);
            }
        }

        index.finish();
        tracing::info!("Loaded {} panel sites total", index.len());

        Ok(index)
    }

    fn load_bcf<P: AsRef<Path>>(path: P) -> Result<Self> {
        use noodles::bcf;

        let file = std::fs::File::open(path.as_ref())
            .with_context(|| format!("failed to open panel {}", path.as_ref().display()))?;

        let mut bcf_reader = bcf::io::Reader::new(BufReader::new(file));
        let header = bcf_reader.read_header()?;

        let chrom_order: Vec<String> = header.contigs().keys().map(|k| k.to_string()).collect();

        let mut index = Self {
            chrom_order,
            ..Self::default()
        };
        let mut record_count = 0u64;

        let string_maps = vcf::header::StringMaps::try_from(&header)?;
        let mut record = bcf::Record::default();
        loop {
            if bcf_reader.read_record(&mut record)? == 0 {
                break;
            }
            record_count += 1;

            let chrom = record.reference_sequence_name(&string_maps)?.to_string();
            let pos = record
                .variant_start()
                .transpose()?
                .map(|p| usize::from(p) as u64)
                .unwrap_or(0);

            let ref_allele = std::str::from_utf8(record.reference_bases().as_ref())?.to_string();
            let alt_alleles: Vec<String> = record
                .alternate_bases()
                .iter()
                .map(|a| a.map(|s| s.to_string()))
                .collect::<std::io::Result<Vec<_>>>()?;

            let (alt_counts, allele_number) =
                allele_counts(&record.info(), &header, alt_alleles.len());

            index.insert(PanelSite {
                chrom,
                pos,
                id: None,
                ref_allele,
                alt_alleles,
                alt_counts,
                allele_number,
            });

            if record_count.is_multiple_of(1_000_000) {
                tracing::info!("Loaded {} panel sites...", record_count);
            }
        }

        index.finish();
        tracing::info!("Loaded {} panel sites total", index.len());

        Ok(index)
    }

    /// Build an index from records (tests and library callers).
    pub fn from_sites(sites: impl IntoIterator<Item = PanelSite>) -> Self {
        let mut index = Self::default();
        for site in sites {
            if !index.chrom_order.contains(&site.chrom) {
                index.chrom_order.push(site.chrom.clone());
            }
            index.insert(site);
        }
        index.finish();
        index
    }

    fn insert(&mut self, site: PanelSite) {
        self.record_count += 1;
        let compact = CompactRecord::from_site(&site);
        match self.records.get_mut(site.chrom.as_str()) {
            Some(records) => records.push(compact),
            None => {
                self.records.insert(site.chrom, vec![compact]);
            }
        }
    }

    fn finish(&mut self) {
        for records in self.records.values_mut() {
            // Stable: records at one position keep their file order.
            records.sort_by_key(|r| r.pos);
            records.shrink_to_fit();
        }
    }

    /// The panel's own name for `chrom`, trying it with and without a `chr` prefix.
    pub fn resolve_chrom(&self, chrom: &str) -> Option<&str> {
        if let Some((name, _)) = self.records.get_key_value(chrom) {
            return Some(name.as_str());
        }
        let alt_chrom = alternate_chrom_name(chrom);
        self.records
            .get_key_value(alt_chrom.as_str())
            .map(|(name, _)| name.as_str())
    }

    fn chrom_records(&self, chrom: &str) -> Option<(&str, &[CompactRecord])> {
        let name = self.resolve_chrom(chrom)?;
        self.records.get(name).map(|r| (name, r.as_slice()))
    }

    /// Every panel record at a position, in panel file order (empty if none).
    pub fn get_all(&self, chrom: &str, pos: u64) -> Vec<PanelSite> {
        let Some((name, records)) = self.chrom_records(chrom) else {
            return Vec::new();
        };
        let start = records.partition_point(|r| r.pos < pos);
        records[start..]
            .iter()
            .take_while(|r| r.pos == pos)
            .map(|r| r.to_site(name))
            .collect()
    }

    /// The first panel record at a position, if any.
    pub fn get(&self, chrom: &str, pos: u64) -> Option<PanelSite> {
        self.get_all(chrom, pos).into_iter().next()
    }

    /// Check if a site exists at the given position.
    pub fn contains(&self, chrom: &str, pos: u64) -> bool {
        self.chrom_records(chrom).is_some_and(|(_, records)| {
            let i = records.partition_point(|r| r.pos < pos);
            records.get(i).is_some_and(|r| r.pos == pos)
        })
    }

    /// Every record on a chromosome in position order (file order within a position).
    pub fn sites_on(&self, chrom: &str) -> impl Iterator<Item = PanelSite> + '_ {
        let (name, records) = self.chrom_records(chrom).unwrap_or(("", &[]));
        records.iter().map(move |r| r.to_site(name))
    }

    /// Number of panel records in the index (split records count separately).
    pub fn len(&self) -> usize {
        self.record_count
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool {
        self.record_count == 0
    }

    /// Get chromosome order from the panel header.
    pub fn chrom_order(&self) -> &[String] {
        &self.chrom_order
    }
}

/// Per-site outcome counts of harmonizing an input against the panel: the
/// Michigan-style match / strand-flip / mismatch / palindrome tally, plus what
/// the absent-genotype fill wrote. Serialized as the report's `panel_qc` block.
#[derive(Debug, Default, Clone, serde::Serialize)]
pub struct PanelQc {
    /// Strand-flip prior the palindrome policy used (from the input header).
    pub strand_flip_prior: f64,
    /// True when absent panel sites were written as homozygous reference.
    pub absent_genotypes_hom_ref: bool,
    pub snv_match: usize,
    pub snv_flip_resolved: usize,
    pub snv_flip_unresolved: usize,
    pub snv_allele_mismatch: usize,
    pub snv_reference_mismatch: usize,
    pub palindrome_het: usize,
    pub palindrome_kept: usize,
    pub palindrome_flipped: usize,
    pub palindrome_ambiguous: usize,
    pub palindrome_unresolved: usize,
    pub multiallelic_match: usize,
    pub multiallelic_mismatch: usize,
    /// A heterozygote whose two ALTs live on different split panel records.
    pub split_record_het: usize,
    pub no_call: usize,
    pub indel_exact: usize,
    pub indel_unmatched: usize,
    pub absent_from_panel: usize,
    pub fill_hom_ref: usize,
    pub fill_missing_no_call: usize,
    pub fill_missing_overlap: usize,
}

impl PanelQc {
    pub fn count(&mut self, class: crate::harmonize::SiteClass) {
        use crate::harmonize::SiteClass as C;
        let slot = match class {
            C::Match => &mut self.snv_match,
            C::FlipResolved => &mut self.snv_flip_resolved,
            C::FlipUnresolved => &mut self.snv_flip_unresolved,
            C::AlleleMismatch => &mut self.snv_allele_mismatch,
            C::PalindromeHet => &mut self.palindrome_het,
            C::PalindromeKept => &mut self.palindrome_kept,
            C::PalindromeFlipped => &mut self.palindrome_flipped,
            C::PalindromeAmbiguous => &mut self.palindrome_ambiguous,
            C::PalindromeUnresolved => &mut self.palindrome_unresolved,
            C::MultiAllelicMatch => &mut self.multiallelic_match,
            C::MultiAllelicMismatch => &mut self.multiallelic_mismatch,
        };
        *slot += 1;
    }
}

/// What the input said about one chromosome, for writing the panel sites it did
/// not mention. Positions are 1-based; spans are inclusive.
#[derive(Debug, Default)]
pub struct ChromClaims {
    /// Contig name the output uses for this chromosome (the input's, after
    /// standardization), so filled records sort and index with the rest.
    pub output_name: String,
    /// Positions with an SNV call: the called alleles after harmonization, or
    /// `None` when the call is missing or could not be placed on the panel.
    pub snv_calls: HashMap<u64, Option<Vec<String>>>,
    /// Panel records the input already wrote: (pos, REF, ALTs joined by ',').
    pub emitted: std::collections::HashSet<(u64, String, String)>,
    /// Bases altered by a called non-SNV variant (deleted or substituted).
    pub altered_spans: Vec<(u64, u64)>,
    /// No-call blocks (gVCF-style records with END and a missing genotype).
    pub no_call_spans: Vec<(u64, u64)>,
}

/// Tracks modifications to the panel (novel alleles and sites).
pub struct PaddedPanel {
    /// Original panel index (immutable reference)
    original: PanelIndex,
    /// Additional ALTs added to existing sites: (chrom, pos) -> new ALTs to append
    added_alts: HashMap<(String, u64), Vec<String>>,
    /// Entirely novel sites not in original panel
    novel_sites: Vec<PanelSite>,
    /// Harmonization outcome counts.
    pub qc: PanelQc,
    /// Per-chromosome claims, keyed by the panel's chromosome name.
    pub claims: HashMap<String, ChromClaims>,
}

impl PaddedPanel {
    /// Create a new padded panel from an original panel index.
    pub fn new(original: PanelIndex) -> Self {
        Self {
            original,
            added_alts: HashMap::new(),
            novel_sites: Vec::new(),
            qc: PanelQc::default(),
            claims: HashMap::new(),
        }
    }

    /// The claims for `chrom` (input naming), created on first use. `None` when
    /// the panel does not carry the chromosome.
    pub fn claims_for(&mut self, chrom: &str, output_name: &str) -> Option<&mut ChromClaims> {
        let panel_chrom = self.original.resolve_chrom(chrom)?.to_string();
        let claims = self.claims.entry(panel_chrom).or_default();
        if claims.output_name.is_empty() {
            claims.output_name = output_name.to_string();
        }
        Some(claims)
    }

    /// Get the original panel site, if any.
    pub fn get_original(&self, chrom: &str, pos: u64) -> Option<PanelSite> {
        self.original.get(chrom, pos)
    }

    /// Get allele index for a base at a position.
    /// If the allele is novel, adds it to the padded panel and returns the new index.
    pub fn get_or_add_allele_index(
        &mut self,
        chrom: &str,
        pos: u64,
        allele: &str,
        ref_base: &str,
    ) -> usize {
        if let Some(site) = self.original.get(chrom, pos) {
            // Site exists in panel - USE CANONICAL CHROM NAME FROM SITE
            // This ensures added_alts is keyed by "chr1" even if user passed "1"
            let canonical_key = (site.chrom.clone(), pos);

            if let Some(idx) = site.allele_index(allele) {
                return idx;
            }

            // Allele not in panel - add it
            let added = self.added_alts.entry(canonical_key).or_default();

            // Check if we already added this allele
            for (i, added_alt) in added.iter().enumerate() {
                if added_alt.eq_ignore_ascii_case(allele) {
                    return site.alt_alleles.len() + 1 + i;
                }
            }

            // Add new allele
            added.push(allele.to_string());
            tracing::debug!(
                "Adding novel ALT {} to panel site {}:{} (index {})",
                allele,
                chrom,
                pos,
                site.alt_alleles.len() + added.len()
            );
            site.alt_alleles.len() + added.len()
        } else {
            // Site not in panel - create novel site
            // Check if we already have this novel site
            for site in &mut self.novel_sites {
                if site.chrom == chrom && site.pos == pos {
                    if let Some(idx) = site.allele_index(allele) {
                        return idx;
                    }
                    // Add allele to existing novel site
                    site.alt_alleles.push(allele.to_string());
                    return site.alt_alleles.len();
                }
            }

            // Create entirely new site
            let alt_alleles = if allele.eq_ignore_ascii_case(ref_base) {
                vec![]
            } else {
                vec![allele.to_string()]
            };

            let site = PanelSite {
                chrom: chrom.to_string(),
                pos,
                id: None,
                ref_allele: ref_base.to_string(),
                alt_alleles,
                alt_counts: Vec::new(),
                allele_number: None,
            };

            tracing::debug!(
                "Adding novel site {}:{} REF={} ALT={:?}",
                chrom,
                pos,
                ref_base,
                site.alt_alleles
            );

            let idx = if allele.eq_ignore_ascii_case(ref_base) {
                0
            } else {
                1
            };
            self.novel_sites.push(site);
            idx
        }
    }

    /// Get count of sites with added ALTs.
    pub fn modified_site_count(&self) -> usize {
        self.added_alts.len()
    }

    /// Get count of entirely novel sites.
    pub fn novel_site_count(&self) -> usize {
        self.novel_sites.len()
    }

    /// Get total number of sites (original + novel).
    pub fn total_site_count(&self) -> usize {
        self.original.len() + self.novel_sites.len()
    }

    /// Get reference to original panel.
    pub fn original(&self) -> &PanelIndex {
        &self.original
    }

    /// The original panel and the QC counters, borrowed together.
    pub fn original_and_qc(&mut self) -> (&PanelIndex, &mut PanelQc) {
        (&self.original, &mut self.qc)
    }

    /// Get the added ALTs for a site.
    /// Performs chr-prefix normalization to handle mismatched naming conventions.
    pub fn added_alts(&self, chrom: &str, pos: u64) -> Option<&Vec<String>> {
        // Try exact match first
        if let Some(v) = self.added_alts.get(&(chrom.to_string(), pos)) {
            return Some(v);
        }
        // Try with/without 'chr' prefix
        self.added_alts.get(&(alternate_chrom_name(chrom), pos))
    }

    /// Iterate over all novel sites.
    pub fn novel_sites(&self) -> impl Iterator<Item = &PanelSite> {
        self.novel_sites.iter()
    }
}

/// `chr1` <-> `1`.
fn alternate_chrom_name(chrom: &str) -> String {
    match chrom.strip_prefix("chr") {
        Some(bare) => bare.to_string(),
        None => format!("chr{chrom}"),
    }
}

/// Read INFO/AC (one per ALT) and INFO/AN from a panel record, when present.
///
/// Missing or malformed values yield `None` rather than an error: allele
/// frequency only informs the palindromic-SNP policy, and a panel without it
/// is still a valid panel.
fn allele_counts(
    info: &dyn noodles::vcf::variant::record::Info,
    header: &vcf::Header,
    n_alts: usize,
) -> (Vec<Option<u32>>, Option<u32>) {
    use noodles::vcf::variant::record::info::field::Value;
    use noodles::vcf::variant::record::info::field::value::Array;

    let to_u32 = |n: i32| u32::try_from(n).ok();
    let allele_number = match info.get(header, "AN") {
        Some(Ok(Some(Value::Integer(n)))) => to_u32(n),
        _ => None,
    };
    let mut alt_counts = vec![None; n_alts];
    match info.get(header, "AC") {
        Some(Ok(Some(Value::Integer(n)))) => {
            if let Some(slot) = alt_counts.first_mut() {
                *slot = to_u32(n);
            }
        }
        Some(Ok(Some(Value::Array(Array::Integer(values))))) => {
            for (slot, value) in alt_counts.iter_mut().zip(values.iter()) {
                *slot = value.ok().flatten().and_then(to_u32);
            }
        }
        _ => {}
    }
    (alt_counts, allele_number)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panel_site_allele_index() {
        let site = PanelSite::new("1", 1000, "A", &["G", "T"]);

        assert_eq!(site.allele_index("A"), Some(0));
        assert_eq!(site.allele_index("a"), Some(0)); // Case insensitive
        assert_eq!(site.allele_index("G"), Some(1));
        assert_eq!(site.allele_index("T"), Some(2));
        assert_eq!(site.allele_index("C"), None);
    }

    #[test]
    fn panel_site_has_allele() {
        let site = PanelSite::new("1", 1000, "A", &["G"]);

        assert!(site.has_allele("A"));
        assert!(site.has_allele("G"));
        assert!(!site.has_allele("T"));
    }

    #[test]
    fn split_multiallelic_records_are_all_kept() {
        // bcftools norm -m- writes A>C and A>G as two records at one POS. Both
        // must survive loading; the last one used to overwrite the first.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("panel.vcf");
        std::fs::write(
            &path,
            concat!(
                "##fileformat=VCFv4.3\n",
                "##contig=<ID=chr1,length=10000>\n",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ALT count\">\n",
                "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "chr1\t1500\t.\tA\tC\t.\tPASS\tAC=30;AN=100\n",
                "chr1\t1500\t.\tA\tG\t.\tPASS\tAC=5;AN=100\n",
                "chr1\t1600\t.\tT\tA\t.\tPASS\t.\n",
            ),
        )
        .unwrap();

        let index = PanelIndex::load(&path).unwrap();
        assert_eq!(index.len(), 3);
        let at_1500 = index.get_all("1", 1500);
        assert_eq!(at_1500.len(), 2, "both split records are kept");
        assert_eq!(at_1500[0].alt_alleles, vec!["C".to_string()]);
        assert_eq!(at_1500[1].alt_alleles, vec!["G".to_string()]);
        assert_eq!(at_1500[0].alt_frequency(0), Some(0.3));
        assert_eq!(at_1500[1].alt_frequency(0), Some(0.05));
        assert_eq!(index.get_all("chr1", 1600)[0].alt_frequency(0), None);
        let positions: Vec<u64> = index.sites_on("1").map(|s| s.pos).collect();
        assert_eq!(positions, vec![1500, 1500, 1600]);
    }
}
