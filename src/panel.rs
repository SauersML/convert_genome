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
}

impl PanelSite {
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

/// Index of panel sites for fast lookup by position.
#[derive(Default)]
pub struct PanelIndex {
    sites: HashMap<(String, u64), PanelSite>,
    /// Chromosome order from the panel header for sorting output
    chrom_order: Vec<String>,
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

        let mut sites = HashMap::with_capacity(estimate_panel_capacity(path));
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

            let site = PanelSite {
                chrom: chrom.clone(),
                pos,
                id: None,
                ref_allele,
                alt_alleles,
            };

            sites.insert((chrom, pos), site);

            if record_count.is_multiple_of(1_000_000) {
                tracing::info!("Loaded {} panel sites...", record_count);
            }
        }

        tracing::info!("Loaded {} panel sites total", sites.len());

        Ok(Self { sites, chrom_order })
    }

    fn load_bcf<P: AsRef<Path>>(path: P) -> Result<Self> {
        use noodles::bcf;

        let file = std::fs::File::open(path.as_ref())
            .with_context(|| format!("failed to open panel {}", path.as_ref().display()))?;

        let mut bcf_reader = bcf::io::Reader::new(BufReader::new(file));
        let header = bcf_reader.read_header()?;

        let chrom_order: Vec<String> = header.contigs().keys().map(|k| k.to_string()).collect();

        let mut sites = HashMap::with_capacity(estimate_panel_capacity(path.as_ref()));
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

            let ref_allele = std::str::from_utf8(record.reference_bases().as_ref())?
                .to_string();
            let alt_alleles: Vec<String> = record
                .alternate_bases()
                .iter()
                .map(|a| a.map(|s| s.to_string()))
                .collect::<std::io::Result<Vec<_>>>()?;

            let site = PanelSite {
                chrom: chrom.clone(),
                pos,
                id: None,
                ref_allele,
                alt_alleles,
            };

            sites.insert((chrom, pos), site);

            if record_count.is_multiple_of(1_000_000) {
                tracing::info!("Loaded {} panel sites...", record_count);
            }
        }

        tracing::info!("Loaded {} panel sites total", sites.len());

        Ok(Self { sites, chrom_order })
    }

    /// Look up a site by chromosome and position.
    pub fn get(&self, chrom: &str, pos: u64) -> Option<&PanelSite> {
        // Try exact match first
        if let Some(site) = self.sites.get(&(chrom.to_string(), pos)) {
            return Some(site);
        }
        // Try with/without 'chr' prefix
        let alt_chrom = if chrom.starts_with("chr") {
            chrom.strip_prefix("chr").unwrap().to_string()
        } else {
            format!("chr{}", chrom)
        };
        self.sites.get(&(alt_chrom, pos))
    }

    /// Check if a site exists at the given position.
    pub fn contains(&self, chrom: &str, pos: u64) -> bool {
        self.get(chrom, pos).is_some()
    }

    /// Get the number of sites in the index.
    pub fn len(&self) -> usize {
        self.sites.len()
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }

    /// Get chromosome order from the panel header.
    pub fn chrom_order(&self) -> &[String] {
        &self.chrom_order
    }
}

/// Tracks modifications to the panel (novel alleles and sites).
pub struct PaddedPanel {
    /// Original panel index (immutable reference)
    original: PanelIndex,
    /// Additional ALTs added to existing sites: (chrom, pos) -> new ALTs to append
    added_alts: HashMap<(String, u64), Vec<String>>,
    /// Entirely novel sites not in original panel
    novel_sites: Vec<PanelSite>,
}

impl PaddedPanel {
    /// Create a new padded panel from an original panel index.
    pub fn new(original: PanelIndex) -> Self {
        Self {
            original,
            added_alts: HashMap::new(),
            novel_sites: Vec::new(),
        }
    }

    /// Get the original panel site, if any.
    pub fn get_original(&self, chrom: &str, pos: u64) -> Option<&PanelSite> {
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

    /// Get the added ALTs for a site.
    /// Performs chr-prefix normalization to handle mismatched naming conventions.
    pub fn added_alts(&self, chrom: &str, pos: u64) -> Option<&Vec<String>> {
        // Try exact match first
        if let Some(v) = self.added_alts.get(&(chrom.to_string(), pos)) {
            return Some(v);
        }
        // Try with/without 'chr' prefix
        let alt_chrom = if chrom.starts_with("chr") {
            chrom.strip_prefix("chr").unwrap().to_string()
        } else {
            format!("chr{}", chrom)
        };
        self.added_alts.get(&(alt_chrom, pos))
    }

    /// Iterate over all novel sites.
    pub fn novel_sites(&self) -> impl Iterator<Item = &PanelSite> {
        self.novel_sites.iter()
    }
}

fn estimate_panel_capacity(path: &Path) -> usize {
    let default_capacity = 1_000_000usize;
    let Ok(meta) = std::fs::metadata(path) else {
        return default_capacity;
    };
    let len = meta.len();
    if len == 0 {
        return default_capacity;
    }
    let estimated = (len / 80) as usize;
    estimated.clamp(100_000, 8_000_000)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panel_site_allele_index() {
        let site = PanelSite {
            chrom: "1".to_string(),
            pos: 1000,
            id: None,
            ref_allele: "A".to_string(),
            alt_alleles: vec!["G".to_string(), "T".to_string()],
        };

        assert_eq!(site.allele_index("A"), Some(0));
        assert_eq!(site.allele_index("a"), Some(0)); // Case insensitive
        assert_eq!(site.allele_index("G"), Some(1));
        assert_eq!(site.allele_index("T"), Some(2));
        assert_eq!(site.allele_index("C"), None);
    }

    #[test]
    fn panel_site_has_allele() {
        let site = PanelSite {
            chrom: "1".to_string(),
            pos: 1000,
            id: None,
            ref_allele: "A".to_string(),
            alt_alleles: vec!["G".to_string()],
        };

        assert!(site.has_allele("A"));
        assert!(site.has_allele("G"));
        assert!(!site.has_allele("T"));
    }
}
