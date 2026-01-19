use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Result, anyhow, bail};
use noodles::vcf::variant::record_buf::RecordBuf;
use rust_lapper::{Interval, Lapper};
use url::Url;

use crate::input::VariantSource;
use crate::reference::ReferenceGenome;
use crate::remote::fetch_remote_resource;
use crate::vcf_utils::remap_sample_genotypes;

/// Represents the strand of a genomic region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    fn from_char(c: char) -> Result<Self> {
        match c {
            '+' => Ok(Strand::Forward),
            '-' => Ok(Strand::Reverse),
            _ => Err(anyhow!("Invalid strand character: {}", c)),
        }
    }
}

/// Errors that can occur during liftover coordinate mapping.
/// These enable fail-closed behavior with explicit rejection reasons.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LiftoverError {
    /// No chain interval found for this position
    Unmapped,
    /// Multiple chain intervals found (ambiguous mapping)
    Ambiguous,
    /// Alleles are incompatible with target reference
    IncompatibleAlleles,
    /// Target contig not found in reference
    ContigNotFound,
    /// Indel/SV spans multiple chain intervals
    Straddled,
}

/// A mapped segment in the destination genome.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ChainMapping {
    /// Destination chromosome ID (internal)
    dest_chrom_id: u32,
    /// Destination start (qStart)
    dest_start: u64,
    /// Destination end (qEnd)
    dest_end: u64,
    /// Destination strand (qStrand)
    dest_strand: Strand,
    /// Destination size (qSize) - needed for reverse strand coordinate correction
    dest_size: u64,
    /// Source start (tStart) - needed for offset calculation
    source_start: u64,
}

/// A map of genomic intervals from Source to Destination.
///
/// In UCSC chain terminology:
/// - Source = Target (tName)
/// - Destination = Query (qName)
pub struct ChainMap {
    /// Map from Source Chromosome -> IntervalTree of Mappings
    map: HashMap<String, Lapper<u64, ChainMapping>>,
    /// Intern table for target chromosomes
    target_chroms: Vec<String>,
}

impl ChainMap {
    /// Load a chain file from a path.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Self::from_reader(reader)
    }

    /// Parse chain file content.
    pub fn from_reader<R: BufRead>(mut reader: R) -> Result<Self> {
        let mut map: HashMap<String, Vec<Interval<u64, ChainMapping>>> = HashMap::new();
        let mut target_chroms: Vec<String> = Vec::new();
        let mut target_chrom_indices: HashMap<String, u32> = HashMap::new();
        let mut line = String::new();

        // Current header info
        struct Header {
            t_name: String,
            // t_size: u64,
            // t_strand: char,
            t_start: u64,
            // t_end: u64, // Unused
            q_name: String,
            q_size: u64,
            q_strand: Strand,
            q_start: u64,
            // q_end: u64,
            // id: String,
        }
        let mut current_header: Option<Header> = None;

        while reader.read_line(&mut line)? > 0 {
            if line.trim().is_empty() || line.starts_with('#') {
                line.clear();
                continue;
            }

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.is_empty() {
                line.clear();
                continue;
            }

            if parts[0] == "chain" {
                // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                if parts.len() < 12 {
                    continue; // Invalid header?
                }
                let t_name = parts[2].to_string();
                let t_start: u64 = parts[5].parse()?;
                // let t_end: u64 = parts[6].parse()?; // Unused

                let q_name = parts[7].to_string();
                let q_size: u64 = parts[8].parse()?;
                let q_strand = Strand::from_char(parts[9].chars().next().unwrap())?;
                let q_start: u64 = parts[10].parse()?;

                current_header = Some(Header {
                    t_name,
                    t_start,
                    // t_end,
                    q_name,
                    q_size,
                    q_strand,
                    q_start,
                });
            } else {
                // Data line: size [dt dq]
                // size: length of ungapped alignment
                // dt: gap in target (source) to next block
                // dq: gap in query (dest) to next block
                if let Some(ref mut header) = current_header {
                    let size: u64 = parts[0].parse()?;

                    // Create interval for this block
                    let t_block_start = header.t_start;
                    let t_block_end = t_block_start + size;

                    let q_block_start = header.q_start;
                    let q_block_end = q_block_start + size;

                    // Intern target chromosome
                    let dest_chrom_id = *target_chrom_indices
                        .entry(header.q_name.clone())
                        .or_insert_with(|| {
                            let id = target_chroms.len() as u32;
                            target_chroms.push(header.q_name.clone());
                            id
                        });

                    let mapping = ChainMapping {
                        dest_chrom_id,
                        dest_start: q_block_start,
                        dest_end: q_block_end,
                        dest_strand: header.q_strand,
                        dest_size: header.q_size,
                        source_start: t_block_start,
                    };

                    map.entry(header.t_name.clone())
                        .or_default()
                        .push(Interval {
                            start: t_block_start,
                            stop: t_block_end,
                            val: mapping,
                        });

                    // Update header current position
                    header.t_start += size;
                    header.q_start += size;

                    if parts.len() == 3 {
                        let dt: u64 = parts[1].parse()?;
                        let dq: u64 = parts[2].parse()?;
                        header.t_start += dt;
                        header.q_start += dq;
                    }
                }
            }
            line.clear();
        }

        // Build Lappers
        let mut final_map = HashMap::new();
        for (chrom, intervals) in map {
            final_map.insert(chrom, Lapper::new(intervals));
        }

        Ok(Self {
            map: final_map,
            target_chroms,
        })
    }

    /// Lift a single coordinate (0-based) from source to destination.
    /// Returns (new_chrom, new_pos, strand) or an error if unmapped/ambiguous.
    ///
    /// Fails closed: ambiguous mappings (multiple chain hits) are rejected.
    pub fn lift(&self, chrom: &str, pos: u64) -> Result<(String, u64, Strand), LiftoverError> {
        // Try normalized chrom names
        let intervals = self
            .map
            .get(chrom)
            .or_else(|| self.map.get(chrom.trim_start_matches("chr")))
            .or_else(|| self.map.get(&format!("chr{}", chrom)))
            .ok_or(LiftoverError::Unmapped)?;

        // Find intersecting interval - use next() instead of collect() for performance
        let mut iter = intervals.find(pos, pos + 1);
        let first = iter.next().ok_or(LiftoverError::Unmapped)?;

        // Check for ambiguous multi-mapping (fail-closed)
        if iter.next().is_some() {
            return Err(LiftoverError::Ambiguous);
        }

        let mapping = &first.val;
        let offset = pos - mapping.source_start;

        let new_pos = if mapping.dest_strand == Strand::Forward {
            mapping.dest_start + offset
        } else {
            // Reverse strand logic
            let q_rev_pos = mapping.dest_start + offset;
            mapping
                .dest_size
                .checked_sub(1 + q_rev_pos)
                .ok_or(LiftoverError::Unmapped)?
        };

        let dest_chrom = self
            .target_chroms
            .get(mapping.dest_chrom_id as usize)
            .cloned()
            .ok_or(LiftoverError::Unmapped)?;

        Ok((dest_chrom, new_pos, mapping.dest_strand))
    }
}

/// Registry for managing chain files.
pub struct ChainRegistry {
    cache_dir: PathBuf,
}

impl ChainRegistry {
    pub fn new() -> Result<Self> {
        let dirs = directories::ProjectDirs::from("com", "convert_genome", "convert_genome")
            .ok_or_else(|| anyhow!("Failed to determine cache directory"))?;
        let cache_dir = dirs.cache_dir().join("chains");
        std::fs::create_dir_all(&cache_dir)?;
        Ok(Self { cache_dir })
    }

    pub fn get_chain(&self, source: &str, target: &str) -> Result<ChainMap> {
        let (filename, url) = match (source, target) {
            ("GRCh37", "GRCh38") | ("hg19", "hg38") => (
                "hg19ToHg38.over.chain.gz",
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
            ),
            ("GRCh38", "GRCh37") | ("hg38", "hg19") => (
                "hg38ToHg19.over.chain.gz",
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
            ),
            ("NCBI36", "GRCh38") | ("hg18", "hg38") => (
                "hg18ToHg38.over.chain.gz",
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz",
            ),
            _ => bail!("Unsupported liftover pair: {} -> {}", source, target),
        };

        let cache_path = self.cache_dir.join(filename);
        if !cache_path.exists() {
            tracing::info!("Chain file not found locally. Downloading from {}", url);
            // Use remote resource logic, but specific to keeping it in cache
            let url = Url::parse(url)?;
            let resource = fetch_remote_resource(&url)?;
            // Move from temp to cache
            // The resource unzips/ungzips. UCSC chains are .gz.
            // fetch_remote_resource returns a path to the decompressed file.
            std::fs::copy(resource.local_path(), &cache_path)?;
        }

        tracing::info!("Loading chain file: {}", cache_path.display());
        ChainMap::load(&cache_path)
    }
}

/// Helper to reverse complement a DNA string.
fn reverse_complement_string(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => c,
        })
        .collect()
}

/// Adapter that lifts variants from Source to Target.
pub struct LiftoverAdapter<S> {
    source: S,
    chain: Arc<ChainMap>,
    target_reference: ReferenceGenome,
}

impl<S> LiftoverAdapter<S> {
    pub fn new(source: S, chain: Arc<ChainMap>, target_reference: ReferenceGenome) -> Self {
        Self {
            source,
            chain,
            target_reference,
        }
    }
}

impl<S: VariantSource> VariantSource for LiftoverAdapter<S> {
    fn next_variant(
        &mut self,
        summary: &mut crate::ConversionSummary,
    ) -> Option<io::Result<RecordBuf>> {
        loop {
            // Get next record from source
            let mut record = match self.source.next_variant(summary)? {
                Ok(r) => r,
                Err(e) => return Some(Err(e)),
            };

            // Original coords
            let chrom = record.reference_sequence_name().to_string();
            let pos: usize = record.variant_start().map(|p| p.into()).unwrap_or(0);

            // Check for indels/SVs - reject them (not supported yet)
            let ref_bases = record.reference_bases();
            let alt_bases = record.alternate_bases();
            let is_snp = ref_bases.len() == 1
                && alt_bases
                    .as_ref()
                    .iter()
                    .all(|a| a.len() == 1 && !a.starts_with('<'));

            if !is_snp {
                // Reject indels and SVs - liftover not safe for them
                summary.liftover_straddled += 1;
                tracing::debug!(
                    chrom = %chrom,
                    pos = pos,
                    "Rejecting non-SNP variant (indel/SV liftover not supported)"
                );
                continue;
            }

            // 0-based for liftover
            let pos_0 = (pos as u64).saturating_sub(1);

            // Attempt lift with explicit error handling
            let lift_result = self.chain.lift(&chrom, pos_0);
            let (new_chrom, new_pos_0, strand) = match lift_result {
                Ok(result) => result,
                Err(LiftoverError::Unmapped) => {
                    summary.liftover_unmapped += 1;
                    continue;
                }
                Err(LiftoverError::Ambiguous) => {
                    summary.liftover_ambiguous += 1;
                    tracing::debug!(
                        chrom = %chrom,
                        pos = pos,
                        "Rejecting variant with ambiguous multi-mapping"
                    );
                    continue;
                }
                Err(_) => {
                    summary.liftover_unmapped += 1;
                    continue;
                }
            };

            // Canonicalize contig name using target reference
            let canonical_chrom = match self.target_reference.resolve_contig_name(&new_chrom) {
                Some(name) => name.to_string(),
                None => {
                    // Try with chr prefix stripped/added
                    let alt_name = if new_chrom.starts_with("chr") {
                        new_chrom.trim_start_matches("chr").to_string()
                    } else {
                        format!("chr{}", new_chrom)
                    };
                    match self.target_reference.resolve_contig_name(&alt_name) {
                        Some(name) => name.to_string(),
                        None => {
                            summary.liftover_contig_missing += 1;
                            tracing::debug!(
                                new_chrom = %new_chrom,
                                "Target contig not found in reference"
                            );
                            continue;
                        }
                    }
                }
            };

            // Update Coords (1-based)
            let new_pos = (new_pos_0 + 1) as usize;

            // Convert to Position
            let new_pos_obj = match noodles::core::Position::new(new_pos) {
                Some(p) => p,
                None => continue,
            };

            // Handle Strand - reverse complement if needed
            if strand == Strand::Reverse {
                let new_ref = reverse_complement_string(ref_bases);
                *record.reference_bases_mut() = new_ref;

                let new_alts: Vec<String> = alt_bases
                    .as_ref()
                    .iter()
                    .map(|a| {
                        if a.starts_with('<') {
                            a.clone()
                        } else {
                            reverse_complement_string(a)
                        }
                    })
                    .collect();
                *record.alternate_bases_mut() = new_alts.into();
            }

            // Get target reference base for compatibility check
            let target_base = match self.target_reference.base(&canonical_chrom, new_pos as u64) {
                Ok(b) => b.to_ascii_uppercase(),
                Err(_) => {
                    summary.reference_failures += 1;
                    continue;
                }
            };

            // Get current alleles after potential strand flip
            let rec_ref = record.reference_bases().to_ascii_uppercase();
            let rec_alts: Vec<String> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|a| a.to_ascii_uppercase())
                .collect();

            // Build set of all called alleles
            let mut all_alleles: Vec<String> = vec![rec_ref.clone()];
            all_alleles.extend(rec_alts.iter().cloned());

            // ALLELE COMPATIBILITY CHECK (fail-closed)
            // For DTC N-ref mode: we need the target base to be one of the alleles
            // For VCF: we need either REF matches target, or we can do a ref swap
            let target_base_str = target_base.to_string();

            if rec_ref == "N" {
                // DTC N-ref mode: target base must be in the alt alleles
                // We'll set REF = target and remove it from ALTs
                if !rec_alts.iter().any(|a| a.eq_ignore_ascii_case(&target_base_str)) {
                    // Target base not in alleles - incompatible
                    summary.liftover_incompatible += 1;
                    tracing::debug!(
                        chrom = %chrom,
                        pos = pos,
                        target_base = %target_base,
                        alleles = ?rec_alts,
                        "Rejecting: target reference base not in called alleles"
                    );
                    continue;
                }

                // Set REF to target base
                *record.reference_bases_mut() = target_base_str.clone();

                // Remove target base from ALTs and remap genotypes
                let mut new_alts = Vec::new();
                let mut allele_mapping = HashMap::new();
                allele_mapping.insert(0, 0); // Old REF (N, unused) -> New REF

                for (i, alt) in rec_alts.iter().enumerate() {
                    let old_idx = i + 1;
                    if alt.eq_ignore_ascii_case(&target_base_str) {
                        // This ALT is now REF
                        allele_mapping.insert(old_idx, 0);
                    } else {
                        new_alts.push(alt.clone());
                        allele_mapping.insert(old_idx, new_alts.len());
                    }
                }

                *record.alternate_bases_mut() = new_alts.into();
                *record.samples_mut() =
                    remap_sample_genotypes(record.samples(), &allele_mapping);
            } else {
                // VCF mode: REF should match target, or we reject
                // (Reference swap is complex and error-prone - fail closed)
                if !rec_ref.eq_ignore_ascii_case(&target_base_str) {
                    // Check if target is in ALTs for possible ref swap
                    if rec_alts.iter().any(|a| a.eq_ignore_ascii_case(&target_base_str)) {
                        // Could do ref swap, but that's complex and risky
                        // For fail-closed, we reject
                        summary.liftover_incompatible += 1;
                        tracing::debug!(
                            chrom = %chrom,
                            pos = pos,
                            rec_ref = %rec_ref,
                            target_base = %target_base,
                            "Rejecting: REF mismatch (ref swap not implemented)"
                        );
                        continue;
                    } else {
                        // Neither REF nor ALTs match target - definitely incompatible
                        summary.liftover_incompatible += 1;
                        tracing::debug!(
                            chrom = %chrom,
                            pos = pos,
                            rec_ref = %rec_ref,
                            target_base = %target_base,
                            "Rejecting: alleles incompatible with target reference"
                        );
                        continue;
                    }
                }
                // REF matches target - good to go
            }

            // Update record with canonical chrom and new position
            *record.reference_sequence_name_mut() = canonical_chrom;
            *record.variant_start_mut() = Some(new_pos_obj);

            return Some(Ok(record));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use assert_fs::TempDir;
    use assert_fs::prelude::*;

    // Helper to create a dummy chain file
    fn create_dummy_chain(dir: &TempDir, filename: &str, content: &str) -> PathBuf {
        let path = dir.child(filename);
        path.write_str(content).unwrap();
        path.path().to_path_buf()
    }

    #[test]
    fn test_chain_map_loading() {
        let dir = TempDir::new().unwrap();
        let chain_content = "chain 100 chr1 1000 + 0 1000 chr1 1000 + 0 1000 1\n100 0 0\n";
        let chain_path = create_dummy_chain(&dir, "test.chain", chain_content);

        let chain_map = ChainMap::load(chain_path).unwrap();
        assert!(chain_map.map.contains_key("chr1"));
        let lapper = chain_map.map.get("chr1").unwrap();
        assert_eq!(lapper.len(), 1);
        assert_eq!(chain_map.target_chroms.len(), 1);
        assert_eq!(chain_map.target_chroms[0], "chr1");
    }

    #[test]
    fn test_lift_point_forward() {
        let dir = TempDir::new().unwrap();
        // chr1:100-200 maps to chr1:200-300 on + strand
        let chain_content = "chain 100 chr1 1000 + 100 200 chr1 1000 + 200 300 1\n100 0 0\n";
        let chain_path = create_dummy_chain(&dir, "fwd.chain", chain_content);
        let chain_map = ChainMap::load(chain_path).unwrap();

        // 0-based input. 100 -> 200
        let (chrom, pos, strand) = chain_map.lift("chr1", 100).unwrap();
        assert_eq!(chrom, "chr1");
        assert_eq!(pos, 200);
        assert_eq!(strand, Strand::Forward);

        // 150 -> 250
        let (_, pos, _) = chain_map.lift("chr1", 150).unwrap();
        assert_eq!(pos, 250);
    }

    #[test]
    fn test_lift_point_reverse() {
        let dir = TempDir::new().unwrap();
        // chr1:100-200 maps to chr1:200-300 on - strand
        // qStart=200, qEnd=300, qSize=1000.
        // If qStrand is -, coords in chain are reversed.

        let chain_content = "chain 100 chr1 1000 + 100 200 chr1 1000 - 100 200 1\n100 0 0\n";
        let chain_path = create_dummy_chain(&dir, "rev.chain", chain_content);
        let chain_map = ChainMap::load(chain_path).unwrap();

        // Input 100 (start of block).
        // Offset = 0.
        // q_rev_pos = 100 + 0 = 100.
        // physical_pos = size - 1 - q_rev_pos = 1000 - 1 - 100 = 899.
        let (_, pos, strand) = chain_map.lift("chr1", 100).unwrap();
        assert_eq!(strand, Strand::Reverse);
        assert_eq!(pos, 899);

        // Input 199 (end of block - 1).
        // Offset = 99.
        // q_rev_pos = 100 + 99 = 199.
        // physical_pos = 1000 - 1 - 199 = 800.
        let (_, pos, _) = chain_map.lift("chr1", 199).unwrap();
        assert_eq!(pos, 800);
    }

    #[test]
    fn test_liftover_adapter_round_trip() {
        // Round trip requires two chains: A->B and B->A.
        // A: chr1 size 1000.
        // B: chr1 size 1000.
        // Mapping: 100..200 -> 200..300 (Offset +100).

        let dir = TempDir::new().unwrap();
        let chain_ab_content = "chain 100 chr1 1000 + 100 200 chr1 1000 + 200 300 1\n100 0 0\n";
        let chain_ba_content = "chain 100 chr1 1000 + 200 300 chr1 1000 + 100 200 1\n100 0 0\n";

        let chain_ab_path = create_dummy_chain(&dir, "ab.chain", chain_ab_content);
        let chain_ba_path = create_dummy_chain(&dir, "ba.chain", chain_ba_content);

        let map_ab = Arc::new(ChainMap::load(chain_ab_path).unwrap());
        let map_ba = Arc::new(ChainMap::load(chain_ba_path).unwrap());

        // Create a dummy source record at 150.
        // Should lift to 250.
        // Then lift back to 150.

        let (c2, p2, _) = map_ab.lift("chr1", 150).unwrap();
        assert_eq!(p2, 250);

        let (_, p3, _) = map_ba.lift(&c2, p2).unwrap();
        assert_eq!(p3, 150);
    }

    #[test]
    fn test_lift_unmapped_returns_error() {
        let dir = TempDir::new().unwrap();
        let chain_content = "chain 100 chr1 1000 + 100 200 chr1 1000 + 200 300 1\n100 0 0\n";
        let chain_path = create_dummy_chain(&dir, "test.chain", chain_content);
        let chain_map = ChainMap::load(chain_path).unwrap();

        // Position outside the mapped interval should return Unmapped
        let result = chain_map.lift("chr1", 50);
        assert_eq!(result, Err(LiftoverError::Unmapped));

        // Unknown chromosome should return Unmapped
        let result = chain_map.lift("chr99", 150);
        assert_eq!(result, Err(LiftoverError::Unmapped));
    }

    #[test]
    fn test_lift_ambiguous_returns_error() {
        let dir = TempDir::new().unwrap();
        // Create overlapping intervals - same source maps to two destinations
        // This simulates a segmental duplication scenario
        let chain_content = concat!(
            "chain 100 chr1 1000 + 100 200 chr1 1000 + 200 300 1\n",
            "100 0 0\n",
            "\n",
            "chain 90 chr1 1000 + 100 200 chr2 1000 + 400 500 2\n",
            "100 0 0\n",
        );
        let chain_path = create_dummy_chain(&dir, "ambig.chain", chain_content);
        let chain_map = ChainMap::load(chain_path).unwrap();

        // Position 150 maps to both chains - should be ambiguous
        let result = chain_map.lift("chr1", 150);
        assert_eq!(result, Err(LiftoverError::Ambiguous));
    }
}
