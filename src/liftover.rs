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
    /// Chain ID (from chain header)
    chain_id: u32,
    /// Chain score (from chain header)
    chain_score: u64,
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
    map: HashMap<String, Arc<Lapper<u64, ChainMapping>>>,
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
            score: u64,
            t_name: String,
            // t_size: u64,
            // t_strand: char,
            t_start: u64,
            // t_end: u64, // Unused
            q_name: String,
            q_size: u64,
            q_strand: Strand,
            q_start: u64,
            chain_id: u32,
            // q_end: u64,
            // id: String,
        }
        let mut current_header: Option<Header> = None;

        while reader.read_line(&mut line)? > 0 {
            if line.trim().is_empty() || line.starts_with('#') {
                line.clear();
                continue;
            }

            let mut it = line.split_whitespace();
            let Some(first) = it.next() else {
                line.clear();
                continue;
            };

            if first == "chain" {
                // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                let Some(score_s) = it.next() else {
                    line.clear();
                    continue;
                };
                let Some(t_name_s) = it.next() else {
                    line.clear();
                    continue;
                };
                // Skip tSize, tStrand
                let _ = it.next();
                let _ = it.next();
                let Some(t_start_s) = it.next() else {
                    line.clear();
                    continue;
                };
                // let t_end: u64 = parts[6].parse()?; // Unused

                // Skip tEnd
                let _ = it.next();

                let Some(q_name_s) = it.next() else {
                    line.clear();
                    continue;
                };
                let Some(q_size_s) = it.next() else {
                    line.clear();
                    continue;
                };
                let Some(q_strand_s) = it.next() else {
                    line.clear();
                    continue;
                };
                let Some(q_start_s) = it.next() else {
                    line.clear();
                    continue;
                };

                // Skip qEnd
                let _ = it.next();

                let Some(chain_id_s) = it.next() else {
                    line.clear();
                    continue;
                };

                let score: u64 = score_s.parse()?;
                let t_name = t_name_s.to_string();
                let t_start: u64 = t_start_s.parse()?;

                let q_name = q_name_s.to_string();
                let q_size: u64 = q_size_s.parse()?;
                let q_strand = Strand::from_char(q_strand_s.chars().next().unwrap())?;
                let q_start: u64 = q_start_s.parse()?;
                let chain_id: u32 = chain_id_s.parse()?;

                current_header = Some(Header {
                    score,
                    t_name,
                    t_start,
                    // t_end,
                    q_name,
                    q_size,
                    q_strand,
                    q_start,
                    chain_id,
                });
            } else {
                // Data line: size [dt dq]
                // size: length of ungapped alignment
                // dt: gap in target (source) to next block
                // dq: gap in query (dest) to next block
                if let Some(ref mut header) = current_header {
                    let size_s = first;
                    let size: u64 = size_s.parse()?;

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
                        chain_id: header.chain_id,
                        chain_score: header.score,
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

                    let dt_opt = it.next();
                    let dq_opt = it.next();
                    if let (Some(dt_s), Some(dq_s)) = (dt_opt, dq_opt) {
                        let dt: u64 = dt_s.parse()?;
                        let dq: u64 = dq_s.parse()?;
                        header.t_start += dt;
                        header.q_start += dq;
                    }
                }
            }
            line.clear();
        }

        // Build Lappers and insert both chr and non-chr aliases to avoid allocations in lift().
        let mut final_map: HashMap<String, Arc<Lapper<u64, ChainMapping>>> = HashMap::new();
        for (chrom, intervals) in map {
            let lapper = Arc::new(Lapper::new(intervals));

            // Insert exact key
            final_map.insert(chrom.clone(), Arc::clone(&lapper));

            // Insert alias key with/without "chr" prefix.
            if let Some(stripped) = chrom.strip_prefix("chr") {
                final_map.insert(stripped.to_string(), Arc::clone(&lapper));
            } else {
                final_map.insert(format!("chr{}", chrom), Arc::clone(&lapper));
            }
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
        // Try normalized chrom names (no allocations in hot loop).
        let intervals = self
            .map
            .get(chrom)
            .or_else(|| self.map.get(chrom.trim_start_matches("chr")))
            .ok_or(LiftoverError::Unmapped)?;

        self.lift_with_intervals(intervals, pos)
    }

    fn lift_with_intervals(
        &self,
        intervals: &Lapper<u64, ChainMapping>,
        pos: u64,
    ) -> Result<(String, u64, Strand), LiftoverError> {
        // Resolve overlaps by selecting the highest-score chain.
        // If there is a score tie across different chain IDs, fail-closed as ambiguous.
        let mut iter = intervals.find(pos, pos + 1);
        let mut best: Option<&Interval<u64, ChainMapping>> = None;
        let mut tied = false;
        while let Some(hit) = iter.next() {
            match best {
                None => {
                    best = Some(hit);
                    tied = false;
                }
                Some(current) => {
                    let a = &current.val;
                    let b = &hit.val;
                    if b.chain_score > a.chain_score {
                        best = Some(hit);
                        tied = false;
                    } else if b.chain_score == a.chain_score && b.chain_id != a.chain_id {
                        tied = true;
                    }
                }
            }
        }

        let best = best.ok_or(LiftoverError::Unmapped)?;
        if tied {
            return Err(LiftoverError::Ambiguous);
        }

        let mapping = &best.val;
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
                "hg19ToHg38.over.chain",
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
            ),
            ("GRCh38", "GRCh37") | ("hg38", "hg19") => (
                "hg38ToHg19.over.chain",
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
            ),
            ("NCBI36", "GRCh38") | ("hg18", "hg38") => (
                "hg18ToHg38.over.chain",
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
    source_reference: Option<ReferenceGenome>,
    last_chrom_key: Option<String>,
    last_intervals: Option<Arc<Lapper<u64, ChainMapping>>>,
}

impl<S> LiftoverAdapter<S> {
    pub fn new(
        source: S,
        chain: Arc<ChainMap>,
        target_reference: ReferenceGenome,
        source_reference: Option<ReferenceGenome>,
    ) -> Self {
        Self {
            source,
            chain,
            target_reference,
            source_reference,
            last_chrom_key: None,
            last_intervals: None,
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
            let chrom_str = chrom.as_str();
            let pos: usize = record.variant_start().map(|p| p.into()).unwrap_or(0);

            // VCF coordinates are 1-based. ChainMap::lift operates on 0-based coordinates.
            // If position is missing/invalid, fail-closed by skipping the record.
            if pos == 0 {
                summary.parse_errors += 1;
                continue;
            }

            // If the record is a DTC placeholder (REF=N) and we have a source reference,
            // pre-fill REF from the source build and optionally complement alleles if they
            // only match the source reference base after complementation.
            if record.reference_bases().eq_ignore_ascii_case("N") {
                if let Some(ref src_ref) = self.source_reference {
                    if let Ok(src_base) = src_ref.base(chrom_str, pos as u64) {
                        let src_base = src_base.to_ascii_uppercase();
                        if matches!(src_base, 'A' | 'C' | 'G' | 'T') {
                            let src_base_str = src_base.to_string();
                            let comp_src_base_str = reverse_complement_string(&src_base_str);

                            // Only attempt per-record complement rescue for simple SNP-like alleles.
                            let alts: Vec<String> = record
                                .alternate_bases()
                                .as_ref()
                                .iter()
                                .map(|s| s.to_ascii_uppercase())
                                .collect();

                            let is_simple_alleles = !alts.is_empty()
                                && alts.iter().all(|a| a.len() == 1)
                                && !alts.iter().any(|a| a.starts_with('<'));

                            if is_simple_alleles {
                                let has_src = alts.iter().any(|a| a == &src_base_str);
                                let has_comp_src = alts.iter().any(|a| a == &comp_src_base_str);

                                if !has_src && has_comp_src {
                                    // Complement all alleles in-place so source ref base appears.
                                    let new_alts: Vec<String> = alts
                                        .iter()
                                        .map(|a| reverse_complement_string(a))
                                        .collect();
                                    *record.alternate_bases_mut() = new_alts.into();
                                }
                            }

                            // Always set REF to the true source base for downstream lifting.
                            *record.reference_bases_mut() = src_base_str;
                        }
                    }
                }
            }

            // Check for indels/SVs - reject them (not supported yet)
            // Check for indels/SVs - reject them (not supported yet)
            // Clone to avoid borrow conflicts later
            let ref_bases = record.reference_bases().to_string();
            let alt_bases: Vec<String> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|s| s.to_string())
                .collect();
            let is_symbolic = alt_bases.iter().any(|a| a.starts_with('<'));
            let is_snp = ref_bases.len() == 1 && alt_bases.iter().all(|a| a.len() == 1) && !is_symbolic;

            // Basic indel endpoint liftover (non-symbolic only). SV/symbolic remain unsupported.
            let is_simple_indel = !is_snp
                && !is_symbolic
                && ref_bases.len() >= 1
                && alt_bases.iter().all(|a| !a.is_empty());

            // Convert to 0-based for liftover (safe because pos > 0).
            let pos_0 = (pos as u64) - 1;

            // Cache chain interval structure per chromosome to avoid repeated HashMap hashing.
            // Inputs are typically sorted by chromosome, so this is a hot-path win.
            let intervals = if let (Some(k), Some(ref intervals)) =
                (self.last_chrom_key.as_deref(), self.last_intervals.as_ref())
                && k == chrom_str
            {
                Arc::clone(intervals)
            } else {
                let trimmed = chrom_str.trim_start_matches("chr");
                let (key_used, intervals) = match self.chain.map.get(chrom_str) {
                    Some(v) => (chrom_str, Arc::clone(v)),
                    None => match self.chain.map.get(trimmed) {
                        Some(v) => (trimmed, Arc::clone(v)),
                        None => {
                            summary.liftover_unmapped += 1;
                            continue;
                        }
                    },
                };
                self.last_chrom_key = Some(key_used.to_string());
                self.last_intervals = Some(Arc::clone(&intervals));
                intervals
            };

            let (new_chrom, new_pos_0, strand) = if is_snp {
                // Attempt lift with explicit error handling
                let lift_result = self.chain.lift_with_intervals(&intervals, pos_0);
                match lift_result {
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
                }
            } else if is_simple_indel {
                // Endpoint liftover: lift start and end of the REF span.
                let ref_len = ref_bases.len() as u64;
                let source_dist = ref_len.saturating_sub(1);
                let end_pos_0 = pos_0 + source_dist;

                let (c1, p1, s1) = match self.chain.lift_with_intervals(&intervals, pos_0) {
                    Ok(r) => r,
                    Err(LiftoverError::Unmapped) => {
                        summary.liftover_unmapped += 1;
                        continue;
                    }
                    Err(LiftoverError::Ambiguous) => {
                        summary.liftover_ambiguous += 1;
                        continue;
                    }
                    Err(_) => {
                        summary.liftover_unmapped += 1;
                        continue;
                    }
                };

                let (c2, p2, s2) = match self.chain.lift_with_intervals(&intervals, end_pos_0) {
                    Ok(r) => r,
                    Err(_) => {
                        summary.liftover_straddled += 1;
                        continue;
                    }
                };

                if c1 != c2 || s1 != s2 {
                    summary.liftover_straddled += 1;
                    continue;
                }

                let target_dist = if s1 == Strand::Forward {
                    p2.saturating_sub(p1)
                } else {
                    p1.saturating_sub(p2)
                };

                if target_dist != source_dist {
                    summary.liftover_straddled += 1;
                    continue;
                }

                // Ensure VCF POS in target is leftmost (forward coordinate).
                let new_start = if s1 == Strand::Forward { p1 } else { p2 };
                (c1, new_start, s1)
            } else {
                summary.liftover_straddled += 1;
                tracing::debug!(
                    chrom = %chrom,
                    pos = pos,
                    "Rejecting unsupported non-SNP variant (symbolic/SV or complex indel)"
                );
                continue;
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
                let new_ref = reverse_complement_string(&ref_bases);
                *record.reference_bases_mut() = new_ref;

                let new_alts: Vec<String> = alt_bases
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
                // DTC N-ref mode: input REF is a placeholder.
                // Always set REF = target base and treat any called allele that equals the
                // target base as a REF allele (GT index 0). Other called alleles become ALTs.
                *record.reference_bases_mut() = target_base_str.clone();

                let mut new_alts: Vec<String> = Vec::new();
                let mut allele_to_new_idx: HashMap<String, usize> = HashMap::new();
                let mut allele_mapping: HashMap<usize, usize> = HashMap::new();

                // Old REF (N placeholder) -> New REF
                allele_mapping.insert(0, 0);

                for (i, alt) in rec_alts.iter().enumerate() {
                    let old_idx = i + 1;

                    if alt.eq_ignore_ascii_case(&target_base_str) {
                        // This called allele matches the target reference base.
                        allele_mapping.insert(old_idx, 0);
                        continue;
                    }

                    // De-duplicate ALTs while maintaining deterministic indices.
                    let key = alt.to_ascii_uppercase();
                    let new_idx = if let Some(idx) = allele_to_new_idx.get(&key) {
                        *idx
                    } else {
                        new_alts.push(alt.clone());
                        let idx = new_alts.len();
                        allele_to_new_idx.insert(key, idx);
                        idx
                    };

                    allele_mapping.insert(old_idx, new_idx);
                }

                *record.alternate_bases_mut() = new_alts.into();
                *record.samples_mut() = remap_sample_genotypes(record.samples(), &allele_mapping);
            } else {
                // VCF mode: REF should match target, or we check for swap
                if !rec_ref.eq_ignore_ascii_case(&target_base_str) {
                    // Check if target is in ALTs for possible ref swap
                    if let Some(target_idx) = rec_alts
                        .iter()
                        .position(|a| a.eq_ignore_ascii_case(&target_base_str))
                    {
                        // Ref swap logic
                        // 1. New REF = Target Base
                        // 2. Old REF becomes an ALT
                        // 3. New ALTs = [Old REF, Other ALTs...]
                        
                        let mut new_alts = Vec::new();
                        let mut allele_mapping = HashMap::new();
                        
                        // Map Old REF (0) -> New position in ALTs (1-based index)
                        // By convention, we'll put the old REF as the first ALT
                        new_alts.push(rec_ref.clone());
                        allele_mapping.insert(0, 1);
                        
                        // Map the Old ALT that matches Target -> New REF (0)
                        // target_idx is 0-based index into rec_alts. So VCF index is target_idx + 1.
                        allele_mapping.insert(target_idx + 1, 0);
                        
                        // Handle other ALTs
                        for (i, alt) in rec_alts.iter().enumerate() {
                            if i == target_idx { continue; } // Already handled (became REF)
                            new_alts.push(alt.clone());
                            // Mapped to new position: existing length of new_alts (1-based)
                            allele_mapping.insert(i + 1, new_alts.len()); 
                        }
                        
                        // Apply changes
                        *record.reference_bases_mut() = target_base_str.clone();
                         *record.alternate_bases_mut() = new_alts.into();
                        *record.samples_mut() = remap_sample_genotypes(record.samples(), &allele_mapping);
                        
                        tracing::debug!(
                            chrom = %chrom,
                            pos = pos,
                            "Performed Reference Swap: OldREF({}) became ALT, OldALT({}) became REF",
                            rec_ref, target_base_str
                        );
                        
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

        // Position 150 maps to both chains - highest score chain should win
        let (chrom, pos, strand) = chain_map.lift("chr1", 150).unwrap();
        assert_eq!(chrom, "chr1");
        assert_eq!(pos, 250);
        assert_eq!(strand, Strand::Forward);
    }

    #[test]
    fn test_lift_overlap_score_tie_is_ambiguous() {
        let dir = TempDir::new().unwrap();
        // Two overlapping chains with the same score but different IDs.
        // This should fail-closed as ambiguous.
        let chain_content = concat!(
            "chain 100 chr1 1000 + 100 200 chr1 1000 + 200 300 1\n",
            "100 0 0\n",
            "\n",
            "chain 100 chr1 1000 + 100 200 chr2 1000 + 400 500 2\n",
            "100 0 0\n",
        );
        let chain_path = create_dummy_chain(&dir, "tie.chain", chain_content);
        let chain_map = ChainMap::load(chain_path).unwrap();

        let result = chain_map.lift("chr1", 150);
        assert_eq!(result, Err(LiftoverError::Ambiguous));
    }

    #[test]
    fn test_ref_swap_logic() {
        use noodles::vcf;
        use noodles::core::Position;
        use noodles::vcf::variant::record_buf::Samples;
        use std::io::Write;
        use crate::input::VariantSource; 

        // 1. Create Target Reference (chr2) with 'G' at position 100
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.child("target.fa");
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            writeln!(f, ">chr2").unwrap();
            let mut seq = vec!['N'; 200];
            seq[100] = 'G'; 
            let s: String = seq.into_iter().collect();
            writeln!(f, "{}", s).unwrap();
        }
        let target_ref = ReferenceGenome::open(fasta_path.path(), None).unwrap();

        // 2. Chain Map (chr1:100 -> chr2:100)
        // ChainMap maps FROM tName/tStart TO qName/qStart.
        // So we need tName=chr1.
        let chain_content = "chain 1000 chr1 200 + 100 101 chr2 200 + 100 101 1\n1 0 0\n";
        let chain_path = create_dummy_chain(&dir, "swap.chain", chain_content);
        let chain_map = ChainMap::load(chain_path).unwrap();

        // 3. Mock Source
        struct MockSource {
            record: Option<vcf::variant::RecordBuf>,
        }
        impl VariantSource for MockSource {
            fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<std::io::Result<vcf::variant::RecordBuf>> {
                // Silence unused warning by accessing it
                summary.parse_errors += 0; 
                self.record.take().map(Ok)
            }
        }

        let input_record = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::new(101).unwrap())
            .set_reference_bases("A")
            .set_alternate_bases(vec![String::from("G")].into())
            .set_samples(Samples::new(
                 vec![String::from("GT")].into_iter().collect(),
                 vec![vec![Some(vcf::variant::record_buf::samples::sample::Value::String("0/1".into()))]]
            ))
            .build();

        let source = MockSource { record: Some(input_record) };

        let mut adapter = LiftoverAdapter {
            source: Box::new(source),
            chain: Arc::new(chain_map),
            target_reference: target_ref,
            source_reference: None,
            last_chrom_key: None,
            last_intervals: None,
        };

        let mut summary = crate::ConversionSummary::default();
        let result = adapter.next_variant(&mut summary).unwrap().unwrap();

        // 4. Verify Swap
        assert_eq!(result.reference_bases().to_string(), "G");
        assert_eq!(result.alternate_bases().as_ref()[0], "A");
        
        let samples = result.samples();
        let sample = samples.values().next().unwrap();
        let genotype = sample.values().iter().next().unwrap().as_ref().unwrap();
        
        if let vcf::variant::record_buf::samples::sample::Value::String(gt) = genotype {
            assert!(gt == "1/0" || gt == "1|0", "Genotype should be swapped to 1/0, got {}", gt);
        } else {
            panic!("Unexpected genotype format");
        }
    }
}
