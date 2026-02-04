use anyhow::{Context, Result, anyhow, bail};
use flate2::read::GzDecoder;
use url::Url;

use crate::dtc::{self, Record as DtcRecord};
use crate::reference::ReferenceGenome;
use crate::remote::{self};
use std::fs;
use std::fs::OpenOptions;
use std::io;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Strand orientation mode for the input file.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StrandMode {
    #[default]
    Auto, // Infer from data (default)
    Forward, // Force forward strand (no flipping)
    Reverse, // Force reverse strand (flip everything)
}

fn xdg_cache_dir() -> Result<PathBuf> {
    if let Ok(dir) = std::env::var("XDG_CACHE_HOME") {
        return Ok(PathBuf::from(dir));
    }
    let home = std::env::var("HOME").context("HOME not set; cannot determine cache directory")?;
    Ok(PathBuf::from(home).join(".cache"))
}

fn check_build_cache_dir() -> Result<PathBuf> {
    Ok(xdg_cache_dir()?.join("check_build"))
}

fn candidate_refs_dirs() -> Vec<PathBuf> {
    let mut candidates = Vec::new();

    if let Some(dirs) = directories::ProjectDirs::from("com", "convert_genome", "convert_genome") {
        candidates.push(dirs.cache_dir().join("refs"));
    }

    if let Ok(dir) = std::env::var("XDG_CACHE_HOME") {
        candidates.push(PathBuf::from(dir).join("convert_genome").join("refs"));
    }

    if let Ok(home) = std::env::var("HOME") {
        candidates.push(
            PathBuf::from(home)
                .join(".cache")
                .join("convert_genome")
                .join("refs"),
        );
    }

    candidates.push(std::env::temp_dir().join("convert_genome").join("refs"));

    let mut deduped = Vec::new();
    for candidate in candidates {
        if !deduped.contains(&candidate) {
            deduped.push(candidate);
        }
    }

    deduped
}

fn ensure_writable_dir(path: &Path) -> Result<()> {
    fs::create_dir_all(path)
        .with_context(|| format!("failed to create {}", path.display()))?;

    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let test_path = path.join(format!(".writetest.{}.{}", std::process::id(), nanos));

    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&test_path)
        .with_context(|| format!("failed to create {}", test_path.display()))?;
    file.write_all(b"")?;
    file.flush()?;
    fs::remove_file(&test_path)
        .with_context(|| format!("failed to clean up {}", test_path.display()))?;

    Ok(())
}

pub fn convert_genome_refs_dir() -> Result<PathBuf> {
    let candidates = candidate_refs_dirs();
    let mut errors = Vec::new();

    for candidate in candidates {
        match ensure_writable_dir(&candidate) {
            Ok(()) => return Ok(candidate),
            Err(err) => errors.push(format!("{} ({})", candidate.display(), err)),
        }
    }

    Err(anyhow!(
        "Failed to find a writable cache directory. Tried: {}",
        errors.join("; ")
    ))
}

fn build_to_check_build_filename(build: &str) -> Option<&'static str> {
    let normalized = build.to_lowercase();
    if normalized.contains("37") || normalized.contains("hg19") {
        Some("hg19.fa.gz")
    } else if normalized.contains("38") || normalized.contains("hg38") {
        Some("hg38.fa.gz")
    } else {
        None
    }
}

fn build_to_uncompressed_name(build: &str) -> Option<&'static str> {
    let normalized = build.to_lowercase();
    if normalized.contains("37") || normalized.contains("hg19") {
        Some("hg19.fa")
    } else if normalized.contains("38") || normalized.contains("hg38") {
        Some("hg38.fa")
    } else {
        None
    }
}

fn decompress_gzip_to_path(gz_path: &Path, out_path: &Path) -> Result<()> {
    let input = fs::File::open(gz_path)
        .with_context(|| format!("failed to open gz reference {}", gz_path.display()))?;
    let mut decoder = GzDecoder::new(input);

    let parent = out_path
        .parent()
        .ok_or_else(|| anyhow!("Invalid output path"))?;
    fs::create_dir_all(parent)?;

    // Create a unique temporary file in the same directory to avoid race conditions
    // and ensure atomic rename.
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let suffix = format!(".tmp.{}.{}", std::process::id(), nanos);

    let file_name = out_path
        .file_name()
        .ok_or_else(|| anyhow!("Invalid filename"))?;
    let mut tmp_name = file_name.to_os_string();
    tmp_name.push(suffix);
    let tmp_path = parent.join(tmp_name);

    let output = fs::File::create(&tmp_path)
        .with_context(|| format!("failed to create output {}", tmp_path.display()))?;
    let mut writer = io::BufWriter::new(output);
    io::copy(&mut decoder, &mut writer)
        .with_context(|| format!("failed to decompress {}", gz_path.display()))?;
    writer.flush()?;
    // Ensure data is on disk before rename
    writer.get_ref().sync_all()?;

    fs::rename(&tmp_path, out_path)
        .with_context(|| format!("failed to finalize {}", out_path.display()))?;
    Ok(())
}

/// Inferred strand state of the file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InferredStrand {
    Forward,
    Reverse,
}

/// Registry of known source reference genomes.
pub struct SourceReferenceRegistry;

impl SourceReferenceRegistry {
    /// Get the download URL for a specific genome build.
    pub fn get_url(build: &str) -> Option<Url> {
        let normalized = build.to_lowercase();
        // URLs for basic FASTA files (no alts/patches for speed/simplicity)
        if normalized.contains("37") || normalized.contains("hg19") {
            // Using a reliable source for hg19 fasta (e.g., from UCSC or similar)
            // Use a known public UCSC URL.
            // Using GRC assembly from NCBI or similar is best.
            // Let's use the one from the project's S3 or similar if available, otherwise a public one.
            // User didn't specify, so I'll suggest a standard UCSC hg19 URL.
            // Warning: These files are large (900MB+ gzipped).
            // Using a smaller masked version or chrom-by-chrom might be better, but we need random access.
            // Let's assume the user has a way to get these or we use a standard one.
            // Since I can't browse the web for a URL, I will use a standard easy-to-guess one or rely on the user providing it?
            // "Source reference loading" implies we fetch it.
            // Let's use the 1000 Genomes hs37d5 or standard hg19.
            Some(
                Url::parse("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz")
                    .unwrap(),
            )
        } else if normalized.contains("38") || normalized.contains("hg38") {
            Some(
                Url::parse("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")
                    .unwrap(),
            )
        } else {
            None
        }
    }
}

/// Download and open a source reference genome.
pub fn load_source_reference(build: &str) -> Result<ReferenceGenome> {
    let refs_dir = convert_genome_refs_dir()?;
    fs::create_dir_all(&refs_dir)?;

    let uncompressed_name = build_to_uncompressed_name(build).ok_or_else(|| {
        anyhow!(
            "Unknown build '{}', cannot map to reference filename",
            build
        )
    })?;
    let reference_path = refs_dir.join(uncompressed_name);

    // If we already have the .fa cached, use it (generate .fai if needed)
    let fai_path = reference_path.with_extension("fa.fai");
    if reference_path.exists() {
        if fai_path.exists() {
            // Both exist, use directly
            return ReferenceGenome::open(&reference_path, Some(fai_path)).with_context(|| {
                format!(
                    "Failed to open cached source reference at {}",
                    reference_path.display()
                )
            });
        } else {
            // .fa exists but .fai missing - generate it
            ReferenceGenome::open(&reference_path, None).with_context(|| {
                format!(
                    "Failed to open cached source reference at {} (generating index)",
                    reference_path.display()
                )
            })?;
            // Now use the generated index
            return ReferenceGenome::open(&reference_path, Some(fai_path)).with_context(|| {
                format!(
                    "Failed to open cached source reference at {}",
                    reference_path.display()
                )
            });
        }
    }
    // Try to reuse check_build's cached .fa.gz.
    if let Some(check_build_filename) = build_to_check_build_filename(build) {
        let check_build_gz = check_build_cache_dir()?.join(check_build_filename);
        if check_build_gz.exists() {
            tracing::info!(
                gz = %check_build_gz.display(),
                out = %reference_path.display(),
                "expanding check_build cached reference"
            );
            decompress_gzip_to_path(&check_build_gz, &reference_path)?;

            // Open once to generate .fai index, which noodles will cache next to the .fa file
            ReferenceGenome::open(&reference_path, None).with_context(|| {
                format!(
                    "Failed to open source reference at {} (index generation)",
                    reference_path.display()
                )
            })?;

            // Now return with the cached index
            let fai_path = reference_path.with_extension("fa.fai");
            return ReferenceGenome::open(&reference_path, Some(fai_path)).with_context(|| {
                format!(
                    "Failed to open source reference at {}",
                    reference_path.display()
                )
            });
        }
    }

    // Fallback: download ourselves if check_build cache is missing.
    let url = SourceReferenceRegistry::get_url(build).ok_or_else(|| {
        anyhow!(
            "Unknown build '{}', cannot download source reference",
            build
        )
    })?;

    tracing::info!("Downloading source reference from {}", url);
    let resource = remote::fetch_remote_resource(&url)?;

    // Copy to unique temp file then rename to avoid race conditions
    let parent = reference_path.parent().unwrap(); // We know this exists from create_dir_all above
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let suffix = format!(".tmp.{}.{}", std::process::id(), nanos);

    let file_name = reference_path
        .file_name()
        .ok_or_else(|| anyhow!("Invalid filename"))?;
    let mut tmp_name = file_name.to_os_string();
    tmp_name.push(suffix);
    let tmp_path = parent.join(tmp_name);

    fs::copy(resource.local_path(), &tmp_path)?;
    fs::rename(&tmp_path, &reference_path)?;

    // Open once to generate .fai index, which noodles will cache next to the .fa file
    ReferenceGenome::open(&reference_path, None).with_context(|| {
        format!(
            "Failed to open source reference at {} (index generation)",
            reference_path.display()
        )
    })?;

    // Now return with the cached index
    let fai_path = reference_path.with_extension("fa.fai");
    ReferenceGenome::open(&reference_path, Some(fai_path)).with_context(|| {
        format!(
            "Failed to open source reference at {}",
            reference_path.display()
        )
    })
}

/// Infer the strand of the input file by comparing against source reference.
pub fn infer_strand_lock(
    records: &[DtcRecord],
    reference: &ReferenceGenome,
) -> Result<InferredStrand> {
    let mut matching_plus = 0;
    let mut matching_minus = 0;
    let mut n_tested = 0;

    for rec in records {
        if n_tested >= 2000 {
            break;
        }

        // Skip ambiguous SNPs
        if is_transversion(rec) {
            // Get reference base
            let chrom = normalize_chrom(&rec.chromosome);
            // Try resolving contig
            let ref_base_res = reference.base(&chrom, rec.position);

            // If failed, try adding/removing chr
            let ref_base = match ref_base_res {
                Ok(b) => b,
                Err(_) => {
                    let alt = if chrom.starts_with("chr") {
                        chrom.trim_start_matches("chr").to_string()
                    } else {
                        format!("chr{}", chrom)
                    };
                    match reference.base(&alt, rec.position) {
                        Ok(b) => b,
                        Err(_) => continue, // Ref not found
                    }
                }
            };

            let ref_base_char = ref_base.to_ascii_uppercase();

            // Parse alleles
            let mut allele_chars = Vec::new();
            if let Ok(alleles) = rec.parse_alleles() {
                for a in alleles {
                    if let dtc::Allele::Base(s) = a {
                        allele_chars.push(s.chars().next().unwrap().to_ascii_uppercase());
                    }
                }
            }

            if allele_chars.is_empty() {
                continue;
            }

            // Check Plus: is ref_base in allele_chars?
            let is_plus = allele_chars.contains(&ref_base_char);

            // Check Minus: is complement(ref_base) in allele_chars?
            let comp_ref = complement(ref_base_char);
            let is_minus = allele_chars.contains(&comp_ref);

            if is_plus && !is_minus {
                matching_plus += 1;
                n_tested += 1;
            } else if is_minus && !is_plus {
                matching_minus += 1;
                n_tested += 1;
            }
        }
    }

    if n_tested < 50 {
        bail!(
            "Too few informative SNPs ({}) for strand inference. Refusing to guess strand.",
            n_tested
        );
    }

    let plus_frac = matching_plus as f64 / n_tested as f64;
    let minus_frac = matching_minus as f64 / n_tested as f64;

    tracing::info!(
        "Strand Inference: {} tested, {:.1}% match Plus, {:.1}% match Minus",
        n_tested,
        plus_frac * 100.0,
        minus_frac * 100.0
    );

    if plus_frac > 0.90 {
        Ok(InferredStrand::Forward)
    } else if minus_frac > 0.90 {
        Ok(InferredStrand::Reverse)
    } else {
        bail!(
            "Inconsistent strand orientation: {:.1}% Plus vs {:.1}% Minus. File may be mixed or corrupted.",
            plus_frac * 100.0,
            minus_frac * 100.0
        )
    }
}

/// Helper to check if record is a transversion (safe for strand inference)
fn is_transversion(rec: &DtcRecord) -> bool {
    // We want to avoid A/T and C/G which are ambiguous if flipped.
    // Also avoid Indels for this check.
    // Safe ones: A/C, A/G, C/A, C/T, G/A, G/T, T/C, T/G
    // Actually, just check if alleles are complements of each other?
    // If Alleles are {A, T} -> Bad. {C, G} -> Bad.

    // Just parse from string for speed
    let genotype = rec.genotype.to_ascii_uppercase();

    // Exclude non-SNP encodings commonly seen in DTC files.
    // We only want clear A/C/G/T alleles for strand inference.
    if genotype.contains('I')
        || genotype.contains('D')
        || genotype.contains('-')
        || genotype.contains('<')
        || genotype.contains('>')
        || genotype.contains('.')
    {
        return false;
    }

    // Quick check logic
    // If it contains A and T -> Bad
    // If it contains C and G -> Bad

    let has_a = genotype.contains('A');
    let has_t = genotype.contains('T');
    let has_c = genotype.contains('C');
    let has_g = genotype.contains('G');

    if has_a && has_t {
        return false;
    }
    if has_c && has_g {
        return false;
    }

    // Must have at least some bases, and avoid N/other ambiguity.
    (has_a || has_t || has_c || has_g) && !genotype.contains('N')
}

fn normalize_chrom(c: &str) -> String {
    if c == "23" {
        "X".to_string()
    } else if c == "24" {
        "Y".to_string()
    } else if c == "25" {
        "MT".to_string()
    }
    // PAR/MT conventions vary
    else {
        c.to_string()
    }
}

fn complement(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'N' => 'N',
        _ => c,
    }
}
