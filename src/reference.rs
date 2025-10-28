use std::{
    collections::HashMap,
    ffi::OsStr,
    fs, io,
    path::{Path, PathBuf},
};

use noodles::{
    core::{Position, Region},
    fasta::{self, fai},
};
use reqwest::StatusCode;
use reqwest::blocking::Client;
use std::str::Utf8Error;
use tempfile::TempDir;
use thiserror::Error;
use zip::ZipArchive;

use flate2::read::MultiGzDecoder;

#[derive(Debug, Clone)]
pub struct ReferenceContig {
    pub name: String,
    pub length: u64,
}

pub struct ReferenceGenome {
    path: PathBuf,
    reader: fasta::io::IndexedReader<fasta::io::BufReader<fs::File>>,
    contigs: Vec<ReferenceContig>,
    alias_to_index: HashMap<String, usize>,
    #[allow(dead_code)]
    temp_dir: Option<TempDir>,
}

#[derive(Debug, Error)]
pub enum ReferenceError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("failed to download reference: {0}")]
    Download(#[from] reqwest::Error),
    #[error("remote server returned status {status}")]
    HttpStatus { status: StatusCode },
    #[error("invalid UTF-8 contig name: {0}")]
    InvalidContigName(#[from] Utf8Error),
    #[error("unknown contig: {query}")]
    UnknownContig { query: String },
    #[error("position {position} is outside contig {contig} length {length}")]
    PositionOutOfBounds {
        contig: String,
        position: u64,
        length: u64,
    },
    #[error("invalid genomic position: {0}")]
    InvalidPosition(#[from] noodles::core::position::TryFromIntError),
    #[error("invalid reference URL: {0}")]
    InvalidUrl(String),
    #[error("archive extraction error: {0}")]
    Archive(#[from] zip::result::ZipError),
    #[error("no FASTA file found inside archive")]
    MissingFastaInArchive,
}

impl ReferenceGenome {
    pub fn open<P: AsRef<Path>>(
        path: P,
        fai_path: Option<PathBuf>,
    ) -> Result<Self, ReferenceError> {
        let path = path.as_ref();
        let (canonical, temp_dir) = if is_remote_path(path) {
            let url = path
                .to_str()
                .ok_or_else(|| ReferenceError::InvalidUrl(path.to_string_lossy().into_owned()))?;
            let downloaded = fetch_remote_reference(url)?;
            let canonical = fs::canonicalize(&downloaded.fasta_path)?;
            (canonical, Some(downloaded.temp_dir))
        } else {
            (fs::canonicalize(path)?, None)
        };

        let index_path = fai_path.unwrap_or_else(|| default_index_path(&canonical));
        let index = if index_path.exists() {
            fai::fs::read(&index_path)?
        } else {
            let index = fasta::fs::index(&canonical)?;
            fai::fs::write(&index_path, &index)?;
            index
        };

        let reader = fasta::io::indexed_reader::Builder::default()
            .set_index(index.clone())
            .build_from_path(&canonical)?;

        let contigs = index
            .as_ref()
            .iter()
            .map(|record| -> Result<ReferenceContig, ReferenceError> {
                let name = std::str::from_utf8(record.name().as_ref())?.to_string();
                Ok(ReferenceContig {
                    name,
                    length: record.length(),
                })
            })
            .collect::<Result<Vec<_>, _>>()?;

        let alias_to_index = build_alias_map(&contigs);

        Ok(Self {
            path: canonical,
            reader,
            contigs,
            alias_to_index,
            temp_dir,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn contigs(&self) -> &[ReferenceContig] {
        &self.contigs
    }

    pub fn resolve_contig_name(&self, query: &str) -> Option<&str> {
        let key = canonical_key(query);
        self.alias_to_index
            .get(&key)
            .map(|&idx| self.contigs[idx].name.as_str())
    }

    pub fn contig(&self, query: &str) -> Option<&ReferenceContig> {
        let key = canonical_key(query);
        self.alias_to_index.get(&key).map(|&idx| &self.contigs[idx])
    }

    pub fn base(&mut self, query: &str, position: u64) -> Result<char, ReferenceError> {
        let contig = self
            .contig(query)
            .ok_or_else(|| ReferenceError::UnknownContig {
                query: query.to_string(),
            })?;

        if position == 0 || position > contig.length {
            return Err(ReferenceError::PositionOutOfBounds {
                contig: contig.name.clone(),
                position,
                length: contig.length,
            });
        }

        let pos = usize::try_from(position).map_err(|_| ReferenceError::PositionOutOfBounds {
            contig: contig.name.clone(),
            position,
            length: contig.length,
        })?;

        let start = Position::try_from(pos)?;
        let region = Region::new(contig.name.clone(), start..=start);
        let record = self.reader.query(&region)?;
        let seq = record.sequence();
        let base = seq.as_ref().first().copied().unwrap_or(b'N');
        Ok(char::from(base).to_ascii_uppercase())
    }
}

fn default_index_path(path: &Path) -> PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(".fai");
    PathBuf::from(s)
}

fn canonical_key(raw: &str) -> String {
    let trimmed = raw.trim();
    let trimmed = trimmed.strip_prefix("chr").unwrap_or(trimmed);
    let upper = trimmed.to_ascii_uppercase();
    match upper.as_str() {
        "M" => "MT".to_string(),
        _ => upper,
    }
}

fn build_alias_map(contigs: &[ReferenceContig]) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for (idx, contig) in contigs.iter().enumerate() {
        let name = contig.name.as_str();
        let canonical = canonical_key(name);
        map.entry(canonical.clone()).or_insert(idx);
        map.entry(name.to_ascii_uppercase()).or_insert(idx);
        if let Some(stripped) = name.strip_prefix("chr") {
            map.entry(canonical_key(stripped)).or_insert(idx);
        }
        if name.eq_ignore_ascii_case("chrM") || name.eq_ignore_ascii_case("MT") {
            map.entry("MT".into()).or_insert(idx);
            map.entry("M".into()).or_insert(idx);
        }
    }
    map
}

struct DownloadedReference {
    temp_dir: TempDir,
    fasta_path: PathBuf,
}

fn fetch_remote_reference(url: &str) -> Result<DownloadedReference, ReferenceError> {
    let parsed =
        reqwest::Url::parse(url).map_err(|err| ReferenceError::InvalidUrl(err.to_string()))?;
    let filename = parsed
        .path_segments()
        .and_then(|segments| segments.last())
        .filter(|segment| !segment.is_empty())
        .unwrap_or("reference.fa");

    let client = Client::builder().build()?;
    let mut response = client.get(parsed.clone()).send()?;
    if !response.status().is_success() {
        return Err(ReferenceError::HttpStatus {
            status: response.status(),
        });
    }

    let temp_dir = tempfile::tempdir()?;
    let download_path = temp_dir.path().join(filename);
    {
        let mut file = fs::File::create(&download_path)?;
        io::copy(&mut response, &mut file)?;
    }

    let fasta_path = prepare_downloaded_file(&download_path, temp_dir.path())?;

    Ok(DownloadedReference {
        temp_dir,
        fasta_path,
    })
}

fn prepare_downloaded_file(
    download_path: &Path,
    temp_dir: &Path,
) -> Result<PathBuf, ReferenceError> {
    match download_path
        .extension()
        .and_then(OsStr::to_str)
        .map(|s| s.to_ascii_lowercase())
    {
        Some(ext) if ext == "gz" => decompress_gzip(download_path, temp_dir),
        Some(ext) if ext == "zip" => decompress_zip(download_path, temp_dir),
        _ => Ok(download_path.to_path_buf()),
    }
}

fn decompress_gzip(path: &Path, temp_dir: &Path) -> Result<PathBuf, ReferenceError> {
    let stem = path
        .file_stem()
        .and_then(OsStr::to_str)
        .map(|s| s.to_string())
        .unwrap_or_else(|| String::from("reference"));
    let output_path = temp_dir.join(stem);
    let file = fs::File::open(path)?;
    let mut decoder = MultiGzDecoder::new(file);
    let mut output = fs::File::create(&output_path)?;
    io::copy(&mut decoder, &mut output)?;
    Ok(output_path)
}

fn decompress_zip(path: &Path, temp_dir: &Path) -> Result<PathBuf, ReferenceError> {
    let file = fs::File::open(path)?;
    let mut archive = ZipArchive::new(file)?;

    for idx in 0..archive.len() {
        let mut entry = archive.by_index(idx)?;
        if entry.is_dir() {
            continue;
        }

        let entry_path = Path::new(entry.name());
        let extension = entry_path
            .extension()
            .and_then(OsStr::to_str)
            .map(|s| s.to_ascii_lowercase());

        if !matches_fasta_extension(extension.as_deref()) {
            continue;
        }

        let filename = entry_path
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or("reference.fa");
        let mut output_path = temp_dir.join(filename);
        {
            let mut output = fs::File::create(&output_path)?;
            io::copy(&mut entry, &mut output)?;
        }

        if extension.as_deref() == Some("gz") {
            output_path = decompress_gzip(&output_path, temp_dir)?;
        }

        return Ok(output_path);
    }

    Err(ReferenceError::MissingFastaInArchive)
}

fn matches_fasta_extension(ext: Option<&str>) -> bool {
    match ext {
        Some("fa") | Some("fasta") | Some("fna") => true,
        Some("gz") => true,
        _ => false,
    }
}

fn is_remote_path(path: &Path) -> bool {
    path.to_str()
        .map(|s| {
            let lower = s.to_ascii_lowercase();
            lower.starts_with("http://") || lower.starts_with("https://")
        })
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;
    use zip::write::FileOptions;

    #[test]
    fn reference_fetches_base() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        let mut file = fs::File::create(&fasta_path).unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();
        drop(file);

        let mut reference = ReferenceGenome::open(&fasta_path, None).unwrap();
        assert_eq!(reference.base("1", 2).unwrap(), 'C');
        assert!(reference.base("1", 0).is_err());
    }

    #[test]
    fn decompresses_gzip_reference() {
        let dir = tempfile::tempdir().unwrap();
        let gz_path = dir.path().join("ref.fa.gz");
        {
            let file = fs::File::create(&gz_path).unwrap();
            let mut encoder = GzEncoder::new(file, Compression::default());
            writeln!(encoder, ">chr1").unwrap();
            writeln!(encoder, "ACGT").unwrap();
            encoder.finish().unwrap();
        }

        let extracted = decompress_gzip(&gz_path, dir.path()).unwrap();
        let contents = fs::read_to_string(&extracted).unwrap();
        assert!(contents.contains("chr1"));
    }

    #[test]
    fn decompresses_zip_reference() {
        let dir = tempfile::tempdir().unwrap();
        let zip_path = dir.path().join("ref.zip");
        {
            let file = fs::File::create(&zip_path).unwrap();
            let mut zip = zip::ZipWriter::new(file);
            let options = FileOptions::default();
            zip.start_file("nested/chr1.fa", options).unwrap();
            writeln!(zip, ">chr1").unwrap();
            writeln!(zip, "ACGT").unwrap();
            zip.finish().unwrap();
        }

        let extracted = decompress_zip(&zip_path, dir.path()).unwrap();
        let contents = fs::read_to_string(&extracted).unwrap();
        assert!(contents.contains("ACGT"));
    }
}
