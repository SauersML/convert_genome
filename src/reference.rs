use std::{
    collections::HashMap,
    fs,
    io::{self, BufReader, Write},
    path::{Path, PathBuf},
};

use curl::easy::{Easy, WriteError};
use flate2::bufread::MultiGzDecoder;
use noodles::{
    core::{Position, Region},
    fasta::{self, fai},
};
use std::str::Utf8Error;
use tempfile::TempDir;
use thiserror::Error;
use url::Url;
use zip::ZipArchive;

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
    _temp_dirs: Vec<TempDir>,
}

#[derive(Debug, Error)]
pub enum ReferenceError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
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
    #[error("reference resource not found: {0}")]
    NotFound(String),
    #[error("failed to parse URL {url}: {source}")]
    UrlParse {
        url: String,
        source: url::ParseError,
    },
    #[error("failed to download {url}: {source}")]
    Download { url: String, source: curl::Error },
    #[error("unsupported archive contents in {0}")]
    UnsupportedArchive(String),
    #[error("archive error: {0}")]
    Archive(#[from] zip::result::ZipError),
}

impl ReferenceGenome {
    pub fn open<P: AsRef<Path>>(
        path: P,
        fai_path: Option<PathBuf>,
    ) -> Result<Self, ReferenceError> {
        let prepared = prepare_reference(path.as_ref(), fai_path)?;
        let canonical = fs::canonicalize(&prepared.fasta_path)?;

        let index_path = prepared
            .fai_path
            .unwrap_or_else(|| default_index_path(&canonical));
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
            _temp_dirs: prepared.temp_dirs,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DownloadKind {
    Plain,
    Gzip,
    Zip,
}

#[derive(Debug)]
struct PreparedReference {
    fasta_path: PathBuf,
    fai_path: Option<PathBuf>,
    temp_dirs: Vec<TempDir>,
}

fn prepare_reference(
    path: &Path,
    fai_path: Option<PathBuf>,
) -> Result<PreparedReference, ReferenceError> {
    if path.exists() {
        let canonical = fs::canonicalize(path)?;
        let fai = match fai_path {
            Some(fai) => Some(canonicalize_or_error(&fai)?),
            None => None,
        };
        return Ok(PreparedReference {
            fasta_path: canonical,
            fai_path: fai,
            temp_dirs: Vec::new(),
        });
    }

    let path_str = path
        .to_str()
        .ok_or_else(|| ReferenceError::NotFound(format!("{}", path.display())))?;

    if !looks_like_url(path_str) {
        return Err(ReferenceError::NotFound(path_str.to_string()));
    }

    let url = Url::parse(path_str).map_err(|source| ReferenceError::UrlParse {
        url: path_str.to_string(),
        source,
    })?;

    let mut prepared = download_reference(&url)?;

    if let Some(fai) = fai_path {
        let aux = prepare_auxiliary_resource(&fai)?;
        prepared.fai_path = Some(canonicalize_or_error(&aux.path)?);
        prepared.temp_dirs.extend(aux.temp_dirs);
    } else if prepared.fai_path.is_some() {
        // already handled during download
    }

    Ok(prepared)
}

#[derive(Debug)]
struct PreparedAuxiliary {
    path: PathBuf,
    temp_dirs: Vec<TempDir>,
}

fn prepare_auxiliary_resource(path: &Path) -> Result<PreparedAuxiliary, ReferenceError> {
    if path.exists() {
        return Ok(PreparedAuxiliary {
            path: canonicalize_or_error(path)?,
            temp_dirs: Vec::new(),
        });
    }

    let path_str = path
        .to_str()
        .ok_or_else(|| ReferenceError::NotFound(format!("{}", path.display())))?;

    if !looks_like_url(path_str) {
        return Err(ReferenceError::NotFound(path_str.to_string()));
    }

    let url = Url::parse(path_str).map_err(|source| ReferenceError::UrlParse {
        url: path_str.to_string(),
        source,
    })?;

    let temp_dir = TempDir::new()?;
    let file_name = infer_filename(&url).unwrap_or_else(|| "reference.fai".to_string());
    let target_path = temp_dir.path().join(file_name);
    download_url_to_file(&url, &target_path)?;

    Ok(PreparedAuxiliary {
        path: target_path,
        temp_dirs: vec![temp_dir],
    })
}

fn canonicalize_or_error(path: &Path) -> Result<PathBuf, ReferenceError> {
    fs::canonicalize(path).map_err(ReferenceError::from)
}

fn download_reference(url: &Url) -> Result<PreparedReference, ReferenceError> {
    let temp_dir = TempDir::new()?;
    let file_name = infer_filename(url).unwrap_or_else(|| "reference.fa".to_string());
    let download_path = temp_dir.path().join(&file_name);
    download_url_to_file(url, &download_path)?;

    let kind = infer_download_kind(url);
    let (fasta_path, fai_path) = finalize_downloaded_file(temp_dir.path(), &download_path, kind)?;

    Ok(PreparedReference {
        fasta_path,
        fai_path,
        temp_dirs: vec![temp_dir],
    })
}

fn download_url_to_file(url: &Url, target: &Path) -> Result<(), ReferenceError> {
    let mut easy = Easy::new();
    easy.url(url.as_str())
        .map_err(|source| ReferenceError::Download {
            url: url.to_string(),
            source,
        })?;
    easy.follow_location(true)
        .map_err(|source| ReferenceError::Download {
            url: url.to_string(),
            source,
        })?;

    let mut file = fs::File::create(target)?;
    let mut writer_error: Option<io::Error> = None;
    {
        let mut transfer = easy.transfer();
        transfer
            .write_function(|data| {
                if let Err(err) = file.write_all(data) {
                    writer_error = Some(err);
                    return Err(WriteError::Pause);
                }
                Ok(data.len())
            })
            .map_err(|source| ReferenceError::Download {
                url: url.to_string(),
                source,
            })?;

        transfer
            .perform()
            .map_err(|source| ReferenceError::Download {
                url: url.to_string(),
                source,
            })?;
    }

    if let Some(err) = writer_error {
        return Err(ReferenceError::from(err));
    }

    Ok(())
}

fn infer_filename(url: &Url) -> Option<String> {
    url.path_segments()
        .and_then(|segments| segments.last())
        .map(|name| name.to_string())
        .filter(|name| !name.is_empty())
}

fn looks_like_url(input: &str) -> bool {
    input.contains("://")
}

fn infer_download_kind(url: &Url) -> DownloadKind {
    let path = url.path().to_ascii_lowercase();
    if path.ends_with(".zip") {
        DownloadKind::Zip
    } else if path.ends_with(".gz") {
        DownloadKind::Gzip
    } else {
        DownloadKind::Plain
    }
}

fn finalize_downloaded_file(
    base_dir: &Path,
    download_path: &Path,
    kind: DownloadKind,
) -> Result<(PathBuf, Option<PathBuf>), ReferenceError> {
    match kind {
        DownloadKind::Plain => {
            let canonical = fs::canonicalize(download_path)?;
            Ok((canonical, None))
        }
        DownloadKind::Gzip => finalize_gzip_download(base_dir, download_path),
        DownloadKind::Zip => finalize_zip_download(base_dir, download_path),
    }
}

fn finalize_gzip_download(
    base_dir: &Path,
    download_path: &Path,
) -> Result<(PathBuf, Option<PathBuf>), ReferenceError> {
    let file_name = download_path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or_else(|| ReferenceError::UnsupportedArchive(download_path.display().to_string()))?;
    let target_path = base_dir.join(file_name);
    let reader = fs::File::open(download_path)?;
    let mut decoder = MultiGzDecoder::new(BufReader::new(reader));
    let mut output = fs::File::create(&target_path)?;
    io::copy(&mut decoder, &mut output)?;
    let canonical = fs::canonicalize(&target_path)?;
    Ok((canonical, None))
}

fn finalize_zip_download(
    base_dir: &Path,
    download_path: &Path,
) -> Result<(PathBuf, Option<PathBuf>), ReferenceError> {
    let file = fs::File::open(download_path)?;
    let mut archive = ZipArchive::new(file)?;
    let mut fasta_path: Option<PathBuf> = None;
    let mut fai_path: Option<PathBuf> = None;

    for i in 0..archive.len() {
        let mut file = archive.by_index(i)?;
        let Some(enclosed_name) = file.enclosed_name().map(|p| p.to_owned()) else {
            continue;
        };
        let out_path = base_dir.join(&enclosed_name);

        if file.is_dir() {
            fs::create_dir_all(&out_path)?;
            continue;
        }

        if let Some(parent) = out_path.parent() {
            if !parent.exists() {
                fs::create_dir_all(parent)?;
            }
        }

        let mut outfile = fs::File::create(&out_path)?;
        io::copy(&mut file, &mut outfile)?;

        if fasta_path.is_none() && is_fasta_path(&out_path) {
            fasta_path = Some(fs::canonicalize(&out_path)?);
        }

        if fai_path.is_none() && out_path.extension().is_some_and(|ext| ext == "fai") {
            fai_path = Some(fs::canonicalize(&out_path)?);
        }
    }

    let fasta = fasta_path
        .ok_or_else(|| ReferenceError::UnsupportedArchive(download_path.display().to_string()))?;

    Ok((fasta, fai_path))
}

fn is_fasta_path(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| matches!(ext.to_ascii_lowercase().as_str(), "fa" | "fasta" | "fna"))
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;

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
    fn detects_known_remote_kinds() {
        let urls = [
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz",
            "https://webdata.illumina.com/downloads/productfiles/microarray-analytics-array/GRCh38_genome.zip",
            "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
            "https://webdata.illumina.com/downloads/productfiles/microarray-analytics-array/GRCh37_genome.zip",
        ];

        for url in urls {
            let parsed = Url::parse(url).unwrap();
            let kind = infer_download_kind(&parsed);
            if url.ends_with(".zip") {
                assert_eq!(kind, DownloadKind::Zip);
            } else {
                assert_eq!(kind, DownloadKind::Gzip);
            }
        }
    }

    #[test]
    fn finalize_gzip_creates_plain_fasta() {
        let temp = TempDir::new().unwrap();
        let gz_path = temp.path().join("test.fa.gz");
        let mut encoder =
            GzEncoder::new(fs::File::create(&gz_path).unwrap(), Compression::default());
        writeln!(encoder, ">chr1").unwrap();
        writeln!(encoder, "ACGT").unwrap();
        encoder.finish().unwrap();

        let (fasta_path, fai_path) =
            finalize_downloaded_file(temp.path(), &gz_path, DownloadKind::Gzip).unwrap();
        assert!(fai_path.is_none());
        let contents = fs::read_to_string(&fasta_path).unwrap();
        assert!(contents.contains("chr1"));
    }

    #[test]
    fn finalize_zip_extracts_fasta_and_index() {
        let temp = TempDir::new().unwrap();
        let zip_path = temp.path().join("archive.zip");
        {
            let file = fs::File::create(&zip_path).unwrap();
            let mut zip = zip::ZipWriter::new(file);
            let options = zip::write::FileOptions::default();
            zip.add_directory("genome", options).unwrap();
            zip.start_file("genome/ref.fa", options).unwrap();
            writeln!(zip, ">chr1").unwrap();
            writeln!(zip, "ACGT").unwrap();
            zip.start_file("genome/ref.fa.fai", options).unwrap();
            writeln!(zip, "ref\t4\t6\t4\t5").unwrap();
            zip.finish().unwrap();
        }

        let (fasta_path, fai_path) =
            finalize_downloaded_file(temp.path(), &zip_path, DownloadKind::Zip).unwrap();
        assert!(fasta_path.ends_with("ref.fa"));
        let fai_path = fai_path.expect("expected fai");
        assert!(fai_path.ends_with("ref.fa.fai"));
    }
}
