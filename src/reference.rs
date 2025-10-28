use std::{
    collections::HashMap,
    fs, io,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    sync::Arc,
};

use lru::LruCache;
use noodles::{
    core::{Position, Region},
    fasta::{self, fai},
};
use parking_lot::Mutex;
use std::str::Utf8Error;
use thiserror::Error;

#[derive(Debug, Clone)]
pub struct ReferenceContig {
    pub name: String,
    pub length: u64,
}

pub struct ReferenceGenome {
    path: PathBuf,
    reader: Arc<Mutex<fasta::io::IndexedReader<fasta::io::BufReader<fs::File>>>>,
    contigs: Vec<ReferenceContig>,
    alias_to_index: HashMap<String, usize>,
    cache: Arc<Mutex<LruCache<(String, u64), u8>>>,
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
}

impl ReferenceGenome {
    pub fn open<P: AsRef<Path>>(
        path: P,
        fai_path: Option<PathBuf>,
    ) -> Result<Self, ReferenceError> {
        let path = path.as_ref();
        let canonical = fs::canonicalize(path)?;

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

        let cache_capacity = NonZeroUsize::new(128 * 1024).expect("non-zero cache capacity");

        Ok(Self {
            path: canonical,
            reader: Arc::new(Mutex::new(reader)),
            contigs,
            alias_to_index,
            cache: Arc::new(Mutex::new(LruCache::new(cache_capacity))),
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

    #[doc(hidden)]
    pub fn cache_len(&self) -> usize {
        self.cache.lock().len()
    }

    pub fn base(&self, query: &str, position: u64) -> Result<char, ReferenceError> {
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
        let key = (contig.name.clone(), position);

        if let Some(&cached) = self.cache.lock().get(&key) {
            return Ok(char::from(cached).to_ascii_uppercase());
        }

        let base = {
            let mut reader = self.reader.lock();
            let record = reader.query(&region)?;
            let seq = record.sequence();
            seq.as_ref().first().copied().unwrap_or(b'N')
        };

        self.cache.lock().put(key, base);

        Ok(char::from(base).to_ascii_uppercase())
    }
}

impl Clone for ReferenceGenome {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            reader: Arc::clone(&self.reader),
            contigs: self.contigs.clone(),
            alias_to_index: self.alias_to_index.clone(),
            cache: Arc::clone(&self.cache),
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn reference_fetches_base() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("ref.fa");
        let mut file = fs::File::create(&fasta_path).unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();
        drop(file);

        let reference = ReferenceGenome::open(&fasta_path, None).unwrap();
        assert_eq!(reference.base("1", 2).unwrap(), 'C');
        assert!(reference.base("1", 0).is_err());

        // Subsequent lookups are served from the cache.
        assert_eq!(reference.base("1", 2).unwrap(), 'C');
        assert!(reference.cache_len() >= 1);
    }
}
