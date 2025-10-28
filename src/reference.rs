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
    contigs: Arc<Vec<ReferenceContig>>,
    alias_to_index: Arc<HashMap<String, usize>>,
    cache: Arc<Mutex<LruCache<(String, u64), u8>>>,
}

impl Clone for ReferenceGenome {
    fn clone(&self) -> Self {
        Self {
            path: self.path.clone(),
            reader: Arc::clone(&self.reader),
            contigs: Arc::clone(&self.contigs),
            alias_to_index: Arc::clone(&self.alias_to_index),
            cache: Arc::clone(&self.cache),
        }
    }
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
            contigs: Arc::new(contigs),
            alias_to_index: Arc::new(alias_to_index),
            cache: Arc::new(Mutex::new(LruCache::new(cache_capacity))),
        })
    }

    pub fn path(&self) -> &Path {
        self.path.as_path()
    }

    pub fn contigs(&self) -> &[ReferenceContig] {
        self.contigs.as_slice()
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

        let cache_key = (contig.name.clone(), position);
        if let Some(base) = self.cache.lock().get(&cache_key).copied() {
            return Ok(char::from(base));
        }

        let start = Position::try_from(pos)?;
        let region = Region::new(contig.name.clone(), start..=start);
        let mut reader = self.reader.lock();
        let record = reader.query(&region)?;
        let seq = record.sequence();
        let base = seq
            .as_ref()
            .first()
            .copied()
            .unwrap_or(b'N')
            .to_ascii_uppercase();
        self.cache.lock().put(cache_key, base);
        Ok(char::from(base))
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
    }
}
