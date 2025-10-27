use std::{
    fs::File,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{Context, Result, anyhow};
use lru::LruCache;
use noodles::core::{Position, Region};
use noodles::fasta::{self as fasta, fai};
use parking_lot::Mutex;

pub struct ReferenceGenome {
    reader: Arc<Mutex<fasta::io::IndexedReader<fasta::io::BufReader<File>>>>,
    index: Arc<fai::Index>,
    cache: Arc<Mutex<LruCache<(String, u64), u8>>>,
    source: PathBuf,
}

impl Clone for ReferenceGenome {
    fn clone(&self) -> Self {
        Self {
            reader: Arc::clone(&self.reader),
            index: Arc::clone(&self.index),
            cache: Arc::clone(&self.cache),
            source: self.source.clone(),
        }
    }
}

impl ReferenceGenome {
    pub fn open<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref();
        ensure_index_exists(path)?;

        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("failed to open reference FASTA at {}", path.display()))?;

        let index = reader.index().clone();
        let reader = Arc::new(Mutex::new(reader));
        let cache_capacity = NonZeroUsize::new(131072).unwrap();

        Ok(Self {
            reader,
            index: Arc::new(index),
            cache: Arc::new(Mutex::new(LruCache::new(cache_capacity))),
            source: path.to_path_buf(),
        })
    }

    pub fn source(&self) -> &Path {
        &self.source
    }

    pub fn index(&self) -> &fai::Index {
        &self.index
    }

    pub fn contigs(&self) -> impl Iterator<Item = fai::Record> + '_ {
        self.index.as_ref().as_ref().iter().cloned()
    }

    pub fn base(&self, chrom: &str, position: u64) -> Result<char> {
        if position == 0 {
            return Err(anyhow!("reference positions must be 1-based"));
        }

        if let Some(base) = {
            let mut cache = self.cache.lock();
            cache.get(&(chrom.to_string(), position)).copied()
        } {
            return Ok(base as char);
        }

        let pos_usize = usize::try_from(position)
            .map_err(|_| anyhow!("reference position {position} exceeds supported range"))?;
        let pos = Position::try_from(pos_usize)
            .map_err(|e| anyhow!("invalid position {position}: {e}"))?;
        let region = Region::new(chrom, pos..=pos);

        let base = {
            let mut reader = self.reader.lock();
            let record = reader
                .query(&region)
                .with_context(|| format!("failed to query reference for {chrom}:{position}"))?;
            let sequence = record.sequence();
            sequence.as_ref().first().copied().ok_or_else(|| {
                anyhow!("reference returned empty sequence for {chrom}:{position}")
            })?
        };

        let upper = base.to_ascii_uppercase();
        {
            let mut cache = self.cache.lock();
            cache.put((chrom.to_string(), position), upper);
        }

        Ok(char::from(upper))
    }
}

fn ensure_index_exists(path: &Path) -> Result<()> {
    let index_path = build_index_path(path);
    if index_path.exists() {
        return Ok(());
    }

    let index = fasta::fs::index(path)
        .with_context(|| format!("failed to build FAI index for {}", path.display()))?;
    fasta::fai::fs::write(&index_path, &index)
        .with_context(|| format!("failed to write index {}", index_path.display()))?;
    Ok(())
}

fn build_index_path(path: &Path) -> PathBuf {
    let mut index_path = path.as_os_str().to_os_string();
    index_path.push(".fai");
    PathBuf::from(index_path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn base_lookup() -> Result<()> {
        let dir = tempdir()?;
        let fasta_path = dir.path().join("ref.fa");
        {
            let mut file = File::create(&fasta_path)?;
            writeln!(file, ">chr1")?;
            writeln!(file, "ACGT")?;
        }

        let reference = ReferenceGenome::open(&fasta_path)?;
        assert_eq!(reference.base("chr1", 1)?, 'A');
        assert_eq!(reference.base("chr1", 2)?, 'C');
        assert_eq!(reference.base("chr1", 4)?, 'T');
        Ok(())
    }
}
