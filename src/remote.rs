use std::{
    fs::File,
    io::{self, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, anyhow, bail};
use curl::easy::{Easy, WriteError};
use flate2::read::MultiGzDecoder;
use tempfile::TempDir;
use url::Url;
use zip::ZipArchive;

const FASTA_EXTENSIONS: &[&str] = &["fa", "fasta", "fna", "fa.gz", "fasta.gz", "fna.gz"];

pub struct RemoteResource {
    _temp_dir: TempDir,
    local_path: PathBuf,
}

impl RemoteResource {
    pub fn local_path(&self) -> &Path {
        &self.local_path
    }
}

pub fn fetch_remote_resource(url: &Url) -> Result<RemoteResource> {
    tracing::info!(target: "convert_genome", source = %url, "downloading remote resource");

    let temp_dir = TempDir::new().context("failed to create temporary directory for download")?;
    let filename = url
        .path_segments()
        .and_then(|segments| segments.last())
        .filter(|segment| !segment.is_empty())
        .map(|segment| segment.to_string())
        .unwrap_or_else(|| String::from("download"));

    let download_path = temp_dir.path().join(filename);
    download_to_path(url, &download_path).with_context(|| format!("failed to download {url}"))?;

    let prepared = prepare_downloaded_file(temp_dir.path(), &download_path)?;

    Ok(RemoteResource {
        _temp_dir: temp_dir,
        local_path: prepared,
    })
}

fn download_to_path(url: &Url, destination: &Path) -> Result<()> {
    let mut easy = Easy::new();
    easy.url(url.as_str())?;
    easy.follow_location(true)?;
    easy.accept_encoding("identity")?;
    easy.fail_on_error(true)?;

    let mut file = File::create(destination)
        .with_context(|| format!("failed to create {}", destination.display()))?;
    let mut write_error: Option<io::Error> = None;

    {
        let mut transfer = easy.transfer();
        transfer.write_function(|data| {
            if let Err(err) = file.write_all(data) {
                write_error = Some(err);
                return Err(WriteError::Pause);
            }
            Ok(data.len())
        })?;
        transfer.perform()?;
    }

    if let Some(err) = write_error {
        return Err(anyhow!(
            "failed to write downloaded data to {}: {err}",
            destination.display()
        ));
    }

    Ok(())
}

fn prepare_downloaded_file(base_dir: &Path, path: &Path) -> Result<PathBuf> {
    if let Some(extension) = path.extension().and_then(|ext| ext.to_str()) {
        match extension.to_ascii_lowercase().as_str() {
            "gz" => return decompress_gzip(path),
            "zip" => return extract_zip_archive(base_dir, path),
            _ => {}
        }
    }

    Ok(path.to_path_buf())
}

fn decompress_gzip(path: &Path) -> Result<PathBuf> {
    let target = path.with_extension("");
    let input = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut decoder = MultiGzDecoder::new(input);
    let mut output =
        File::create(&target).with_context(|| format!("failed to create {}", target.display()))?;
    io::copy(&mut decoder, &mut output)
        .with_context(|| format!("failed to decompress {}", path.display()))?;
    Ok(target)
}

fn extract_zip_archive(base_dir: &Path, path: &Path) -> Result<PathBuf> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut archive = ZipArchive::new(file)
        .with_context(|| format!("failed to read zip archive {}", path.display()))?;

    let mut selected_index = None;
    let mut selected_name = None;

    for i in 0..archive.len() {
        let entry = archive.by_index(i)?;
        if entry.is_dir() {
            continue;
        }

        let Some(name) = entry.enclosed_name() else {
            continue;
        };

        let name_str = name.to_string_lossy();
        let name_lower = name_str.to_ascii_lowercase();
        let is_fasta = FASTA_EXTENSIONS.iter().any(|ext| name_lower.ends_with(ext));
        if selected_index.is_none() || is_fasta {
            selected_index = Some(i);
            selected_name = Some(name.to_owned());
        }

        if is_fasta {
            break;
        }
    }

    let Some(index) = selected_index else {
        bail!("zip archive {} does not contain files", path.display());
    };

    let mut entry = archive.by_index(index)?;
    let output_name = selected_name
        .and_then(|name| name.file_name().map(|os| os.to_owned()))
        .unwrap_or_else(|| entry.name().into());
    let mut output_path = base_dir.to_path_buf();
    output_path.push(output_name);

    {
        let mut output = File::create(&output_path)
            .with_context(|| format!("failed to create {}", output_path.display()))?;
        io::copy(&mut entry, &mut output)
            .with_context(|| format!("failed to extract {}", entry.name()))?;
    }

    prepare_downloaded_file(base_dir, &output_path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use zip::ZipWriter;

    #[test]
    fn decompresses_gzip_files() {
        let dir = TempDir::new().unwrap();
        let source_path = dir.path().join("test.fa.gz");
        let data = b"ACGT\n";
        {
            let file = File::create(&source_path).unwrap();
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            encoder.write_all(data).unwrap();
            encoder.finish().unwrap();
        }

        let decompressed = prepare_downloaded_file(dir.path(), &source_path).unwrap();
        let contents = fs::read_to_string(decompressed).unwrap();
        assert_eq!(contents, "ACGT\n");
    }

    #[test]
    fn extracts_zip_archives() {
        let dir = TempDir::new().unwrap();
        let zip_path = dir.path().join("archive.zip");
        {
            let file = File::create(&zip_path).unwrap();
            let mut zip = ZipWriter::new(file);
            let options = zip::write::FileOptions::default();
            zip.start_file("test.fa", options).unwrap();
            zip.write_all(b">chr1\nACGT\n").unwrap();
            zip.finish().unwrap();
        }

        let extracted = prepare_downloaded_file(dir.path(), &zip_path).unwrap();
        let contents = fs::read_to_string(extracted).unwrap();
        assert!(contents.contains("ACGT"));
    }

    #[test]
    fn extracts_nested_compressed_files() {
        let dir = TempDir::new().unwrap();
        let nested_zip = dir.path().join("nested.zip");

        let inner_gz_path = dir.path().join("inner.fa.gz");
        {
            let file = File::create(&inner_gz_path).unwrap();
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            encoder.write_all(b">chr1\nNNNN\n").unwrap();
            encoder.finish().unwrap();
        }

        {
            let file = File::create(&nested_zip).unwrap();
            let mut zip = ZipWriter::new(file);
            let options = zip::write::FileOptions::default();
            let data = fs::read(&inner_gz_path).unwrap();
            zip.start_file("inner.fa.gz", options).unwrap();
            zip.write_all(&data).unwrap();
            zip.finish().unwrap();
        }

        let extracted = prepare_downloaded_file(dir.path(), &nested_zip).unwrap();
        let contents = fs::read_to_string(extracted).unwrap();
        assert!(contents.contains("NNNN"));
    }
}
