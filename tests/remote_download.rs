use std::fs;

use anyhow::Result;
use convert_genome::remote;
use serial_test::serial;
use url::Url;

/// Downloads a remote reference and validates it exists.
/// Cleans up temp files after validation to conserve disk space.
fn ensure_remote_reference(url: &str) -> Result<()> {
    let parsed = Url::parse(url)?;
    let resource = remote::fetch_remote_resource(&parsed)?;
    let path = resource.local_path().to_path_buf();
    let metadata = fs::metadata(&path)?;
    assert!(metadata.len() > 0, "downloaded file was empty");

    // Clean up to conserve disk space on CI runners
    // The resource will be dropped, but explicitly remove the file too
    if let Some(parent) = path.parent() {
        drop(fs::remove_dir_all(parent));
    }

    Ok(())
}

#[test]
#[serial]
fn download_grch38_primary_assembly_from_ebi() -> Result<()> {
    ensure_remote_reference(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz",
    )
}

#[test]
#[serial]
fn download_grch38_illumina_bundle() -> Result<()> {
    ensure_remote_reference(
        "https://webdata.illumina.com/downloads/productfiles/microarray-analytics-array/GRCh38_genome.zip",
    )
}

#[test]
#[serial]
fn download_grch38_primary_assembly_from_ensembl() -> Result<()> {
    ensure_remote_reference(
        "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    )
}

#[test]
#[serial]
fn download_grch37_illumina_bundle() -> Result<()> {
    ensure_remote_reference(
        "https://webdata.illumina.com/downloads/productfiles/microarray-analytics-array/GRCh37_genome.zip",
    )
}
