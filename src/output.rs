use std::{fs::File, io::BufWriter, path::Path};

use anyhow::{Context, Result, anyhow};
use noodles::bcf;
use noodles::vcf::{
    self as vcf,
    variant::{RecordBuf, io::Write},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Vcf,
    Bcf,
}

impl std::str::FromStr for OutputFormat {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "vcf" => Ok(Self::Vcf),
            "bcf" => Ok(Self::Bcf),
            other => Err(anyhow!("unsupported output format: {other}")),
        }
    }
}

pub fn write_records<P>(
    path: P,
    format: OutputFormat,
    header: &vcf::Header,
    records: &[RecordBuf],
) -> Result<()>
where
    P: AsRef<Path>,
{
    match format {
        OutputFormat::Vcf => write_vcf(path, header, records),
        OutputFormat::Bcf => write_bcf(path, header, records),
    }
}

fn write_vcf<P>(path: P, header: &vcf::Header, records: &[RecordBuf]) -> Result<()>
where
    P: AsRef<Path>,
{
    let writer = File::create(path.as_ref())
        .map(BufWriter::new)
        .with_context(|| format!("failed to create VCF file at {}", path.as_ref().display()))?;
    let mut writer = vcf::io::Writer::new(writer);
    writer
        .write_header(header)
        .context("failed to write VCF header")?;
    for record in records {
        writer
            .write_variant_record(header, record)
            .context("failed to write VCF record")?;
    }
    Ok(())
}

fn write_bcf<P>(path: P, header: &vcf::Header, records: &[RecordBuf]) -> Result<()>
where
    P: AsRef<Path>,
{
    let writer = File::create(path.as_ref())
        .map(BufWriter::new)
        .with_context(|| format!("failed to create BCF file at {}", path.as_ref().display()))?;
    let mut writer = bcf::io::Writer::new(writer);
    writer
        .write_header(header)
        .context("failed to write BCF header")?;
    for record in records {
        writer
            .write_variant_record(header, record)
            .context("failed to write BCF record")?;
    }
    Ok(())
}
