pub mod converter;
pub mod dtc;
pub mod output;
pub mod reference;

use std::{fs::File, path::Path};

use anyhow::{Context, Result};

pub use converter::{ConversionOptions, ConversionResult, convert};
pub use dtc::{DtcGenotype, DtcReader, DtcRecord};
pub use output::{OutputFormat, write_records};
pub use reference::ReferenceGenome;

pub fn load_dtc_records<P>(path: P) -> Result<Vec<DtcRecord>>
where
    P: AsRef<Path>,
{
    let reader = File::open(path.as_ref())
        .with_context(|| format!("failed to open DTC file {}", path.as_ref().display()))?;
    let reader = DtcReader::new(reader);
    let mut records = Vec::new();
    for record in reader.records() {
        records.push(record?);
    }
    Ok(records)
}
