//! Panel writer for outputting padded reference panels.
//!
//! This module handles writing the modified (padded) reference panel
//! that includes novel alleles from user data.

use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use noodles::vcf::{self, variant::io::Write as VariantWrite};

use crate::panel::{PaddedPanel, PanelSite};

/// Write the padded panel to a VCF file.
///
/// This reads the original panel and injects modifications as it streams through.
/// Novel sites are appended at the end.
pub fn write_padded_panel<P: AsRef<Path>>(
    original_panel_path: P,
    padded: &PaddedPanel,
    output_path: P,
) -> Result<()> {
    let original = original_panel_path.as_ref();
    let output = output_path.as_ref();

    // Open original panel for reading
    let file = std::fs::File::open(original)
        .with_context(|| format!("failed to open original panel {}", original.display()))?;

    let reader: Box<dyn BufRead> = if original.to_string_lossy().ends_with(".gz") {
        Box::new(io::BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(io::BufReader::new(file))
    };

    let mut vcf_reader = vcf::io::Reader::new(reader);
    let header = vcf_reader.read_header()?;

    // Create output file
    let out_file = std::fs::File::create(output)
        .with_context(|| format!("failed to create padded panel {}", output.display()))?;

    // Wrap in gzip if needed
    let writer: Box<dyn Write> = if output.to_string_lossy().ends_with(".gz") {
        Box::new(flate2::write::GzEncoder::new(
            BufWriter::new(out_file),
            flate2::Compression::default(),
        ))
    } else {
        Box::new(BufWriter::new(out_file))
    };

    let mut vcf_writer = vcf::io::Writer::new(writer);

    // Write header (same as original)
    vcf_writer.write_header(&header)?;

    let mut records_written = 0u64;

    // Stream through original panel records
    for result in vcf_reader.record_bufs(&header) {
        let mut record = result?;

        let chrom = record.reference_sequence_name().to_string();
        let pos = record.variant_start().map(usize::from).unwrap_or(0) as u64;

        // Check if we need to add alleles to this site
        if let Some(added_alts) = padded.added_alts(&chrom, pos) {
            // Modify the record's ALT alleles
            let mut new_alts: Vec<_> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|a| a.to_string())
                .collect();
            new_alts.extend(added_alts.iter().cloned());

            // Create new AlternateBases from the allele strings
            let alt_bases = vcf::variant::record_buf::AlternateBases::from(new_alts);
            *record.alternate_bases_mut() = alt_bases;

            tracing::debug!(
                "Padded site {}:{} with {} new alleles",
                chrom,
                pos,
                added_alts.len()
            );
        }

        vcf_writer.write_variant_record(&header, &record)?;
        records_written += 1;

        if records_written.is_multiple_of(1_000_000) {
            tracing::info!("Written {} panel records...", records_written);
        }
    }

    // Write novel sites at the end
    let novel_count = padded.novel_site_count();
    if novel_count > 0 {
        tracing::info!("Writing {} novel sites...", novel_count);

        for site in padded.novel_sites() {
            let record = site_to_record(site)?;
            vcf_writer.write_variant_record(&header, &record)?;
            records_written += 1;
        }
    }

    tracing::info!(
        "Padded panel complete: {} records, {} modified sites, {} novel sites",
        records_written,
        padded.modified_site_count(),
        novel_count
    );

    Ok(())
}

/// Convert a PanelSite to a VCF RecordBuf.
fn site_to_record(site: &PanelSite) -> Result<vcf::variant::record_buf::RecordBuf> {
    use noodles::core::Position;
    use vcf::variant::record_buf::{AlternateBases, Ids, RecordBuf};

    let pos = Position::try_from(site.pos as usize)
        .map_err(|e| anyhow::anyhow!("invalid position {}: {}", site.pos, e))?;

    // Create AlternateBases from the allele strings
    let alt_bases = AlternateBases::from(site.alt_alleles.clone());

    let ids: Ids = if let Some(id) = &site.id {
        vec![id.clone()].into_iter().collect()
    } else {
        // Generate synthetic ID: chr:pos:ref:alt
        let alt_str = if site.alt_alleles.is_empty() {
            ".".to_string()
        } else {
            site.alt_alleles.join(",")
        };
        let synthetic_id = format!(
            "{}:{}:{}:{}",
            site.chrom, site.pos, site.ref_allele, alt_str
        );
        vec![synthetic_id].into_iter().collect()
    };

    let record = RecordBuf::builder()
        .set_reference_sequence_name(site.chrom.clone())
        .set_variant_start(pos)
        .set_ids(ids)
        .set_reference_bases(site.ref_allele.clone())
        .set_alternate_bases(alt_bases)
        .build();

    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_site_to_record() {
        let site = PanelSite {
            chrom: "1".to_string(),
            pos: 1000,
            id: Some("rs123".to_string()),
            ref_allele: "A".to_string(),
            alt_alleles: vec!["G".to_string()],
        };

        let record = site_to_record(&site).unwrap();
        assert_eq!(record.reference_sequence_name(), "1");
        assert_eq!(record.reference_bases().to_string(), "A");
    }
}
