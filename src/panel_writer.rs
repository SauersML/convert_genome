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

    let chrom_rank = header
        .contigs()
        .keys()
        .enumerate()
        .map(|(i, c)| (c.to_string(), i))
        .collect::<std::collections::HashMap<_, _>>();

    let mut novel_sites: Vec<&PanelSite> = padded.novel_sites().collect();
    novel_sites.sort_by(|a, b| compare_sites(a, b, &chrom_rank));
    let mut next_novel_idx = 0usize;

    // Stream through original panel records
    for result in vcf_reader.record_bufs(&header) {
        let mut record = result?;

        let chrom = record.reference_sequence_name().to_string();
        let pos = record.variant_start().map(usize::from).unwrap_or(0) as u64;

        // Emit novel sites that should come before this record
        while let Some(site) = novel_sites.get(next_novel_idx) {
            if compare_site_to_record(site, &chrom, pos, &chrom_rank) == std::cmp::Ordering::Less {
                let mut novel_record = site_to_record(site)?;
                if let Some(samples) = build_missing_samples(&header)? {
                    *novel_record.samples_mut() = samples;
                }
                vcf_writer.write_variant_record(&header, &novel_record)?;
                records_written += 1;
                next_novel_idx += 1;
            } else {
                break;
            }
        }

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

        if header.sample_names().len() > 0 && record.samples().values().count() == 0 {
            if let Some(samples) = build_missing_samples(&header)? {
                *record.samples_mut() = samples;
            }
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

        for site in novel_sites.iter().skip(next_novel_idx) {
            let mut record = site_to_record(site)?;
            if let Some(samples) = build_missing_samples(&header)? {
                *record.samples_mut() = samples;
            }
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

fn compare_sites(
    a: &PanelSite,
    b: &PanelSite,
    chrom_rank: &std::collections::HashMap<String, usize>,
) -> std::cmp::Ordering {
    let a_rank = contig_rank(&a.chrom, chrom_rank);
    let b_rank = contig_rank(&b.chrom, chrom_rank);
    match a_rank.cmp(&b_rank) {
        std::cmp::Ordering::Equal => match a.chrom.cmp(&b.chrom) {
            std::cmp::Ordering::Equal => a.pos.cmp(&b.pos),
            other => other,
        },
        other => other,
    }
}

fn compare_site_to_record(
    site: &PanelSite,
    chrom: &str,
    pos: u64,
    chrom_rank: &std::collections::HashMap<String, usize>,
) -> std::cmp::Ordering {
    let site_rank = contig_rank(&site.chrom, chrom_rank);
    let record_rank = contig_rank(chrom, chrom_rank);
    match site_rank.cmp(&record_rank) {
        std::cmp::Ordering::Equal => match site.chrom.as_str().cmp(chrom) {
            std::cmp::Ordering::Equal => site.pos.cmp(&pos),
            other => other,
        },
        other => other,
    }
}

fn contig_rank(chrom: &str, chrom_rank: &std::collections::HashMap<String, usize>) -> usize {
    if let Some(rank) = chrom_rank.get(chrom) {
        return *rank;
    }
    if chrom.starts_with("chr") {
        if let Some(rank) = chrom_rank.get(chrom.trim_start_matches("chr")) {
            return *rank;
        }
    } else {
        let alt = format!("chr{}", chrom);
        if let Some(rank) = chrom_rank.get(&alt) {
            return *rank;
        }
    }
    usize::MAX
}

fn build_missing_samples(header: &vcf::Header) -> Result<Option<vcf::variant::record_buf::Samples>> {
    use vcf::variant::record_buf::samples::{sample::Value, Keys};
    use vcf::variant::record_buf::Samples;

    let sample_count = header.sample_names().len();
    if sample_count == 0 {
        return Ok(None);
    }

    let format_keys: Vec<String> = header.formats().keys().map(|k| k.to_string()).collect();
    if format_keys.is_empty() {
        anyhow::bail!("panel header has samples but no FORMAT definitions");
    }

    let keys: Keys = format_keys.into_iter().collect();
    let keys_len = keys.as_ref().len();
    let empty_sample: Vec<Option<Value>> = vec![None; keys_len];
    let values = vec![empty_sample; sample_count];

    Ok(Some(Samples::new(keys, values)))
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

    #[test]
    fn test_padded_panel_includes_samples_for_novel_sites() {
        let dir = tempfile::tempdir().unwrap();
        let original_path = dir.path().join("panel.vcf");
        let output_path = dir.path().join("panel_out.vcf");

        let vcf = concat!(
            "##fileformat=VCFv4.3\n",
            "##contig=<ID=chr22,length=1000>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n",
            "chr22\t1\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n"
        );
        std::fs::write(&original_path, vcf).unwrap();

        let panel_index = crate::panel::PanelIndex::load(&original_path).unwrap();
        let mut padded = crate::panel::PaddedPanel::new(panel_index);
        padded.get_or_add_allele_index("chr22", 2, "C", "A");

        write_padded_panel(&original_path, &padded, &output_path).unwrap();

        let reader = std::fs::File::open(&output_path).unwrap();
        let mut vcf_reader = vcf::io::Reader::new(std::io::BufReader::new(reader));
        let header = vcf_reader.read_header().unwrap();
        let expected_samples = header.sample_names().len();

        for result in vcf_reader.record_bufs(&header) {
            let record = result.unwrap();
            let sample_count = record.samples().values().count();
            assert_eq!(
                sample_count,
                expected_samples,
                "record missing samples at {}:{}",
                record.reference_sequence_name(),
                record.variant_start().map(usize::from).unwrap_or(0)
            );
        }
    }

    #[test]
    fn test_padded_panel_output_sorted_by_position() {
        let dir = tempfile::tempdir().unwrap();
        let original_path = dir.path().join("panel.vcf");
        let output_path = dir.path().join("panel_out.vcf");

        let vcf = concat!(
            "##fileformat=VCFv4.3\n",
            "##contig=<ID=chr1,length=1000>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n",
            "chr1\t200\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
        );
        std::fs::write(&original_path, vcf).unwrap();

        let panel_index = crate::panel::PanelIndex::load(&original_path).unwrap();
        let mut padded = crate::panel::PaddedPanel::new(panel_index);
        padded.get_or_add_allele_index("chr1", 100, "C", "A");

        write_padded_panel(&original_path, &padded, &output_path).unwrap();

        let reader = std::fs::File::open(&output_path).unwrap();
        let mut vcf_reader = vcf::io::Reader::new(std::io::BufReader::new(reader));
        let header = vcf_reader.read_header().unwrap();

        let mut last_pos = 0usize;
        for result in vcf_reader.record_bufs(&header) {
            let record = result.unwrap();
            let pos = record.variant_start().map(usize::from).unwrap_or(0);
            assert!(
                pos >= last_pos,
                "panel output not sorted: {} then {}",
                last_pos,
                pos
            );
            last_pos = pos;
        }
    }
}
