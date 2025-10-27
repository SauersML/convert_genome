use std::{cmp::Ordering, sync::Arc};

use anyhow::{Context, Result};
use rayon::prelude::*;

use noodles::core::Position;
use noodles::vcf::{
    self as vcf,
    header::{
        FileFormat,
        record::Value as HeaderValue,
        record::{
            key::Other,
            value::{
                Map,
                map::{Contig, Filter, Format},
            },
        },
    },
    variant::{
        RecordBuf,
        record::samples::keys::key,
        record::samples::series::value::genotype::Phasing,
        record_buf::{
            AlternateBases, Filters, Ids, Info, Samples,
            samples::{
                Keys,
                sample::value::{
                    Genotype as GenotypeValue, Value as SampleValue,
                    genotype::Allele as GenotypeAllele,
                },
            },
        },
    },
};

use crate::{dtc::DtcRecord, reference::ReferenceGenome};

#[derive(Debug, Clone)]
pub struct ConversionOptions {
    pub sample_name: String,
    pub skip_reference_calls: bool,
}

impl Default for ConversionOptions {
    fn default() -> Self {
        Self {
            sample_name: String::from("sample"),
            skip_reference_calls: false,
        }
    }
}

pub struct ConversionResult {
    pub header: vcf::Header,
    pub records: Vec<RecordBuf>,
}

pub fn convert(
    records: Vec<DtcRecord>,
    reference: ReferenceGenome,
    options: &ConversionOptions,
) -> Result<ConversionResult> {
    if records.is_empty() {
        let header = build_header(&reference, options)?;
        return Ok(ConversionResult {
            header,
            records: Vec::new(),
        });
    }

    let mut records = records;
    records.sort_by(|a, b| compare_records(a, b));

    let reference = Arc::new(reference);
    let shared_options = Arc::new(options.clone());

    let converted: Vec<RecordBuf> = records
        .into_par_iter()
        .filter_map(|record| {
            let reference = Arc::clone(&reference);
            let options = Arc::clone(&shared_options);
            match convert_record(reference.as_ref(), options.as_ref(), record) {
                Ok(Some(record)) => Some(Ok(record)),
                Ok(None) => None,
                Err(err) => Some(Err(err)),
            }
        })
        .collect::<Result<Vec<_>>>()?;

    let header = build_header(reference.as_ref(), shared_options.as_ref())?;

    Ok(ConversionResult {
        header,
        records: converted,
    })
}

fn compare_records(a: &DtcRecord, b: &DtcRecord) -> Ordering {
    let a_rank = contig_rank(&a.chrom);
    let b_rank = contig_rank(&b.chrom);

    match a_rank.cmp(&b_rank) {
        Ordering::Equal => a.position.cmp(&b.position),
        other => other,
    }
}

fn contig_rank(chrom: &str) -> u32 {
    let raw = chrom.trim_start_matches("chr");
    match raw {
        "X" | "x" => 23,
        "Y" | "y" => 24,
        "M" | "m" | "MT" | "mt" => 25,
        "XY" | "xy" => 26,
        _ => raw.parse::<u32>().unwrap_or(u32::MAX - 1),
    }
}

fn convert_record(
    reference: &ReferenceGenome,
    options: &ConversionOptions,
    record: DtcRecord,
) -> Result<Option<RecordBuf>> {
    let ref_base = reference.base(&record.chrom, record.position)?;

    let mut alt_bases: Vec<char> = Vec::new();
    let allele_chars = record.genotype.allele_chars();
    let mut allele_indices: Vec<Option<usize>> = Vec::new();

    for allele in allele_chars.iter() {
        match allele {
            Some(base) => {
                let base_upper = base.to_ascii_uppercase();
                if base_upper == ref_base {
                    allele_indices.push(Some(0));
                } else {
                    let alt_index = match alt_bases.iter().position(|&a| a == base_upper) {
                        Some(i) => i + 1,
                        None => {
                            alt_bases.push(base_upper);
                            alt_bases.len()
                        }
                    };
                    allele_indices.push(Some(alt_index));
                }
            }
            None => allele_indices.push(None),
        }
    }

    if options.skip_reference_calls
        && allele_indices
            .iter()
            .all(|allele| matches!(allele, Some(0)))
    {
        return Ok(None);
    }

    let alt_strings: Vec<String> = alt_bases.iter().map(|&c| c.to_string()).collect();

    let mut genotype = GenotypeValue::default();
    {
        let genotype_alleles = genotype.as_mut();
        for allele in &allele_indices {
            let allele_value = match allele {
                Some(index) => GenotypeAllele::new(Some(*index), Phasing::Unphased),
                None => GenotypeAllele::new(None, Phasing::Unphased),
            };
            genotype_alleles.push(allele_value);
        }
    }
    let sample_value = SampleValue::from(genotype);

    let keys: Keys = [String::from(key::GENOTYPE)].into_iter().collect();
    let samples = Samples::new(keys, vec![vec![Some(sample_value)]]);

    let ids: Ids = if record.id == "." {
        Ids::default()
    } else {
        [record.id].into_iter().collect()
    };

    let pos = Position::try_from(usize::try_from(record.position).context("position overflow")?)
        .with_context(|| format!("invalid position for {}", record.chrom))?;

    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(record.chrom)
        .set_variant_start(pos)
        .set_reference_bases(ref_base.to_string())
        .set_samples(samples)
        .set_filters(Filters::pass())
        .set_info(Info::default());

    if !ids.as_ref().is_empty() {
        builder = builder.set_ids(ids);
    }

    if !alt_strings.is_empty() {
        builder = builder.set_alternate_bases(AlternateBases::from(alt_strings));
    }

    let record = builder.build();

    Ok(Some(record))
}

#[cfg(test)]
fn format_genotype(indices: &[Option<usize>]) -> String {
    if indices.is_empty() {
        return String::from("./.");
    }

    let mut parts = Vec::with_capacity(indices.len());

    for index in indices {
        match index {
            Some(i) => parts.push(i.to_string()),
            None => parts.push(String::from(".")),
        }
    }

    if parts.iter().all(|part| part == ".") {
        if parts.len() == 1 {
            String::from(".")
        } else {
            parts.join("/")
        }
    } else {
        parts.join("/")
    }
}

fn build_header(reference: &ReferenceGenome, options: &ConversionOptions) -> Result<vcf::Header> {
    let mut builder = vcf::Header::builder()
        .set_file_format(FileFormat::default())
        .add_filter("PASS", Map::<Filter>::pass())
        .add_format(key::GENOTYPE, Map::<Format>::from(key::GENOTYPE))
        .add_sample_name(options.sample_name.clone());

    for contig in reference.index().as_ref() {
        let mut map = Map::<Contig>::new();
        if let Ok(length) = usize::try_from(contig.length()) {
            *map.length_mut() = Some(length);
        }
        builder = builder.add_contig(contig.name().to_string(), map);
    }

    let source_value = HeaderValue::from("convert_genome");
    let source_key = "source".parse::<Other>()?;
    let reference_key = "reference".parse::<Other>()?;
    builder = builder.insert(source_key, source_value)?.insert(
        reference_key,
        HeaderValue::from(reference.source().display().to_string()),
    )?;

    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dtc::DtcGenotype;

    #[test]
    fn genotype_formatting() {
        assert_eq!(format_genotype(&[Some(0), Some(1)]), "0/1");
        assert_eq!(format_genotype(&[Some(1), Some(1)]), "1/1");
        assert_eq!(format_genotype(&[None, None]), "./.");
    }

    #[test]
    fn convert_produces_reference_genotype() -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Write;
        use tempfile::tempdir;

        let dir = tempdir()?;
        let fasta_path = dir.path().join("ref.fa");
        {
            let mut file = std::fs::File::create(&fasta_path)?;
            writeln!(file, ">chr1")?;
            writeln!(file, "T")?;
        }

        let reference = ReferenceGenome::open(&fasta_path)?;
        let records = vec![DtcRecord {
            id: "rsTest".into(),
            chrom: "chr1".into(),
            position: 1,
            genotype: DtcGenotype::from_call("TT"),
        }];

        let result = super::convert(records, reference, &ConversionOptions::default())?;
        assert_eq!(result.records.len(), 1);

        let record = &result.records[0];
        let sample = record.samples().values().next().expect("missing sample");
        let genotype_value = sample
            .values()
            .get(0)
            .and_then(|value| value.as_ref())
            .expect("missing genotype value");

        if let SampleValue::Genotype(genotype) = genotype_value {
            let positions: Vec<_> = genotype
                .as_ref()
                .iter()
                .map(|allele| allele.position())
                .collect();
            assert_eq!(positions, vec![Some(0), Some(0)]);
        } else {
            panic!("expected genotype value");
        }

        Ok(())
    }

    #[test]
    fn sort_order() {
        let a = DtcRecord {
            id: "rs1".into(),
            chrom: "chr1".into(),
            position: 100,
            genotype: DtcGenotype::missing(),
        };
        let b = DtcRecord {
            id: "rs2".into(),
            chrom: "chr2".into(),
            position: 50,
            genotype: DtcGenotype::missing(),
        };

        assert!(compare_records(&a, &b) == Ordering::Less);
    }
}
