use std::io;
use noodles::vcf::variant::record_buf::{RecordBuf, AlternateBases, Samples, samples::Keys};
use noodles::vcf::variant::record_buf::samples::sample::Value;
use noodles::core::Position;

use crate::dtc::{self, Record as DtcRecord, Allele as DtcAllele};
use crate::reference::ReferenceGenome;
use crate::conversion::{ConversionConfig, Ploidy, determine_ploidy, format_genotype};
use clap::ValueEnum;

#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum InputFormat {
    /// Direct-to-consumer text output (23andMe, etc.)
    Dtc,
    /// Variant Call Format
    Vcf,
    /// Binary Call Format
    Bcf,
    /// Detect format automatically
    Auto,
}

impl InputFormat {
    pub fn detect(path: &std::path::Path) -> Self {
        if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
            match ext.to_lowercase().as_str() {
                "vcf" | "vcf.gz" => return Self::Vcf,
                "bcf" => return Self::Bcf,
                "txt" | "csv" | "tsv" => return Self::Dtc,
                _ => {}
            }
        }
        // TODO: magic bytes check
        Self::Dtc // Default fallback
    }
}

/// Trait for a source of genomic variants.
/// Yields `RecordBuf` (generic VCF record wrapper).
pub trait VariantSource {
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>>;
}

impl<T: VariantSource + ?Sized> VariantSource for Box<T> {
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>> {
        (**self).next_variant(summary)
    }
}

/// Adapter for reading DTC (23andMe) files.
pub struct DtcSource<R> {
    reader: dtc::Reader<R>,
    reference: ReferenceGenome,
    config: ConversionConfig,
    header_keys: Keys, 
}

impl<R: std::io::BufRead> DtcSource<R> {
    pub fn new(reader: dtc::Reader<R>, reference: ReferenceGenome, config: ConversionConfig) -> Self {
        let keys: Keys = vec![String::from("GT")].into_iter().collect();
        Self {
            reader,
            reference,
            config,
            header_keys: keys,
        }
    }

    fn convert_record(&self, record: DtcRecord, summary: &mut crate::ConversionSummary) -> io::Result<Option<RecordBuf>> {
        // Validation and normalization logic moved from conversion.rs
        
        let canonical_name = if let Some(name) = self.reference.resolve_contig_name(&record.chromosome) {
             name.to_string()
        } else {
             // Skipping unknown chromosome
             summary.unknown_chromosomes += 1;
             return Ok(None);
        };
        
        // Validate position
        let position = record.position;
        
        // Use reference.base instead of get_context
        let reference_base = match self.reference.base(&canonical_name, position) {
            Ok(base) => base,
            Err(_) => {
                summary.reference_failures += 1;
                return Ok(None); 
            }
        };
        
        // Parse alleles
        let allele_states = match record.parse_alleles() {
            Ok(states) => states,
            Err(_) => {
                summary.invalid_genotypes += 1;
                return Ok(None);
            }
        };
        
        // Biological Ploidy
        let ploidy = determine_ploidy(&canonical_name, position, self.config.sex, self.config.par_boundaries.as_ref());
        
        // Normalize
        let normalized_alleles = match ploidy {
            Ploidy::Zero => {
                summary.missing_genotype_records += 1;
                return Ok(None);
            }
            Ploidy::Haploid => {
                 let unique: Vec<_> = allele_states
                     .iter()
                     .filter(|&a| *a != DtcAllele::Missing)
                     .cloned()
                     .collect();
                 
                 let mut uniq_vals = unique.clone();
                 uniq_vals.dedup();
                 
                 if uniq_vals.is_empty() {
                     vec![DtcAllele::Missing]
                 } else if uniq_vals.len() == 1 {
                     vec![uniq_vals[0].clone()]
                 } else {
                     vec![DtcAllele::Missing]
                 }
             },
             Ploidy::Diploid => allele_states,
        };
        
        // Build Alt Bases
        let mut alt_bases = Vec::new();
        for allele in normalized_alleles.iter() {
             match allele {
                 DtcAllele::Base(base) => {
                     let b = base.to_uppercase();
                     let r = reference_base.to_string().to_uppercase();
                     if b != r && !alt_bases.contains(&b) {
                         alt_bases.push(b);
                     }
                 },
                 DtcAllele::Deletion => {
                     let del = String::from("DEL");
                     if !alt_bases.contains(&del) { alt_bases.push(del); }
                 },
                 DtcAllele::Insertion => {
                     let ins = String::from("INS");
                     if !alt_bases.contains(&ins) { alt_bases.push(ins); }
                 },
                 DtcAllele::Missing => {}
             }
        }
        
        if alt_bases.is_empty() && !self.config.include_reference_sites {
             summary.skipped_reference_sites += 1;
             return Ok(None);
        }
        
        // Formatted Alts (add < > if needed)
        let formatted_alts: Vec<String> = alt_bases.iter().map(|a| {
            if a == "DEL" || a == "INS" { format!("<{}>", a) } else { a.clone() }
        }).collect();
        
        let alternate_bases = AlternateBases::from(formatted_alts);
        
        // Genotype string
        let genotype_string = format_genotype(&normalized_alleles, reference_base, &alt_bases)
             .map_err(|e| {
                 summary.invalid_genotypes += 1;
                 io::Error::new(io::ErrorKind::InvalidData, e.to_string())
             })?;
             
        // Build RecordBuf
        // IDs
        let ids = record.id.as_ref().map(|id| [id.clone()].into_iter().collect::<noodles::vcf::variant::record_buf::Ids>());
        
        let samples_buf = Samples::new(
            self.header_keys.clone(),
            vec![vec![Some(Value::String(genotype_string))]]
        );
        
        let mut builder = RecordBuf::builder()
            .set_reference_sequence_name(canonical_name)
            .set_variant_start(Position::try_from(position as usize).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?)
            .set_reference_bases(reference_base.to_string())
            .set_alternate_bases(alternate_bases)
            .set_samples(samples_buf);
            
        if let Some(i) = ids {
            builder = builder.set_ids(i);
        }
        
        Ok(Some(builder.build()))
    }
}

impl<R: std::io::BufRead> VariantSource for DtcSource<R> {
    fn next_variant(&mut self, summary: &mut crate::ConversionSummary) -> Option<io::Result<RecordBuf>> {
        loop {
            // Read next generic record
            let record = match self.reader.next() {
                Some(Ok(r)) => r,
                Some(Err(e)) => return Some(Err(io::Error::new(io::ErrorKind::InvalidData, e.to_string()))),
                None => return None,
            };
            
            match self.convert_record(record, summary) {
                Ok(Some(buf)) => return Some(Ok(buf)),
                Ok(None) => continue, // Filtered, try next
                Err(e) => return Some(Err(e)),
            }
        }
    }
}
