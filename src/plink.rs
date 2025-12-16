use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record_buf::RecordBuf;
use noodles::vcf::variant::record_buf::samples::sample::Value;
// Import traits for RecordBuf methods
use noodles::vcf::variant::record::{AlternateBases, Ids};

use crate::cli::Sex;

/// PLINK 1.9 Writer for .bed, .bim, .fam files
pub struct PlinkWriter {
    bed: BufWriter<File>,
    bim: BufWriter<File>,
    fam: BufWriter<File>,
    samples_written: bool,
}

impl PlinkWriter {
    /// Create a new PlinkWriter.
    /// `prefix` is the base path (e.g., "output"). Extensions (.bed, .bim, .fam) will be appended.
    pub fn new<P: AsRef<Path>>(prefix: P) -> io::Result<Self> {
        let prefix = prefix.as_ref();
        // Strip extension if provided (simple handling)
        let stem = if prefix.extension().map_or(false, |e| e == "vcf" || e == "bcf") {
            prefix.with_extension("")
        } else {
            prefix.to_path_buf()
        };

        let bed_path = stem.with_extension("bed");
        let bim_path = stem.with_extension("bim");
        let fam_path = stem.with_extension("fam");

        let mut bed = BufWriter::new(File::create(bed_path)?);
        let bim = BufWriter::new(File::create(bim_path)?);
        let fam = BufWriter::new(File::create(fam_path)?);

        // Write Magic Bytes for .bed (SNP-major mode)
        // 0x6C 0x1B 0x01
        bed.write_all(&[0x6C, 0x1B, 0x01])?;

        Ok(Self {
            bed,
            bim,
            fam,
            samples_written: false,
        })
    }

    /// Write the .fam file using sample ID and Sex.
    /// This should be called once before writing variants.
    pub fn write_fam(&mut self, sample_id: &str, sex: Sex) -> io::Result<()> {
        if self.samples_written {
            return Ok(());
        }
        // FAM format: FID IID PID MID Sex Phenotype
        // Sex: 1=Male, 2=Female, 0=Unknown
        let sex_code = match sex {
            Sex::Male => "1",
            Sex::Female => "2",
        };
        // Phenotype -9 = missing
        writeln!(
            self.fam,
            "{0}\t{0}\t0\t0\t{1}\t-9",
            sample_id, sex_code
        )?;
        self.samples_written = true;
        Ok(())
    }

    pub fn write_variant(&mut self, record: &RecordBuf) -> io::Result<()> {
        // Enforce Biallelic: Ref + 1 Alt
        if record.alternate_bases().len() > 1 {
             // RecordBuf usually uses variant_start() or alignment_start(). 
             let pos = record.variant_start().map(|p| usize::from(p)).unwrap_or(0);
             tracing::warn!("Skipping multiallelic site at {}:{}", record.reference_sequence_name(), pos);
             return Ok(());
        }

        // Write BIM line
        // Chrom Ident cm Pos Ref Alt
        let chrom = record.reference_sequence_name();
        
        // Ids: join if multiple
        let id = if record.ids().is_empty() {
             ".".to_string()
        } else {
             record.ids().iter().map(|s| s.to_string()).collect::<Vec<_>>().join(";")
        };
        
        // Pos
        let pos = record.variant_start().map(|p| usize::from(p)).unwrap_or(0);
        
        let ref_base = record.reference_bases();
        
        // alternate_bases
        let alt_base_str = if record.alternate_bases().is_empty() {
             ".".to_string()
        } else {
             if let Some(res) = record.alternate_bases().as_ref().iter().next() {
                 res.to_string()
             } else {
                 ".".to_string()
             }
        };

        // Write to BIM
        writeln!(self.bim, "{}\t{}\t0\t{}\t{}\t{}", chrom, id, pos, ref_base, alt_base_str)?;

        // Write BED Genotypes
        // User mapping (Variant-Major):
        // 00 (0) Homozygous for first allele (Ref)
        // 01 (1) Missing genotype
        // 10 (2) Heterozygous
        // 11 (3) Homozygous for second allele (Alt)
        
        let mut byte = 0u8;
        let mut count = 0;
        
        for sample in record.samples().values() {
             let genotype_val = sample.get(key::GENOTYPE).expect("Genotype not found");
             
             let bits = match genotype_val {
                 Some(Value::String(s)) => {
                     match s.as_str() {
                         "0/0" | "0|0" | "0" => 0, // HomRef (00)
                         "0/1" | "0|1" | "1/0" | "1|0" | "1" => 2, // Het (10)
                         "1/1" | "1|1" => 3, // HomAlt (11)
                         "./." | "." | ".|." => 1, // Missing (01)
                         _ => 1, // Treat unknown as missing
                     }
                 },
                 _ => 1, // Missing
             };
             
             // Check haploid '1' case mapping logic
             let corrected_bits = if let Some(Value::String(s)) = genotype_val {
                  if s == "1" { 
                      3 // HomAlt (11)
                  } else if s == "0" {
                      0 // HomRef (00)
                  } else {
                      bits
                  }
             } else { bits };

             // Pack bits: sample 0 in low bits (0-1), etc.
             // byte |= bits << (2 * count)
             byte |= corrected_bits << (2 * count);
             count += 1;
             
             if count == 4 {
                 self.bed.write_all(&[byte])?;
                 byte = 0;
                 count = 0;
             }
        }
        
        // Pad and flush last byte if needed
        if count > 0 {
             self.bed.write_all(&[byte])?;
        }
        
        Ok(())
    }
}
