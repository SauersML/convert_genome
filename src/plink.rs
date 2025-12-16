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
        let stem = if prefix.extension().is_some_and(|e| e == "vcf" || e == "bcf") {
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
        writeln!(self.fam, "{0}\t{0}\t0\t0\t{1}\t-9", sample_id, sex_code)?;
        self.samples_written = true;
        Ok(())
    }

    pub fn write_variant(&mut self, record: &RecordBuf) -> io::Result<()> {
        let alt_bases = record.alternate_bases();

        // If no ALTs, just skip or write as monomorphic?
        // PLINK usually wants at least one ALT or forced dummy.
        // Existing logic handled 1 ALT.
        // Let's handle 0 ALTs by writing a monomorphic site (ALT=".")
        if alt_bases.is_empty() {
            return self.write_biallelic_variant(record, 0, None, ".");
        }

        // Iterate over all ALT alleles (1-based index)
        for (i, alt_allele) in alt_bases.as_ref().iter().enumerate() {
            let alt_idx = i + 1; // 1-based index of the ALT allele we are targeting

            // If multiallelic (more than 1 ALT), append suffix to ID
            let suffix = if alt_bases.len() > 1 {
                Some(format!("_ALT{}", alt_idx))
            } else {
                None
            };

            self.write_biallelic_variant(record, alt_idx, suffix, alt_allele)?;
        }

        Ok(())
    }

    /// Write a single biallelic variant line to BIM and BED
    fn write_biallelic_variant(
        &mut self,
        record: &RecordBuf,
        target_alt_idx: usize,
        id_suffix: Option<String>,
        alt_allele_str: &str,
    ) -> io::Result<()> {
        // Write BIM line
        let chrom = record.reference_sequence_name();

        // Ids
        let mut id = if record.ids().is_empty() {
            ".".to_string()
        } else {
            record
                .ids()
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>()
                .join(";")
        };

        if let Some(suffix) = id_suffix {
            if id == "." {
                // Generate synthetic ID if missing
                let pos = record.variant_start().map(usize::from).unwrap_or(0);
                id = format!("{}:{}{}", chrom, pos, suffix);
            } else {
                id.push_str(&suffix);
            }
        }

        let pos = record.variant_start().map(usize::from).unwrap_or(0);
        let ref_base = record.reference_bases();

        writeln!(
            self.bim,
            "{}\t{}\t0\t{}\t{}\t{}",
            chrom, id, pos, ref_base, alt_allele_str
        )?;

        // Write BED Genotypes
        let mut byte = 0u8;
        let mut count = 0;

        for sample in record.samples().values() {
            let genotype_val = sample.get(key::GENOTYPE);

            // Map genotype bits logic:
            // 00 (0) = HomRef (0 copies of target alt)
            // 01 (1) = Missing
            // 10 (2) = Het (1 copy of target alt)
            // 11 (3) = HomAlt (2 copies of target alt)

            // Map genotype string to PLINK bed bits
            // 00 (0) = HomRef, 01 (1) = Missing, 10 (2) = Het, 11 (3) = HomAlt
            let bits: u8 = match genotype_val {
                Some(Some(Value::String(s))) => {
                    let (a1, a2) = parse_gt_indices(s);
                    // Count how many match target_alt_idx
                    let mut copies = 0;
                    let mut missing = false;

                    // Helper to check allele
                    let check = |idx_opt: Option<usize>| {
                        if let Some(idx) = idx_opt {
                            if idx == target_alt_idx { 1 } else { 0 }
                        } else {
                            // If index is missing (None), is the whole GT missing?
                            // Usually in VCF "./." means both missing.
                            // If "0/." -> one missing. PLINK doesn't support half-calls well.
                            // Treat any missing as Missing GT.
                            2 // Special marker
                        }
                    };

                    let c1 = check(a1);
                    let c2 = check(a2);

                    if c1 == 2 || c2 == 2 {
                        missing = true;
                    } else {
                        copies = c1 + c2;
                    }

                    if missing {
                        1 // Missing (01)
                    } else {
                        match copies {
                            0 => 0, // HomRef (00)
                            1 => 2, // Het (10)
                            2 => 3, // HomAlt (11)
                            _ => 1, // Should not happen for diploid
                        }
                    }
                }
                _ => 1, // Missing (None, None, or not a String)
            };

            byte |= bits << (2 * count);
            count += 1;

            if count == 4 {
                self.bed.write_all(&[byte])?;
                byte = 0;
                count = 0;
            }
        }

        if count > 0 {
            self.bed.write_all(&[byte])?;
        }

        Ok(())
    }
}

/// Helper to parse GT string into two indices.
/// Returns (Option<usize>, Option<usize>). None means missing ('.').
fn parse_gt_indices(gt: &str) -> (Option<usize>, Option<usize>) {
    // Split by '/' or '|'
    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();

    // Handle haploid (1 allele) or diploid (2 alleles)
    match alleles.len() {
        1 => {
            let a1 = alleles[0].parse::<usize>().ok();
            // Create homozygous diploid equivalent for PLINK
            (a1, a1)
        }
        2 => {
            let a1 = alleles[0].parse::<usize>().ok();
            let a2 = alleles[1].parse::<usize>().ok();
            (a1, a2)
        }
        _ => (None, None), // Empty or complex >2 ploidy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::vcf::variant::record_buf::{RecordBuf, Samples};

    #[test]
    fn test_multiallelic_split() -> io::Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let prefix = temp_dir.path().join("test");
        let mut writer = PlinkWriter::new(&prefix)?;

        // Create a record with 2 ALTs: REF=A, ALT=C,G
        // Sample 1: 1/2 (C/G) -> should be Het for site 1 (A/C) and Het for site 2 (A/G)
        // Sample 2: 0/1 (A/C) -> Het for site 1 (A/C), HomRef for site 2 (A/G)
        // Sample 3: 2/2 (G/G) -> HomRef for site 1 (A/C), HomAlt for site 2 (A/G)

        use noodles::vcf::variant::record_buf::AlternateBases;

        // Mock samples
        // noodles 0.x Samples::new takes (Keys, Vec<Vec<Option<Value>>>)
        use noodles::vcf::variant::record_buf::samples::Keys;
        let keys: Keys = vec![String::from("GT")].into_iter().collect();
        let samples = vec![
            vec![Some(Value::String("1/2".into()))],
            vec![Some(Value::String("0/1".into()))],
            vec![Some(Value::String("2/2".into()))],
        ];
        let samples = Samples::new(keys, samples);

        let record = RecordBuf::builder()
            .set_reference_sequence_name("1")
            .set_variant_start(noodles::core::Position::try_from(100).unwrap())
            .set_reference_bases("A")
            .set_alternate_bases(AlternateBases::from(vec!["C".into(), "G".into()]))
            .set_samples(samples)
            .build();

        writer.write_variant(&record)?;

        // Verify BIM output
        drop(writer); // Flush

        let bim_content = std::fs::read_to_string(prefix.with_extension("bim"))?;
        let lines: Vec<&str> = bim_content.lines().collect();
        assert_eq!(lines.len(), 2);

        assert!(lines[0].contains("1:100_ALT1"));
        assert!(lines[0].ends_with("A\tC"));

        assert!(lines[1].contains("1:100_ALT2"));
        assert!(lines[1].ends_with("A\tG"));

        Ok(())
    }
}
