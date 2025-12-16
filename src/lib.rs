pub mod cli;
pub mod conversion;
pub mod dtc;
pub mod input;
pub mod panel;
pub mod plink;
pub mod reference;
pub mod remote;

pub use conversion::{ConversionConfig, OutputFormat, convert_dtc_file};

/// Summary of the conversion process logic.
#[derive(Debug, Default, Clone)]
pub struct ConversionSummary {
    pub total_records: usize,
    pub emitted_records: usize,
    pub variant_records: usize,
    pub reference_records: usize,
    pub missing_genotype_records: usize,
    pub skipped_reference_sites: usize,
    pub unknown_chromosomes: usize,
    pub reference_failures: usize,
    pub invalid_genotypes: usize,
    pub symbolic_allele_records: usize,
    pub parse_errors: usize,
}

impl ConversionSummary {
    pub fn record_emission(&mut self, has_alt: bool) {
        self.emitted_records += 1;
        if has_alt {
            self.variant_records += 1;
        } else {
            self.reference_records += 1;
        }
    }
}
