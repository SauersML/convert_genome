pub mod cli;
pub mod conversion;
pub mod dtc;
pub mod harmonize;
pub mod inference;
pub mod input;
pub mod liftover;
pub mod panel;
pub mod panel_writer;
pub mod plink;
pub mod reference;
pub mod remote;
pub mod report;
pub mod smart_reader;
pub mod vcf_utils;

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
    // Liftover stats - specific rejection reasons
    pub liftover_unmapped: usize,       // No chain interval found
    pub liftover_ambiguous: usize,      // Multiple chain intervals (rejected)
    pub liftover_incompatible: usize,   // Alleles incompatible with target reference
    pub liftover_straddled: usize,      // Indel spans multiple intervals
    pub liftover_contig_missing: usize, // Target contig not in reference
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
