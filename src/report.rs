//! Structured run report for downstream tool consumption.
//!
//! Writes a JSON file alongside the output containing all metadata
//! about the conversion run: inferred sex, detected build, statistics, etc.

use serde::Serialize;
use std::path::Path;

use crate::ConversionSummary;
use crate::cli::Sex;
use crate::conversion::OutputFormat;
use crate::input::InputFormat;

/// Complete report of a conversion run.
/// Serialized to JSON alongside the output file.
#[derive(Debug, Clone, Serialize)]
pub struct RunReport {
    /// Tool version
    pub version: String,
    /// Timestamp of run (ISO 8601)
    pub timestamp: String,

    /// Input/output configuration
    pub input: InputInfo,
    pub output: OutputInfo,
    pub reference: ReferenceInfo,

    /// Whether allele standardization was enabled
    pub standardize: bool,

    /// Reference panel if used
    #[serde(skip_serializing_if = "Option::is_none")]
    pub panel: Option<PanelInfo>,

    /// Inferred or provided sample metadata
    pub sample: SampleInfo,

    /// Build detection results
    #[serde(skip_serializing_if = "Option::is_none")]
    pub build_detection: Option<BuildDetection>,

    /// Conversion statistics
    pub statistics: Statistics,
}

#[derive(Debug, Clone, Serialize)]
pub struct InputInfo {
    pub path: String,
    pub format: String,
    pub origin: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct OutputInfo {
    pub path: String,
    pub format: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReferenceInfo {
    pub path: String,
    pub origin: String,
    pub assembly: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct PanelInfo {
    pub path: String,
    pub total_sites: usize,
    pub modified_sites: usize,
    pub novel_sites: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct SampleInfo {
    pub id: String,
    pub sex: String,
    pub sex_inferred: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct BuildDetection {
    pub detected_build: String,
    pub hg19_match_rate: f64,
    pub hg38_match_rate: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct Statistics {
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

impl From<&ConversionSummary> for Statistics {
    fn from(s: &ConversionSummary) -> Self {
        Statistics {
            total_records: s.total_records,
            emitted_records: s.emitted_records,
            variant_records: s.variant_records,
            reference_records: s.reference_records,
            missing_genotype_records: s.missing_genotype_records,
            skipped_reference_sites: s.skipped_reference_sites,
            unknown_chromosomes: s.unknown_chromosomes,
            reference_failures: s.reference_failures,
            invalid_genotypes: s.invalid_genotypes,
            symbolic_allele_records: s.symbolic_allele_records,
            parse_errors: s.parse_errors,
        }
    }
}

impl RunReport {
    /// Write the report as JSON to a file alongside the output.
    /// For output.vcf, writes output_report.json
    pub fn write(&self, output_path: &Path) -> std::io::Result<()> {
        let stem = output_path
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy();
        let report_name = format!("{}_report.json", stem);
        let report_path = output_path.with_file_name(report_name);

        let json = serde_json::to_string_pretty(self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        std::fs::write(&report_path, json)?;
        tracing::info!("Wrote run report to {}", report_path.display());

        Ok(())
    }
}

/// Builder for constructing a RunReport during conversion.
#[derive(Debug, Default)]
pub struct RunReportBuilder {
    pub input_path: String,
    pub input_format: Option<InputFormat>,
    pub input_origin: String,
    pub output_path: String,
    pub output_format: Option<OutputFormat>,
    pub reference_path: String,
    pub reference_origin: String,
    pub assembly: String,
    pub standardize: bool,
    pub panel: Option<PanelInfo>,
    pub sample_id: String,
    pub sex: Option<Sex>,
    pub sex_inferred: bool,
    pub build_detection: Option<BuildDetection>,
}

impl RunReportBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn build(self, summary: &ConversionSummary) -> RunReport {
        let now = time::OffsetDateTime::now_utc();
        let timestamp = now
            .format(&time::format_description::well_known::Rfc3339)
            .unwrap_or_else(|_| "unknown".to_string());

        RunReport {
            version: env!("CARGO_PKG_VERSION").to_string(),
            timestamp,
            input: InputInfo {
                path: self.input_path,
                format: format_name(self.input_format),
                origin: self.input_origin,
            },
            output: OutputInfo {
                path: self.output_path,
                format: output_format_name(self.output_format),
            },
            reference: ReferenceInfo {
                path: self.reference_path,
                origin: self.reference_origin,
                assembly: self.assembly,
            },
            standardize: self.standardize,
            panel: self.panel,
            sample: SampleInfo {
                id: self.sample_id,
                sex: sex_name(self.sex),
                sex_inferred: self.sex_inferred,
            },
            build_detection: self.build_detection,
            statistics: Statistics::from(summary),
        }
    }
}

fn format_name(f: Option<InputFormat>) -> String {
    match f {
        Some(InputFormat::Dtc) => "dtc".to_string(),
        Some(InputFormat::Vcf) => "vcf".to_string(),
        Some(InputFormat::Bcf) => "bcf".to_string(),
        Some(InputFormat::Auto) | None => "auto".to_string(),
    }
}

fn output_format_name(f: Option<OutputFormat>) -> String {
    match f {
        Some(OutputFormat::Vcf) => "vcf".to_string(),
        Some(OutputFormat::Bcf) => "bcf".to_string(),
        Some(OutputFormat::Plink) => "plink".to_string(),
        None => "vcf".to_string(),
    }
}

fn sex_name(s: Option<Sex>) -> String {
    match s {
        Some(Sex::Male) => "male".to_string(),
        Some(Sex::Female) | None => "female".to_string(),
    }
}
