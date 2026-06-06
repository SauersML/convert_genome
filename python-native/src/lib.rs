//! SPIKE: native, in-process pyo3 bindings for `convert_genome`.
//!
//! This crate wraps the existing `convert_genome` library functions and
//! exposes them to Python as native functions that run *in-process*
//! (no fork/exec of the CLI binary, no `_report.json` sidecar round-trip).
//!
//! It deliberately does NOT reimplement any pipeline logic: every function
//! here is a thin adapter over a `pub` function in the `convert_genome`
//! crate. The point of the spike is to measure, not to fork the codebase.

use std::path::PathBuf;

use convert_genome::{
    ConversionConfig, OutputFormat, convert_dtc_file,
    conversion::prescan_dtc_records,
    inference::{
        detect_build_from_dtc, detect_build_from_variant_file,
        infer_sex_detailed_from_variant_file, infer_sex_from_records,
    },
    input::InputFormat,
};

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

fn to_py_err(e: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(format!("{e:#}"))
}

fn parse_output_format(s: &str) -> PyResult<OutputFormat> {
    match s.to_lowercase().as_str() {
        "vcf" => Ok(OutputFormat::Vcf),
        "bcf" => Ok(OutputFormat::Bcf),
        "plink" => Ok(OutputFormat::Plink),
        other => Err(PyValueError::new_err(format!("unknown output format: {other}"))),
    }
}

fn resolve_input_format(path: &std::path::Path, hint: Option<&str>) -> PyResult<InputFormat> {
    let fmt = match hint.map(|h| h.to_lowercase()).as_deref() {
        None | Some("auto") => InputFormat::detect(path),
        Some("dtc") => InputFormat::Dtc,
        Some("vcf") => InputFormat::Vcf,
        Some("bcf") => InputFormat::Bcf,
        Some(other) => {
            return Err(PyValueError::new_err(format!("unknown input format: {other}")));
        }
    };
    Ok(fmt)
}

fn sex_to_str(sex: convert_genome::cli::Sex) -> &'static str {
    match sex {
        convert_genome::cli::Sex::Male => "male",
        convert_genome::cli::Sex::Female => "female",
        convert_genome::cli::Sex::Unknown => "unknown",
    }
}

/// Run a full DTC->VCF/BCF/PLINK conversion *in process* and return the
/// summary statistics as a dict. Unlike the subprocess wrapper, this does
/// not write or parse a `_report.json` sidecar.
#[pyfunction]
#[pyo3(signature = (input, output, *, assembly="GRCh38", input_format=None,
                    output_format="vcf", sample=None, reference=None,
                    reference_fai=None, panel=None, input_build=None,
                    variants_only=false, standardize=false))]
#[allow(clippy::too_many_arguments)]
fn convert(
    py: Python<'_>,
    input: PathBuf,
    output: PathBuf,
    assembly: &str,
    input_format: Option<&str>,
    output_format: &str,
    sample: Option<String>,
    reference: Option<PathBuf>,
    reference_fai: Option<PathBuf>,
    panel: Option<PathBuf>,
    input_build: Option<String>,
    variants_only: bool,
    standardize: bool,
) -> PyResult<Py<PyDict>> {
    let in_fmt = resolve_input_format(&input, input_format)?;
    let out_fmt = parse_output_format(output_format)?;
    let sample_id = sample.unwrap_or_else(|| "sample".to_string());

    let config = ConversionConfig {
        input: input.clone(),
        input_format: in_fmt,
        input_origin: input.to_string_lossy().to_string(),
        reference_fasta: reference.clone(),
        reference_origin: reference.as_ref().map(|p| p.to_string_lossy().to_string()),
        reference_fai: reference_fai.clone(),
        reference_fai_origin: reference_fai.as_ref().map(|p| p.to_string_lossy().to_string()),
        output,
        output_dir: None,
        output_format: out_fmt,
        sample_id,
        assembly: assembly.to_string(),
        include_reference_sites: !variants_only,
        sex: None,
        par_boundaries: convert_genome::reference::ParBoundaries::new(assembly),
        standardize,
        panel,
        input_build,
    };

    // Heavy work: release the GIL so other Python threads can run.
    let summary = py
        .detach(|| convert_dtc_file(config))
        .map_err(to_py_err)?;

    let d = PyDict::new(py);
    d.set_item("total_records", summary.total_records)?;
    d.set_item("emitted_records", summary.emitted_records)?;
    d.set_item("variant_records", summary.variant_records)?;
    d.set_item("reference_records", summary.reference_records)?;
    d.set_item("missing_genotype_records", summary.missing_genotype_records)?;
    d.set_item("skipped_reference_sites", summary.skipped_reference_sites)?;
    d.set_item("unknown_chromosomes", summary.unknown_chromosomes)?;
    d.set_item("reference_failures", summary.reference_failures)?;
    d.set_item("invalid_genotypes", summary.invalid_genotypes)?;
    d.set_item("symbolic_allele_records", summary.symbolic_allele_records)?;
    d.set_item("parse_errors", summary.parse_errors)?;
    Ok(d.into())
}

/// Infer biological sex from a genome file, returning the call plus the
/// underlying evidence metrics (which the subprocess wrapper cannot expose
/// — it never surfaces these from a standalone call).
#[pyfunction]
#[pyo3(signature = (input, *, input_format=None, build="GRCh38"))]
fn infer_sex(
    py: Python<'_>,
    input: PathBuf,
    input_format: Option<&str>,
    build: &str,
) -> PyResult<Py<PyDict>> {
    let in_fmt = resolve_input_format(&input, input_format)?;
    let d = PyDict::new(py);

    match in_fmt {
        InputFormat::Vcf | InputFormat::Bcf => {
            let detail = py
                .detach(|| infer_sex_detailed_from_variant_file(&input, in_fmt, build))
                .map_err(to_py_err)?;
            d.set_item("sex", sex_to_str(detail.sex))?;
            let metrics = PyDict::new(py);
            metrics.set_item("y_genome_density", detail.y_genome_density)?;
            metrics.set_item("x_autosome_het_ratio", detail.x_autosome_het_ratio)?;
            metrics.set_item("composite_sex_index", detail.composite_sex_index)?;
            metrics.set_item("n_autosomes", detail.n_autosomes)?;
            metrics.set_item("n_y_nonpar", detail.n_y_nonpar)?;
            // The algorithm exposes a composite index rather than a single
            // [0,1] "confidence"; surface it as the confidence proxy.
            d.set_item("confidence", detail.composite_sex_index)?;
            d.set_item("metrics", metrics)?;
        }
        InputFormat::Dtc => {
            let records = py
                .detach(|| prescan_dtc_records(&input))
                .map_err(to_py_err)?;
            let sex = infer_sex_from_records(&records, build).map_err(to_py_err)?;
            d.set_item("sex", sex_to_str(sex))?;
            // DTC path in the library returns only the call (no metrics).
            d.set_item("confidence", py.None())?;
            d.set_item("metrics", PyDict::new(py))?;
        }
        InputFormat::Auto => {
            return Err(PyValueError::new_err("could not auto-detect input format"));
        }
    }
    Ok(d.into())
}

/// Detect the genome build (GRCh37/hg19 vs GRCh38/hg38) from a genome file,
/// returning the detected build plus the per-build match rates.
#[pyfunction]
#[pyo3(signature = (input, *, input_format=None))]
fn detect_build(
    py: Python<'_>,
    input: PathBuf,
    input_format: Option<&str>,
) -> PyResult<Py<PyDict>> {
    let in_fmt = resolve_input_format(&input, input_format)?;
    let d = PyDict::new(py);

    let result = match in_fmt {
        InputFormat::Vcf | InputFormat::Bcf => py
            .detach(|| detect_build_from_variant_file(&input, in_fmt))
            .map_err(to_py_err)?,
        InputFormat::Dtc => {
            let records = py
                .detach(|| prescan_dtc_records(&input))
                .map_err(to_py_err)?;
            Some(detect_build_from_dtc(&records).map_err(to_py_err)?)
        }
        InputFormat::Auto => {
            return Err(PyValueError::new_err("could not auto-detect input format"));
        }
    };

    match result {
        Some(r) => {
            d.set_item("build", r.detected_build)?;
            d.set_item("hg19_match_rate", r.hg19_match_rate)?;
            d.set_item("hg38_match_rate", r.hg38_match_rate)?;
            // Confidence: margin of the winning build over the other.
            let conf = (r.hg19_match_rate - r.hg38_match_rate).abs();
            d.set_item("confidence", conf)?;
        }
        None => {
            d.set_item("build", py.None())?;
            d.set_item("confidence", py.None())?;
        }
    }
    Ok(d.into())
}

#[pymodule]
#[pyo3(name = "_native")]
fn native_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(convert, m)?)?;
    m.add_function(wrap_pyfunction!(infer_sex, m)?)?;
    m.add_function(wrap_pyfunction!(detect_build, m)?)?;
    Ok(())
}
