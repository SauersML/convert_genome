use anyhow::{Context, Result, anyhow}; // Added anyhow
use assert_fs::{TempDir, prelude::*};
use convert_genome::{ConversionConfig, convert_dtc_file};
use noodles::vcf::variant::record::samples::Sample; // Needed for iter() on record_buf::Sample
use std::io::BufRead;
use std::fs;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;

fn get_data_files() -> Vec<PathBuf> {
    let mut files = Vec::new();
    let data_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("data");

    for entry in WalkDir::new(data_dir).into_iter().filter_map(|e| e.ok()) {
        let path = entry.path();
        if path.is_file() && path.extension().is_some_and(|ext| ext == "txt") {
            files.push(path.to_path_buf());
        }
    }
    files.sort();
    files
}

fn base_config(input: PathBuf, output: PathBuf, assembly: &str) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_format: convert_genome::input::InputFormat::Dtc,
        input_origin: input.to_string_lossy().to_string(),
        reference_fasta: None, // Will be set in tests if needed
        reference_origin: None,
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_dir: None,
        output_format: convert_genome::OutputFormat::Vcf,
        sample_id: "SAMPLE".into(),
        assembly: assembly.into(),
        include_reference_sites: true,
        sex: None, // Auto-detect
        par_boundaries: None,
        standardize: false,
        panel: None,
    }
}

#[test]
fn test_build_detection_hg19() -> Result<()> {
    let temp = TempDir::new()?;
    let files = get_data_files();
    assert!(!files.is_empty(), "No data files found in data/ directory");

    for file in &files {
        let file_name = file.file_name().unwrap().to_string_lossy();
        println!("Testing build detection for {}", file_name);

        let output = temp.child(format!("{}.vcf", file_name));

        // Ensure we have a reference. For this test, we can use hg19 (GRCh37).
        let hg19_ref = convert_genome::source_ref::load_source_reference("GRCh37")
            .with_context(|| "Failed to load hg19 reference for test")?;

        let mut config = base_config(file.clone(), output.path().to_path_buf(), "GRCh37");
        config.reference_fasta = Some(PathBuf::from(hg19_ref.path()));

        convert_dtc_file(config)?;

        // Read the report
        let report_path = output.path().with_file_name(format!(
            "{}_report.json",
            output.path().file_stem().unwrap().to_string_lossy()
        ));
        let report_json = fs::read_to_string(&report_path)
            .with_context(|| format!("Failed to read report at {}", report_path.display()))?;
        let report: serde_json::Value = serde_json::from_str(&report_json)?;

        let detected_build = report
            .get("build_detection")
            .and_then(|bd| bd.get("detected_build"))
            .and_then(|s| s.as_str())
            .ok_or_else(|| anyhow!("Report missing build_detection.detected_build"))?;

        println!("  Detected build: {}", detected_build);

        let is_hg19 = detected_build.contains("37") || detected_build.contains("hg19");
        assert!(
            is_hg19,
            "File {} detected as {}, expected hg19/GRCh37",
            file_name, detected_build
        );
    }

    Ok(())
}

#[test]
fn test_liftover_fidelity_hg38() -> Result<()> {
    let temp = TempDir::new()?;
    let files = get_data_files();
    assert!(!files.is_empty(), "No data files found in data/ directory");

    // Load hg38 reference for verification final check
    let hg38_ref = convert_genome::source_ref::load_source_reference("GRCh38")
        .with_context(|| "Failed to load hg38 reference for validation")?;

    for file in &files {
        let file_name = file.file_name().unwrap().to_string_lossy();
        println!("Testing liftover fidelity for {}", file_name);

        let output = temp.child(format!("{}_hg38.vcf", file_name));

        let mut config = base_config(file.clone(), output.path().to_path_buf(), "GRCh38");
        config.reference_fasta = Some(PathBuf::from(hg38_ref.path()));

        let summary = convert_dtc_file(config)?;

        // 1. Retention Check: > 80% of variants retained
        let retention = summary.emitted_records as f64 / summary.total_records as f64;
        println!(
            "  Retention: {:.1}% ({}/{})",
            retention * 100.0,
            summary.emitted_records,
            summary.total_records
        );

        assert!(
            retention > 0.80,
            "Retention {:.1}% is below 80% threshold for {}",
            retention * 100.0,
            file_name
        );

        // 2. Reference Match Check: > 75% of output records match hg38 reference base

        use noodles::vcf;
        let mut reader = vcf::io::reader::Builder::default().build_from_path(output.path())?;
        let header = reader.read_header()?;

        let mut matching_genotypes = 0;
        let mut total_checks = 0;
        let mut ref_base_matches = 0;

        // Use record_bufs for easier sample access
        let mut last_seen: Option<(String, usize)> = None;
        let output_path = output.path().to_path_buf();
        for (idx, result) in reader.record_bufs(&header).enumerate() {
            let record = match result {
                Ok(record) => record,
                Err(e) => {
                    println!(
                        "  Record parse failed at index {} (last seen {:?}): {}",
                        idx + 1,
                        last_seen.as_ref().map(|(c, p)| format!("{c}:{p}")),
                        e
                    );
                    println!("  Output file: {}", output_path.display());
                    dump_record_line(&output_path, idx + 1);
                    return Err(e.into());
                }
            };

            // Check REF base integrity (sanity check)
            let chrom = record.reference_sequence_name().to_string();
            // Position conversion: record.variant_start() returns Option<Position>
            // Position can be converted to usize (via From).
            let pos_raw = record
                .variant_start()
                .ok_or_else(|| anyhow!("missing pos"))?;
            let pos = usize::from(pos_raw);
            last_seen = Some((chrom.clone(), pos));

            // hg38_ref lookup
            if let Ok(ref_base) = hg38_ref.base(&chrom, pos as u64)
                && record.reference_bases().to_string().to_uppercase()
                    == ref_base.to_string().to_uppercase()
            {
                ref_base_matches += 1;
            }

            // Check Genotype
            let samples = record.samples();
            // RecordBuf Samples: `values()` returns iterator over samples' values.
            // Each item is `&Vec<Option<Value>>`.
            // We want values for the FIRST sample.

            if let Some(sample_values) = samples.values().next() {
                // sample_values is `&Vec<Option<Value>>` (wrapped in Sample struct).
                // First field usually GT.
                // iter(&header) returns Iterator<Item = io::Result<(&str, Option<Value>)>>

                if let Some(Ok((_, Some(value)))) = sample_values.iter(&header).next() {
                    match value {
                        noodles::vcf::variant::record::samples::series::Value::String(s) => {
                            if s == "0/0" || s == "0|0" || s == "0" {
                                matching_genotypes += 1;
                            }
                        }
                        noodles::vcf::variant::record::samples::series::Value::Genotype(gt) => {
                            // Check if all alleles are 0 (Ref) and we have at least one allele
                            let is_hom_ref = gt.iter().count() > 0
                                && gt.iter().all(|res| match res {
                                    Ok(allele) => allele.0 == Some(0),
                                    Err(_) => false,
                                });

                            if is_hom_ref {
                                matching_genotypes += 1;
                            }
                        }
                        _ => {}
                    }
                }
            }

            total_checks += 1;
        }

        let genotype_concordance = matching_genotypes as f64 / total_checks as f64;
        let ref_validity = ref_base_matches as f64 / total_checks as f64;

        println!(
            "  Genotype Concordance (Match Reference): {:.1}%",
            genotype_concordance * 100.0
        );
        println!("  Ref Base Validity: {:.1}%", ref_validity * 100.0);

        assert!(
            genotype_concordance > 0.50,
            "Genotype concordance {:.1}% is below 50% threshold for {}",
            genotype_concordance * 100.0,
            file_name
        );
    }

    Ok(())
}

fn dump_record_line(path: &Path, record_index: usize) {
    let mut prev = None;
    let mut curr = None;
    let mut next = None;
    let mut seen = 0usize;

    if let Ok(file) = fs::File::open(path) {
        let reader = std::io::BufReader::new(file);
        for line in reader.lines().flatten() {
            if line.starts_with('#') {
                continue;
            }
            seen += 1;
            if seen == record_index.saturating_sub(1) {
                prev = Some(line);
            } else if seen == record_index {
                curr = Some(line);
            } else if seen == record_index + 1 {
                next = Some(line);
                break;
            }
        }
    }

    println!("  Record context for index {}:", record_index);
    if let Some(line) = prev {
        println!("    prev: {}", line);
    }
    if let Some(line) = curr {
        let cols: Vec<&str> = line.split('\t').collect();
        println!("    curr: {}", line);
        println!("    curr columns: {}", cols.len());
        if cols.len() >= 9 {
            println!("    curr FORMAT: {:?}", cols[8]);
            if cols.len() >= 10 {
                println!("    curr SAMPLE: {:?}", cols[9]);
                println!("    curr SAMPLE bytes: {:?}", cols[9].as_bytes());
            }
        }
    } else {
        println!("    curr: <not found>");
    }
    if let Some(line) = next {
        println!("    next: {}", line);
    }
}
