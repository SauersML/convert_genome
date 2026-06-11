use std::{fs, path::PathBuf};

use anyhow::Result;
use assert_fs::{TempDir, prelude::*};
use convert_genome::input::InputFormat;
use convert_genome::reference::ReferenceGenome;
use convert_genome::{ConversionConfig, ConversionSummary, OutputFormat, convert_dtc_file};
use noodles::bcf::io::reader::Builder as BcfReaderBuilder;
use noodles::vcf;
use noodles::vcf::variant::record::samples::Sample;
use rayon::ThreadPoolBuilder;

fn decode_gt_indices_from_value(
    gt: vcf::variant::record::samples::series::Value,
) -> Vec<Option<usize>> {
    match gt {
        vcf::variant::record::samples::series::Value::String(s) => s
            .split(|c| c == '/' || c == '|')
            .map(|tok| {
                if tok == "." || tok.is_empty() {
                    None
                } else {
                    tok.parse::<usize>().ok()
                }
            })
            .collect(),
        vcf::variant::record::samples::series::Value::Genotype(geno) => {
            geno.iter().map(|a| a.ok().and_then(|a| a.0)).collect()
        }
        _ => vec![],
    }
}

fn gt_alleles_for_record(record: &vcf::variant::RecordBuf, header: &vcf::Header) -> Vec<String> {
    let ref_base = record.reference_bases().to_string().to_uppercase();
    let alts: Vec<String> = record
        .alternate_bases()
        .as_ref()
        .iter()
        .map(|s| s.to_string().to_uppercase())
        .collect();

    let samples = record.samples();
    let sample_values = samples.values().next().expect("missing sample");
    let gt_value = sample_values
        .iter(header)
        .next()
        .expect("missing GT")
        .expect("invalid sample value")
        .1
        .expect("missing GT value");

    let indices = decode_gt_indices_from_value(gt_value);
    indices
        .into_iter()
        .filter_map(|idx| match idx {
            None => None,
            Some(0) => Some(ref_base.clone()),
            Some(n) => alts.get(n.saturating_sub(1)).cloned(),
        })
        .collect()
}

fn write_reference(dir: &TempDir) -> Result<PathBuf> {
    let fasta = dir.child("ref.fa");
    fasta.write_str(">chr1\nACGTACGT\n>chr2\nTTTTCCCC\n")?;
    Ok(fasta.path().to_path_buf())
}

fn write_dtc(dir: &TempDir, contents: &str) -> Result<PathBuf> {
    let input = dir.child("input.txt");
    input.write_str(contents)?;
    Ok(input.path().to_path_buf())
}

fn write_genomestudio(dir: &TempDir, data_rows: &str) -> Result<PathBuf> {
    let input = dir.child("report.txt");
    input.write_str(&format!(
        "[Header]\nGSGT Version,2.0.5\nNum Samples,1\n[Data]\n\
         Sample ID,SNP Name,Chr,Position,Allele1 - Plus,Allele2 - Plus\n{data_rows}"
    ))?;
    Ok(input.path().to_path_buf())
}

fn write_panel_vcf(dir: &TempDir) -> Result<PathBuf> {
    let vcf_path = dir.child("panel.vcf");
    vcf_path.write_str(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\tP3\tP4\tP5\n1\t1500\t.\tT\tTA\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n",
    )?;
    Ok(vcf_path.path().to_path_buf())
}

fn write_panel_vcf_single_alt(dir: &TempDir) -> Result<PathBuf> {
    let vcf_path = dir.child("panel_single_alt.vcf");
    vcf_path.write_str(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\tP3\tP4\tP5\n1\t1500\t.\tT\tA\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\n",
    )?;
    Ok(vcf_path.path().to_path_buf())
}

fn write_reference_2000(dir: &TempDir) -> Result<PathBuf> {
    let fasta = dir.child("ref.fa");
    let seq: String = std::iter::repeat('T').take(2000).collect();
    fasta.write_str(&format!(">1\n{}\n", seq))?;
    Ok(fasta.path().to_path_buf())
}

fn write_reference_with_base(dir: &TempDir, base: char, len: usize) -> Result<PathBuf> {
    let fasta = dir.child("ref.fa");
    let seq: String = std::iter::repeat(base).take(len).collect();
    fasta.write_str(&format!(">1\n{}\n", seq))?;
    Ok(fasta.path().to_path_buf())
}

fn base_config(input: PathBuf, reference: PathBuf, output: PathBuf) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_format: convert_genome::input::InputFormat::Dtc,
        input_origin: input.to_string_lossy().to_string(),
        reference_fasta: Some(reference.clone()),
        reference_origin: Some(reference.to_string_lossy().to_string()),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_dir: None,
        output_format: convert_genome::OutputFormat::Vcf,
        sample_id: "SAMPLE".into(),
        assembly: "GRCh38".into(),
        include_reference_sites: true,
        sex: Some(convert_genome::cli::Sex::Female),
        par_boundaries: None,
        standardize: false,
        panel: None,
        input_build: None,
    }
}

fn run_conversion_with_threads(
    config: ConversionConfig,
    threads: usize,
) -> Result<(ConversionSummary, Vec<u8>)> {
    let output = config.output.clone();
    let pool = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let summary = pool.install(|| convert_dtc_file(config))?;
    let data = fs::read(output)?;
    Ok((summary, data))
}

#[test]
fn converts_to_vcf_and_bcf() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    // Input modified to be concordant with reference (>chr1\nACGTACGT\n>chr2\nTTTTCCCC\n)
    // chr1:2 is C. Input CC -> Match.
    // chr1:3 is G. Input AG -> Match (G).
    // chr2:4 is T. Input TT -> Match.
    // Concordance 3/3 = 1.0. Avoids liftover trigger (which fails in sandbox).
    let input = write_dtc(&temp, "rs1\t1\t2\tCC\nrs2\t1\t3\tAG\nrs3\t2\t4\tTT\n")?;

    let vcf_path = temp.child("out.vcf");
    let config = base_config(
        input.clone(),
        reference.clone(),
        vcf_path.path().to_path_buf(),
    );
    let summary = convert_dtc_file(config.clone())?;

    assert_eq!(summary.emitted_records, 3);
    // Updated expectation due to concordant input change:
    // rs1 (CC) is Ref match (was Variant AA)
    // rs2 (AG) is Variant
    // rs3 (TT) is Ref match
    // So 1 Variant, 2 Reference.
    assert_eq!(summary.variant_records, 1);
    assert_eq!(summary.reference_records, 2);

    let vcf_data = fs::read_to_string(vcf_path.path())?;
    assert!(vcf_data.contains("#CHROM\tPOS"));
    assert!(vcf_data.contains("1\t2"));

    let bcf_path = temp.child("out.bcf");
    let mut bcf_config = config;
    bcf_config.output = bcf_path.path().to_path_buf();
    bcf_config.output_format = OutputFormat::Bcf;
    let bcf_summary = convert_dtc_file(bcf_config)?;
    assert_eq!(bcf_summary.emitted_records, 3);

    let mut reader = BcfReaderBuilder::default().build_from_path(bcf_path.path())?;
    reader.read_header()?;
    let mut records = reader.records();
    assert!(records.next().transpose()?.is_some());

    Ok(())
}

#[test]
fn malformed_variant_id_is_skipped_not_crash() -> Result<()> {
    // Regression test for the "failed to spill sorted records / invalid ID"
    // crash: a DTC record whose ID contains a ';' (forbidden in the VCF ID
    // field) used to survive parsing and abort the entire conversion deep in
    // the external-sort spill writer with a bare, unactionable "invalid ID".
    // It must now be rejected at parse time and skipped with a counted parse
    // error, so the surrounding good records still convert cleanly.
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    // rs2 carries an ID with a ';' -> must be skipped; rs1 and rs3 are good.
    let input = write_dtc(
        &temp,
        "rs1\t1\t2\tCC\nrs2;bad\t1\t3\tAG\nrs3\t2\t4\tTT\n",
    )?;

    let vcf_path = temp.child("malformed_id.vcf");
    let config = base_config(input, reference, vcf_path.path().to_path_buf());
    // Must NOT error (previously: Err "failed to spill sorted records").
    let summary = convert_dtc_file(config)?;

    // The two good records are emitted; the malformed one is counted+skipped.
    assert_eq!(summary.emitted_records, 2);
    assert_eq!(summary.parse_errors, 1);

    let vcf_data = fs::read_to_string(vcf_path.path())?;
    assert!(vcf_data.contains("#CHROM\tPOS"));
    // The forbidden ID never reaches the output.
    assert!(!vcf_data.contains("rs2;bad"));

    Ok(())
}

#[test]
fn genomestudio_malformed_snp_name_is_skipped_not_crash() -> Result<()> {
    // Same robustness guarantee as the DTC case, for Illumina GenomeStudio
    // Final Reports: a row whose SNP Name carries a ';' (forbidden in the VCF
    // ID field) must be skipped with a counted parse error rather than aborting
    // the whole conversion in the serializer. Declaring input_build == target
    // keeps this fast (no liftover/network path).
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    // rs1 (chr1:2 ref C) and rs3 (chr2:4 ref T) are good; rs2;bad is rejected.
    let input = write_genomestudio(
        &temp,
        "S1,rs1,1,2,C,C\nS1,rs2;bad,1,3,A,G\nS1,rs3,2,4,T,T\n",
    )?;

    let vcf_path = temp.child("gs_malformed_id.vcf");
    let mut config = base_config(input, reference, vcf_path.path().to_path_buf());
    config.input_format = InputFormat::GenomeStudio;
    config.input_build = Some("GRCh38".to_string());

    let summary = convert_dtc_file(config)?;

    assert_eq!(summary.emitted_records, 2);
    assert_eq!(summary.parse_errors, 1);

    let vcf_data = fs::read_to_string(vcf_path.path())?;
    assert!(!vcf_data.contains("rs2;bad"));

    Ok(())
}

#[test]
fn handles_empty_input() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    let input = write_dtc(&temp, "")?;

    let vcf_path = temp.child("empty.vcf");
    let config = base_config(input, reference, vcf_path.path().to_path_buf());
    let summary = convert_dtc_file(config)?;
    assert_eq!(summary.total_records, 0);
    let vcf_data = fs::read_to_string(vcf_path.path())?;
    assert!(vcf_data.contains("#CHROM"));
    Ok(())
}

#[test]
fn reports_unknown_chromosomes() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    // Add valid records to pass concordance check (>0.7)
    // 3 valid records (concordant), 1 unknown.
    // Concordance = 3/3 = 1.0 (unknown is skipped in check).
    let input = write_dtc(
        &temp,
        "rs1\tUn\t1\tAA\nrs2\t1\t1\tAA\nrs3\t1\t2\tCC\nrs4\t1\t3\tGG\n",
    )?;

    let vcf_path = temp.child("unknown.vcf");
    let config = base_config(input, reference, vcf_path.path().to_path_buf());
    let summary = convert_dtc_file(config)?;
    assert_eq!(summary.unknown_chromosomes, 1);
    assert_eq!(summary.emitted_records, 3);
    Ok(())
}

#[test]
fn parallel_matches_single_thread() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference(&temp)?;
    let input = write_dtc(
        &temp,
        "rs1\t1\t2\tCC\nrs2\t1\t3\tAG\nrs3\t2\t4\tTT\nrs4\t2\t5\tCT\n",
    )?;

    let single_output = temp.child("single.vcf");
    let parallel_output = temp.child("parallel.vcf");

    let mut single_config = base_config(
        input.clone(),
        reference.clone(),
        single_output.path().to_path_buf(),
    );
    // DTC input with default CLI settings runs hg19/hg38 build detection which
    // touches the network. Declare the build to short-circuit it for this test.
    single_config.input_build = Some("GRCh38".into());
    let mut parallel_config = base_config(input, reference, parallel_output.path().to_path_buf());
    parallel_config.input_build = Some("GRCh38".into());

    let (single_summary, single_bytes) = run_conversion_with_threads(single_config, 1)?;
    let (parallel_summary, parallel_bytes) = run_conversion_with_threads(parallel_config, 4)?;

    assert_eq!(
        single_summary.emitted_records,
        parallel_summary.emitted_records
    );
    assert_eq!(single_bytes, parallel_bytes);
    Ok(())
}

/// Memory-bounded streaming regression for the `--standardize` path.
///
/// The dense-WGS OOM (NA12878, ~4M variants, SIGKILL at 14 GiB) came from the
/// old parallel path draining the entire input into a `HashMap<chrom,
/// Vec<RecordBuf>>` in RAM before transforming. The replacement transforms in
/// bounded batches and spills to disk via the external sorter, so live record
/// residency is a constant independent of genome size.
///
/// This test feeds far more records than the external sorter's per-spill
/// `RECORDBUF_CHUNK_SIZE` (50_000) and more than several `TRANSFORM_BATCH_SIZE`
/// (16_384) parallel batches, forcing the spill+merge and multi-batch paths to
/// run. It asserts that the parallel (4-thread) output is byte-identical to the
/// single-threaded output. Byte-identity across the spill/merge boundary is the
/// observable proof that the bounded-memory rewrite preserves output exactly.
#[test]
fn large_standardize_parallel_matches_single_thread_with_spills() -> Result<()> {
    // 120_000 records > 2x RECORDBUF_CHUNK_SIZE (so the sorter spills to
    // multiple temp files and merges) and > 7x TRANSFORM_BATCH_SIZE (so the
    // parallel batch loop iterates many times). Two chromosomes interleaved so
    // the external sort genuinely reorders across chunk boundaries.
    const N: u64 = 120_000;
    let temp = TempDir::new()?;
    // Two-contig reference (chr1 + chr2), each long enough for all positions.
    // All-`A` so every `AA` call standardizes to a reference site deterministically.
    let seq: String = std::iter::repeat('A').take((N as usize) + 16).collect();
    let reference_child = temp.child("ref_large.fa");
    reference_child.write_str(&format!(">1\n{seq}\n>2\n{seq}\n"))?;
    let reference = reference_child.path().to_path_buf();

    let mut dtc = String::new();
    // Emit chr1 and chr2 in descending position order so the input is NOT
    // pre-sorted — this exercises the external sort's reorder path, which is
    // exactly what `--standardize` relies on.
    for pos in (1..=N).rev() {
        dtc.push_str(&format!("rs1_{pos}\t1\t{pos}\tAA\n"));
        dtc.push_str(&format!("rs2_{pos}\t2\t{pos}\tAA\n"));
    }
    let input = write_dtc(&temp, &dtc)?;

    let single_output = temp.child("single_large.vcf");
    let parallel_output = temp.child("parallel_large.vcf");

    let mut single_config = base_config(
        input.clone(),
        reference.clone(),
        single_output.path().to_path_buf(),
    );
    single_config.standardize = true;
    single_config.input_build = Some("GRCh38".into());

    let mut parallel_config = base_config(input, reference, parallel_output.path().to_path_buf());
    parallel_config.standardize = true;
    parallel_config.input_build = Some("GRCh38".into());

    let (single_summary, single_bytes) = run_conversion_with_threads(single_config, 1)?;
    let (parallel_summary, parallel_bytes) = run_conversion_with_threads(parallel_config, 4)?;

    assert_eq!(
        single_summary.emitted_records, parallel_summary.emitted_records,
        "emitted record count diverged between serial and parallel"
    );
    assert!(
        single_summary.emitted_records as u64 >= N,
        "expected at least {N} emitted records, got {}",
        single_summary.emitted_records
    );
    assert_eq!(
        single_bytes, parallel_bytes,
        "parallel bounded-batch output is not byte-identical to serial across spill/merge"
    );
    Ok(())
}

#[test]
fn reference_cache_populates() -> Result<()> {
    let temp = TempDir::new()?;
    let fasta = write_reference(&temp)?;
    let reference = ReferenceGenome::open(&fasta, None)?;
    assert_eq!(reference.cache_len(), 0);
    assert_eq!(reference.base("1", 1)?, 'A');
    assert!(reference.cache_len() >= 1);
    Ok(())
}

#[test]
fn preserves_private_multiallelic_site_with_panel() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference_2000(&temp)?;
    let panel = write_panel_vcf(&temp)?;

    let mut dtc = String::new();
    for pos in 1..=2000u64 {
        if pos == 1500 {
            dtc.push_str(&format!("rs{pos}\t1\t{pos}\tAG\n"));
        } else {
            dtc.push_str(&format!("rs{pos}\t1\t{pos}\tTT\n"));
        }
    }
    let input = write_dtc(&temp, &dtc)?;

    let vcf_path = temp.child("out.vcf");
    let mut config = base_config(input, reference, vcf_path.path().to_path_buf());
    config.standardize = true;
    config.panel = Some(panel);
    convert_dtc_file(config)?;

    let mut reader = vcf::io::reader::Builder::default().build_from_path(vcf_path.path())?;
    let header = reader.read_header()?;

    let mut found = false;
    for result in reader.record_bufs(&header) {
        let record = result?;
        let chrom = record.reference_sequence_name().to_string();
        let pos_raw = record.variant_start().expect("missing pos");
        let pos = usize::from(pos_raw) as u64;
        if chrom == "1" && pos == 1500 {
            found = true;

            assert_eq!(record.reference_bases().to_string().to_uppercase(), "T");

            let alts: Vec<String> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|s| s.to_string().to_uppercase())
                .collect();
            assert_eq!(alts.len(), 2);
            assert!(alts.contains(&"A".to_string()));
            assert!(alts.contains(&"G".to_string()));

            let mut gt_alleles = gt_alleles_for_record(&record, &header);
            gt_alleles.sort();
            assert_eq!(gt_alleles, vec!["A".to_string(), "G".to_string()]);

            break;
        }
    }

    assert!(found, "did not find chr1:1500 in output VCF");
    Ok(())
}

#[test]
fn preserves_one_panel_alt_and_one_private_alt_with_panel() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference_2000(&temp)?;
    let panel = write_panel_vcf_single_alt(&temp)?;

    let mut dtc = String::new();
    for pos in 1..=2000u64 {
        if pos == 1500 {
            dtc.push_str(&format!("rs{pos}\t1\t{pos}\tAG\n"));
        } else {
            dtc.push_str(&format!("rs{pos}\t1\t{pos}\tTT\n"));
        }
    }
    let input = write_dtc(&temp, &dtc)?;

    let vcf_path = temp.child("out_one_private.vcf");
    let mut config = base_config(input, reference, vcf_path.path().to_path_buf());
    config.standardize = true;
    config.panel = Some(panel);
    convert_dtc_file(config)?;

    let mut reader = vcf::io::reader::Builder::default().build_from_path(vcf_path.path())?;
    let header = reader.read_header()?;

    let mut found = false;
    for result in reader.record_bufs(&header) {
        let record = result?;
        let chrom = record.reference_sequence_name().to_string();
        let pos_raw = record.variant_start().expect("missing pos");
        let pos = usize::from(pos_raw) as u64;
        if chrom == "1" && pos == 1500 {
            found = true;

            assert_eq!(record.reference_bases().to_string().to_uppercase(), "T");

            let alts: Vec<String> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|s| s.to_string().to_uppercase())
                .collect();
            assert_eq!(alts.len(), 2);
            assert!(alts.contains(&"A".to_string()));
            assert!(alts.contains(&"G".to_string()));

            let mut gt_alleles = gt_alleles_for_record(&record, &header);
            gt_alleles.sort();
            assert_eq!(gt_alleles, vec!["A".to_string(), "G".to_string()]);

            break;
        }
    }

    assert!(found, "did not find chr1:1500 in output VCF");
    Ok(())
}

#[test]
fn standardize_outputs_sorted_positions() -> Result<()> {
    let temp = TempDir::new()?;
    let reference = write_reference_with_base(&temp, 'A', 2000)?;
    let input = write_dtc(
        &temp,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=2000>\n\
##contig=<ID=chr1,length=2000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
1\t20\t.\tA\tG\t.\t.\t.\tGT\t0/1\n\
chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\n",
    )?;

    let vcf_path = temp.child("out_standardize_sorted.vcf");
    let mut config = base_config(input, reference, vcf_path.path().to_path_buf());
    config.input_format = InputFormat::Vcf;
    config.standardize = true;
    convert_dtc_file(config)?;

    let mut reader = vcf::io::reader::Builder::default().build_from_path(vcf_path.path())?;
    let header = reader.read_header()?;

    let mut last_pos: Option<usize> = None;
    for result in reader.record_bufs(&header) {
        let record = result?;
        assert_eq!(record.reference_sequence_name(), "1");
        let pos = usize::from(record.variant_start().expect("missing pos"));
        if let Some(prev) = last_pos {
            assert!(
                pos >= prev,
                "positions not sorted after standardize: {} then {}",
                prev,
                pos
            );
        }
        last_pos = Some(pos);
    }

    Ok(())
}
