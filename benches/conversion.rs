use std::{fs, io::Cursor, path::PathBuf};

use convert_genome::{
    conversion::{ConversionConfig, OutputFormat, convert_dtc_file},
    dtc,
    reference::ReferenceGenome,
};
use criterion::{BatchSize, BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use rayon::ThreadPoolBuilder;
use tempfile::{NamedTempFile, tempdir};

fn create_reference(sequence: &str) -> (tempfile::TempDir, PathBuf) {
    let dir = tempdir().unwrap();
    let path = dir.path().join("ref.fa");
    fs::write(&path, format!(">chr1\n{}\n", sequence)).unwrap();
    (dir, path)
}

fn create_dtc_file(dir: &tempfile::TempDir, records: usize) -> PathBuf {
    let path = dir.path().join("input.txt");
    let mut content = String::new();
    for i in 1..=records {
        content.push_str(&format!("rs{0}\t1\t{0}\tAG\n", i));
    }
    fs::write(&path, content).unwrap();
    path
}

fn base_config(
    input: PathBuf,
    reference: &PathBuf,
    output: PathBuf,
    format: OutputFormat,
) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_origin: input.display().to_string(),
        reference_fasta: reference.clone(),
        reference_origin: reference.display().to_string(),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_format: format,
        sample_id: String::from("sample"),
        assembly: String::from(""),
        include_reference_sites: true,
    }
}

fn bench_reference_lookup(c: &mut Criterion) {
    let sequence = "ACGT".repeat(256);
    let (_dir, reference_path) = create_reference(&sequence);
    let positions: Vec<u64> = (1..=512).collect();

    c.bench_function("reference_lookup_uncached", |b| {
        b.iter_batched(
            || ReferenceGenome::open(&reference_path, None).unwrap(),
            |reference| {
                for &pos in &positions {
                    black_box(reference.base("1", pos).unwrap());
                }
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("reference_lookup_cached", |b| {
        let reference = ReferenceGenome::open(&reference_path, None).unwrap();
        for &pos in &positions {
            let _ = reference.base("1", pos).unwrap();
        }
        b.iter(|| {
            for &pos in &positions {
                black_box(reference.base("1", pos).unwrap());
            }
        });
    });
}

fn bench_dtc_parsing(c: &mut Criterion) {
    let mut content = String::new();
    for i in 0..1000 {
        content.push_str(&format!("rs{0}\t1\t{0}\tAG\n", i));
    }
    let data = content.into_bytes();

    c.bench_function("dtc_parsing", |b| {
        b.iter(|| {
            let cursor = Cursor::new(&data);
            let reader = dtc::Reader::new(cursor);
            for result in reader {
                black_box(&result);
            }
        });
    });
}

fn bench_conversion_pipeline(c: &mut Criterion) {
    let (dir, reference_path) = create_reference(&"ACGT".repeat(512));
    let input_path = create_dtc_file(&dir, 1000);
    let dir_path = dir.path().to_path_buf();

    let mut group = c.benchmark_group("conversion_pipeline");
    group.bench_function(BenchmarkId::new("vcf_parallel", 1000), |b| {
        b.iter_batched(
            || {
                let output = NamedTempFile::new_in(&dir_path).unwrap();
                let config = base_config(
                    input_path.clone(),
                    &reference_path,
                    output.path().to_path_buf(),
                    OutputFormat::Vcf,
                );
                (output, config)
            },
            |(output, config)| {
                convert_dtc_file(config).expect("conversion");
                output.close().unwrap();
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function(BenchmarkId::new("bcf_parallel", 1000), |b| {
        b.iter_batched(
            || {
                let output = NamedTempFile::new_in(&dir_path).unwrap();
                let config = base_config(
                    input_path.clone(),
                    &reference_path,
                    output.path().to_path_buf(),
                    OutputFormat::Bcf,
                );
                (output, config)
            },
            |(output, config)| {
                convert_dtc_file(config).expect("conversion");
                output.close().unwrap();
            },
            BatchSize::SmallInput,
        );
    });
    group.finish();
}

fn bench_parallel_vs_sequential(c: &mut Criterion) {
    let (dir, reference_path) = create_reference(&"ACGT".repeat(512));
    let input_path = create_dtc_file(&dir, 2000);
    let dir_path = dir.path().to_path_buf();
    let sequential_pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let parallel_pool = ThreadPoolBuilder::new().num_threads(4).build().unwrap();

    let mut group = c.benchmark_group("parallel_vs_sequential");
    group.bench_function("sequential", |b| {
        b.iter_batched(
            || NamedTempFile::new_in(&dir_path).unwrap(),
            |output| {
                let config = base_config(
                    input_path.clone(),
                    &reference_path,
                    output.path().to_path_buf(),
                    OutputFormat::Vcf,
                );
                sequential_pool
                    .install(|| convert_dtc_file(config).expect("sequential conversion"));
                output.close().unwrap();
            },
            BatchSize::SmallInput,
        );
    });

    group.bench_function("parallel", |b| {
        b.iter_batched(
            || NamedTempFile::new_in(&dir_path).unwrap(),
            |output| {
                let config = base_config(
                    input_path.clone(),
                    &reference_path,
                    output.path().to_path_buf(),
                    OutputFormat::Vcf,
                );
                parallel_pool.install(|| convert_dtc_file(config).expect("parallel conversion"));
                output.close().unwrap();
            },
            BatchSize::SmallInput,
        );
    });
    group.finish();
}

criterion_group!(
    conversion_benches,
    bench_reference_lookup,
    bench_dtc_parsing,
    bench_conversion_pipeline,
    bench_parallel_vs_sequential
);
criterion_main!(conversion_benches);
