use std::{
    fs::File,
    io::{BufReader, Write},
    path::{Path, PathBuf},
};

use convert_genome::{
    ConversionConfig, OutputFormat, convert_dtc_file, dtc, reference::ReferenceGenome,
};
use criterion::{BatchSize, Criterion, black_box, criterion_group, criterion_main};
use tempfile::tempdir;

fn write_reference(dir: &Path) -> PathBuf {
    let fasta_path = dir.join("bench.fa");
    let mut file = File::create(&fasta_path).expect("create reference");
    writeln!(file, ">1").unwrap();
    writeln!(file, "ACGTACGTACGTACGT").unwrap();
    writeln!(file, ">2").unwrap();
    writeln!(file, "TTTTCCCCAAAAGGGG").unwrap();
    fasta_path
}

fn write_dtc_records(dir: &Path, count: usize) -> PathBuf {
    let dtc_path = dir.join("input.txt");
    let mut file = File::create(&dtc_path).expect("create dtc");
    for i in 0..count {
        writeln!(file, "rs{0}\t1\t{1}\tAC", i + 1, (i % 8) + 1).unwrap();
    }
    dtc_path
}

fn build_config(
    input: PathBuf,
    reference: PathBuf,
    output: PathBuf,
    format: OutputFormat,
) -> ConversionConfig {
    ConversionConfig {
        input: input.clone(),
        input_origin: input.to_string_lossy().into_owned(),
        reference_fasta: reference.clone(),
        reference_origin: reference.to_string_lossy().into_owned(),
        reference_fai: None,
        reference_fai_origin: None,
        output,
        output_format: format,
        sample_id: String::from("sample"),
        assembly: String::from("GRCh38"),
        include_reference_sites: true,
    }
}

fn bench_reference(c: &mut Criterion) {
    let temp = tempdir().unwrap();
    let fasta = write_reference(temp.path());

    c.bench_function("reference_lookup_cold", |b| {
        b.iter_batched(
            || ReferenceGenome::open(&fasta, None).expect("open reference"),
            |reference| {
                black_box(reference.base("1", 4).unwrap());
            },
            BatchSize::SmallInput,
        );
    });

    let reference = ReferenceGenome::open(&fasta, None).expect("open reference");
    c.bench_function("reference_lookup_cached", |b| {
        b.iter(|| {
            black_box(reference.base("1", 4).unwrap());
        });
    });
}

fn bench_dtc_reader(c: &mut Criterion) {
    let temp = tempdir().unwrap();
    let dtc_path = write_dtc_records(temp.path(), 10_000);

    c.bench_function("dtc_reader_throughput", |b| {
        b.iter(|| {
            let file = File::open(&dtc_path).unwrap();
            let reader = BufReader::new(file);
            let mut dtc_reader = dtc::Reader::new(reader);
            let mut count = 0;
            while let Some(record) = dtc_reader.next() {
                let _ = record.unwrap();
                count += 1;
            }
            black_box(count);
        });
    });
}

fn bench_conversion(c: &mut Criterion) {
    c.bench_function("conversion_parallel", |b| {
        b.iter_batched(
            || {
                let temp = tempdir().unwrap();
                let fasta = write_reference(temp.path());
                let input = write_dtc_records(temp.path(), 5_000);
                let output = temp.path().join("out.vcf");
                let config = build_config(input, fasta, output, OutputFormat::Vcf);
                (temp, config)
            },
            |(temp, config)| {
                unsafe {
                    std::env::set_var("RAYON_NUM_THREADS", "4");
                }
                let summary = convert_dtc_file(config).unwrap();
                black_box(summary.emitted_records);
                unsafe {
                    std::env::remove_var("RAYON_NUM_THREADS");
                }
                drop(temp);
            },
            BatchSize::SmallInput,
        );
    });

    c.bench_function("conversion_single_thread", |b| {
        b.iter_batched(
            || {
                let temp = tempdir().unwrap();
                let fasta = write_reference(temp.path());
                let input = write_dtc_records(temp.path(), 5_000);
                let output = temp.path().join("out.vcf");
                let config = build_config(input, fasta, output, OutputFormat::Vcf);
                (temp, config)
            },
            |(temp, config)| {
                unsafe {
                    std::env::set_var("RAYON_NUM_THREADS", "1");
                }
                let summary = convert_dtc_file(config).unwrap();
                black_box(summary.emitted_records);
                unsafe {
                    std::env::remove_var("RAYON_NUM_THREADS");
                }
                drop(temp);
            },
            BatchSize::SmallInput,
        );
    });
}

criterion_group!(benches, bench_reference, bench_dtc_reader, bench_conversion);
criterion_main!(benches);
