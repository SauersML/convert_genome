use std::{fs, io::Write, path::PathBuf};

use convert_genome::dtc;
use convert_genome::reference::ReferenceGenome;
use convert_genome::{ConversionConfig, OutputFormat, convert_dtc_file};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rayon::ThreadPoolBuilder;
use std::hint::black_box;

struct Fixtures {
    dir: tempfile::TempDir,
    reference: PathBuf,
    input: PathBuf,
}

impl Fixtures {
    fn new(record_count: usize) -> Self {
        let dir = tempfile::tempdir().expect("temp dir");
        let reference = dir.path().join("ref.fa");
        let sequence = "ACGT".repeat(4_096);
        fs::write(&reference, format!(">chr1\n{}\n", sequence)).expect("write reference");

        let input = dir.path().join("input.txt");
        let mut file = fs::File::create(&input).expect("create input");
        for i in 0..record_count {
            writeln!(file, "rs{}\t1\t{}\tAC", i, i + 1).expect("write record");
        }

        Self {
            dir,
            reference,
            input,
        }
    }

    fn config(&self, output: PathBuf, format: OutputFormat) -> ConversionConfig {
        ConversionConfig {
            input: self.input.clone(),
            input_format: convert_genome::input::InputFormat::Dtc,
            input_origin: self.input.display().to_string(),
            reference_fasta: self.reference.clone(),
            reference_origin: self.reference.display().to_string(),
            reference_fai: None,
            reference_fai_origin: None,
            output,
            output_dir: None,
            output_format: format,
            sample_id: "BENCH".into(),
            assembly: "GRCh38".into(),
            include_reference_sites: true,
            sex: Some(convert_genome::cli::Sex::Female),
            par_boundaries: None,
            standardize: false,
            panel: None,
        }
    }
}

fn benchmark_reference_lookup(c: &mut Criterion) {
    let fixtures = Fixtures::new(1);
    let cached_reference = ReferenceGenome::open(&fixtures.reference, None).expect("reference");

    c.bench_function("reference_lookup_cached", |b| {
        b.iter(|| cached_reference.base("1", black_box(1000)).unwrap())
    });

    c.bench_function("reference_lookup_uncached", |b| {
        b.iter(|| {
            let reference = ReferenceGenome::open(&fixtures.reference, None).unwrap();
            reference.base("1", black_box(1000)).unwrap()
        })
    });
}

fn benchmark_dtc_parsing(c: &mut Criterion) {
    let fixtures = Fixtures::new(10_000);
    let data = fs::read(&fixtures.input).expect("read input");

    c.bench_function("dtc_parsing", |b| {
        b.iter(|| {
            let cursor = std::io::Cursor::new(&data);
            let reader = dtc::Reader::new(cursor);
            let mut count = 0;
            for record in reader {
                drop(record);
                count += 1;
            }
            black_box(count)
        })
    });
}

fn benchmark_conversion(c: &mut Criterion) {
    let fixtures = Fixtures::new(25_000);
    let mut group = c.benchmark_group("conversion_pipeline");
    group.throughput(Throughput::Elements(25_000_u64));

    group.bench_function(BenchmarkId::new("parallel", 25_000), |b| {
        b.iter(|| {
            let output = fixtures.dir.path().join("parallel.vcf");
            let config = fixtures.config(output.clone(), OutputFormat::Vcf);
            convert_dtc_file(config).expect("convert");
            let _ = fs::remove_file(&output);
        })
    });

    group.bench_function(BenchmarkId::new("single_thread", 25_000), |b| {
        b.iter(|| {
            let output = fixtures.dir.path().join("single.vcf");
            let config = fixtures.config(output.clone(), OutputFormat::Vcf);
            let pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
            pool.install(|| convert_dtc_file(config)).expect("convert");
            let _ = fs::remove_file(&output);
        })
    });

    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(20);
    targets = benchmark_reference_lookup, benchmark_dtc_parsing, benchmark_conversion
);
criterion_main!(benches);
