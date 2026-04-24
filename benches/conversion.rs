use std::{fs, io::Write, path::PathBuf};

use convert_genome::dtc;
use convert_genome::input::InputFormat;
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
            reference_fasta: Some(self.reference.clone()),
            reference_origin: Some(self.reference.display().to_string()),
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
            input_build: None,
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
            black_box(reader.count())
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
            let config = fixtures.config(output, OutputFormat::Vcf);
            convert_dtc_file(config).expect("convert")
        })
    });

    group.bench_function(BenchmarkId::new("single_thread", 25_000), |b| {
        b.iter(|| {
            let output = fixtures.dir.path().join("single.vcf");
            let config = fixtures.config(output, OutputFormat::Vcf);
            let pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
            pool.install(|| convert_dtc_file(config)).expect("convert")
        })
    });

    group.finish();
}

struct MultiChromVcfFixtures {
    dir: tempfile::TempDir,
    reference: PathBuf,
    input: PathBuf,
    total_variants: u64,
}

impl MultiChromVcfFixtures {
    /// Produces a deterministic multi-chromosome VCF matching the parity
    /// fixture's shape: repeating ACGT reference, positions spaced 4bp apart,
    /// a mix of 0/0, 0/1, 1/1 genotypes. This is the target of the
    /// parallelization rewrite — every record in a chromosome is independent.
    fn new(chroms: &[(&str, usize)]) -> Self {
        let dir = tempfile::tempdir().expect("temp dir");
        let max_pos = chroms.iter().map(|(_, n)| *n).max().unwrap_or(0);
        let seq_len = (max_pos * 4 + 8).max(32);

        let reference = dir.path().join("ref.fa");
        {
            let mut file = fs::File::create(&reference).expect("create ref");
            let pattern = b"ACGT";
            let mut buf = Vec::with_capacity(seq_len);
            for i in 0..seq_len {
                buf.push(pattern[i % pattern.len()]);
            }
            for (chrom, _) in chroms {
                writeln!(file, ">{}", chrom).expect("write header");
                file.write_all(&buf).expect("write seq");
                writeln!(file).expect("newline");
            }
        }

        let input = dir.path().join("input.vcf");
        let mut total = 0u64;
        {
            let mut file = fs::File::create(&input).expect("create input");
            writeln!(file, "##fileformat=VCFv4.2").expect("header");
            for (chrom, count) in chroms {
                writeln!(
                    file,
                    "##contig=<ID={},length={}>",
                    chrom,
                    count * 4 + 8
                )
                .expect("contig");
            }
            writeln!(
                file,
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
            )
            .expect("cols");
            let alts = ['C', 'G', 'T'];
            for (chrom, count) in chroms {
                for i in 0..*count {
                    let pos = i * 4 + 1;
                    let alt = alts[i % alts.len()];
                    let gt = if i % 7 == 0 {
                        "0/0"
                    } else if i % 5 == 0 {
                        "0/1"
                    } else {
                        "1/1"
                    };
                    writeln!(
                        file,
                        "{}\t{}\t.\tA\t{}\t.\t.\t.\tGT\t{}",
                        chrom, pos, alt, gt
                    )
                    .expect("record");
                    total += 1;
                }
            }
        }

        Self {
            dir,
            reference,
            input,
            total_variants: total,
        }
    }

    fn config(&self, output: PathBuf) -> ConversionConfig {
        ConversionConfig {
            input: self.input.clone(),
            input_format: InputFormat::Vcf,
            input_origin: self.input.display().to_string(),
            reference_fasta: Some(self.reference.clone()),
            reference_origin: Some(self.reference.display().to_string()),
            reference_fai: None,
            reference_fai_origin: None,
            output,
            output_dir: None,
            output_format: OutputFormat::Vcf,
            sample_id: "BENCH".into(),
            assembly: "GRCh38".into(),
            include_reference_sites: true,
            sex: Some(convert_genome::cli::Sex::Female),
            par_boundaries: convert_genome::reference::ParBoundaries::new("GRCh38"),
            standardize: false,
            panel: None,
            // Declared build skips network-based detection and liftover.
            input_build: Some("GRCh38".into()),
        }
    }
}

fn benchmark_multichrom_vcf_throughput(c: &mut Criterion) {
    // ~14k variants across 6 chromosomes (including chrX) — same shape as
    // the parity fixture. Small enough to bench in seconds, large enough
    // that per-chromosome parallelism has something to chew on.
    let chroms: &[(&str, usize)] = &[
        ("1", 4_000),
        ("2", 3_000),
        ("3", 2_500),
        ("4", 2_000),
        ("5", 1_500),
        ("X", 1_000),
    ];
    let fixtures = MultiChromVcfFixtures::new(chroms);

    let mut group = c.benchmark_group("multichrom_vcf");
    group.throughput(Throughput::Elements(fixtures.total_variants));
    group.sample_size(15);

    group.bench_function(BenchmarkId::new("default_pool", fixtures.total_variants), |b| {
        b.iter(|| {
            let output = fixtures.dir.path().join("default.vcf");
            let config = fixtures.config(output);
            convert_dtc_file(config).expect("convert")
        })
    });

    group.bench_function(BenchmarkId::new("single_thread", fixtures.total_variants), |b| {
        b.iter(|| {
            let output = fixtures.dir.path().join("single.vcf");
            let config = fixtures.config(output);
            let pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
            pool.install(|| convert_dtc_file(config)).expect("convert")
        })
    });

    group.finish();
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(20);
    targets =
        benchmark_reference_lookup,
        benchmark_dtc_parsing,
        benchmark_conversion,
        benchmark_multichrom_vcf_throughput
);
criterion_main!(benches);
