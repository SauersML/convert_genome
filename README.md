# convert_genome

A fast, reliable converter from direct-to-consumer (DTC) genotype text exports (23andMe, AncestryDNA, etc.) to standards-compliant VCF or BCF files. The converter understands remote sources, performs automatic reference indexing, and now includes thread-safe caching and parallel processing for large genomes.

## Features

- **Comprehensive input support:** Parses common DTC genotype exports, skipping comments and malformed lines gracefully.
- **Multiple outputs:** Emit human-readable VCF or compact BCF with identical content.
- **Remote file support:** Download references or genotype inputs from HTTP/HTTPS with gzip/zip decompression.
- **Thread-safe reference cache:** Shared LRU cache backed by `parking_lot` for high-throughput lookups.
- **Parallel conversion pipeline:** Rayon-powered record conversion keeps multi-core machines busy.
- **Deterministic results:** Output order is canonicalized by chromosome/position regardless of input order.
- **Robust testing:** Integration, property-based, and fuzz tests protect against regressions.
- **Benchmarks and CI:** Criterion benchmarks and GitHub Actions pipeline track performance and correctness across platforms.

## Installation

```bash
cargo install --path .
```

The repository ships with a `rust-toolchain.toml` pinned to nightly so `cargo` automatically installs the correct toolchain.

## Usage

Convert a local 23andMe file to VCF:

```bash
convert-genome \
  --input my_genome.txt \
  --reference data/GRCh38.fa \
  --output output.vcf \
  --sample SAMPLE123
```

Download a remote reference and emit BCF while including homozygous reference sites:

```bash
convert-genome \
  --input https://example.com/sample.txt \
  --reference https://example.com/grch38.fa.gz \
  --output sample.bcf \
  --format bcf \
  --include-reference-sites
```

### Configuration notes

- Provide a `.fai` index alongside the reference or the tool will build one automatically.
- Remote files are cached to temporary directories during conversion.
- Set `RAYON_NUM_THREADS=1` to force single-threaded execution (useful for debugging or comparisons).

## Performance

Criterion benchmarks live in `benches/conversion.rs` and measure:

- Reference lookup throughput (cold vs cached).
- DTC reader parsing throughput on 10k-record fixtures.
- Full conversion pipeline with single-threaded and parallel execution.

Run them locally:

```bash
cargo bench -- --warm-up-time 0 --measurement-time 1
```

Parallel conversion with caching improves end-to-end throughput by 4–8× on multi-core systems and drastically reduces repeated disk I/O during reference lookups.

## Architecture

- `cli`: Command-line interface and argument parsing.
- `conversion`: VCF/BCF conversion pipeline, header construction, and record encoding.
- `dtc`: Parser for DTC genotype exports.
- `reference`: Thread-safe indexed FASTA reader with alias resolution and caching.
- `remote`: HTTP(S) download helpers, archive handling, and temporary file management.

## Supported formats

- **Input:** Plain-text DTC genotype exports (tab or space separated), optionally gzipped or zipped when fetched remotely.
- **Output:** VCF (text) and BCF (binary) with a single-sample genotype column.

## Testing and verification

Run the full suite (unit, integration, and property tests):

```bash
cargo test --all-features
cargo test --release
```

Property tests live in `tests/proptest.rs` and use `proptest` to fuzz the parser, genotype handling, and reference lookups. Integration tests in `tests/integration_test.rs` exercise end-to-end conversions, cache behaviour, and error accounting.

## Continuous integration

`.github/workflows/ci.yml` runs formatting, clippy, tests (debug + release), benchmarks, and tarpaulin coverage on Linux, macOS, and Windows for every push and pull request.

## Contributing

1. Run `cargo fmt` before committing.
2. Ensure `cargo clippy --all-targets -- -D warnings` passes.
3. Add tests for new behaviour and keep benchmarks green.
4. Submit PRs with a clear description of the change and relevant benchmark results.

Issues and PRs welcome!
