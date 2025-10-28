# convert_genome

`convert_genome` converts direct-to-consumer (DTC) genotype exports (23andMe, AncestryDNA, etc.) into standard [VCF](https://samtools.github.io/hts-specs/VCFv4.5.pdf) or [BCF](https://samtools.github.io/hts-specs/BCFv2_qref.pdf) files. The converter understands remote references, streams compressed archives, and now includes high-performance, parallel processing suitable for multi-million record datasets.

## Features

- **Parallel conversion pipeline** powered by [`rayon`](https://docs.rs/rayon), scaling across all cores.
- **Thread-safe reference genome access** with a shared [`LRU`](https://docs.rs/lru) cache for rapid base lookups.
- **Remote reference support** for `http://` and `https://` URLs with transparent decompression of `.gz` and `.zip` archives.
- **Robust parsing** with property-based and integration tests covering malformed input, missing fields, and concurrent access.
- **Benchmark suite** (Criterion) to track performance of parsing, reference lookups, and pipeline throughput.
- **Comprehensive CI/CD** across Linux, macOS, and Windows with formatting, linting, testing, coverage, and benchmarks.

## Installation

The project targets Rust **nightly** (see `rust-toolchain.toml`). Install the converter directly from the repository:

```bash
cargo install --path .
```

Alternatively, build the binary without installing:

```bash
cargo build --release
```

The resulting executable lives at `target/release/convert_genome`.

## Usage

The CLI accepts both local files and remote resources. A minimal invocation converts a DTC file to VCF:

```bash
convert_genome \
  --input data/genotypes.txt \
  --reference GRCh38.fa \
  --output genotypes.vcf \
  --sample SAMPLE_ID
```

Generate a BCF with explicit assembly metadata and keep homozygous reference calls:

```bash
convert_genome \
  --input https://example.org/sample.txt.gz \
  --reference https://example.org/GRCh38.fa.gz \
  --output sample.bcf \
  --output-format bcf \
  --assembly GRCh38 \
  --sample SAMPLE_ID \
  --include-reference-sites
```

If a `.fai` index is not provided the converter will generate one next to the FASTA automatically.

## Performance

Reference lookups use a shared, thread-safe LRU cache sized for 128k entries, dramatically reducing random I/O. The conversion pipeline collects DTC records, sorts them for cache locality, and processes them in parallel; results are written sequentially to keep deterministic ordering.

The Criterion benchmarks can be executed with:

```bash
cargo bench
```

Benchmarks cover:

- Cached vs. uncached reference lookups.
- DTC parsing throughput.
- Full conversion pipeline comparisons (parallel vs. single-threaded execution).

## Testing

Unit, integration, and property-based tests ensure correctness across a wide surface area:

```bash
cargo test              # Debug builds, all tests
cargo test --release    # Property tests under release optimizations
```

Ignored integration tests in `tests/remote_download.rs` exercise real-world genome downloads; run them manually as needed.

## Continuous Integration

See [`.github/workflows/ci.yml`](.github/workflows/ci.yml). The workflow performs:

- Formatting (`cargo fmt --check`)
- Linting (`cargo clippy --all-targets -- -D warnings`)
- Cross-platform builds
- Test suites (debug + release/property)
- Benchmarks (`cargo bench --no-fail-fast`)
- Coverage reporting via `cargo tarpaulin` on Linux

## Project Architecture

- `src/cli.rs` – Argument parsing and top-level command dispatch.
- `src/conversion.rs` – Conversion pipeline, header construction, and record translation.
- `src/dtc.rs` – Streaming parser for DTC genotype exports.
- `src/reference.rs` – Reference genome loader, contig metadata, and cached base access.
- `src/remote.rs` – Remote fetching with HTTP(S) support and archive extraction.

Additional resources:

- `tests/` – Integration and property-based test suites.
- `benches/` – Criterion benchmarks for core subsystems.

## Contributing

1. Install the nightly toolchain (`rustup toolchain install nightly`).
2. Run formatting and linting before submitting: `cargo fmt` and `cargo clippy --all-targets -- -D warnings`.
3. Execute the full test suite (debug + release) and benchmarks.
4. For large datasets or new reference assemblies, add integration tests with representative fixtures.

Issues and pull requests are welcome! Please include benchmark results when proposing performance-sensitive changes.
