# convert_genome

`convert_genome` transforms direct-to-consumer (DTC) genotype exports (23andMe, AncestryDNA, etc.) into standards-compliant VCF or BCF files. The tool is optimized for large genomes, supports remote downloads, and embraces modern Rust features for safety and performance.

## Features
- **VCF/BCF conversion** – ingest whitespace-delimited DTC genotype text files and emit VCF or BCF outputs with rich metadata.
- **Parallel pipeline** – Rayon-powered record processing accelerates large datasets on multi-core machines.
- **Thread-safe reference cache** – shared LRU cache avoids repeated FASTA I/O and enables efficient concurrent lookups.
- **Remote input support** – download HTTP/HTTPS inputs and references with transparent gzip/zip decompression.
- **Robust validation** – integration tests, property-based fuzzing, and coverage guard against regressions.
- **Comprehensive tooling** – Criterion benchmarks, Clippy/format enforcement, and GitHub Actions CI/CD keep quality high.

## Installation
The project targets Rust 2024 on the nightly toolchain (configured via `rust-toolchain.toml`). Install directly from the repository:

```bash
cargo install --path .
```

To build the CLI without installing:

```bash
cargo build --release
```

## Usage
Run the CLI to convert a local genotype file:

```bash
convert-genome \
  --input path/to/genome.txt \
  --reference path/to/reference.fa \
  --output sample.vcf \
  --format vcf \
  --sample SAMPLE_ID \
  --assembly GRCh38
```

Remote inputs are supported transparently. For example:

```bash
convert-genome \
  --input https://example.org/genome.zip \
  --reference https://example.org/reference.fa.gz \
  --output sample.bcf \
  --format bcf
```

The converter downloads archives, decompresses gzip/zip content, and streams data without manual staging.

## Performance
Benchmarks (via Criterion) cover reference lookups, DTC parsing throughput, full conversion pipelines, and sequential-versus-parallel comparisons. On typical four-core systems you can expect:

- **4–8× faster** record processing compared to the sequential baseline.
- **Order-of-magnitude improvements** for repeated reference lookups thanks to the shared 128k-entry LRU cache.
- **Reduced allocations** in hot paths by reusing `ReferenceGenome` clones and avoiding redundant string copies.

Run the full benchmark suite locally:

```bash
cargo bench --no-fail-fast
```

## Architecture
The crate is organized into focused modules:

- `cli` – command-line argument parsing and execution wiring.
- `conversion` – core DTC → VCF/BCF pipeline, parallel processing, and reporting.
- `dtc` – streaming parser for genotype text exports.
- `reference` – thread-safe FASTA reader with shared LRU cache and alias handling.
- `remote` – HTTP(S) download and archive extraction utilities.
- `tests` & `benches` – integration tests, property-based fuzzing, and Criterion benchmarks.

`src/lib.rs` re-exports the primary conversion API for programmatic use.

## Supported formats
- **Inputs:** whitespace-delimited genotype text exports with optional rsID columns (23andMe, AncestryDNA, MyHeritage, etc.).
- **Outputs:** VCF (text) and BCF (binary) files compatible with standard tooling. Reference FASTA files are indexed on demand (FAI).

## Testing & Quality Gates
Quality checks are automated in CI (`.github/workflows/ci.yml`) and reproducible locally:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --all-features
cargo test --release --test proptest
cargo bench --no-fail-fast
cargo tarpaulin --out Xml   # Linux coverage (optional)
```

Integration tests cover end-to-end conversions (VCF/BCF, malformed inputs, missing chromosomes) while property-based tests fuzz the DTC parser, genotype handling, and concurrent reference access.

## Contributing
1. Fork and clone the repository.
2. Install the nightly toolchain (automatically handled via `rustup` when entering the project).
3. Run the quality checks listed above before submitting a pull request.
4. Open issues or discussions for substantial changes—performance, new formats, or architectural tweaks.

Contributions that improve performance, broaden format support, or enhance tooling are welcome! Please include relevant tests or benchmarks with every change.
