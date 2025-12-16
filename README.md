# convert_genome

`convert_genome` converts direct-to-consumer (DTC) genotype exports (23andMe, AncestryDNA, MyHeritage, etc.) into standard [VCF](https://samtools.github.io/hts-specs/VCFv4.5.pdf), [BCF](https://samtools.github.io/hts-specs/BCFv2_qref.pdf), or [PLINK](https://www.cog-genomics.org/plink/1.9/formats) binary formats. The converter supports remote references via HTTP(S) and handles compressed `.gz` and `.zip` archives.

## Features

- **Multiple output formats**: VCF text, BCF binary, and PLINK 1.9 binary (.bed/.bim/.fam).
- **Flexible input parsing**: Supports 23andMe, AncestryDNA (5-column), MyHeritage (CSV), deCODEme (6-column with strand flipping), and standard VCF/BCF.
- **Remote reference support**: Fetch `http://` and `https://` URLs with transparent decompression of `.gz` and `.zip` archives.
- **Sex chromosome handling**: Correct ploidy enforcement for X, Y, and MT chromosomes with PAR region awareness.
- **Allele polarization**: Optional `--standardize` mode to normalize alleles against the reference genome.
- **Thread-safe reference access**: Uses an [`LRU`](https://docs.rs/lru) cache for efficient base lookups.
- **Comprehensive CI/CD**: Formatting, linting, testing, and coverage across Linux (with macOS/Windows support).

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

The CLI accepts three positional arguments: `INPUT`, `REFERENCE`, and `OUTPUT`. The `--sex` flag is required for correct sex chromosome ploidy handling.

### Basic conversion to VCF

```bash
convert_genome \
  data/genotypes.txt \
  GRCh38.fa \
  genotypes.vcf \
  --sex male
```

### BCF output with explicit sample ID

```bash
convert_genome \
  https://example.org/sample.txt.gz \
  https://example.org/GRCh38.fa.gz \
  sample.bcf \
  --format bcf \
  --sex female \
  --sample SAMPLE_ID \
  --assembly GRCh38
```

### PLINK output (variants only)

```bash
convert_genome \
  input.txt \
  reference.fa \
  output \
  --format plink \
  --sex male \
  --variants-only
```

If a `.fai` index is not provided, the converter will generate one next to the FASTA automatically.

### Command-line options

| Option | Description |
|--------|-------------|
| `--format <vcf\|bcf\|plink>` | Output format (default: vcf) |
| `--sex <male\|female>` | **Required.** Sample sex for X/Y ploidy |
| `--sample <NAME>` | Sample identifier for VCF header |
| `--assembly <NAME>` | Assembly label for metadata (default: GRCh38) |
| `--variants-only` | Omit reference-only sites from output |
| `--standardize` | Normalize chromosome names and polarize alleles |
| `--input-format <dtc\|vcf\|bcf\|auto>` | Input format (default: auto-detect) |
| `--reference-fai <PATH>` | Explicit FASTA index path |
| `--log-level <LEVEL>` | Logging verbosity (default: info) |

## Performance

Reference lookups use a shared, thread-safe LRU cache. Input records are sorted by genomic position before processing to maximize cache locality.

Run benchmarks with:

```bash
cargo bench
```

Benchmarks cover:

- Cached vs. uncached reference lookups.
- DTC parsing throughput.
- Full conversion pipeline comparisons.

## Testing

Unit, integration, and property-based tests ensure correctness:

```bash
cargo test                # Debug builds, all tests
cargo test --release      # Property tests under release optimizations
```


The `tests/remote_download.rs` file contains integration tests that download real genome reference files from external servers (EBI, Ensembl, Illumina).


## Continuous Integration

See [`.github/workflows/ci.yml`](.github/workflows/ci.yml). The workflow performs:

- Formatting (`cargo fmt --check`)
- Linting (`cargo clippy --all-targets -- -D warnings`)
- Cross-platform builds
- Test suites (debug + release/property)
- Benchmarks (`cargo bench --no-fail-fast`)
- Coverage reporting via `cargo tarpaulin` on Linux

## Project Architecture

### Core modules

- `src/cli.rs` – Argument parsing and top-level command dispatch.
- `src/conversion.rs` – Conversion pipeline, header construction, and record translation.
- `src/dtc.rs` – Parser for DTC genotype exports (23andMe, AncestryDNA, etc.).
- `src/input.rs` – Input source abstraction (DTC, VCF, BCF) with sorting for cache locality.
- `src/reference.rs` – Reference genome loader, contig metadata, and cached base access.
- `src/remote.rs` – Remote fetching with HTTP(S) support and archive extraction.
- `src/plink.rs` – PLINK 1.9 binary format writer (.bed/.bim/.fam).

### Additional resources

- `tests/` – Integration and property-based test suites.
- `benches/` – Criterion benchmarks for core subsystems.

## Contributing

1. Install the nightly toolchain (`rustup toolchain install nightly`).
2. Run formatting and linting before submitting: `cargo fmt` and `cargo clippy --all-targets -- -D warnings`.
3. Execute the full test suite (debug + release) and benchmarks.
4. For large datasets or new reference assemblies, add integration tests with representative fixtures.

Issues and pull requests are welcome! Please include benchmark results when proposing performance-sensitive changes.
