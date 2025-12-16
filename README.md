# convert_genome

`convert_genome` converts direct-to-consumer (DTC) genotype exports (23andMe, AncestryDNA, MyHeritage, etc.) into standard [VCF](https://samtools.github.io/hts-specs/VCFv4.5.pdf), [BCF](https://samtools.github.io/hts-specs/BCFv2_qref.pdf), or [PLINK](https://www.cog-genomics.org/plink/1.9/formats) binary formats. The converter supports remote references via HTTP(S) and handles compressed `.gz` and `.zip` archives.

## Features

- **Multiple output formats**: VCF text, BCF binary, and PLINK 1.9 binary (.bed/.bim/.fam).
- **Flexible input parsing**: Supports 23andMe, AncestryDNA (5-column), MyHeritage (CSV), deCODEme (6-column with strand flipping), and standard VCF/BCF.
- **Auto-detection**:
  - **Genome Build**: Automatically detects `GRCh37` vs `GRCh38` vs `hg19`.
  - **Biological Sex**: Infers sex from X/Y chromosome heterozygosity and density.
- **Reference Panel Support**: Harmonize your data against a reference VCF (e.g., 1000 Genomes) to correct strand flips and align alleles.
- **Structured Reporting**: Automatically generates a `<output>_report.json` file with full run metadata, stats, and inference results.
- **Remote reference support**: Fetch `http://` and `https://` URLs with transparent decompression of `.gz` and `.zip` archives.
- **Sex chromosome handling**: Correct ploidy enforcement for X, Y, and MT chromosomes with PAR region awareness.
- **Allele polarization**: Optional `--standardize` mode to normalize alleles against the reference genome.
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

### Basic conversion (auto-detect everything)

```bash
convert_genome \
  --input data/genotypes.txt \
  --reference GRCh38.fa \
  --output genotypes.vcf
```

This will:
1. Detect the genome build (and warn if it differs from the default GRCh38).
2. Infer biological sex from the data.
3. Convert to VCF.
4. Produce `genotypes.vcf` and `genotypes_report.json`.

### Advanced pipeline: Standardize, Harmonize, and Convert

Perform a full imputation-ready conversion in a single pass:

```bash
convert_genome \
  --input sample.txt \
  --reference hg38.fa \
  --standardize \        # Standardize alleles to reference forward strand
  --panel panel.vcf \    # Harmonize against a reference panel (corrects flips)
  --output sample.bcf \  # Output compact binary BCF
  --format bcf
```

### PLINK output

```bash
convert_genome \
  --input input.txt \
  --reference reference.fa \
  --output output \
  --format plink \
  --sex male \          # Override auto-detection
  --variants-only       # Drop reference-matching sites
```

### Command-line options

| Option | Description |
|--------|-------------|
| `--input <PATH>` | **Required.** Input genotype file (DTC, VCF, or BCF) |
| `--reference <PATH>` | **Required.** Reference genome FASTA (local or URL) |
| `--output <PATH>` | **Required.** Output file path (or prefix for PLINK) |
| `--format <vcf\|bcf\|plink>` | Output format (default: vcf) |
| `--sex <male\|female>` | Explicitly set sex (disables auto-inference) |
| `--sample <NAME>` | Sample identifier for VCF header |
| `--assembly <NAME>` | Assembly label for metadata (default: GRCh38) |
| `--panel <PATH>` | VCF panel for harmonization (e.g., 1000G sites) |
| `--standardize` | Standardize alleles to reference forward strand |
| `--variants-only` | Omit reference-only sites from output |
| `--input-format <dtc\|vcf\|bcf\|auto>` | Input format (default: auto-detect) |
| `--reference-fai <PATH>` | Explicit FASTA index path |
| `--log-level <LEVEL>` | Logging verbosity (default: info) |

## Output Report

Every run produces a JSON report alongside the output file (e.g., `sample_report.json`) containing:

```json
{
  "version": "0.1.0",
  "timestamp": "2025-12-16T22:44:00Z",
  "input": { "format": "dtc", ... },
  "output": { "format": "vcf", ... },
  "sample": { "id": "SAMPLE", "sex": "male", "sex_inferred": true },
  "build_detection": { "detected_build": "GRCh38" },
  "standardize": true,
  "panel": {
    "total_sites": 120000,
    "modified_sites": 150,  // Sites flipped/swapped to match panel
    "novel_sites": 25       // Sites not in panel
  },
  "statistics": {
    "total_records": 638234,
    "emitted_records": 612847,
    ...
  }
}
```

## Project Architecture

### Core modules

- `src/cli.rs` – Argument parsing and top-level command dispatch.
- `src/conversion.rs` – Conversion pipeline, report generation, and record translation.
- `src/dtc.rs` – Parser for DTC genotype exports (23andMe, AncestryDNA, etc.).
- `src/imputation.rs` – Logic for sex inference and build detection.
- `src/inference.rs` – Logic for sex inference and build detection.
- `src/harmonize.rs` – Allele harmonization against reference panels.
- `src/panel.rs` – Reference panel loading and management.
- `src/report.rs` – JSON run report generation.
- `src/reference.rs` – Reference genome loader, contig metadata, and cached base access.
- `src/remote.rs` – Remote fetching with HTTP(S) support and archive extraction.
- `src/plink.rs` – PLINK 1.9 binary format writer (.bed/.bim/.fam).

## Contributing

1. Install the nightly toolchain (`rustup toolchain install nightly`).
2. Run formatting and linting before submitting: `cargo fmt` and `cargo clippy --all-targets -- -D warnings`.
3. Execute the full test suite (debug + release) and benchmarks.
