# convert_genome

[![Crates.io](https://img.shields.io/crates/v/convert_genome.svg)](https://crates.io/crates/convert_genome)
[![docs.rs](https://img.shields.io/docsrs/convert_genome)](https://docs.rs/convert_genome)
[![CI](https://github.com/SauersML/convert_genome/actions/workflows/ci.yml/badge.svg)](https://github.com/SauersML/convert_genome/actions/workflows/ci.yml)

`convert_genome` converts direct-to-consumer (DTC) genotype exports (23andMe, AncestryDNA, MyHeritage, etc.) into standard [VCF](https://samtools.github.io/hts-specs/VCFv4.5.pdf), [BCF](https://samtools.github.io/hts-specs/BCFv2_qref.pdf), or [PLINK](https://www.cog-genomics.org/plink/1.9/formats) binary formats. The converter supports remote references via HTTP(S) and handles compressed `.gz` and `.zip` archives.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
  - [Automatic Install (Recommended)](#automatic-install-recommended)
  - [Manual installation](#manual-installation)
- [CLI Usage](#cli-usage)
  - [Supported Inputs](#supported-inputs)
  - [Workflow Examples](#workflow-examples)
    - [Simple Format Conversion](#simple-format-conversion)
    - [Preparing for Imputation](#preparing-for-imputation)
    - [Handling Ancient/Old Data](#handling-ancientold-data)
  - [Flag Deep Dives](#flag-deep-dives)
    - [`--standardize`](#--standardize)
    - [`--panel`](#--panel)
    - [`--sex`](#--sex)
    - [`--assembly`](#--assembly)
- [Using as a Rust Library](#using-as-a-rust-library)
  - [Add the Dependency](#add-the-dependency)
  - [Core Types and Entry Points](#core-types-and-entry-points)
  - [Minimal Example](#minimal-example)
- [Genome Build Detection \& Liftover](#genome-build-detection--liftover)
  - [How does it know?](#how-does-it-know)
  - [Liftover Details](#liftover-details)
- [How It Works (Data Processing Logic)](#how-it-works-data-processing-logic)
  - [Strand Flipping Logic](#strand-flipping-logic)
  - [Ploidy Enforcement](#ploidy-enforcement)
  - [Build Detection](#build-detection)
- [Output Report](#output-report)
  - [Why the report matters](#why-the-report-matters)
  - [Key field definitions](#key-field-definitions)
- [Project Architecture](#project-architecture)
- [Contributing](#contributing)

## Features

- **Multiple output formats**: VCF text, BCF binary, and PLINK 1.9 binary (.bed/.bim/.fam).
- **Flexible input parsing**: Supports 23andMe, AncestryDNA (5-column), MyHeritage (CSV), deCODEme (6-column with strand flipping), and standard VCF/BCF.
- **Auto-detection**:
  - **Genome Build**: Automatically detects `GRCh37` vs `GRCh38` vs `hg19`.
  - **Biological Sex**: Infers sex from X/Y chromosome heterozygosity and density.
- **Automatic Liftover**: Seamlessly converts data between genome builds (e.g., `hg19` to `GRCh38`) by automatically downloading UCSC chain files when a build mismatch is detected.
- **Reference Panel Support**: Harmonize your data against a reference VCF (e.g., 1000 Genomes) to correct strand flips and align alleles.
- **Structured Reporting**: Automatically generates a `<output>_report.json` file with full run metadata, stats, and inference results.
- **Remote reference support**: Fetch `http://` and `https://` URLs with transparent decompression of `.gz` and `.zip` archives.
- **Sex chromosome handling**: Correct ploidy enforcement for X, Y, and MT chromosomes with PAR region awareness.
- **Allele polarization**: Optional `--standardize` mode to normalize alleles against the reference genome.
- **Comprehensive CI/CD**: Formatting, linting, testing, and coverage across Linux (with macOS/Windows support).

## Installation

### Automatic Install (Recommended)
Installs the latest binary for your platform (macOS/Linux/Windows):
```bash
# macOS / Linux / Windows (Git Bash)
curl -fsSL https://raw.githubusercontent.com/SauersML/convert_genome/main/install.sh | bash
```

## Manual installation

The project targets Rust **nightly** (see `rust-toolchain.toml`). Install the converter directly from the repository:

```bash
cargo install --path .
```

Alternatively, build the binary without installing:

```bash
cargo build --release
```

The resulting executable lives at `target/release/convert_genome`.

## CLI Usage

### Supported Inputs

- **DTC genotype exports**
  - **23andMe**
  - **AncestryDNA** (5-column)
  - **MyHeritage** (CSV)
  - **FTDNA** (CSV / whitespace-delimited, depending on export)
  - **deCODEme** (6-column with strand indicator)
- **VCF / BCF** (for format conversion or downstream harmonization)

Compressed inputs are supported directly:

- **Gzip**: `.gz`
- **Zip**: `.zip`

### Workflow Examples

#### Simple Format Conversion

Convert a DTC genotype export to VCF without extra harmonization steps:

```bash
convert_genome \
  --input data/genotypes.txt \
  --reference GRCh38.fa \
  --output genotypes.vcf
```

This will:

1. Detect the genome build.
2. Infer biological sex from the data.
3. Convert to VCF.
4. Produce `genotypes.vcf` and `genotypes_report.json`.

#### Preparing for Imputation

Imputation tools (Beagle, Shapeit, Eagle, Minimac, etc.) tend to be strict about:

- reference/alternate allele consistency
- strand/polarity
- expected ploidy on sex chromosomes
- compactness and indexing (BCF is usually preferred)

An “imputation-ready” run typically includes `--standardize` plus `--panel`:

```bash
convert_genome \
  --input sample.txt \
  --reference hg38.fa \
  --standardize \
  --panel panel.vcf \
  --output sample.bcf \
  --format bcf
```

#### Handling Ancient/Old Data

`--assembly` defines the **target** output build label and drives liftover decisions.

If your input is detected as `hg19` / `GRCh37`, but you request `GRCh38` output, the tool will automatically trigger liftover:

```bash
convert_genome \
  --input ancient.txt \
  --reference GRCh38.fa \
  --assembly GRCh38 \
  --output ancient.vcf
```

Bi-directional liftover is supported (up-lift to `GRCh38` or down-lift to `GRCh37`/`hg19`), and the chain registry includes older conversions (e.g., `NCBI36`).

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

### Flag Deep Dives

#### `--standardize`

Enforces reference-matching alleles against the provided reference FASTA.

- For SNPs, this can include allele polarization (swapping `REF/ALT` and remapping `GT`) when the reference base is found among the alleles.
- Applies ploidy rules for sex chromosomes (with PAR region awareness) when sex is known or inferred.

#### `--panel`

Harmonizes alleles against a reference VCF panel (e.g., 1000 Genomes sites).

- This is especially important for strand-ambiguous SNPs (A/T and C/G), where naive complement checks cannot uniquely determine orientation.
- The panel provides an additional anchor for allele alignment and padding.

#### `--sex`

If omitted, `convert_genome` infers sex using X/Y heterozygosity and variant density.

You can override inference with:

- `--sex male`
- `--sex female`

This influences X/Y ploidy enforcement (e.g., male non-PAR X/Y is treated as haploid).

#### `--assembly`

Defines the **target** output build label (default: `GRCh38`).

If the input build is detected as `GRCh37`/`hg19` but you specify `--assembly GRCh38`, the tool will automatically trigger the liftover engine.

## Using as a Rust Library

### Add the Dependency

Add this to your `Cargo.toml`:

```toml
[dependencies]
convert_genome = "0.1.2"
```

### Core Types and Entry Points

- `ConversionConfig`: configuration struct for controlling input/output formats, paths, inference toggles, and harmonization.
- `convert_dtc_file`: primary entry point for running a conversion and getting a `ConversionSummary` back.

The library uses `anyhow::Result` for top-level errors and also has internal error types for record-level failures.

### Minimal Example

```rust
use convert_genome::{convert_dtc_file, ConversionConfig, OutputFormat};
use convert_genome::cli::Sex;
use convert_genome::input::InputFormat;
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    let summary = convert_dtc_file(ConversionConfig {
        input: PathBuf::from("data/genotypes.txt"),
        input_format: InputFormat::Auto,
        input_origin: String::from("local"),
        reference_fasta: Some(PathBuf::from("GRCh38.fa")),
        reference_origin: Some(String::from("local")),
        reference_fai: None,
        reference_fai_origin: None,
        output: PathBuf::from("out.vcf"),
        output_dir: None,
        output_format: OutputFormat::Vcf,
        sample_id: String::from("SAMPLE"),
        assembly: String::from("GRCh38"),
        include_reference_sites: true,
        sex: Some(Sex::Female),
        par_boundaries: None,
        standardize: true,
        panel: None,
    })?;

    eprintln!("emitted_records={}", summary.emitted_records);
    Ok(())
}
```

## Genome Build Detection & Liftover

### How does it know?

The tool samples variants from your input and checks concordance against expected reference alleles to distinguish common builds (e.g., `GRCh37`/`hg19` vs `GRCh38`). This prevents accidentally emitting a file labeled as one build while containing coordinates from another.

### Liftover Details

- **Trigger:** `--assembly` defines the **target** build.
  If the input is detected as `hg19`/`GRCh37` but you request `GRCh38`, liftover is automatically enabled.
- **Automatic downloads:** required UCSC chain files (e.g., `hg19ToHg38.over.chain.gz`) are fetched on demand and cached locally. The first run may require an internet connection.
- **Fail-safe behavior:** variants that cannot be mapped (deleted regions, gaps, ambiguous mappings) are safely filtered and counted in the output report.

## How It Works (Data Processing Logic)

### Strand Flipping Logic

At a trusted SNP site, one allele should match the reference base for the assumed build.

- If alleles do not match the reference, the tool checks the complement.
- If the complement matches, the allele representation is flipped.
- Strand-ambiguous SNPs (A/T and C/G) require additional context and are best handled with `--panel`.

### Ploidy Enforcement

The converter enforces expected ploidy on sex chromosomes:

- Male non-PAR X/Y: haploid
- Female X: diploid
- Female Y: dropped

PAR boundaries are assembly-specific.

### Build Detection

Build detection samples input sites and estimates concordance against expected reference alleles to infer `GRCh37`/`hg19` vs `GRCh38`.

## Output Report

### Why the report matters

Every run produces a JSON report alongside the output file (e.g., `sample_report.json`). Treat this as an audit trail:

- what the tool inferred (sex, build)
- whether liftover was applied
- how many sites were standardized or harmonized
- how many variants were filtered or failed verification

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

### Key field definitions

- `total_records`: number of input rows processed from the source (after basic parsing).
- `emitted_records`: number of records written to the output.
- `variant_records`: emitted records with at least one ALT allele.
- `reference_records`: emitted reference-matching records (ALT empty).

Liftover-specific counters are also included:

- `liftover_unmapped`: no chain interval found
- `liftover_ambiguous`: multiple chains/intervals eligible (rejected)
- `liftover_incompatible`: alleles incompatible with target reference checks
- `liftover_straddled`: endpoints do not lift consistently (e.g., indel spans blocks)
- `liftover_contig_missing`: lifted contig does not exist in target reference

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
