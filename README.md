# convert_genome

[![Crates.io](https://img.shields.io/crates/v/convert_genome.svg)](https://crates.io/crates/convert_genome)
[![docs.rs](https://img.shields.io/docsrs/convert_genome)](https://docs.rs/convert_genome)
[![CI](https://github.com/SauersML/convert_genome/actions/workflows/ci.yml/badge.svg)](https://github.com/SauersML/convert_genome/actions/workflows/ci.yml)

`convert_genome` converts genomes in one format to another. This includes, for example, direct-to-consumer (DTC) genotype exports (23andMe, AncestryDNA, MyHeritage, etc.) into standard [VCF](https://samtools.github.io/hts-specs/VCFv4.5.pdf), [BCF](https://samtools.github.io/hts-specs/BCFv2_qref.pdf), or [PLINK](https://www.cog-genomics.org/plink/1.9/formats) binary formats. The converter supports remote references via HTTP(S) and handles compressed `.gz` and `.zip` archives. You can also convert from one genome build from another. It automatically detects genome builds, infers biological sex, handles coordinate liftovers, and harmonizes alleles against reference panels.

## Installation

### Automatic Install (Recommended)
Installs the latest binary for your platform (macOS/Linux/Windows):
```bash
# macOS / Linux / Windows (Git Bash)
curl -fsSL https://raw.githubusercontent.com/SauersML/convert_genome/main/install.sh | bash
```

## Core Capabilities

1.  **Format Conversion**: Transforms generic CSV/TSV genotype tables into compliant VCF v4.5, BCF v2.2, or PLINK 1.9 binary formats.
2.  **Smart Input Handling**: Transparently handles plain text, GZIP-compressed files, and ZIP archives.
3.  **Automatic Inference**:
    *   **Genome Build**: Detects if input is GRCh37/hg19 or GRCh38/hg19 based on coordinate/allele concordance.
    *   **Biological Sex**: Infers sex based on X-chromosome heterozygosity and Y-chromosome variant density (needed to determine ploidy for VCF).
    *   **Strand Orientation**: Detects if the input file is reported on the forward or reverse strand relative to the reference.
4.  **Automatic Liftover**: If the detected input build differs from the requested output assembly (e.g., input is hg19, output is GRCh38), the tool automatically downloads UCSC chain files and lifts coordinates over.
5.  **Standardization & Harmonization**:
    *   **Allele Polarization**: Swaps REF/ALT alleles to match the provided reference genome.
    *   **Ploidy Enforcement**: Enforces correct haploid/diploid states for sex chromosomes based on sex and Pseudoautosomal Regions (PAR).
    *   **Panel Alignment**: Aligns alleles against an external reference VCF (e.g., 1000 Genomes) to resolve strand ambiguity (A/T and C/G SNPs).

---

## Usage Logic & Behavior

### Input Parsing
The tool automatically detects the input format. It supports:
*   **DTC Formats**: 23andMe, AncestryDNA (5-column), MyHeritage (CSV), deCODEme (6-column with strand), and generic whitespace-delimited formats.
*   **Standard Formats**: VCF and BCF inputs are also supported for normalization/liftover tasks.

**Compression**: Inputs can be `.gz`, `.zip`, bgzipped, or some combination thereof.

### Genome Build Detection & Liftover
The conversion pipeline uses a specific logic flow to ensure coordinates are correct:

1.  **Detection**: The tool samples variants from the input and checks concordance against known reference alleles for both `GRCh37` and `GRCh38`.
2.  **Trigger**: The user specifies a *target* `--assembly` (default: `GRCh38`).
3.  **Action**:
    *   If **Input Build == Target Assembly**, conversion proceeds directly.
    *   If **Input Build != Target Assembly**, the "Liftover Engine" is engaged.
4.  **Liftover Execution**:
    *   Required chain files (e.g., `hg19ToHg38.over.chain.gz`) are automatically downloaded from UCSC to a local cache.
    *   Coordinates are remapped.
    *   **Fail-Closed Behavior**: Variants that map to multiple locations (ambiguous), do not map at all, or "straddle" chain boundaries (split indels) are discarded to preserve data integrity.

### Sex Inference & Ploidy
Unless explicitly overridden via flags, the tool infers sex to ensure correct VCF representation:
*   **Female**: X is Diploid, Y is absent (or filtered).
*   **Male**: X and Y are Haploid, *except* in Pseudoautosomal Regions (PAR1/PAR2), where they remain Diploid.
*   **Mitochondrial (MT)**: Always treated as Haploid.

### Allele Standardization (Polarization)
When `--standardize` is enabled, the tool ensures the `REF` allele in the output VCF matches the FASTA reference provided.
*   If the input reports "A" but the Reference genome says "G":
    *   The tool checks if the "A" is a valid ALT.
    *   It swaps the alleles (REF becomes G, ALT becomes A) and updates the Genotype (GT) indices (e.g., `0/0` -> `1/1`).
*   Synthetic IDs are generated for variants lacking identifiers (formatted as `chrom:pos:ref:alt`).

### Panel Harmonization
For preparing data for imputation (e.g., Beagle), strand ambiguity is a major issue (e.g., an A/T SNP is indistinguishable from a T/A SNP if the strand is unknown).
*   By providing a `--panel` (VCF/BCF), the tool checks if the input alleles match the panel's alleles.
*   It attempts to align the input to the panel, flipping strands if necessary.
*   It outputs a "padded" panel containing any novel alleles found in the user data, ensuring the reference panel and target VCF are perfectly compatible for imputation tools.

---

## Output Formats

### VCF (Variant Call Format)
*   **Version**: 4.5.
*   **Encoding**: Standard GT (Genotype) field.
*   **Metadata**: Headers include assembly, conversion software version, and date.
*   **Symbolic Alleles**: Large deletions/insertions are normalized to `<DEL>` or `<INS>` symbolic alleles with `SVTYPE` info fields.

### BCF (Binary Call Format)
*   **Version**: 2.2.
*   **Behavior**: Functionally identical to VCF but highly compressed and indexed. Recommended for large-scale pipelines.

### PLINK 1.9 (Binary)
Produces a file trio using the output filename as a prefix:
1.  **`.bed`**: Primary binary genotype matrix.
2.  **`.bim`**: Variant information (Chromosome, SNP ID, cM, Position, Allele 1, Allele 2).
3.  **`.fam`**: Sample information (FID, IID, Paternal ID, Maternal ID, Sex, Phenotype).
    *   *Note*: Sex is encoded as 1 (Male) or 2 (Female). Phenotype is set to -9 (Missing).

### Run Report (`_report.json`)
Every execution produces a sidecar JSON file containing a comprehensive audit trail:
*   **Inference Results**: What sex and build were detected.
*   **Statistics**: Total records, valid records, filtered variants.
*   **Liftover Details**: Specific counts for unmapped variants, ambiguous mappings, or reference mismatches.
*   **Panel Stats**: How many sites were harmonized vs. novel.

---

## CLI Options Overview

The tool is controlled via a unified command-line interface.

**Required:**
*   `--input`: Path to the genotype file.
*   `--reference`: Path to the reference FASTA (or a URL).
*   `--output` (or `--output-dir`): Destination for the converted data.

**Key Flags:**
*   `--format`: `vcf`, `bcf`, or `plink`.
*   `--assembly`: The **target** assembly (e.g., `GRCh38`). Drive liftover logic.
*   `--sex`: Override automatic sex inference (`male` or `female`).
*   `--standardize`: Enable reference-based allele polarization and normalization.
*   `--panel`: Path to a VCF/BCF reference panel for harmonization.
*   `--variants-only`: Output only sites where the sample differs from the reference.

---

## Library Usage (API Concepts)

For Rust developers, the core logic is exposed via the `conversion` module. The primary entry point is the `convert_dtc_file` function, driven by a configuration struct.

### `ConversionConfig`
The configuration object controls the entire pipeline state:
*   **Paths**: Input, output, reference FASTA, and optional panel paths.
*   **Format Enums**: `InputFormat` (DTC/VCF/BCF) and `OutputFormat`.
*   **Biological Context**: `Sex` enum, `ParBoundaries` (defined ranges for X/Y recombination), and `assembly` string.
*   **Flags**: `standardize` (bool), `include_reference_sites` (bool).

### `convert_dtc_file(config: ConversionConfig) -> Result<ConversionSummary>`
This function executes the pipeline:
1.  **Pre-scan**: Reads a subset of the file to run the Inference Engine (Sex/Build).
2.  **Resource Loading**: Fetches/Loads reference genomes, chain files, and panels.
3.  **Source Iterator**: Wraps the input in a smart iterator that handles parsing, sorting, and initial validation.
4.  **Transformation Stream**:
    *   **Liftover Adapter**: If chains are loaded, maps coordinates on-the-fly.
    *   **Standardizer**: Polarizes alleles against the loaded ReferenceGenome.
    *   **Harmonizer**: Aligns against the PaddedPanel.
5.  **Writing**: Streams processed records to the specific output writer (VCF/BCF/PLINK).

### `ConversionSummary`
The return object provides precise metrics on the run, useful for quality control or integration testing:
*   `total_records` vs `emitted_records`.
*   `liftover_unmapped`: Count of variants lost due to missing chain mappings.
*   `invalid_genotypes`: Count of malformed input lines.
*   `reference_failures`: Count of sites where reference lookup failed.
