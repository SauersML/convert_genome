# convert_genome (Python)

Python wrapper for the [`SauersML/convert_genome`](https://github.com/SauersML/convert_genome)
CLI: convert direct-to-consumer genotype dumps (23andMe, AncestryDNA, ...)
and standard VCF/BCF inputs into compliant VCF, BCF, or PLINK 1.9 binary.

## Install

```bash
pip install convert_genome
# you also need the Rust binary:
cargo install convert_genome
```

The wrapper finds the binary on `PATH`, via `$CONVERT_GENOME_BIN`, or
through the `binary=` keyword argument.

## Quick start

```python
from convert_genome import convert, OutputFormat

result = convert(
    input="23andme.txt",
    output="out.vcf",
    format=OutputFormat.VCF,
    assembly="hg38",
    standardize=True,
)

print(result.statistics.emitted_records)
print(result.sample.sex_inferred)
print(result.build_detection.detected_build)
print(result.report_path)
```

The Rust tool writes `<stem>_report.json` next to each output. The
wrapper loads that JSON and returns a typed `ConversionResult` whose
`statistics`, `sample`, `build_detection`, and friends are immutable
dataclasses.

## Shortcuts: avoid downloads and inference

If you already have the reference, panel, target build, or sex on hand,
pass them in — convert_genome will skip every corresponding
auto-discovery step.

```python
convert(
    input="raw.txt",
    output="out.vcf",
    reference="/cache/hg38.fa",         # skip FASTA download
    reference_fai="/cache/hg38.fa.fai", # skip .fai indexing
    input_build="hg19",                  # skip build detection
    assembly="GRCh38",                   # target build (still does liftover)
    panel="/cache/1kg_panel.vcf",        # supply harmonisation panel
    sex="female",                        # skip sex inference
    standardize=True,
)
```

Every flag the upstream CLI accepts (`--input-format`, `--variants-only`,
`--log-level`, ...) is also reachable as a kwarg or via the `Converter`
builder's `with_*` methods.

## Builder

```python
from convert_genome import Converter, Sex, OutputFormat

plan = (
    Converter(input="raw.txt", output_dir="out/", format=OutputFormat.PLINK)
        .with_assembly("GRCh38")
        .with_reference("/cache/hg38.fa", "/cache/hg38.fa.fai")
        .with_panel("/data/1kg_panel.vcf.gz")
        .with_standardize()
        .with_sex(Sex.MALE)
)

result = plan.run()
print(plan.argv())  # everything the wrapper would pass to the CLI
```

`Converter` is a frozen dataclass; every `with_*` method returns a new
instance, so branching off of a partially-built plan is safe.

## Errors

* `ConvertGenomeBinaryNotFound` — CLI not installed / not on PATH.
* `InvalidConfig` — argument combination rejected before launching.
* `ConvertGenomeFailed` — CLI exited non-zero (stdout/stderr available).
* `ReportNotFound` — CLI ran clean but didn't produce a JSON report.

All subclass `ConvertGenomeError`.
