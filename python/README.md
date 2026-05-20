# convert_genome (Python)

Python wrapper for the
[`SauersML/convert_genome`](https://github.com/SauersML/convert_genome) CLI.
Convert direct-to-consumer dumps (23andMe, AncestryDNA, MyHeritage,
deCODEme) and standard VCF/BCF into compliant VCF, BCF, or PLINK 1.9
binary — with build detection, sex inference, liftover, and panel
harmonisation, all controllable from kwargs.

```python
from convert_genome import convert, OutputFormat

result = convert(
    input="23andme.txt",
    output="out.vcf",
    format=OutputFormat.VCF,
    assembly="hg38",
    standardize=True,
)

result.statistics.emitted_records       # int
result.sample.sex_inferred              # bool
result.build_detection.detected_build   # 'GRCh37' / 'GRCh38' / ...
result.report_path                      # path to <stem>_report.json
result.output_paths                     # files that actually exist on disk
result.yield_rate                       # emitted / total
```

The wrapper runs the Rust binary, parses the sidecar
`<stem>_report.json` into typed frozen dataclasses, and returns a
single `ConversionResult`.

## Install

```bash
pip install convert_genome
# the Rust binary:
cargo install convert_genome
```

Binary located via `binary=`, `$CONVERT_GENOME_BIN`, or `PATH` (in that
order). Missing binary → `ConvertGenomeBinaryNotFound` with the
suggested install command.

## Shortcuts: skip every auto-discovery step

The CLI will download/auto-detect things it doesn't need to. Pass them
in directly:

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

`sex` is lenient: passing `"unknown"` or `"indeterminate"` (e.g. when
chaining out of `infer_sex`) silently omits the `--sex` flag and lets
the CLI run its own inference.

## Builder

`Converter` is a frozen dataclass; every `with_*` returns a new
instance, so branching is safe.

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

print(plan.argv())   # exact argv that would be passed to the CLI
result = plan.run()
```

## Enums

```python
InputFormat.AUTO / .DTC / .VCF / .BCF
OutputFormat.VCF / .BCF / .PLINK
Sex.MALE / .FEMALE
Assembly.GRCH37 / .GRCH38     # plus a `.parse()` classmethod that
                              # accepts 'hg19' / 'hg38' / 'build38' / ...
```

## Output

The Rust tool writes `<stem>_report.json` alongside the main output.
The wrapper loads it into `ConversionResult`, with sub-dataclasses for
each section:

```python
result.input         # InputInfo (path, format, origin)
result.output        # OutputInfo (path, format)
result.reference     # ReferenceInfo (path, origin, assembly)
result.panel         # PanelInfo | None
result.sample        # SampleInfo (id, sex, sex_inferred)
result.build_detection  # BuildDetection | None (detected_build, match rates)
result.statistics    # Statistics (total / emitted / variant / ... records)
result.report_path   # path to the JSON sidecar
result.output_paths  # tuple[Path] — files that actually exist on disk
```

For PLINK output, `output_paths` includes the `.bed/.bim/.fam` trio. For
`output_dir` with a panel, it includes `panel.vcf`. Non-existent paths
are filtered out automatically.

## Errors

* `ConvertGenomeBinaryNotFound` — CLI not installed / not on PATH.
* `InvalidConfig` — argument combination rejected before launching
  (e.g. missing input file, conflicting output/output_dir).
* `ConvertGenomeFailed` — CLI exited non-zero. The exception carries
  `stdout`, `stderr`, `returncode`.
* `ReportNotFound` — CLI ran clean but didn't write a JSON sidecar.

All subclass `ConvertGenomeError`.
