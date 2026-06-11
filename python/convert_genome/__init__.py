"""convert_genome — Python bindings for the SauersML/convert_genome CLI.

Convert direct-to-consumer (23andMe, AncestryDNA, ...) and standard
VCF/BCF inputs into compliant VCF, BCF, or PLINK 1.9 binary, with build
detection, sex inference, liftover, and panel harmonisation.

This package shells out to the `convert_genome` Rust binary and parses
its sidecar `_report.json` into typed dataclasses. The Python API is
typed kwargs end-to-end; you never need to remember CLI flag names.

Quick start
-----------

>>> from convert_genome import convert, OutputFormat
>>> result = convert(
...     input="23andme.txt",
...     output="out.vcf",
...     format=OutputFormat.VCF,
...     assembly="GRCh38",
...     standardize=True,
... )
>>> result.sample.sex_inferred
True
>>> result.statistics.emitted_records
612_345
>>> result.build_detection.detected_build
'GRCh37'
"""

from ._api import (
    convert,
    Converter,
    ConversionResult,
    InputInfo,
    OutputInfo,
    ReferenceInfo,
    PanelInfo,
    SampleInfo,
    BuildDetection,
    Statistics,
    InputFormat,
    OutputFormat,
    Sex,
    Assembly,
    ConvertGenomeError,
    ConvertGenomeBinaryNotFound,
    ConvertGenomeFailed,
    InvalidConfig,
    ReportNotFound,
    locate_binary,
)

__all__ = [
    "convert",
    "Converter",
    "ConversionResult",
    "InputInfo",
    "OutputInfo",
    "ReferenceInfo",
    "PanelInfo",
    "SampleInfo",
    "BuildDetection",
    "Statistics",
    "InputFormat",
    "OutputFormat",
    "Sex",
    "Assembly",
    "ConvertGenomeError",
    "ConvertGenomeBinaryNotFound",
    "ConvertGenomeFailed",
    "InvalidConfig",
    "ReportNotFound",
    "locate_binary",
]

__version__ = "0.2.0"
