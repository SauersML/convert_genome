"""convert_genome — Python bindings for the SauersML/convert_genome engine.

Convert direct-to-consumer (23andMe, AncestryDNA, ...) and standard
VCF/BCF inputs into compliant VCF, BCF, or PLINK 1.9 binary, with build
detection, sex inference, liftover, and panel harmonisation.

This package ships a **native, in-process** binding to the Rust core
(built with maturin/PyO3): the compiled extension
``convert_genome._convert_genome`` is loaded directly into the Python
process, so the core runs with no subprocess to a ``convert_genome``
binary and no ``_report.json`` sidecar round-trip.

In-process core (native, the default)
-------------------------------------

When the native extension is present (a normal wheel install or
``maturin develop``), these call the Rust engine directly, in-process:

  * ``convert_genome.convert(input, output, ...)``  – run a conversion,
    returning the summary statistics as a dict
  * ``convert_genome.infer_sex(input, ...)``        – sex call + metrics
  * ``convert_genome.detect_build(input, ...)``     – build call + rates

The legacy typed subprocess wrapper (``Converter`` / the dataclass-
returning ``convert``, driving the ``convert_genome`` binary) remains
importable from ``convert_genome._api`` as a fallback for environments
without the compiled extension; its result/exception types are
re-exported below. ``HAVE_NATIVE`` reports which path is active, and the
subprocess entry point is always available as ``subprocess_convert``.

Quick start (native)
--------------------

>>> import convert_genome
>>> convert_genome.HAVE_NATIVE
True
>>> stats = convert_genome.convert("23andme.txt", "out.vcf", assembly="GRCh38")
>>> stats["emitted_records"]
612345
>>> convert_genome.detect_build("23andme.txt")["build"]
'GRCh37'
"""

# --- legacy typed subprocess API + shared result/exception types -------------
# Imported unconditionally: these define the public dataclasses/exceptions and
# provide a fallback path when the native extension is unavailable.
from ._api import (
    convert as _subprocess_convert,
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

# --- native, in-process core (preferred) -------------------------------------
try:
    from . import _convert_genome as _native
except ImportError:  # pragma: no cover - only when the extension isn't built
    _native = None

if _native is not None:
    # The native extension is the default in-process implementation.
    convert = _native.convert
    infer_sex = _native.infer_sex
    detect_build = _native.detect_build
    HAVE_NATIVE = True
else:
    # Fall back to the typed subprocess wrapper. ``infer_sex`` / ``detect_build``
    # are native-only conveniences and are not provided on this path.
    convert = _subprocess_convert
    infer_sex = None
    detect_build = None
    HAVE_NATIVE = False

# Always expose the subprocess wrapper explicitly for callers that want the
# typed dataclass result regardless of which path ``convert`` takes.
subprocess_convert = _subprocess_convert

__all__ = [
    # in-process core (native when available)
    "convert",
    "infer_sex",
    "detect_build",
    "HAVE_NATIVE",
    # typed subprocess wrapper + helpers
    "subprocess_convert",
    "Converter",
    "locate_binary",
    # result types
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
    # exceptions
    "ConvertGenomeError",
    "ConvertGenomeBinaryNotFound",
    "ConvertGenomeFailed",
    "InvalidConfig",
    "ReportNotFound",
]

__version__ = "0.2.0"
