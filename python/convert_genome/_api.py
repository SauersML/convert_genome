"""Pythonic wrapper around the convert_genome CLI.

Design
------
The CLI already writes a structured ``<stem>_report.json`` next to each
output, so we don't parse stdout — we run the binary, wait for it to
finish, then load that JSON into typed frozen dataclasses.

The ``Converter`` class is an immutable builder. The top-level
``convert(...)`` is the one-shot convenience.

We deliberately don't shadow the CLI's auto-detection logic: pass
``input_format=InputFormat.AUTO`` (the default) and let the Rust tool
sniff. Where we *do* validate eagerly is on parameter combinations that
the CLI rejects late and noisily — e.g. ``--output`` xor ``--output-dir``,
``--reference-fai`` requiring ``--reference``.
"""

from __future__ import annotations

import enum
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Iterable, List, Mapping, Optional, Tuple, Union

PathLike = Union[str, os.PathLike]


# ---------------------------------------------------------------------------
# Enums (mirror src/cli.rs)
# ---------------------------------------------------------------------------


class InputFormat(str, enum.Enum):
    AUTO = "auto"
    DTC = "dtc"
    VCF = "vcf"
    BCF = "bcf"


class OutputFormat(str, enum.Enum):
    VCF = "vcf"
    BCF = "bcf"
    PLINK = "plink"


class Sex(str, enum.Enum):
    MALE = "male"
    FEMALE = "female"


class Assembly(str, enum.Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"

    @classmethod
    def parse(cls, value: Union[str, "Assembly", None]) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, cls):
            return value.value
        norm = str(value).strip()
        low = norm.lower().replace("-", "").replace("_", "")
        aliases = {
            "grch37": "GRCh37",
            "hg19": "GRCh37",
            "build37": "GRCh37",
            "grch38": "GRCh38",
            "hg38": "GRCh38",
            "build38": "GRCh38",
        }
        return aliases.get(low, norm)


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class ConvertGenomeError(Exception):
    """Base class for all convert_genome wrapper errors."""


class ConvertGenomeBinaryNotFound(ConvertGenomeError, FileNotFoundError):
    """The `convert_genome` binary could not be located."""


class InvalidConfig(ConvertGenomeError, ValueError):
    """A combination of arguments is mutually exclusive or incomplete."""


class ConvertGenomeFailed(ConvertGenomeError, RuntimeError):
    """The binary ran but returned a non-zero exit code."""

    def __init__(self, message: str, *, stdout: str = "", stderr: str = "", returncode: int = 0):
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode


class ReportNotFound(ConvertGenomeError, FileNotFoundError):
    """The binary exited 0 but produced no sidecar `_report.json`."""


# ---------------------------------------------------------------------------
# Result dataclasses — mirror src/report.rs
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class InputInfo:
    path: str
    format: str
    origin: str


@dataclass(frozen=True)
class OutputInfo:
    path: str
    format: str


@dataclass(frozen=True)
class ReferenceInfo:
    path: str
    origin: str
    assembly: str


@dataclass(frozen=True)
class PanelInfo:
    path: str
    total_sites: int
    modified_sites: int
    novel_sites: int


@dataclass(frozen=True)
class SampleInfo:
    id: str
    sex: str
    sex_inferred: bool


@dataclass(frozen=True)
class BuildDetection:
    detected_build: str
    hg19_match_rate: float
    hg38_match_rate: float


@dataclass(frozen=True)
class Statistics:
    total_records: int
    emitted_records: int
    variant_records: int
    reference_records: int
    missing_genotype_records: int
    skipped_reference_sites: int
    unknown_chromosomes: int
    reference_failures: int
    invalid_genotypes: int
    symbolic_allele_records: int
    parse_errors: int


@dataclass(frozen=True)
class ConversionResult:
    """The full run report parsed from ``<stem>_report.json``."""

    version: str
    timestamp: str
    input: InputInfo
    output: OutputInfo
    reference: ReferenceInfo
    standardize: bool
    sample: SampleInfo
    statistics: Statistics
    panel: Optional[PanelInfo] = None
    build_detection: Optional[BuildDetection] = None
    report_path: Optional[Path] = None
    output_paths: Tuple[Path, ...] = field(default_factory=tuple)
    stdout: str = ""
    stderr: str = ""

    @property
    def main_output(self) -> Path:
        return Path(self.output.path)

    @property
    def emitted_records(self) -> int:
        return self.statistics.emitted_records

    @property
    def total_records(self) -> int:
        return self.statistics.total_records

    @property
    def yield_rate(self) -> float:
        if self.statistics.total_records == 0:
            return 0.0
        return self.statistics.emitted_records / self.statistics.total_records


# ---------------------------------------------------------------------------
# Binary location
# ---------------------------------------------------------------------------


def locate_binary(override: Optional[PathLike] = None) -> Path:
    """Locate `convert_genome` or raise `ConvertGenomeBinaryNotFound`.

    Resolution: explicit ``override`` → ``convert_genome`` on PATH.
    No environment-variable indirection.
    """
    if override is not None:
        p = Path(override)
        if not p.exists():
            raise ConvertGenomeBinaryNotFound(f"convert_genome binary not at {p}")
        return p
    which = shutil.which("convert_genome")
    if which:
        return Path(which)
    raise ConvertGenomeBinaryNotFound(
        "convert_genome not found. Install with: cargo install convert_genome, "
        "or pass binary=... explicitly."
    )


# ---------------------------------------------------------------------------
# Converter
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Converter:
    """Immutable conversion plan. Call ``.run()`` to execute."""

    input: Path
    output: Optional[Path] = None
    output_dir: Optional[Path] = None
    format: OutputFormat = OutputFormat.VCF
    input_format: InputFormat = InputFormat.AUTO
    assembly: str = "GRCh38"
    input_build: Optional[str] = None
    reference: Optional[Path] = None
    reference_fai: Optional[Path] = None
    panel: Optional[Path] = None
    sample: Optional[str] = None
    sex: Optional[Sex] = None
    standardize: bool = False
    variants_only: bool = False
    log_level: str = "info"
    binary: Optional[Path] = None
    timeout: Optional[float] = None
    extra_args: Tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.output is None and self.output_dir is None:
            raise InvalidConfig("Either output= or output_dir= must be provided.")
        if self.output is not None and self.output_dir is not None:
            raise InvalidConfig("output= and output_dir= are mutually exclusive.")
        if self.reference_fai is not None and self.reference is None:
            raise InvalidConfig("reference_fai= requires reference=")
        if not Path(self.input).exists():
            raise InvalidConfig(f"Input file does not exist: {self.input}")

    # Builder helpers — return a new Converter with the field replaced.

    def with_output(self, output: PathLike) -> "Converter":
        return replace(self, output=Path(output), output_dir=None)

    def with_output_dir(self, output_dir: PathLike) -> "Converter":
        return replace(self, output_dir=Path(output_dir), output=None)

    def with_reference(self, ref: PathLike, fai: Optional[PathLike] = None) -> "Converter":
        return replace(self, reference=Path(ref), reference_fai=Path(fai) if fai else None)

    def with_panel(self, panel: PathLike) -> "Converter":
        return replace(self, panel=Path(panel))

    def with_sex(self, sex: Sex) -> "Converter":
        return replace(self, sex=sex)

    def with_sample(self, sample: str) -> "Converter":
        return replace(self, sample=sample)

    def with_standardize(self, on: bool = True) -> "Converter":
        return replace(self, standardize=on)

    def with_variants_only(self, on: bool = True) -> "Converter":
        return replace(self, variants_only=on)

    def with_assembly(self, assembly: str) -> "Converter":
        return replace(self, assembly=Assembly.parse(assembly) or "GRCh38")

    def with_input_build(self, build: Optional[str]) -> "Converter":
        return replace(self, input_build=Assembly.parse(build))

    def with_binary(self, path: PathLike) -> "Converter":
        return replace(self, binary=Path(path))

    def with_timeout(self, seconds: Optional[float]) -> "Converter":
        return replace(self, timeout=seconds)

    def with_log_level(self, level: str) -> "Converter":
        return replace(self, log_level=level)

    def with_extra_args(self, args: Iterable[str]) -> "Converter":
        return replace(self, extra_args=tuple(args))

    # --- Execution ---------------------------------------------------------

    def argv(self) -> List[str]:
        """Compute the argv that would be invoked. Useful for tests / dry-runs."""
        binary = locate_binary(self.binary)
        argv: List[str] = [str(binary)]

        argv += ["--input-format", self.input_format.value]
        argv += ["--format", self.format.value]
        argv += ["--output-build", Assembly.parse(self.assembly) or self.assembly]

        if self.reference is not None:
            argv += ["--reference", str(self.reference)]
        if self.reference_fai is not None:
            argv += ["--reference-fai", str(self.reference_fai)]
        if self.panel is not None:
            argv += ["--panel", str(self.panel)]
        if self.sample is not None:
            argv += ["--sample", self.sample]
        if self.input_build is not None:
            argv += ["--input-build", self.input_build]
        if self.sex is not None:
            argv += ["--sex", self.sex.value]
        if self.standardize:
            argv += ["--standardize"]
        if self.variants_only:
            argv += ["--variants-only"]
        if self.log_level and self.log_level != "info":
            argv += ["--log-level", self.log_level]

        if self.output_dir is not None:
            argv += ["--output-dir", str(self.output_dir)]

        argv += list(self.extra_args)

        # positional: INPUT first, then OUTPUT if not --output-dir
        argv.append(str(self.input))
        if self.output is not None:
            argv.append(str(self.output))
        return argv

    def run(self, *, capture: bool = True) -> ConversionResult:
        argv = self.argv()
        try:
            completed = subprocess.run(
                argv,
                capture_output=capture,
                text=True,
                timeout=self.timeout,
                check=False,
            )
        except FileNotFoundError as e:
            raise ConvertGenomeBinaryNotFound(str(e)) from e

        if completed.returncode != 0:
            raise ConvertGenomeFailed(
                f"convert_genome exited with status {completed.returncode}",
                stdout=completed.stdout or "",
                stderr=completed.stderr or "",
                returncode=completed.returncode,
            )

        report_path = self._resolve_report_path(completed.stdout or "", completed.stderr or "")
        if not report_path.exists():
            raise ReportNotFound(f"Expected report at {report_path} but it is missing.")

        with open(report_path) as f:
            data = json.load(f)
        return _result_from_report(
            data,
            report_path=report_path,
            output_paths=self._resolve_outputs(),
            stdout=completed.stdout or "",
            stderr=completed.stderr or "",
        )

    # --- helpers -----------------------------------------------------------

    def _resolve_outputs(self) -> Tuple[Path, ...]:
        outs: List[Path] = []
        if self.output is not None:
            outs.append(self.output)
            if self.format is OutputFormat.PLINK:
                base = self.output.with_suffix("")
                outs += [base.with_suffix(s) for s in (".bed", ".bim", ".fam")]
        if self.output_dir is not None:
            d = self.output_dir
            fname = {
                OutputFormat.VCF: "genotypes.vcf",
                OutputFormat.BCF: "genotypes.bcf",
                OutputFormat.PLINK: "genotypes",
            }[self.format]
            primary = d / fname
            outs.append(primary)
            if self.format is OutputFormat.PLINK:
                outs += [d / f"genotypes{s}" for s in (".bed", ".bim", ".fam")]
            if self.panel is not None:
                outs.append(d / "panel.vcf")
        return tuple(p for p in outs if p.exists())

    def _resolve_report_path(self, stdout: str, stderr: str) -> Path:
        m = re.search(r"Wrote run report to ([^\r\n]+)", stdout + "\n" + stderr)
        if m:
            return Path(m.group(1).strip())

        if self.output is not None:
            stem = self.output.stem
            return self.output.with_name(f"{stem}_report.json")
        # output_dir case (validated in __post_init__)
        d = self.output_dir
        assert d is not None
        return d / "genotypes_report.json"


# ---------------------------------------------------------------------------
# Top-level convenience
# ---------------------------------------------------------------------------


def convert(
    *,
    input: PathLike,
    output: Optional[PathLike] = None,
    output_dir: Optional[PathLike] = None,
    format: Union[OutputFormat, str] = OutputFormat.VCF,
    input_format: Union[InputFormat, str] = InputFormat.AUTO,
    assembly: str = "GRCh38",
    input_build: Optional[str] = None,
    reference: Optional[PathLike] = None,
    reference_fai: Optional[PathLike] = None,
    panel: Optional[PathLike] = None,
    sample: Optional[str] = None,
    sex: Optional[Union[Sex, str]] = None,
    standardize: bool = False,
    variants_only: bool = False,
    log_level: str = "info",
    binary: Optional[PathLike] = None,
    timeout: Optional[float] = None,
    extra_args: Optional[Iterable[str]] = None,
    capture: bool = True,
) -> ConversionResult:
    """Run one conversion. Returns the parsed run report."""
    # Lenient sex coercion: callers chaining results out of infer_sex
    # may pass `InferredSex.INDETERMINATE.value` ("indeterminate") or
    # `gnomon`'s "unknown". Neither is a valid convert_genome --sex
    # value, but the most useful behaviour is "no override — let the
    # CLI run its own inference", which is the same as `sex=None`.
    sex_coerced: Optional[Sex]
    if sex is None:
        sex_coerced = None
    elif isinstance(sex, str) and sex.strip().lower() in {"unknown", "indeterminate", ""}:
        sex_coerced = None
    elif isinstance(sex, Sex):
        sex_coerced = sex
    else:
        sex_coerced = _coerce_enum(sex, Sex)

    converter = Converter(
        input=Path(input),
        output=Path(output) if output else None,
        output_dir=Path(output_dir) if output_dir else None,
        format=_coerce_enum(format, OutputFormat),
        input_format=_coerce_enum(input_format, InputFormat),
        assembly=Assembly.parse(assembly) or "GRCh38",
        input_build=Assembly.parse(input_build),
        reference=Path(reference) if reference else None,
        reference_fai=Path(reference_fai) if reference_fai else None,
        panel=Path(panel) if panel else None,
        sample=sample,
        sex=sex_coerced,
        standardize=standardize,
        variants_only=variants_only,
        log_level=log_level,
        binary=Path(binary) if binary else None,
        timeout=timeout,
        extra_args=tuple(extra_args) if extra_args else (),
    )
    return converter.run(capture=capture)


def _coerce_enum(value, enum_cls):
    if isinstance(value, enum_cls):
        return value
    if isinstance(value, str):
        try:
            return enum_cls(value.lower())
        except ValueError:
            pass
        for member in enum_cls:
            if member.name.lower() == value.lower():
                return member
    raise InvalidConfig(f"Cannot coerce {value!r} to {enum_cls.__name__}")


# ---------------------------------------------------------------------------
# Report parsing
# ---------------------------------------------------------------------------


def _result_from_report(
    data: Mapping[str, Any],
    *,
    report_path: Path,
    output_paths: Tuple[Path, ...],
    stdout: str,
    stderr: str,
) -> ConversionResult:
    try:
        stats = Statistics(**data["statistics"])
        input_info = InputInfo(**data["input"])
        output_info = OutputInfo(**data["output"])
        reference_info = ReferenceInfo(**data["reference"])
        sample_info = SampleInfo(**data["sample"])
        panel_info = PanelInfo(**data["panel"]) if data.get("panel") else None
        build = BuildDetection(**data["build_detection"]) if data.get("build_detection") else None
        return ConversionResult(
            version=data["version"],
            timestamp=data["timestamp"],
            input=input_info,
            output=output_info,
            reference=reference_info,
            standardize=data["standardize"],
            sample=sample_info,
            statistics=stats,
            panel=panel_info,
            build_detection=build,
            report_path=report_path,
            output_paths=output_paths,
            stdout=stdout,
            stderr=stderr,
        )
    except KeyError as e:
        raise ConvertGenomeFailed(
            f"Report {report_path} missing expected field: {e}",
            stdout=stdout,
            stderr=stderr,
        ) from e
    except TypeError as e:
        raise ConvertGenomeFailed(
            f"Report {report_path} has an unexpected schema: {e}",
            stdout=stdout,
            stderr=stderr,
        ) from e
