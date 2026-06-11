"""Tests for the convert_genome Python wrapper.

We don't run the real Rust binary in CI. Tests build a fake binary as a
small Python script that:

  1. Echoes argv (so we can assert flag mapping).
  2. Writes a minimal but valid run-report JSON wherever the wrapper
     expects to find it.
  3. Exits 0 (or fails on demand, to test error paths).
"""

from __future__ import annotations

import json
import stat
import textwrap
from pathlib import Path

import pytest

from convert_genome import (
    Assembly,
    ConversionResult,
    Converter,
    ConvertGenomeBinaryNotFound,
    ConvertGenomeError,
    ConvertGenomeFailed,
    InputFormat,
    InvalidConfig,
    OutputFormat,
    ReportNotFound,
    Sex,
    convert,
    locate_binary,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_input(tmp_path: Path, name: str = "in.txt") -> Path:
    p = tmp_path / name
    p.write_text("# rsid\tchromosome\tposition\tgenotype\nrs1\t1\t12345\tAG\n")
    return p


_FAKE_FAILING = textwrap.dedent(
    """\
    #!/usr/bin/env python3
    import sys
    sys.stderr.write('ERROR: synthetic failure\\n')
    sys.exit(2)
    """
)


_FAKE_NOREPORT = textwrap.dedent(
    """\
    #!/usr/bin/env python3
    print('Did some work, but forgot to write a report.')
    """
)


def _good_fake_body(report: dict, *, log_line: str = "") -> str:
    """A fake CLI that writes ``report`` to the wrapper's expected path."""
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env python3
        import json, sys, pathlib

        argv = sys.argv[1:]
        log_argv = pathlib.Path(sys.argv[0]).parent / 'argv.json'
        log_argv.write_text(json.dumps(argv))

        # Figure out where the wrapper expects the report.
        report = {json.dumps(report)!r}
        if '--output-dir' in argv:
            d = pathlib.Path(argv[argv.index('--output-dir') + 1])
            d.mkdir(parents=True, exist_ok=True)
            report_path = d / 'genotypes_report.json'
        else:
            # Positional layout: ... INPUT OUTPUT
            output_path = pathlib.Path(argv[-1])
            stem = output_path.stem
            report_path = output_path.with_name(stem + '_report.json')

        # Also write a sentinel output file so output_paths picks it up.
        if '--output-dir' in argv:
            (d / 'genotypes.vcf').write_text('##fake\\n')
        else:
            output_path.write_text('##fake\\n')

        report_path.write_text(report)
        sys.stdout.write({log_line!r})
        """
    )


def _make_fake(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "fake_convert_genome"
    p.write_text(body)
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return p


def _minimal_report() -> dict:
    return {
        "version": "0.1.2",
        "timestamp": "2026-05-19T00:00:00Z",
        "input": {"path": "/in.txt", "format": "DTC", "origin": "local"},
        "output": {"path": "/out.vcf", "format": "VCF"},
        "reference": {"path": "/ref.fa", "origin": "downloaded", "assembly": "GRCh38"},
        "standardize": True,
        "sample": {"id": "S1", "sex": "Female", "sex_inferred": True},
        "build_detection": {
            "detected_build": "GRCh38",
            "hg19_match_rate": 12.3,
            "hg38_match_rate": 98.7,
        },
        "statistics": {
            "total_records": 1000,
            "emitted_records": 990,
            "variant_records": 800,
            "reference_records": 190,
            "missing_genotype_records": 10,
            "skipped_reference_sites": 0,
            "unknown_chromosomes": 0,
            "reference_failures": 0,
            "invalid_genotypes": 0,
            "symbolic_allele_records": 0,
            "parse_errors": 0,
        },
    }


# ---------------------------------------------------------------------------
# Locator
# ---------------------------------------------------------------------------


def test_locate_binary_not_on_path(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", str(tmp_path))  # empty PATH
    with pytest.raises(ConvertGenomeBinaryNotFound):
        locate_binary()


def test_locate_binary_override(tmp_path):
    fake = _make_fake(tmp_path, "#!/usr/bin/env python3\nprint('x')\n")
    assert locate_binary(fake) == fake


def test_locate_binary_override_missing_raises(tmp_path):
    with pytest.raises(ConvertGenomeBinaryNotFound):
        locate_binary(tmp_path / "no-such-binary")


# ---------------------------------------------------------------------------
# Eager validation
# ---------------------------------------------------------------------------


def test_requires_output_or_output_dir(tmp_path):
    in_ = _make_input(tmp_path)
    with pytest.raises(InvalidConfig):
        Converter(input=in_)


def test_output_xor_output_dir(tmp_path):
    in_ = _make_input(tmp_path)
    with pytest.raises(InvalidConfig):
        Converter(input=in_, output=tmp_path / "out.vcf", output_dir=tmp_path / "d")


def test_reference_fai_requires_reference(tmp_path):
    in_ = _make_input(tmp_path)
    fai = tmp_path / "x.fai"
    fai.write_text("")
    with pytest.raises(InvalidConfig):
        Converter(input=in_, output=tmp_path / "o.vcf", reference_fai=fai)


def test_input_must_exist(tmp_path):
    with pytest.raises(InvalidConfig):
        Converter(input=tmp_path / "nope.txt", output=tmp_path / "o.vcf")


# ---------------------------------------------------------------------------
# Assembly alias normalisation
# ---------------------------------------------------------------------------


def test_assembly_parses_aliases():
    assert Assembly.parse("hg19") == "GRCh37"
    assert Assembly.parse("GRCh38") == "GRCh38"
    assert Assembly.parse("build38") == "GRCh38"
    assert Assembly.parse("Hg-38") == "GRCh38"
    # Unknown strings pass through (CLI will decide).
    assert Assembly.parse("CHM13") == "CHM13"


# ---------------------------------------------------------------------------
# argv generation
# ---------------------------------------------------------------------------


def test_argv_for_simple_vcf(tmp_path):
    in_ = _make_input(tmp_path)
    fake = _make_fake(tmp_path, "")
    argv = (
        Converter(input=in_, output=tmp_path / "o.vcf", binary=fake)
        .with_assembly("hg38")
        .argv()
    )
    assert argv[0] == str(fake)
    assert "--format" in argv and argv[argv.index("--format") + 1] == "vcf"
    assert argv[argv.index("--output-build") + 1] == "GRCh38"
    assert argv[-2] == str(in_)
    assert argv[-1] == str(tmp_path / "o.vcf")
    # No --output-dir.
    assert "--output-dir" not in argv


def test_argv_for_plink_output_dir(tmp_path):
    in_ = _make_input(tmp_path)
    fake = _make_fake(tmp_path, "")
    out_dir = tmp_path / "d"
    argv = (
        Converter(
            input=in_,
            output_dir=out_dir,
            format=OutputFormat.PLINK,
            binary=fake,
        )
        .with_standardize()
        .with_variants_only()
        .with_sex(Sex.MALE)
        .argv()
    )
    assert "--standardize" in argv
    assert "--variants-only" in argv
    assert argv[argv.index("--sex") + 1] == "male"
    assert argv[argv.index("--format") + 1] == "plink"
    assert argv[argv.index("--output-dir") + 1] == str(out_dir)
    # Last positional is INPUT, no trailing OUTPUT since --output-dir was set.
    assert argv[-1] == str(in_)


def test_argv_includes_reference_and_panel(tmp_path):
    in_ = _make_input(tmp_path)
    ref = tmp_path / "ref.fa"
    fai = tmp_path / "ref.fa.fai"
    panel = tmp_path / "panel.vcf"
    for p in (ref, fai, panel):
        p.write_text("")
    fake = _make_fake(tmp_path, "")
    argv = (
        Converter(input=in_, output=tmp_path / "o.vcf", binary=fake)
        .with_reference(ref, fai)
        .with_panel(panel)
        .with_sample("ID1")
        .with_input_build("hg19")
        .argv()
    )
    assert argv[argv.index("--reference") + 1] == str(ref)
    assert argv[argv.index("--reference-fai") + 1] == str(fai)
    assert argv[argv.index("--panel") + 1] == str(panel)
    assert argv[argv.index("--sample") + 1] == "ID1"
    assert argv[argv.index("--input-build") + 1] == "GRCh37"


def test_argv_for_genome_studio_input(tmp_path):
    in_ = _make_input(tmp_path)
    fake = _make_fake(tmp_path, "")
    argv = (
        Converter(
            input=in_,
            output=tmp_path / "o.vcf",
            input_format=InputFormat.GENOME_STUDIO,
            binary=fake,
        ).argv()
    )
    assert argv[argv.index("--input-format") + 1] == "genome-studio"


def test_genome_studio_input_format_coercion():
    # convert(...) coerces the string form, including the kebab-case value and
    # the enum member name, to InputFormat.GENOME_STUDIO.
    from convert_genome._api import _coerce_enum

    assert _coerce_enum("genome-studio", InputFormat) is InputFormat.GENOME_STUDIO
    assert _coerce_enum("GENOME_STUDIO", InputFormat) is InputFormat.GENOME_STUDIO
    assert InputFormat.GENOME_STUDIO.value == "genome-studio"


def test_run_parses_genome_studio_report(tmp_path):
    # The report's input.format is surfaced verbatim, including "genome-studio".
    in_ = _make_input(tmp_path)
    out = tmp_path / "out.vcf"
    report = _minimal_report()
    report["input"] = {"path": "/in.csv", "format": "genome-studio", "origin": "local"}
    fake = _make_fake(tmp_path, _good_fake_body(report))

    result = convert(input=in_, output=out, binary=fake)
    assert result.input.format == "genome-studio"


# ---------------------------------------------------------------------------
# Run + JSON parsing
# ---------------------------------------------------------------------------


def test_run_parses_report(tmp_path):
    in_ = _make_input(tmp_path)
    out = tmp_path / "out.vcf"
    fake = _make_fake(tmp_path, _good_fake_body(_minimal_report()))

    result = convert(
        input=in_,
        output=out,
        binary=fake,
        standardize=True,
        assembly="hg38",
        sex="female",
    )
    assert isinstance(result, ConversionResult)
    assert result.statistics.emitted_records == 990
    assert result.statistics.total_records == 1000
    assert result.yield_rate == pytest.approx(0.99)
    assert result.sample.sex_inferred is True
    assert result.build_detection is not None
    assert result.build_detection.detected_build == "GRCh38"
    assert result.report_path is not None
    assert result.report_path.exists()
    assert any(p.suffix == ".vcf" for p in result.output_paths)


def test_run_with_output_dir(tmp_path):
    in_ = _make_input(tmp_path)
    out_dir = tmp_path / "outdir"
    fake = _make_fake(tmp_path, _good_fake_body(_minimal_report()))

    result = convert(
        input=in_,
        output_dir=out_dir,
        format=OutputFormat.VCF,
        binary=fake,
    )
    assert result.report_path == out_dir / "genotypes_report.json"
    assert (out_dir / "genotypes.vcf") in result.output_paths


def test_run_locates_report_via_log_line(tmp_path):
    """The CLI logs 'Wrote run report to <path>' — prefer that over guessing."""
    in_ = _make_input(tmp_path)
    out = tmp_path / "out.vcf"
    # Move the report to an unusual location, then advertise it via log.
    alt = tmp_path / "side_report.json"
    body = textwrap.dedent(
        f"""\
        #!/usr/bin/env python3
        import json, pathlib, sys
        pathlib.Path({str(alt)!r}).write_text(json.dumps({_minimal_report()!r}))
        pathlib.Path({str(out)!r}).write_text('##')
        print('Wrote run report to', {str(alt)!r})
        """
    )
    fake = _make_fake(tmp_path, body)
    result = convert(input=in_, output=out, binary=fake)
    assert result.report_path == alt


def test_failure_exit_code(tmp_path):
    in_ = _make_input(tmp_path)
    fake = _make_fake(tmp_path, _FAKE_FAILING)
    with pytest.raises(ConvertGenomeFailed) as ei:
        convert(input=in_, output=tmp_path / "o.vcf", binary=fake)
    assert ei.value.returncode == 2
    assert "synthetic failure" in ei.value.stderr


def test_missing_report_raises(tmp_path):
    in_ = _make_input(tmp_path)
    fake = _make_fake(tmp_path, _FAKE_NOREPORT)
    with pytest.raises(ReportNotFound):
        convert(input=in_, output=tmp_path / "o.vcf", binary=fake)


def test_unexpected_schema_field_is_tolerated(tmp_path):
    # Forward-compatibility: a report from a newer CLI that grows extra fields
    # is parsed by ignoring the unknown ones (see _api._only_known), not by
    # raising. Known fields still populate correctly.
    in_ = _make_input(tmp_path)
    fwd = _minimal_report()
    fwd["statistics"]["mystery_new_field"] = 42
    fwd["mystery_top_level"] = {"anything": True}
    fake = _make_fake(tmp_path, _good_fake_body(fwd))

    result = convert(input=in_, output=tmp_path / "o.vcf", binary=fake)
    assert result.statistics.emitted_records == 990
    assert not hasattr(result.statistics, "mystery_new_field")


def test_argv_log_actually_reflects_invocation(tmp_path):
    in_ = _make_input(tmp_path)
    fake = _make_fake(tmp_path, _good_fake_body(_minimal_report()))
    convert(
        input=in_,
        output=tmp_path / "o.vcf",
        binary=fake,
        format="bcf",
        standardize=True,
        variants_only=True,
        extra_args=["--log-level", "warn"],
    )
    argv = json.loads((tmp_path / "argv.json").read_text())
    assert argv[argv.index("--format") + 1] == "bcf"
    assert "--standardize" in argv
    assert "--variants-only" in argv
    assert argv[argv.index("--log-level") + 1] == "warn"


def test_error_hierarchy():
    assert issubclass(ConvertGenomeBinaryNotFound, ConvertGenomeError)
    assert issubclass(ConvertGenomeFailed, ConvertGenomeError)
    assert issubclass(InvalidConfig, ConvertGenomeError)
    assert issubclass(ReportNotFound, ConvertGenomeError)


def test_sex_indeterminate_maps_to_no_flag(tmp_path):
    """Regression: callers chaining through infer_sex may pass
    'indeterminate' or 'unknown'. Neither is a valid --sex value, but
    omitting the flag and letting the CLI infer is the right behaviour."""
    in_ = _make_input(tmp_path)
    fake = _make_fake(tmp_path, _good_fake_body(_minimal_report()))
    for value in ("indeterminate", "unknown", "INDETERMINATE", "  unknown  ", ""):
        convert(input=in_, output=tmp_path / "o.vcf", binary=fake, sex=value)
        argv = json.loads((tmp_path / "argv.json").read_text())
        assert "--sex" not in argv, f"sex={value!r} should not produce --sex flag"


def test_converter_is_immutable(tmp_path):
    in_ = _make_input(tmp_path)
    base = Converter(input=in_, output=tmp_path / "o.vcf")
    new = base.with_standardize()
    assert base.standardize is False
    assert new.standardize is True
    assert new is not base
