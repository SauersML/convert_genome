#!/usr/bin/env python3
"""Empirical benchmark: native pyo3 bindings vs the subprocess wrapper.

Run with the spike venv:  .venv/bin/python bench.py
"""
import json
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import convert_genome_native as native
import convert_genome as wrapper

BIN = str(Path(__file__).resolve().parent.parent / "target/release/convert_genome")
DTC = "/tmp/e2e_genomes/genome_23andme_grch37.txt"      # 4000 variants
VCF = "/tmp/e2e_genomes/genome_grch38.vcf.gz"            # 400 variants


def timeit(fn, n):
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return ts


def summ(ts):
    ts_ms = [t * 1000 for t in ts]
    return {
        "n": len(ts),
        "mean_ms": round(statistics.mean(ts_ms), 3),
        "median_ms": round(statistics.median(ts_ms), 3),
        "min_ms": round(min(ts_ms), 3),
        "max_ms": round(max(ts_ms), 3),
    }


results = {}

# ---------------------------------------------------------------------------
# (a) infer_sex latency — many calls. This is where fork/exec would show.
#     The subprocess wrapper has NO standalone infer_sex, so the cheapest
#     subprocess equivalent is spawning the CLI binary at all. We measure:
#       - native infer_sex (in-process)
#       - bare subprocess spawn floor (`--version`): the unavoidable
#         fork/exec + dynamic-link + process-teardown cost per call.
# ---------------------------------------------------------------------------
print("== (a) infer_sex latency (DTC, 4000 variants) ==", file=sys.stderr)

results["native_infer_sex_dtc"] = summ(
    timeit(lambda: native.infer_sex(DTC, input_format="dtc", build="GRCh37"), 50)
)
results["native_infer_sex_vcf"] = summ(
    timeit(lambda: native.infer_sex(VCF, input_format="vcf", build="GRCh38"), 50)
)

# bare process spawn floor
results["subprocess_spawn_floor"] = summ(
    timeit(lambda: subprocess.run([BIN, "--version"], capture_output=True), 50)
)

# ---------------------------------------------------------------------------
# (b) full convert latency: native in-process vs subprocess wrapper.
#     Same input, same output target, declared build (skip detection) so
#     both do identical work.
# ---------------------------------------------------------------------------
print("== (b) full convert latency (DTC, 4000 variants) ==", file=sys.stderr)


def native_convert():
    with tempfile.TemporaryDirectory() as d:
        native.convert(DTC, str(Path(d) / "o.vcf"),
                       assembly="GRCh37", input_format="dtc", input_build="GRCh37")


def wrapper_convert():
    with tempfile.TemporaryDirectory() as d:
        wrapper.convert(input=DTC, output=str(Path(d) / "o.vcf"),
                        assembly="GRCh37", input_format="dtc",
                        input_build="GRCh37", binary=BIN, capture=True)


results["native_convert_dtc"] = summ(timeit(native_convert, 15))
results["wrapper_convert_dtc"] = summ(timeit(wrapper_convert, 15))

# ---------------------------------------------------------------------------
# (c) in-memory return: native returns a dict directly; the wrapper must
#     write + parse a _report.json sidecar. Demonstrate both shapes.
# ---------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as d:
    nat = native.convert(DTC, str(Path(d) / "o.vcf"),
                         assembly="GRCh37", input_format="dtc", input_build="GRCh37")
    results["native_convert_returns"] = {"type": "dict in-memory", "keys": sorted(nat.keys())}
with tempfile.TemporaryDirectory() as d:
    wr = wrapper.convert(input=DTC, output=str(Path(d) / "o.vcf"),
                         assembly="GRCh37", input_format="dtc",
                         input_build="GRCh37", binary=BIN)
    results["wrapper_convert_returns"] = {
        "type": "dataclass parsed from _report.json sidecar",
        "emitted_records": wr.emitted_records,
        "sample_sex": wr.sample.sex,
    }

print(json.dumps(results, indent=2))
