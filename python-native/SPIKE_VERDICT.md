# pyo3 native bindings for `convert_genome` — SPIKE verdict

**Branch:** `spike/pyo3-native-bindings` (SauersML/convert_genome)
**Question:** Is building real pyo3 in-process bindings worth it, vs. the existing
typed *subprocess* wrapper in `python/convert_genome/_api.py`?
**Verdict (one line):** **KEEP THE SUBPROCESS WRAPPER.** Native bindings are
real and work, but the measured win only exists for many-small-calls (not the
real WGS workload), and the deploy/CI cost goes up materially. Cherry-pick ONE
idea from the spike instead (below).

---

## What was built (real, not speculation)

A working pyo3 crate `python-native/` (crate-type `cdylib`, detached from the CLI
build via its own `[workspace]`) that wraps the existing library functions
**without reimplementing any pipeline logic**:

- `convert(input, output, assembly=…, …)` → wraps `convert_dtc_file`, returns the
  `ConversionSummary` as an in-memory dict (no `_report.json` sidecar).
- `infer_sex(input) → {sex, confidence, metrics}` → wraps
  `infer_sex_detailed_from_variant_file` (VCF/BCF) / `infer_sex_from_records` (DTC).
- `detect_build(input) → {build, hg19_match_rate, hg38_match_rate, confidence}` →
  wraps `detect_build_from_variant_file` / `detect_build_from_dtc`.

It builds and installs cleanly: `maturin develop --release` produces a
`cp39-abi3` wheel (one wheel covers all Python ≥ 3.9). Errors surface as Python
exceptions; heavy work releases the GIL (`py.detach`).

Two tiny, behavior-preserving library changes were needed (both legit "reuse"):
- `inference.rs`: added `SexInference` struct + `infer_sex_detailed_from_variant_file`;
  the old `infer_sex_from_variant_file` now delegates to it. (Old fn discarded the
  evidence metrics — the wrapper literally cannot get confidence from a standalone call.)
- `conversion.rs`: made `prescan_dtc_records` `pub` so the DTC inference path is reachable.

---

## Benchmarks (real inputs, this machine: M-series arm64, release builds)

Inputs: `genome_23andme_grch37.txt` (4,000 variants, DTC),
`genome_grch38.vcf.gz` (400 variants), a 199,879-variant slice of the real WGS
`real_genome.vcf.gz`.

### (a) infer_sex latency — where fork/exec would show
| call | median |
|---|---|
| native infer_sex, DTC 4k | **1.5 ms** |
| native infer_sex, VCF 400 | **0.5 ms** |
| bare CLI spawn floor (`--version`, *no work at all*) | 7.7 ms |
| 200× native infer_sex | 312 ms total (1.56 ms/call) |
| 200× bare CLI spawn (no work) | 1505 ms total (7.52 ms/call) |

Native is **~4.8× faster than the fork/exec floor alone** for small DTC files —
and the subprocess wrapper has **no standalone infer_sex at all** (it only runs
inference as a side effect of a full `convert`, then parses it back out of JSON).

**BUT** on a realistic ~200k-variant VCF, native infer_sex is **~1040 ms/call** —
the file parsing dominates and the 7.7 ms spawn cost is **<1%**. The fork/exec
saving is irrelevant at real WGS scale; it only matters if you call inference
thousands of times on tiny DTC files.

### (b) full convert latency (DTC 4k → VCF)
| path | median | min | max |
|---|---|---|---|
| native (in-process) | **25.0 ms** | 24.6 | 501* |
| subprocess wrapper | 35.9 ms | 35.1 | 39.9 |

Native is ~11 ms faster per convert (the spawn floor + sidecar write/parse).
*The 501 ms native outlier = first-call reference cache warmup happening
in-process. For a real WGS convert (seconds-to-minutes of actual work) an 11 ms
difference is noise.

### (c) in-memory return
Native `convert` returns a dict directly (all 11 summary stat fields). The wrapper
must have the CLI **write** `<stem>_report.json`, then **read+parse** it back into
dataclasses. Native *can* return data Python-side without a disk round-trip —
**but it does not return variant records in memory**: the pipeline writes VCF/BCF
to disk regardless (Beagle consumes files), so the only thing native saves is the
small JSON sidecar, not the bulk output. The hypothesis that "the pipeline is
file-oriented so in-memory return adds little" is **confirmed**.

---

## Ergonomics

Honest answer: **basically a wash.** Both are `import x; x.f(...)`.
- Native is marginally nicer for the *standalone* `infer_sex` / `detect_build`
  (the wrapper exposes neither — you'd have to run a full convert and dig the
  result out of the report). This is the one genuine ergonomic gap.
- The wrapper is *richer* for convert: immutable builder (`Converter`,
  `.with_*()`), typed frozen dataclasses, eager arg-combo validation, full
  `RunReport` schema. The native `convert` returns only the bare summary; to reach
  parity you'd re-expose the whole report struct through pyo3 (more wrapping work).

## Deploy / CI (the deciding factor)

The native extension links the **same** C deps as the CLI binary
(`libcurl`, openssl/`Security.framework`, libz — from the reference-download path);
`cargo tree` shows `curl-sys`, `openssl-sys`, `libz-sys`. So the wheel adds **no
new** C-dep burden beyond the binary — but the **packaging mechanism gets much
heavier**:

| | today (wrapper) | with native bindings |
|---|---|---|
| CI to publish | one `python -m build` job, **pure-Python universal wheel**, no Rust toolchain | maturin **per-platform matrix** (manylinux arm64 + x86_64 + macOS), Rust toolchain, **auditwheel** to bundle libcurl/openssl into the manylinux wheel |
| install on arm64-Linux (AWS Batch) | trivial, any Python, seconds | needs the matching prebuilt arm64 wheel, or it tries to compile from sdist (Rust + C dev headers on the box) |
| failure modes | ~none | abi3 tag mismatches, missing manylinux wheel → sdist compile, openssl version skew |

The repo ships the wrapper **today** as exactly that pure-Python wheel
(`.github/workflows/pypi-publish.yml`), and ships the CLI binary separately. The
AWS Batch image already has the binary, so the wrapper has **zero** extra deploy
surface. Native bindings would add an ongoing compiled-wheel matrix to maintain.

## Correctness / maintenance

No correctness advantage either way — both call the identical Rust code paths.
Maintenance *favors the wrapper*: native bindings are a second public API surface
to keep in sync with the report schema, plus a cross-compile/auditwheel pipeline.

---

## Recommendation

**KEEP the subprocess wrapper.** pyo3 is not worth a full migration here because:
1. The latency win (fork/exec + JSON sidecar, ~11 ms/convert, ~7 ms/infer) is
   **negligible at real WGS scale** — confirmed: a real ~200k-variant infer_sex
   is ~1 s, making spawn <1%. It only helps for thousands of tiny-DTC calls,
   which is not the production workload.
2. In-memory record return — the main thing native fundamentally enables — is
   **moot**: the pipeline writes VCF/BCF to disk for Beagle anyway.
3. Ergonomics are a wash; the wrapper is actually richer for `convert`.
4. Deploy/CI cost goes **up** (pure-Python wheel → per-arch compiled-wheel matrix
   + auditwheel for libcurl/openssl), for no production speed benefit.

**Cheap win to cherry-pick instead (no pyo3):** the one real gap is that the
wrapper has **no standalone `infer_sex` / `detect_build`**. Add `--mode infer-sex`
/ `--mode detect-build` subcommands to the CLI that print the JSON the native
functions return, and add thin `infer_sex()` / `detect_build()` helpers to the
*subprocess* wrapper. That closes the only genuine ergonomic gap, keeps the
pure-Python wheel, and ships today.

**If** a future workload ever does need thousands of inference calls per process
(e.g. a long-lived service classifying many small DTC uploads), revisit native
bindings for `infer_sex`/`detect_build` only — the spike code here is ready and
the 4.8× per-call win is real for that specific shape.
