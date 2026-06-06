"""Native in-process bindings for convert_genome (pyo3 spike).

Re-exports the Rust extension's functions at the package top level so
callers can do ``import convert_genome_native as cg; cg.infer_sex(...)``.
"""

from ._native import convert, detect_build, infer_sex

__all__ = ["convert", "detect_build", "infer_sex"]
