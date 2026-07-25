"""
ABCD longitudinal brain-development pipeline.

Subject-level developmental *change* phenotypes from repeated structural MRI,
built to be linked to the AHBA C3 spatial programme and the single-cell PC1
maturation axis already in this repository.

Layout follows the repo convention that computation is separate from plotting:

    config.py     RunConfig - every analysis decision, content-hashed
    paths.py      root resolution; no cwd-relative paths anywhere
    io.py         ReleaseAdapter (5.1 concrete, 7.0 stub) - raw table access
    qc.py         QC as logged, composable predicates
    assemble.py   tidy long-format assembly to Parquet
    phenotype.py  subject-level slope phenotypes + reliability
    spatial.py    parcellation handling, spin tests
    plotting.py   figures only

Model fitting lives in ``R/`` and communicates via Parquet, never rpy2.
"""

from .config import RunConfig, ConfigError  # noqa: F401
from . import paths  # noqa: F401

__all__ = ["RunConfig", "ConfigError", "paths"]
