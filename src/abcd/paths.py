"""
Path resolution for the ABCD pipeline.

Every path derives from two roots, so nothing depends on the working directory
a notebook happened to be launched from (the 5.1 code used ``../`` relatives
and broke whenever a script moved):

``REPO_ROOT``
    This repository.  Code, committed reference data, outputs.

``ABCD_ROOT``
    The ABCD analysis tree holding the raw release directories.  Not in the
    repo (it is large and access-controlled).  Resolved from, in order:
    the ``ABCD_ROOT`` environment variable; ``~/Git/ABCD`` if it exists.
    :func:`abcd_root` raises with a clear message if neither is available,
    rather than failing later with a confusing FileNotFoundError.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"
OUT_DIR = REPO_ROOT / "out"
FIG_DIR = REPO_ROOT / "figures"
CONFIG_DIR = REPO_ROOT / "configs"

_ABCD_FALLBACKS = (Path.home() / "Git" / "ABCD",)


class DataRootError(FileNotFoundError):
    """Raised when the ABCD analysis tree cannot be located."""


def abcd_root() -> Path:
    """Locate the ABCD analysis tree (raw releases, QC lists, GWAS sumstats)."""
    env = os.environ.get("ABCD_ROOT")
    if env:
        p = Path(env).expanduser()
        if not p.exists():
            raise DataRootError(f"ABCD_ROOT={env} does not exist")
        return p
    for p in _ABCD_FALLBACKS:
        if p.exists():
            return p
    raise DataRootError(
        "Could not locate the ABCD analysis tree. Set the ABCD_ROOT "
        f"environment variable, or place it at one of: {list(_ABCD_FALLBACKS)}"
    )


_AHBA_FALLBACKS = (
    Path.home() / "Git" / "AHBA" / "data" / "abagen-data" / "expression",
)


def ahba_expression_path(parcellation: str = "hcp") -> Path:
    """Locate the AHBA gene-by-region expression matrix.

    These matrices are not redistributed with this repo: they are the
    donor-normalised, probe-selected output of the pipeline in Dear et al.
    (2024) and are several MB each.  Set ``ABCD_AHBA_DIR`` to the directory
    holding ``{dk,hcp}_3d_ds5.csv``, or place it at the fallback path.
    """
    env = os.environ.get("ABCD_AHBA_DIR")
    roots = (Path(env).expanduser(),) if env else _AHBA_FALLBACKS
    fname = f"{parcellation}_3d_ds5.csv"
    for r in roots:
        p = r / fname
        if p.exists():
            return p
    raise DataRootError(
        f"Could not locate AHBA expression matrix {fname}. Set ABCD_AHBA_DIR "
        f"to the directory containing it (looked in: {[str(r) for r in roots]})."
    )


def release_dir(release: str) -> Path:
    """Directory of a raw ABCD release, e.g. ``abcd-data-release-5.1/``."""
    d = abcd_root() / f"abcd-data-release-{release}"
    if not d.exists():
        raise DataRootError(
            f"Release {release} not found at {d}. "
            f"Available: {sorted(p.name for p in abcd_root().glob('abcd-data-release-*'))}"
        )
    return d


def run_dir(run_id: str, create: bool = True) -> Path:
    """Output directory for one config, holding its manifest and results."""
    d = OUT_DIR / run_id
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


def interim_dir(create: bool = True) -> Path:
    """Assembled long-format tables, shared across runs of the same release."""
    d = OUT_DIR / "interim"
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d
