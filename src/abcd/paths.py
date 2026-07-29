"""
Path resolution for the ABCD pipeline.

Every path derives from two roots, so nothing depends on the working directory
a notebook happened to be launched from (the 5.1 code used ``../`` relatives
and broke whenever a script moved):

``REPO_ROOT``
    This repository.  Code, committed reference data, outputs.

``ABCD_ROOT``
    A tree holding raw release directories.  Resolved from, in order:
    :data:`REPO_ROOT` itself; the ``ABCD_ROOT`` environment variable;
    ``~/Git/ABCD`` if it exists.  :func:`abcd_root` raises with a clear
    message if none is available, rather than failing later with a confusing
    FileNotFoundError.

Releases are looked up across *all* of those roots rather than under a single
one (:func:`abcd_roots`, :func:`release_dir`).  That is what lets a release
live inside the repo -- ``abcd-data-release-7.0/``, gitignored -- so the common
case needs no environment variable at all, while 5.1 (2.4 GB) stays outside it
under ``~/Git/ABCD``.  A single-root design would force a choice between the
two; searching several means both resolve with nothing exported.
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


class SourceUnavailable(RuntimeError):
    """Raised when a required input is not present in this release.

    Lives here rather than in ``qc.py`` because both ``io.py`` (adapters) and
    ``qc.py`` (predicates) raise it, and ``paths`` is the one module both
    already depend on.  ``qc.SourceUnavailable`` remains a valid alias.
    """


def abcd_roots() -> tuple[Path, ...]:
    """Every tree that may hold a raw release, in search order.

    ``REPO_ROOT`` comes first so a release vendored into the repo wins without
    any environment variable.  ``ABCD_ROOT`` is honoured next -- an explicit
    export still overrides -- then the historical ``~/Git/ABCD`` location.

    Returns the roots that exist, deduplicated and order-preserving.  Never
    raises: callers that need at least one root call :func:`abcd_root`.
    """
    cands = [REPO_ROOT]
    env = os.environ.get("ABCD_ROOT")
    if env:
        p = Path(env).expanduser()
        if not p.exists():
            raise DataRootError(f"ABCD_ROOT={env} does not exist")
        cands.append(p)
    cands.extend(_ABCD_FALLBACKS)

    seen: dict[Path, None] = {}
    for c in cands:
        try:
            r = c.resolve()
        except OSError:  # pragma: no cover - unreadable mount
            continue
        if r.exists():
            seen.setdefault(r, None)
    return tuple(seen)


def abcd_root() -> Path:
    """The first tree that actually contains a release directory.

    Prefers a root holding releases over one that merely exists, so a bare
    ``~/Git/ABCD`` cannot shadow the repo-vendored release.
    """
    roots = abcd_roots()
    for r in roots:
        if any(
            (r / pat.format(r=rel)).exists()
            for rel in ("5.1", "7.0")
            for pat in _RELEASE_DIR_PATTERNS
        ):
            return r
    if roots:
        return roots[0]
    raise DataRootError(
        "Could not locate an ABCD data tree. Place a release directory in the "
        f"repo ({REPO_ROOT}), set the ABCD_ROOT environment variable, or place "
        f"it at one of: {list(_ABCD_FALLBACKS)}"
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


#: Release directory naming is not consistent across ABCD releases: 5.1 ships
#: as ``abcd-data-release-5.1/``, 7.0 as ``abcd-7.0/``.  Candidates are tried in
#: order, so a release can be re-homed by adding a pattern rather than editing
#: call sites.
_RELEASE_DIR_PATTERNS = (
    "abcd-data-release-{r}",
    "abcd-{r}",
    "abcd_data_release_{r}",
    "release-{r}",
)


def release_dir(release: str) -> Path:
    """Directory of a raw ABCD release.

    Handles the several release naming conventions, and searches every root in
    :func:`abcd_roots` -- so 7.0 vendored into the repo and 5.1 kept outside it
    both resolve without an environment variable.
    """
    roots = abcd_roots()
    if not roots:
        abcd_root()  # raises DataRootError with the actionable message
    tried: list[str] = []
    for root in roots:
        for pat in _RELEASE_DIR_PATTERNS:
            d = root / pat.format(r=release)
            tried.append(str(d))
            if d.exists():
                return d
    available = sorted(
        {
            p.name
            for root in roots
            for p in root.glob("*")
            if p.is_dir()
            and ("release" in p.name.lower() or p.name.startswith("abcd-"))
        }
    )
    raise DataRootError(
        f"Release {release} not found in any of {[str(r) for r in roots]}. "
        f"Available release directories: {available}. Tried: {tried}"
    )


def find_in_roots(name: str) -> Path | None:
    """Locate a loose analysis file across every root in :func:`abcd_roots`.

    Auxiliary inputs (static subject lists, GWAS sumstats) sit beside the
    release directories rather than inside them, and need not live in the same
    tree as the release currently being analysed.  Returns ``None`` when absent
    so callers can degrade -- the legacy QC policy, for instance, reports the
    input as unavailable rather than failing the run.
    """
    for root in abcd_roots():
        p = root / name
        if p.exists():
            return p
    return None


def find_table(root: Path, stem: str) -> Path:
    """Locate a release table by filename stem, wherever it sits in the tree.

    Release 7.0's documented layout nests each table in its own directory
    (``y/qc__incl/mr_y_qc__incl.tsv``), but a partial download may place them
    flat (``y/mr_y_qc__incl.tsv``).  Both are valid; searching by stem means the
    adapter does not care which, and a missing table raises with the paths
    actually searched rather than a bare FileNotFoundError.
    """
    for suffix in (".tsv", ".csv", ".parquet"):
        direct = root / f"{stem}{suffix}"
        if direct.exists():
            return direct
    hits = sorted(
        p for p in root.rglob(f"{stem}.*")
        if p.suffix in {".tsv", ".csv", ".parquet"}
    )
    if not hits:
        raise DataRootError(
            f"table {stem!r} not found anywhere under {root}. "
            "Check the release download is complete."
        )
    # prefer .tsv, then shallowest path, for determinism
    hits.sort(key=lambda p: ({".tsv": 0, ".csv": 1, ".parquet": 2}[p.suffix],
                             len(p.relative_to(root).parts)))
    return hits[0]


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
