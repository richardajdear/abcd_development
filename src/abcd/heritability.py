"""Twin-design heritability from pipeline phenotypes.

Falconer's estimator, h2 = 2 * (r_MZ - r_DZ), applied to the subject-level
phenotypes written by ``abcd.phenotype``.  This buys a heritability estimate
without genotypes, a GRM, or a GWAS: relatedness is known from the study design
rather than estimated from data, so a few hundred twin pairs answer what would
otherwise need tens of thousands of unrelated subjects.

Two things about the estimator that have already caused one wrong answer in this
project, both documented at length in ``docs/REPORT_7.0.md``:

1.  It is only valid on runs fitted WITHOUT a family random effect.  A model
    containing ``(1 | family_id)`` absorbs the between-family variance, leaving
    subject BLUPs that are within-family deviations; siblings' deviations then
    anti-correlate by construction.  :func:`falconer` refuses to guess about
    this -- pass phenotypes from a ``family_effect: false`` run.

2.  Because it *differences* two correlations, artefacts that affect both pair
    types partly cancel, so a plausible-looking h2 can come from invalid inputs.
    Every function here therefore returns ``r_MZ`` and ``r_DZ`` alongside the
    estimate.  Always look at them.  A negative or near-zero ``r_DZ`` on a
    familial trait means something is wrong upstream, whatever h2 says.

The DZ class in ABCD is diluted with non-twin full siblings (different ages,
lower phenotypic correlation), which biases h2 upward.  These estimates are a
screening tool -- "is this phenotype worth cluster time?" -- not a final value.
GCTA-GRM/REML is the defensible replacement.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from . import io

__all__ = ["pair_table", "falconer", "falconer_by_region", "MZ", "DZ"]

MZ = "MZ"
DZ = "DZ_or_sib"

#: Genetic sharing difference between the two pair classes.  Falconer's factor
#: of 2 is 1 / (1.0 - 0.5); named here so the assumption is visible.
_SHARING_GAP = 0.5


def _normalise_ids(s: pd.Series) -> pd.Series:
    """Reconcile subject-ID spellings across releases and tables.

    7.0's ``participant_id`` is ``sub-0A4P0LWM`` (no NDARINV token) while the
    imaging tables and the 5.1 relatedness file use ``NDAR_INV...``/``NDARINV...``.
    Joins that skip this silently return zero overlap.
    """
    return (
        s.astype(str)
        .str.replace("sub-", "", regex=False)
        .str.replace("NDAR_INV", "NDARINV", regex=False)
    )


def pair_table(release: str = "7.0") -> pd.DataFrame:
    """Subject -> (family_id, pair_type) using genotype-confirmed zygosity.

    7.0 does not ship pi-hat, so zygosity is taken from the 5.1 genotype file
    and matched on normalised subject ID; subjects without a 5.1 record are
    dropped.  Returns only ``MZ`` and ``DZ_or_sib`` rows -- the classes Falconer
    needs -- indexed by the *release's* subject ID so it joins straight onto
    phenotype tables.
    """
    r51 = io.Release51Adapter().relatedness()
    if "pair_type" not in r51.columns:
        raise KeyError(
            f"5.1 relatedness has no 'pair_type' column (got {list(r51.columns)}); "
            "zygosity cannot be established"
        )
    zyg = r51.assign(_k=_normalise_ids(r51.subject))
    zyg = zyg.loc[zyg.pair_type.isin([MZ, DZ]), ["_k", "pair_type"]]

    adapter = io.Release70Adapter() if release == "7.0" else io.Release51Adapter()
    rel = adapter.relatedness()
    out = (
        rel.assign(_k=_normalise_ids(rel.subject))
        .merge(zyg, on="_k", how="inner", suffixes=("_design", ""))
        .set_index("subject")
    )
    keep = ["family_id", "pair_type"]
    return out[keep].dropna(subset=["family_id"])


def _pair_matrix(y: pd.Series, pairs: pd.DataFrame, pair_type: str) -> np.ndarray:
    """(n_pairs, 2) array of phenotype values for complete pairs of one class."""
    d = pairs.loc[pairs.pair_type == pair_type, ["family_id"]].join(
        y.rename("_y"), how="inner"
    )
    d = d.dropna(subset=["_y"])
    grouped = d.groupby("family_id")["_y"].apply(list)
    return np.array([v for v in grouped if len(v) == 2], dtype=float)


def falconer(
    y: pd.Series,
    pairs: pd.DataFrame | None = None,
    *,
    n_boot: int = 0,
    seed: int = 0,
    min_pairs: int = 20,
) -> dict[str, Any]:
    """Falconer h2 for one subject-level phenotype.

    Parameters
    ----------
    y
        Phenotype indexed by subject ID (one value per subject).
    pairs
        Output of :func:`pair_table`; built on demand if omitted.
    n_boot
        Bootstrap resamples over pairs for a 95% interval.  0 skips it.
    min_pairs
        Refuse to estimate below this many complete pairs in either class.

    Returns
    -------
    dict with ``h2``, ``r_MZ``, ``r_DZ``, ``n_MZ``, ``n_DZ`` and -- if
    ``n_boot`` -- ``ci_lo``/``ci_hi``.  **Read r_MZ and r_DZ, not just h2**: the
    estimator differences them, so it cannot detect an artefact that shifts both.
    """
    if pairs is None:
        pairs = pair_table()
    y = y[~y.index.duplicated()]

    mats = {pt: _pair_matrix(y, pairs, pt) for pt in (MZ, DZ)}
    for pt, A in mats.items():
        if len(A) < min_pairs:
            raise ValueError(
                f"only {len(A)} complete {pt} pairs (need >= {min_pairs}); "
                "check that the phenotype index is subject IDs matching pair_table()"
            )

    r = {pt: float(pearsonr(A[:, 0], A[:, 1]).statistic) for pt, A in mats.items()}
    h2 = (r[MZ] - r[DZ]) / _SHARING_GAP
    out: dict[str, Any] = {
        "h2": h2,
        "r_MZ": r[MZ],
        "r_DZ": r[DZ],
        "n_MZ": len(mats[MZ]),
        "n_DZ": len(mats[DZ]),
    }

    if n_boot:
        rng = np.random.default_rng(seed)
        draws = np.empty(n_boot)
        for b in range(n_boot):
            rb = {}
            for pt, A in mats.items():
                idx = rng.integers(0, len(A), len(A))
                rb[pt] = pearsonr(A[idx, 0], A[idx, 1]).statistic
            draws[b] = (rb[MZ] - rb[DZ]) / _SHARING_GAP
        out["ci_lo"], out["ci_hi"] = (float(v) for v in np.percentile(draws, [2.5, 97.5]))
    return out


def falconer_by_region(
    ph: pd.DataFrame,
    pairs: pd.DataFrame | None = None,
    *,
    phenotype: str = "slope",
    min_pairs: int = 30,
) -> pd.DataFrame:
    """Per-region Falconer h2 from a long-format phenotype table.

    Regions with too few complete pairs get NaN rather than an estimate.  Note
    that per-region estimates are noisy (split-half rho ~0.28 in 7.0) -- use the
    spatial pattern, not individual values.
    """
    if pairs is None:
        pairs = pair_table()
    sub = ph[ph.phenotype == phenotype]
    wide = sub.pivot_table(index="subject", columns="label", values="value")

    rows = []
    for label in wide.columns:
        try:
            res = falconer(wide[label], pairs, min_pairs=min_pairs)
        except ValueError:
            res = {"h2": np.nan, "r_MZ": np.nan, "r_DZ": np.nan, "n_MZ": 0, "n_DZ": 0}
        rows.append({"label": label, **res})
    return pd.DataFrame(rows).set_index("label")
