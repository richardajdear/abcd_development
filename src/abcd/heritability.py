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

__all__ = ["pair_table", "falconer", "falconer_by_region",
           "MZ", "DZ", "DZ_TWIN", "SIB"]

MZ = "MZ"
#: Pooled DZ-twin + non-twin-sibling class.  Only 5.1 forces this pooling; on
#: 7.0 the two are separable and pooling biases h2 upward (see ``falconer``).
DZ = "DZ_or_sib"
DZ_TWIN = "DZ_twin"
SIB = "full_sib"

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
    """Explicit relative pairs with genotype-confirmed zygosity.

    Returns one row per *unordered* pair: ``a``, ``b`` (subject IDs), plus
    ``family_id``, ``pihat`` and ``pair_type``.  On 7.0, ``pair_type`` is one of
    ``MZ``, ``DZ_twin`` or ``full_sib``, read from the release's own pi-hat
    table (:meth:`io.Release70Adapter.genotyped_pairs`).  On 5.1 the source can
    only distinguish ``MZ`` from a pooled ``DZ_or_sib``.

    Two design notes, both of which changed answers in this project:

    * **Pairs are explicit, not inferred from family size.**  An earlier version
      returned subject -> family_id and let the caller treat any family with two
      phenotyped members as a pair.  That drops families holding a twin pair
      plus a third sibling, and it cannot represent a family containing two
      different pair classes.  On the settled 7.0 run, using explicit pairs
      raises the usable MZ count from 260 to 266 and the DZ+sibling count from
      871 to 967.
    * **DZ twins are not pooled with non-twin siblings.**  See :func:`falconer`.
    """
    if release == "7.0":
        return io.Release70Adapter().genotyped_pairs()

    r51 = io.Release51Adapter().relatedness()
    if "pair_type" not in r51.columns:
        raise KeyError(
            f"5.1 relatedness has no 'pair_type' column (got {list(r51.columns)}); "
            "zygosity cannot be established"
        )
    # 5.1 ships subject -> (family, pair_type) rather than explicit pairs, so
    # pairs are reconstructed within family among same-class members.  This is
    # an approximation the 7.0 path does not need.
    z = r51.loc[r51.pair_type.isin([MZ, DZ]), ["subject", "family_id", "pair_type"]]
    rows = []
    for (fam, cls), g in z.groupby(["family_id", "pair_type"]):
        subs = sorted(g.subject)
        for i in range(0, len(subs) - 1, 2):
            rows.append({"a": subs[i], "b": subs[i + 1], "family_id": fam,
                         "pihat": np.nan, "pair_type": cls})
    return pd.DataFrame(rows, columns=["a", "b", "family_id", "pihat", "pair_type"])


def _pair_matrix(y: pd.Series, pairs: pd.DataFrame, pair_type: str | list[str]) -> np.ndarray:
    """(n_pairs, 2) array of phenotype values for complete pairs of one class.

    ``pair_type`` may be a list, which pools those classes into one correlation.
    """
    want = [pair_type] if isinstance(pair_type, str) else list(pair_type)
    d = pairs.loc[pairs.pair_type.isin(want), ["a", "b"]]
    v = np.column_stack([y.reindex(d.a).to_numpy(dtype=float),
                         y.reindex(d.b).to_numpy(dtype=float)])
    return v[~np.isnan(v).any(axis=1)]


def falconer(
    y: pd.Series,
    pairs: pd.DataFrame | None = None,
    *,
    dz_class: str = "twins_only",
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
    dz_class
        Which pairs form the second Falconer class.  ``"twins_only"`` (default)
        uses DZ twins alone.  ``"pooled"`` adds non-twin full siblings, which is
        all the 5.1 source could do -- it roughly doubles the DZ sample but
        biases h2 *upward*, because siblings correlate less than DZ twins at the
        same 0.5 genetic sharing (differing ages, non-shared gestation).  On the
        7.0 slope phenotype the gap is r = 0.279 (DZ twins) vs 0.090 (siblings),
        difference +0.189 with a bootstrap 95% CI of [+0.046, +0.321], so the
        pooling is not a harmless convenience: it moves global-slope h2 from
        0.44 to 0.64.  Use ``"pooled"`` only to reproduce older numbers.
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

    have = set(pairs.pair_type.unique())
    if dz_class not in ("twins_only", "pooled"):
        raise ValueError(f"dz_class must be 'twins_only' or 'pooled', got {dz_class!r}")
    if DZ in have:  # 5.1-style pooled source; the choice does not exist
        dz_spec: str | list[str] = DZ
        dz_name = DZ
    elif dz_class == "twins_only":
        dz_spec, dz_name = DZ_TWIN, DZ_TWIN
    else:
        dz_spec, dz_name = [DZ_TWIN, SIB], DZ
    if dz_spec != DZ and not ({DZ_TWIN} & have):
        raise KeyError(
            f"pair table has classes {sorted(have)}; expected {DZ_TWIN} or {DZ}"
        )

    mats = {MZ: _pair_matrix(y, pairs, MZ), dz_name: _pair_matrix(y, pairs, dz_spec)}
    for pt, A in mats.items():
        if len(A) < min_pairs:
            raise ValueError(
                f"only {len(A)} complete {pt} pairs (need >= {min_pairs}); "
                "check that the phenotype index is subject IDs matching pair_table()"
            )

    r = {pt: float(pearsonr(A[:, 0], A[:, 1]).statistic) for pt, A in mats.items()}
    h2 = (r[MZ] - r[dz_name]) / _SHARING_GAP
    out: dict[str, Any] = {
        "h2": h2,
        "r_MZ": r[MZ],
        "r_DZ": r[dz_name],
        "n_MZ": len(mats[MZ]),
        "n_DZ": len(mats[dz_name]),
        "dz_class": dz_name,
    }
    if SIB in have:  # report the class we did not use, so the choice is auditable
        other = _pair_matrix(y, pairs, SIB if dz_name == DZ_TWIN else [DZ_TWIN, SIB])
        if len(other) >= min_pairs:
            out["r_sib_excluded" if dz_name == DZ_TWIN else "r_pooled"] = float(
                pearsonr(other[:, 0], other[:, 1]).statistic
            )

    if n_boot:
        rng = np.random.default_rng(seed)
        draws = np.empty(n_boot)
        for b in range(n_boot):
            rb = {}
            for pt, A in mats.items():
                idx = rng.integers(0, len(A), len(A))
                rb[pt] = pearsonr(A[idx, 0], A[idx, 1]).statistic
            draws[b] = (rb[MZ] - rb[dz_name]) / _SHARING_GAP
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
