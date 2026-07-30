"""
Structural covariance of developmental slopes, and its principal components.

The structural covariance matrix here is the region-by-region correlation of
subject slopes: entry (i, j) is the correlation across subjects between region
i's rate of change and region j's.  It answers "which regions mature together
across individuals", which is a different question from "which regions mature
fastest" (that is the group map in :mod:`abcd.maps`).

An identity worth stating plainly
---------------------------------
The 5.1 thesis took its "structural covariance PCs" by running PCA on the
68x68 correlation matrix ``R``, treating each region as an observation.  Those
are **the same vectors** as the region loadings of a PCA on the standardised
subject-by-region slope matrix, to within numerical noise (|r| > 0.996 for
PC1-3 on 7.0).  This is algebra, not coincidence: PCA on ``R`` diagonalises
``R'R``, and since ``R = V L V'`` is symmetric, ``R'R = V L^2 V'`` -- the same
eigenvectors, with squared eigenvalues.

So the two analyses are one analysis.  ``sc_pcs`` and ``slope_pcs`` are both
provided because the report compares them, but the comparison is a
verification that the pipeline is self-consistent, not independent evidence.

The squaring is not harmless for reporting.  PCA on ``R`` reports
``L^2 / sum(L^2)`` as its explained-variance ratio, which on 7.0 gives PC1 =
21% -- but this is a share of squared eigenvalues, not of variance.  The
variance PC1 actually explains in the slope data is ``L / sum(L)`` = 7.1%.
Quoting the 21% figure overstates the dominance of the leading component
threefold.  :func:`sc_pcs` therefore returns ``L / sum(L)`` and names it
``var_explained``, with the inflated quantity available as
``var_explained_sq_convention`` for comparison against the thesis text.

Sign convention
---------------
Eigenvector signs are arbitrary.  Rather than hard-code the thesis's
``flip=[1, 1, -1]``, :func:`sc_pcs` orients each component deterministically:
the hemisphere-averaged loading is made to correlate positively with an
anterior-posterior axis, so re-running on new data cannot silently flip a map
and invert a reported correlation.  ``orient=False`` disables this.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .phenotype import load_fits


def slope_matrix(run_dir: str | Path, fits_name: str = "fits",
                 value: str = "re_slope") -> pd.DataFrame:
    """Subject-by-region matrix of slope BLUPs, from one fitted run."""
    bl = load_fits(run_dir, fits_name)["blups"]
    W = bl.pivot(index="subject", columns="label", values=value)
    return W


def sc_matrix(W: pd.DataFrame, min_overlap: int = 100) -> pd.DataFrame:
    """Region-by-region correlation of subject slopes.

    Pairwise-complete, like the thesis code's ``cor(use='p')``.  Region pairs
    sharing fewer than ``min_overlap`` subjects are set to NaN rather than
    given a correlation nobody should trust -- the thesis code would have
    returned a value computed on however many subjects happened to overlap.
    In practice the LMM fits every region on every retained subject, so the
    matrix is complete and this guard never fires; it exists so that changing
    the region set later cannot silently introduce unstable entries.
    """
    R = W.corr(min_periods=min_overlap)
    if R.isna().to_numpy().sum() > 0:
        n_bad = int(R.isna().to_numpy().sum() // 2)
        raise ValueError(
            f"{n_bad} region pairs have fewer than {min_overlap} shared "
            "subjects -- inspect the slope matrix before taking PCs"
        )
    return R


@dataclass
class PCResult:
    """Region loadings plus the variance they explain."""

    loadings: pd.DataFrame          # region x component
    var_explained: np.ndarray       # L / sum(L) -- the honest quantity
    var_explained_sq_convention: np.ndarray | None = None


def _ap_axis(labels) -> np.ndarray:
    """Crude anterior-posterior ordering used only to fix component signs.

    Frontal regions score +1, occipital -1, everything else 0.  This is a sign
    convention, not an analysis: it never enters a correlation or a test.
    """
    ant = {"superiorfrontal", "rostralmiddlefrontal", "caudalmiddlefrontal",
           "parsopercularis", "parstriangularis", "parsorbitalis",
           "lateralorbitofrontal", "medialorbitofrontal", "frontalpole",
           "rostralanteriorcingulate"}
    post = {"lateraloccipital", "lingual", "cuneus", "pericalcarine",
            "superiorparietal", "inferiorparietal", "precuneus"}
    out = []
    for lab in labels:
        reg = lab.split("_", 1)[1] if "_" in lab else lab
        out.append(1.0 if reg in ant else (-1.0 if reg in post else 0.0))
    return np.asarray(out)


def _orient(V: np.ndarray, labels) -> np.ndarray:
    """Flip each column so it correlates positively with the A-P axis."""
    ap = _ap_axis(labels)
    for k in range(V.shape[1]):
        v = V[:, k]
        s = float(np.dot(v - v.mean(), ap - ap.mean()))
        if s < 0:
            V[:, k] = -v
    return V


def sc_pcs(R: pd.DataFrame, n_components: int = 3,
           orient: bool = True) -> PCResult:
    """Leading eigenvectors of the structural covariance matrix.

    Returned loadings are z-scored across regions, matching the thesis's
    ``x / sd(x)`` normalisation, so they are comparable in scale to the AHBA
    component scores.
    """
    A = R.to_numpy()
    L, V = np.linalg.eigh(A)
    order = np.argsort(L)[::-1]
    L, V = L[order], V[:, order]
    Vk = V[:, :n_components].copy()
    if orient:
        Vk = _orient(Vk, R.index)
    Z = (Vk - Vk.mean(0)) / Vk.std(0)
    load = pd.DataFrame(Z, index=R.index,
                        columns=[f"SC{i+1}" for i in range(n_components)])
    # eigenvalues of a correlation matrix sum to its trace = n regions
    return PCResult(
        loadings=load,
        var_explained=(L / L.sum())[:n_components],
        var_explained_sq_convention=(L**2 / (L**2).sum())[:n_components],
    )


def slope_pcs(W: pd.DataFrame, n_components: int = 3,
              orient: bool = True) -> PCResult:
    """Region loadings of a PCA on the standardised subject-by-region slopes.

    Mathematically the same vectors as :func:`sc_pcs` (see module docstring);
    computed independently so the report can demonstrate that rather than
    assert it.
    """
    Wz = (W - W.mean()) / W.std()
    A = Wz.to_numpy()
    A = A - A.mean(0)
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    V = Vt.T[:, :n_components].copy()
    if orient:
        V = _orient(V, W.columns)
    Z = (V - V.mean(0)) / V.std(0)
    load = pd.DataFrame(Z, index=W.columns,
                        columns=[f"PC{i+1}" for i in range(n_components)])
    ev = S**2
    return PCResult(loadings=load, var_explained=(ev / ev.sum())[:n_components])


def subject_scores(W: pd.DataFrame, loadings: pd.DataFrame) -> pd.DataFrame:
    """Project subjects onto region loadings -> candidate GWAS phenotypes.

    Standardises each region first, so a component is not dominated by
    whichever region happens to have the largest slope variance.
    """
    Wz = (W - W.mean()) / W.std()
    common = [c for c in loadings.index if c in Wz.columns]
    S = Wz[common].to_numpy() @ loadings.loc[common].to_numpy()
    return pd.DataFrame(S, index=Wz.index, columns=loadings.columns)
