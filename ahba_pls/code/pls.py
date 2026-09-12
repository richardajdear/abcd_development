"""
PLS-SVD (PLS correlation / PLSC) between regional gene expression and regional
imaging maps, with spin-permutation inference on the components and
region-bootstrap Z-scores on the gene weights.

Conventions
-----------
X : (n_regions, n_genes)  expression, z-scored per gene across regions.
Y : (n_regions, k)        imaging maps, z-scored per map across regions.
R = X'Y / (n-1)           cross-covariance (= correlations, since both z-scored).
R = U S V'                U: gene weights (genes x k), V: imaging saliences (k x k).

Component j:
  gene-side region scores   Lx_j = X @ U[:, j]     ("spatial scores" of the genes)
  imaging-side region scores Ly_j = Y @ V[:, j]
  covariance explained       S_j^2 / sum(S^2)

With a single Y column, U[:, 0] is proportional to the vector of gene-map
Pearson correlations, so option 1 of the design is the exact special case.

Inference
---------
* Spin permutation: Y (complete on all 34 LH DK parcels) is rotated on the
  sphere with a bijective parcel assignment (as in abcd.spatial.spin_indices),
  then restricted to the 33 AHBA-covered regions and re-fitted. The p-value
  for component j compares its singular value to the null distribution of the
  j-th singular value (ordered-component null; conservative for later
  components, standard in PLSC).
* Bootstrap: regions resampled with replacement; each refit is
  Procrustes-aligned (orthogonal rotation + sign) to the observed U so that
  weights are comparable across resamples; bootstrap Z = U / SD_boot(U).
  This is the "bootstrapped weight" convention of Whitaker & Vertes 2016.

Sign convention: each component is oriented so that the salience of the FIRST
Y column (dCT in every option) is positive -- i.e. positive gene weights mean
"higher expression where thinning rate is more positive (less thinning)".
Callers should flip as needed for presentation and say so.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

REPO = Path(__file__).resolve().parents[2]
DK_CENTROIDS = REPO / "data" / "dk_centroids.csv"


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def zscore_cols(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, float)
    mu = a.mean(0, keepdims=True)
    sd = a.std(0, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    return (a - mu) / sd


def _procrustes_rotation(U_ref: np.ndarray, U_new: np.ndarray) -> np.ndarray:
    """Orthogonal Q minimising ||U_new Q - U_ref||_F (includes sign flips)."""
    A = U_new.T @ U_ref
    P, _, Qt = np.linalg.svd(A)
    return P @ Qt


# --------------------------------------------------------------------------
# core fit
# --------------------------------------------------------------------------

@dataclass
class PLSFit:
    U: np.ndarray            # genes x k  (gene weights)
    V: np.ndarray            # k x k      (imaging saliences)
    S: np.ndarray            # k          (singular values)
    Lx: np.ndarray           # regions x k (gene-side scores)
    Ly: np.ndarray           # regions x k (imaging-side scores)
    genes: pd.Index
    ycols: list[str]
    regions: pd.Index
    cov_explained: np.ndarray = field(init=False)
    lx_ly_corr: np.ndarray = field(init=False)

    def __post_init__(self):
        self.cov_explained = self.S ** 2 / (self.S ** 2).sum()
        self.lx_ly_corr = np.array([
            np.corrcoef(self.Lx[:, j], self.Ly[:, j])[0, 1] for j in range(len(self.S))
        ])

    def weights(self) -> pd.DataFrame:
        return pd.DataFrame(self.U, index=self.genes,
                            columns=[f"PLS{j+1}" for j in range(self.U.shape[1])])

    def scores(self) -> pd.DataFrame:
        k = self.U.shape[1]
        out = {}
        for j in range(k):
            out[f"PLS{j+1}_gene_scores"] = self.Lx[:, j]
            out[f"PLS{j+1}_imaging_scores"] = self.Ly[:, j]
        return pd.DataFrame(out, index=self.regions)

    def saliences(self) -> pd.DataFrame:
        return pd.DataFrame(self.V, index=self.ycols,
                            columns=[f"PLS{j+1}" for j in range(self.V.shape[1])])


def _svd_fit(Xz: np.ndarray, Yz: np.ndarray):
    n = Xz.shape[0]
    R = Xz.T @ Yz / (n - 1)
    U, S, Vt = np.linalg.svd(R, full_matrices=False)
    V = Vt.T
    # orient: salience of first Y column positive
    sgn = np.sign(V[0, :])
    sgn[sgn == 0] = 1
    return U * sgn, S, V * sgn


def pls_svd(X: pd.DataFrame, Y: pd.DataFrame) -> PLSFit:
    """PLS-SVD of z-scored X (regions x genes) against z-scored Y (regions x k)."""
    if not X.index.equals(Y.index):
        raise ValueError("X and Y must share an identical region index (same order)")
    Xz = zscore_cols(X.to_numpy(float))
    Yz = zscore_cols(Y.to_numpy(float))
    U, S, V = _svd_fit(Xz, Yz)
    return PLSFit(U=U, V=V, S=S, Lx=Xz @ U, Ly=Yz @ V,
                  genes=X.columns, ycols=list(Y.columns), regions=X.index)


# --------------------------------------------------------------------------
# spin permutation (left hemisphere, bijective assignment)
# --------------------------------------------------------------------------

def load_lh_centroids(labels: pd.Index | list[str] | None = None,
                      path: Path = DK_CENTROIDS) -> pd.DataFrame:
    """Spherical centroids for LH DK parcels, indexed by ggseg label (lh_*)."""
    df = pd.read_csv(path)
    df = df[df.hemi == "lh"].set_index("label")[["x", "y", "z"]]
    if labels is not None:
        missing = [l for l in labels if l not in df.index]
        if missing:
            raise ValueError(f"labels absent from centroid file: {missing[:5]}")
        df = df.loc[list(labels)]
    return df


def spin_permutations(coords: np.ndarray, n_perm: int = 1000, seed: int = 0) -> np.ndarray:
    """(n_perm, n) index arrays: perm[i] gives the source parcel for each target parcel.

    Random 3D rotation of spherical centroids, then a one-to-one parcel
    assignment by minimum total displacement (Hungarian algorithm), so every
    null map is a true permutation of the original values.
    """
    rng = np.random.default_rng(seed)
    n = coords.shape[0]
    out = np.empty((n_perm, n), dtype=int)
    for i in range(n_perm):
        q, r = np.linalg.qr(rng.normal(size=(3, 3)))
        rot = q * np.sign(np.diag(r))
        rotated = coords @ rot.T
        d = ((rotated[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
        _, col = linear_sum_assignment(d)
        out[i] = col
    return out


def spin_test_pls(X: pd.DataFrame, Y_full: pd.DataFrame, n_perm: int = 5000,
                  seed: int = 0) -> dict:
    """Spin-permutation p-values for PLS singular values.

    X is on the covered regions (33); Y_full is on ALL LH parcels (34) so
    that the rotation is performed on the complete map and the uncovered
    parcel is dropped afterwards.
    """
    cov = load_lh_centroids(Y_full.index)
    perms = spin_permutations(cov.to_numpy(float), n_perm=n_perm, seed=seed)
    pos = {l: i for i, l in enumerate(Y_full.index)}
    keep = np.array([pos[l] for l in X.index])
    Yfull = Y_full.to_numpy(float)
    obs = pls_svd(X, Y_full.loc[X.index])
    Xz = zscore_cols(X.to_numpy(float))
    k = len(obs.S)
    null_S = np.empty((n_perm, k))
    null_cov = np.empty((n_perm, k))
    for i, p in enumerate(perms):
        Yp = zscore_cols(Yfull[p][keep])
        _, S, _ = _svd_fit(Xz, Yp)
        null_S[i] = S
        null_cov[i] = S ** 2 / (S ** 2).sum()
    p_S = (1 + (null_S >= obs.S[None, :]).sum(0)) / (1 + n_perm)
    p_cov = (1 + (null_cov >= obs.cov_explained[None, :]).sum(0)) / (1 + n_perm)
    return {"fit": obs, "p_singular": p_S, "p_cov_explained": p_cov,
            "null_S": null_S, "null_cov": null_cov, "perms": perms}


_PERM_CACHE: dict = {}


def cached_spin_permutations(labels, n_perm: int = 5000, seed: int = 0) -> np.ndarray:
    """spin_permutations for a label set, memoised on (labels, n_perm, seed)."""
    key = (tuple(labels), n_perm, seed)
    if key not in _PERM_CACHE:
        cov = load_lh_centroids(list(labels))
        _PERM_CACHE[key] = spin_permutations(cov.to_numpy(float), n_perm=n_perm, seed=seed)
    return _PERM_CACHE[key]


def spin_corr(a: pd.Series, b_full: pd.Series, n_perm: int = 5000, seed: int = 0,
              method: str = "spearman") -> tuple[float, float, np.ndarray]:
    """Spin-test the correlation of map `a` (covered regions) with `b_full`.

    `b_full` should be the more complete map (it is the one rotated); it is
    aligned to `a`'s regions after rotation.  If `b_full` covers only `a`'s
    regions the rotation is done on that subset (a close approximation for
    33/34 parcels).
    """
    from scipy import stats
    b_full = b_full.dropna()
    common = a.index.intersection(b_full.index)
    a = a.loc[common]
    perms = cached_spin_permutations(b_full.index, n_perm=n_perm, seed=seed)
    pos = {l: i for i, l in enumerate(b_full.index)}
    keep = np.array([pos[l] for l in common])
    bv = b_full.to_numpy(float)
    f = (lambda x, y: stats.spearmanr(x, y).statistic) if method == "spearman" \
        else (lambda x, y: stats.pearsonr(x, y).statistic)
    r_obs = float(f(a.to_numpy(float), bv[keep]))
    null = np.array([f(a.to_numpy(float), bv[p][keep]) for p in perms])
    p = (1 + (np.abs(null) >= abs(r_obs)).sum()) / (1 + n_perm)
    return r_obs, float(p), null


# --------------------------------------------------------------------------
# bootstrap gene weights
# --------------------------------------------------------------------------

def bootstrap_weights(X: pd.DataFrame, Y: pd.DataFrame, n_boot: int = 1000,
                      seed: int = 0, align: str = "procrustes") -> dict:
    """Region-bootstrap Z-scores for gene weights and imaging saliences.

    Each resample refits the PLS-SVD, is Procrustes-aligned to the observed
    solution (rotation within the k-dimensional subspace, which for k=1 is a
    sign flip), and the SD across resamples gives Z = U_obs / SD(U_boot).
    Also returns the bootstrap mean weights and the per-component
    reproducibility (mean correlation of resampled with observed weights).
    """
    rng = np.random.default_rng(seed)
    obs = pls_svd(X, Y)
    Xa, Ya = X.to_numpy(float), Y.to_numpy(float)
    n, k = Xa.shape[0], len(obs.S)
    Us = np.empty((n_boot, *obs.U.shape))
    Vs = np.empty((n_boot, *obs.V.shape))
    corr_rep = np.empty((n_boot, k))
    done = 0
    while done < n_boot:
        idx = rng.integers(0, n, n)
        if len(np.unique(idx)) < max(5, k + 2):
            continue
        Xz, Yz = zscore_cols(Xa[idx]), zscore_cols(Ya[idx])
        if not np.all(np.isfinite(Xz)) or not np.all(np.isfinite(Yz)):
            continue
        U, S, V = _svd_fit(Xz, Yz)
        if align == "procrustes":
            Q = _procrustes_rotation(obs.U, U)
            U, V = U @ Q, V @ Q
        else:  # sign only
            sgn = np.sign((U * obs.U).sum(0)); sgn[sgn == 0] = 1
            U, V = U * sgn, V * sgn
        Us[done], Vs[done] = U, V
        corr_rep[done] = [np.corrcoef(U[:, j], obs.U[:, j])[0, 1] for j in range(k)]
        done += 1
    sd = Us.std(0, ddof=1)
    sd[sd == 0] = np.nan
    Z = obs.U / sd
    cols = [f"PLS{j+1}" for j in range(k)]
    return {
        "fit": obs,
        "Z": pd.DataFrame(Z, index=X.columns, columns=cols),
        "U_boot_mean": pd.DataFrame(Us.mean(0), index=X.columns, columns=cols),
        "U_boot_sd": pd.DataFrame(sd, index=X.columns, columns=cols),
        "V_boot_Z": pd.DataFrame(obs.V / Vs.std(0, ddof=1), index=Y.columns, columns=cols),
        "reproducibility": pd.Series(corr_rep.mean(0), index=cols),
    }


# --------------------------------------------------------------------------
# convenience: project reference gene weights to region scores
# --------------------------------------------------------------------------

def project_weights(X: pd.DataFrame, w: pd.Series) -> tuple[pd.Series, int]:
    """Region scores = z(X)[:, shared] @ w[shared]/||w||; returns (scores, n_shared)."""
    shared = X.columns.intersection(w.dropna().index)
    Xz = zscore_cols(X[shared].to_numpy(float))
    wv = w.loc[shared].to_numpy(float)
    wv = wv / np.linalg.norm(wv)
    return pd.Series(Xz @ wv, index=X.index), len(shared)
