"""Shared inputs and decompositions for the c3axis exploratory scripts (README_HPC.md section 8).

Everything here is map-level or reads per-child data that stays in gitignored locations
(out/, abcd-data-release-7.0/). Scripts write summary tables only (rule 16).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path[:0] = [str(REPO / "src"), str(REPO / "ahba_pls" / "code")]
from abcd import covariance as cov  # noqa: E402
import pls  # noqa: E402  (spin_corr, cached_spin_permutations, HCP_CENTROIDS)

RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
AP = REPO / "ahba_pls"
REFS_L = ["CT_lh", "dCT_lh", "PLS1", "PLS2", "C1", "C2", "C3"]


def sp(a, b) -> float:
    return stats.spearmanr(a, b, nan_policy="omit").statistic


def zs(W: pd.DataFrame) -> pd.DataFrame:
    return (W - W.mean()) / W.std()


def load_inputs():
    """Slope BLUP matrices (all 358 parcels; LH 179 relabelled to the map labels) and LH reference maps."""
    bl = pd.read_parquet(RUN / "fits/blups.parquet")
    Wh = bl.pivot(index="subject", columns="label", values="re_slope")
    assert Wh.notna().all().all()
    maps = pd.read_csv(AP / "results/hcp_summary_maps.csv", index_col=0)[["CT", "dCT", "PLS1", "PLS2", "C1", "C3"]]
    ylab = {l.lower(): l for l in maps.index}
    ul = pd.read_csv(AP / "data/reference/hcp_cortices/HCP-MMP1_UniqueRegionList.csv", encoding="utf-8-sig")
    id2 = {int(r.regionID): ylab.get(f"lh_{r.region}".lower()) for r in ul[ul.LR == "L"].itertuples()}
    c = pd.read_csv(REPO / "data/ahba_dme_hcp_top8kgenes_scores.csv")
    c = c[c.id.isin(id2)].assign(label=lambda d: d.id.map(id2)).dropna(subset=["label"]).set_index("label")
    maps["C2"] = c.C2.reindex(maps.index)
    fx = pd.read_parquet(RUN / "fits/fixed.parquet").pivot(index="label", columns="term", values="estimate")
    lhfx = fx[fx.index.str.startswith("lh_")].copy()
    lhfx.index = [ylab.get(i.lower(), i) for i in lhfx.index]
    REFL = maps[["PLS1", "PLS2", "C1", "C2", "C3"]].join(
        pd.DataFrame({"CT_lh": lhfx["(Intercept)"], "dCT_lh": lhfx["age_c"]}))
    WL = Wh[[c for c in Wh.columns if c.startswith("lh_")]].copy()
    WL.columns = [ylab.get(c.lower(), c) for c in WL.columns]
    return Wh, WL, REFL, maps, ylab


def child_covariates(index) -> pd.DataFrame:
    mt = pd.read_parquet(RUN / "model_table.parquet",
                         columns=["subject", "sex", "site", "family_id", "age_first", "age_span", "n_visits"])
    cv = mt.drop_duplicates("subject").set_index("subject").loc[index]
    cv["male"] = (cv.sex == "M").astype(float)
    return cv


# ---------------------------------------------------------------- decompositions
def regress_global(W: pd.DataFrame) -> pd.DataFrame:
    """Each parcel regressed on the child's mean standardised slope, with its own beta."""
    Z = zs(W); g = Z.mean(1).to_numpy()
    B = Z.T.to_numpy() @ g / (g @ g)
    return pd.DataFrame(Z.to_numpy() - np.outer(g, B), index=W.index, columns=W.columns)


def row_center(W: pd.DataFrame) -> pd.DataFrame:
    """Each child's cortex-wide mean of the standardised slopes subtracted uniformly."""
    Z = zs(W); return Z.sub(Z.mean(1), axis=0)


def _svd(Z: pd.DataFrame, k: int):
    A = Z.to_numpy() - Z.to_numpy().mean(0)
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    return pd.DataFrame(Vt.T[:, :k], index=Z.columns, columns=[f"c{i+1}" for i in range(k)]), (S ** 2 / (S ** 2).sum())[:k]


def pca_maps(W: pd.DataFrame, k: int = 10):          # standardise, then PCA
    return _svd(zs(W), k)


def pca_raw(Z: pd.DataFrame, k: int = 10):           # PCA without re-standardising
    return _svd(Z, k)


def dme_maps(W: pd.DataFrame, k: int = 10, alpha: float = 1.0, sparsity: float = 0):
    """Diffusion map embedding with the Dear et al. 2024 settings (negative affinities kept)."""
    import brainspace.gradient.gradient as bg
    import brainspace.gradient.kernels as bk
    bg.compute_affinity = lambda x, kernel=None, sparsity=.9, pre_sparsify=True, non_negative=False, gamma=None: \
        bk.compute_affinity(x, kernel=kernel, sparsity=sparsity, pre_sparsify=pre_sparsify, non_negative=False, gamma=gamma)
    gm = bg.GradientMaps(n_components=k, approach="dm", kernel="normalized_angle", random_state=0)
    gm.fit(zs(W).T.to_numpy(), sparsity=sparsity, alpha=alpha)
    return (pd.DataFrame(gm.gradients_, index=W.columns, columns=[f"c{i+1}" for i in range(k)]),
            np.asarray(gm.lambdas_))


def orient(G: pd.DataFrame) -> pd.DataFrame:
    """Flip each component to correlate positively with the anterior-posterior axis."""
    return pd.DataFrame(cov._orient(G.to_numpy().copy(), G.index), index=G.index, columns=G.columns)


def bilat(df: pd.DataFrame, ylab: dict) -> pd.DataFrame:
    d = df.copy(); d["key"] = d.index.str.replace(r"^(lh|rh)_", "", regex=True)
    out = d.groupby("key").mean(); out.index = [ylab.get(("lh_" + i).lower(), "lh_" + i) for i in out.index]
    return out
