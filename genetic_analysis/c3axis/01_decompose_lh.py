"""01 -- LEFT-hemisphere decompositions of between-child slope covariance (README_HPC.md 8.1, 8.6).

Input: standardised child x parcel slope BLUPs of the HCP-MMP run, 179 LH parcels (AHBA
C1-C3 are left-hemisphere maps). Compares plain PCA, PCA after regressing out the global
mean (per-parcel beta), PCA on row-centred slopes (each child's mean subtracted), and DME
at alpha 0 / 0.5 / 1, against CT, dCT, PLS1, PLS2 and AHBA C1-C3.

Writes (map-level / aggregate only):
  lh_method_comparison.tsv        Spearman rho of components 1-10 of every method with each map
  lh_components_vs_refs.tsv       components 1-8 of PCA, row-centred PCA and DME, with spin p (5,000)
  lh_component_maps.csv           selected oriented component maps + reference maps (group level)
  lh_child_global_coupling.tsv    child-score correlation of each component with the child's mean slope
  lh_child_score_correlations.tsv row-centred scores vs map projections and scan-design covariates
  table_component_split_half.tsv  split-half |r| of loadings (100 halves) + % variance, three decompositions
  table_rc45_plane_stability.tsv  stability of the row-centred 4-5 plane and its C3-aligned direction

Run from the repo root:  python genetic_analysis/c3axis/01_decompose_lh.py   (about 4 min)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (HERE, REFS_L, bilat, child_covariates, dme_maps, load_inputs, orient, pca_maps,
                    pca_raw, pls, regress_global, row_center, sp, zs)

K = 10
Wh, WL, REFL, maps, ylab = load_inputs()
WLr, WLc = regress_global(WL), row_center(WL)
RC = "PCA, row-centred (constant removed)"
methods = {
    "PCA": pca_maps(WL),
    "PCA, global regressed": pca_maps(WLr),
    "DME alpha=0": dme_maps(WL, alpha=0),
    "DME alpha=0.5": dme_maps(WL, alpha=0.5),
    "DME alpha=1": dme_maps(WL, alpha=1),
    "DME alpha=1, global regressed": dme_maps(WLr, alpha=1),
    RC: pca_raw(WLc),
    "PCA, row-centred then re-standardised": pca_maps(WLc),
}

# ---- 1. every method x component x map -------------------------------------------------
LM = pd.DataFrame([dict(method=m, k=j + 1, eig=ev[j], **{t: sp(G[c].reindex(REFL.index), REFL[t]) for t in REFS_L})
                   for m, (G, ev) in methods.items() for j, c in enumerate(G)])
LM.to_csv(HERE / "lh_method_comparison.tsv", sep="\t", index=False, float_format="%.4g")

# ---- 2. spin-tested matrix for PCA, row-centred PCA and DME ------------------------------
Gp = methods["PCA"][0]
pc1 = Gp.c1 * np.sign(Gp.c1.mean())                               # raw PC1, all loadings positive
rows = []
for m in ["PCA", RC, "DME alpha=1"]:
    G = orient(methods[m][0]); G = (G - G.mean()) / G.std()
    if m == "PCA":
        G["c1"] = G.c1 * np.sign(stats.pearsonr(G.c1, pc1).statistic)   # high = strong global loading
    for c in list(G.columns)[:8]:
        for t in REFS_L:
            r, p, _ = pls.spin_corr(REFL[t].dropna(), G[c], n_perm=5000, seed=0, centroids=pls.HCP_CENTROIDS)
            rows.append(dict(method=m, comp=c, ref=t, rho=r, p_spin=p, eig=methods[m][1][int(c[1:]) - 1]))
pd.DataFrame(rows).to_csv(HERE / "lh_components_vs_refs.tsv", sep="\t", index=False, float_format="%.4g")

# ---- 3. group-level component maps for the figure ----------------------------------------
Gso = {m: orient(methods[m][0]) for m in ("PCA", "DME alpha=1")}
F = pd.DataFrame({"PCA_c1": Gso["PCA"].c1 * np.sign(stats.pearsonr(Gso["PCA"].c1, pc1).statistic),
                  "PCA_c2": Gso["PCA"].c2, "PCA_c3": Gso["PCA"].c3,
                  **{f"DME_c{i}": Gso["DME alpha=1"][f"c{i}"] for i in (1, 2, 5, 7)}})
((F - F.mean()) / F.std()).join(REFL).rename_axis("label").to_csv(HERE / "lh_component_maps.csv", float_format="%.6g")

# ---- 4. child-level coupling with the global slope ---------------------------------------
Gc = methods[RC][0]; Gd = methods["DME alpha=1"][0]
ZL = zs(WL); gmean = ZL.mean(1)
sc_rc = WLc.to_numpy() @ Gc.to_numpy()
sc_p = (ZL.to_numpy() - ZL.to_numpy().mean(0)) @ Gp.to_numpy()
pd.DataFrame({"comp": [f"c{i+1}" for i in range(K)],
              "r_global_rowcentred": [stats.pearsonr(gmean, sc_rc[:, j]).statistic for j in range(K)],
              "r_global_PCA": [stats.pearsonr(gmean, sc_p[:, j]).statistic for j in range(K)],
              "r_DME_vs_rowcentredPCA": [abs(stats.pearsonr(Gd[c], Gc[c]).statistic) for c in Gd]}
             ).to_csv(HERE / "lh_child_global_coupling.tsv", sep="\t", index=False, float_format="%.4g")

# ---- 5. row-centred child scores vs map projections and scan design ----------------------
S = pd.DataFrame(sc_rc[:, :8], index=WL.index, columns=[f"c{i+1}" for i in range(8)])
S["global_mean"] = gmean


def proj(m: pd.Series) -> pd.Series:      # slope of the child's row-centred raw BLUP map on the map
    m = m.dropna(); X = WL[m.index].to_numpy(); v = (m - m.mean()).to_numpy()
    return pd.Series((X - X.mean(1, keepdims=True)) @ v / (v @ v), index=WL.index)


S["proj_dCT"], S["proj_C3"], S["proj_PLS2"] = proj(REFL.dCT_lh), proj(REFL.C3), proj(REFL.PLS2)
cv = child_covariates(WL.index)
C = pd.concat([S, cv[["age_first", "age_span", "n_visits", "male"]]], axis=1).corr(method="spearman")
C.loc[[f"c{i}" for i in range(1, 9)] + ["global_mean", "proj_dCT", "proj_C3", "proj_PLS2"],
      ["global_mean", "proj_dCT", "proj_C3", "proj_PLS2", "age_first", "age_span", "n_visits", "male"]
      ].rename_axis("row_centred_component_or_score").to_csv(HERE / "lh_child_score_correlations.tsv", sep="\t",
                                                             float_format="%.3f")

# ---- 6. split-half reliability of the component maps -------------------------------------
def _load(Z, k=10):
    A = Z.to_numpy() - Z.to_numpy().mean(0); return np.linalg.svd(A, full_matrices=False)[2][:k].T


def split_half(W, prep, n_rep=100, k=10, Kout=7, seed=0):
    rng = np.random.default_rng(seed); rows = []
    for rep in range(n_rep):
        idx = rng.permutation(len(W)); h = len(W) // 2
        C_ = np.abs(_load(prep(W.iloc[idx[:h]]), k).T @ _load(prep(W.iloc[idx[h:]]), k))
        ri, ci = linear_sum_assignment(-C_); match = dict(zip(ri, ci))
        rows += [dict(comp=j + 1, same=C_[j, j], matched=C_[j, match[j]], stable=match[j] == j) for j in range(Kout)]
    return pd.DataFrame(rows).groupby("comp").agg(
        same_index_median=("same", "median"), same_index_p05=("same", lambda s: s.quantile(.05)),
        matched_median=("matched", "median"), index_stable_frac=("stable", "mean"))


Zh = zs(Wh); Sh = np.linalg.svd(Zh.to_numpy() - Zh.to_numpy().mean(0), compute_uv=False)
decs = {"HCP both hemispheres, plain PCA": (Wh, zs, Sh ** 2 / (Sh ** 2).sum()),
        "LH only, plain PCA": (WL, zs, methods["PCA"][1]),
        "LH only, row-centred PCA (= DME)": (WL, row_center, methods[RC][1])}
pd.concat({k: split_half(W, f).assign(var_explained_pct=np.asarray(ev)[:7] * 100) for k, (W, f, ev) in decs.items()},
          names=["decomposition"]).reset_index().to_csv(HERE / "table_component_split_half.tsv", sep="\t",
                                                         index=False, float_format="%.4g")

# ---- 7. the row-centred 4-5 plane and its C3-aligned direction ---------------------------
c3v = REFL.C3.reindex(WL.columns); cidx = c3v.notna().to_numpy(); c3z = ((c3v - c3v.mean()) / c3v.std()).to_numpy()


def c3_dir(V2):
    B = V2[cidx] - V2[cidx].mean(0); w = np.linalg.lstsq(B, c3z[cidx], rcond=None)[0]
    v = V2 @ w; return v / np.linalg.norm(v), np.corrcoef(B @ w, c3z[cidx])[0, 1]


rng = np.random.default_rng(3); rr = []
for rep in range(30):
    idx = rng.permutation(len(WL)); h = len(WL) // 2
    VA, VB = _load(row_center(WL.iloc[idx[:h]]), 7), _load(row_center(WL.iloc[idx[h:]]), 7)
    cos = np.linalg.svd(VA[:, 3:5].T @ VB[:, 3:5], compute_uv=False)
    (dA, rA), (dB, rB) = c3_dir(VA[:, 3:5]), c3_dir(VB[:, 3:5])
    rr.append(dict(plane_cos1=cos[0], plane_cos2=cos[1], c3dir_halfA_vs_B=abs(dA @ dB), c3_r_A=rA, c3_r_B=rB))
pd.DataFrame(rr).describe().to_csv(HERE / "table_rc45_plane_stability.tsv", sep="\t", float_format="%.4g")
print("wrote 7 tables to", HERE)
