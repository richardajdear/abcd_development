"""
22_hcp_summary.py -- inputs for the HCP-MMP-only summary slide
(fig5_hcp_summary.R): the option-2 (Y = dCT + CT) fit at 137 parcels, both
components, against AHBA C1 and C3.

PLS1 and PLS2 of that fit:
  PLS1  the static component: CT salience +0.99, 86% of the cross-covariance.
        Oriented so positive = higher baseline thickness.
  PLS2  the thinning signature: dCT salience +0.99. Oriented so positive =
        faster thinning (the thinning orientation used everywhere here).

12_hcp_pls.py saved gene weights for both components but region scores only for
PLS2, so the fit is repeated here (same X, same Y, same code path) and the
refit's gene weights are asserted identical to the saved ones before anything
is used. Nothing else is refitted.

Writes, all under results/hcp_summary_*:
  maps.csv          CT, dCT (all 179 parcels), PLS1, PLS2, C1, C3 (137 covered)
  map_pairs.tsv     Spearman rho + spin p (5,000 rotations, HCP centroids)
  gene_pairs.tsv    Spearman rho between gene vectors (7,973 genes)
  gene_weights.tsv  the four gene vectors, one universe
  sets.tsv          cell-class (Seidlitz 2020) and cortical-layer (Maynard 2021;
                    FDR < 0.05, t > 0, as in Dear et al. 2024) enrichment
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls                                             # noqa: E402

ROOT = HERE.parent; REPO = ROOT.parent
RES = ROOT / "results"
AHBA = Path.home() / "Git" / "AHBA" / "data"
EXPR = AHBA / "abagen-data" / "expression"
N_SPIN, N_PERM, SEED = 5000, 20000, 0

# ---- refit (X and Y exactly as 12_hcp_pls.py) --------------------------------
Y180 = pd.read_csv(RES / "hcp_y_maps_180.csv", index_col=0)
c123w = pd.read_csv(REPO / "data" / "weights.csv", index_col=0)
genes = list(c123w.index)
Xh = pd.read_csv(EXPR / "hcp_3d.csv", index_col=0)[genes]
c123s = pd.read_csv(REPO / "data" / "ahba_dme_hcp_top8kgenes_scores.csv")
id2label = dict(zip(c123s.id, "lh_" + c123s.label.astype(str)))
Xh.index = [id2label.get(i, f"id{i}") for i in Xh.index]
Xh = Xh.loc[[l for l in Xh.index if l in Y180.index]]
fit = pls.pls_svd(Xh, Y180.loc[Xh.index, ["dCT", "CT"]])

saved = pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)
W = fit.weights()
for comp in ("PLS1", "PLS2"):
    d = np.abs(W[comp] - saved[f"hcp_opt2_dCT_CT_{comp}_U"].loc[W.index]).max()
    assert d < 1e-4, f"refit {comp} differs from 12_hcp_pls.py by {d}"
sal = fit.saliences()
assert sal.loc["CT", "PLS1"] > 0.9 and sal.loc["dCT", "PLS2"] > 0.9, sal

sc = fit.scores()
S1 = sc["PLS1_gene_scores"]                 # CT salience > 0: positive = thicker
S2 = -sc["PLS2_gene_scores"]                # dCT salience > 0: flip -> faster thinning
old = pd.read_csv(RES / "hcp_pls_scores.csv", index_col=0)["thinning_score"]
assert np.corrcoef(S2.loc[old.index], old)[0, 1] > 0.9999, "PLS2 scores do not match the saved ones"

Z = pd.DataFrame({"PLS1": saved["hcp_opt2_dCT_CT_PLS1_Z"],
                  "PLS2": -saved["hcp_opt2_dCT_CT_PLS2_Z"],
                  "C1": c123w["C1"], "C3": c123w["C3"]}).dropna()
Z.rename_axis("gene").to_csv(RES / "hcp_summary_gene_weights.tsv", sep="\t", float_format="%.6g")

C = c123s.assign(label="lh_" + c123s.label.astype(str)).set_index("label")[["C1", "C3"]]
maps = Y180[["CT", "dCT"]].join(pd.DataFrame({"PLS1": S1, "PLS2": S2})).join(C)
maps.rename_axis("label").to_csv(RES / "hcp_summary_maps.csv", float_format="%.6g")

# ---- map pairs (spin) --------------------------------------------------------
rows = []
for comp, s in (("PLS1", S1), ("PLS2", S2)):
    for ref in ("CT", "dCT", "C1", "C3"):
        rho, p, _ = pls.spin_corr(s, maps[ref], n_perm=N_SPIN, seed=SEED,
                                  centroids=pls.HCP_CENTROIDS)
        rows.append(dict(x=ref, y=comp, rho=rho, p_spin=p,
                         n=int(s.index.intersection(maps[ref].dropna().index).size)))
MP = pd.DataFrame(rows)
MP.to_csv(RES / "hcp_summary_map_pairs.tsv", sep="\t", index=False, float_format="%.4g")

GP = pd.DataFrame([dict(x=r, y=c, rho=stats.spearmanr(Z[c], Z[r]).statistic, n=len(Z))
                   for c in ("PLS1", "PLS2") for r in ("C1", "C3")])
GP.to_csv(RES / "hcp_summary_gene_pairs.tsv", sep="\t", index=False, float_format="%.4g")

# ---- gene-set enrichment: cell classes and layers ----------------------------
cell = pd.read_csv(AHBA / "seidlitz_cell_genes.csv")
gcols = [c for c in cell.columns if c not in ("Type", "Paper", "Cluster", "Class")]
SETS = {("cell class", k): {x for x in pd.unique(g[gcols].to_numpy().ravel())
                            if isinstance(x, str) and x.strip()}
        for k, g in cell.groupby("Class")}
may = pd.read_csv(AHBA / "maynard_layers.csv", encoding="utf-8-sig")
for lay, col in [("L1", "Layer1"), ("L2", "Layer2"), ("L3", "Layer3"), ("L4", "Layer4"),
                 ("L5", "Layer5"), ("L6", "Layer6"), ("WM", "WM")]:
    SETS[("layer", lay)] = set(may.loc[(may[f"fdr_{col}"] < 0.05) & (may[f"t_stat_{col}"] > 0), "gene"])

rng_seed = SEED
out = []
for (kind, name), gset in SETS.items():
    idx = np.flatnonzero(Z.index.isin(gset))
    if len(idx) < 10:
        continue
    rng = np.random.default_rng(rng_seed)              # same draws for every vector
    draws = np.stack([rng.choice(len(Z), len(idx), replace=False) for _ in range(N_PERM)])
    for v in Z.columns:
        z = ((Z[v] - Z[v].mean()) / Z[v].std()).to_numpy()
        obs, null = z[idx].mean(), z[draws].mean(1)
        out.append(dict(kind=kind, set=name, vector=v, n_genes=len(idx), mean_z=obs,
                        z=(obs - null.mean()) / null.std(),
                        p_perm=(1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (1 + N_PERM)))
ST = pd.DataFrame(out)
ST["q_bh"] = np.nan
for v, g in ST.groupby("vector"):                      # BH within vector, across all sets
    p = g.p_perm.to_numpy(); o = np.argsort(p); m = len(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    ST.loc[g.index[o], "q_bh"] = np.minimum(q, 1)
ST.to_csv(RES / "hcp_summary_sets.tsv", sep="\t", index=False, float_format="%.4g")

pd.set_option("display.width", 200)
print(MP.round(3).to_string(index=False)); print(GP.round(3).to_string(index=False))
print(ST.pivot(index=["kind", "set"], columns="vector", values="z").round(1).to_string())
print("\nset sizes:", ST.drop_duplicates("set").set_index("set").n_genes.to_dict())
