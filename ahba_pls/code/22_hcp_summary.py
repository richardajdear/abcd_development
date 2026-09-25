"""
22_hcp_summary.py -- inputs for the HCP-MMP summary slide (fig5_hcp_summary.R):
the option-2 PLS (Y = dCT + CT) against AHBA expression, both components, vs
AHBA C1 and C3.  Run once per AHBA matrix:

    python code/22_hcp_summary.py 3d_ds5   # hcp_3d_ds5.csv: >=3-donor region filter + DS5 gene filter
                                           #   (137 parcels x 7,973 genes; the matrix C1-C3 were fit on)
    python code/22_hcp_summary.py ds5      # hcp_ds5.csv: no region filter, DS5 gene filter
                                           #   (177 parcels x 7,973 genes, its own DS5 set)
    python code/22_hcp_summary.py 3d       # hcp_3d.csv: >=3-donor region filter, no gene filter
    python code/22_hcp_summary.py base     # hcp_base.csv: neither filter (177 parcels x 15,637 genes)

Components are fitted from scratch for every matrix (spin test of the singular
values, 1,000-sample region bootstrap for gene Z).  For 3d_ds5 the fit is
asserted to reproduce 12_hcp_pls.py's saved weights, Z and PLS2 scores.

Orientation: PLS1 (CT salience ~ +0.99) is aligned with AHBA C1, positive =
thinner baseline cortex; PLS2 (dCT salience ~ +0.99) with faster thinning.
AHBA C1 and C3 are always the published maps and weights (fitted on 3d_ds5).

Writes results/<P>_{maps.csv, map_pairs.tsv, gene_pairs.tsv, gene_weights.tsv,
sets.tsv, magma.tsv, components.tsv}, with P = hcp_summary (3d_ds5) or
hcp_summary_<variant>.
  map_pairs   Spearman rho + spin p (5,000 rotations, HCP centroids)
  gene_pairs  Spearman rho over the genes shared with the C1-C3 set
  sets        cell class (Seidlitz 2020) and layer (Maynard 2021; FDR < 0.05, t > 0)
              enrichment on the SHARED universe (PLS genes x C1-C3 genes)
  magma       gene-property, two-sided, each vector on its own universe, against
              SCZ25_META, MDD_div, ASD and the four global-slope gene analyses
              (whichever are in genetic_analysis/inputs/magma/)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls, magma_utils                                 # noqa: E402

VARIANTS = {"3d_ds5": "hcp_3d_ds5.csv", "ds5": "hcp_ds5.csv", "3d": "hcp_3d.csv", "base": "hcp_base.csv"}
VARIANT = sys.argv[1] if len(sys.argv) > 1 else "3d_ds5"
assert VARIANT in VARIANTS, f"variant must be one of {list(VARIANTS)}"
ROOT = HERE.parent; REPO = ROOT.parent
RES = ROOT / "results"
P = "hcp_summary" if VARIANT == "3d_ds5" else f"hcp_summary_{VARIANT}"
AHBA = Path.home() / "Git" / "AHBA" / "data"
EXPR = AHBA / "abagen-data" / "expression"
N_SPIN, N_BOOT, N_PERM, SEED = 5000, 1000, 20000, 0
out = lambda name: RES / f"{P}_{name}"

# ---- X: parcel ids 1-180 are the left-hemisphere regionIDs of HCP-MMP1 --------
Y180 = pd.read_csv(RES / "hcp_y_maps_180.csv", index_col=0)
ulist = pd.read_csv(ROOT / "data" / "reference" / "hcp_cortices" / "HCP-MMP1_UniqueRegionList.csv",
                    encoding="utf-8-sig")
ylab = {l.lower(): l for l in Y180.index}                # case-insensitive: atlas 7Pl vs ABCD 7PL
id2label = {int(r.regionID): ylab.get(f"lh_{r.region}".lower())
            for r in ulist[ulist.LR == "L"].itertuples()}
X = pd.read_csv(EXPR / VARIANTS[VARIANT], index_col=0)
X = X.loc[[i for i in X.index if i in id2label]]         # left hemisphere only (hcp_ds5 has 360 rows)
X = X.dropna(how="all")
n_partial = int(X.isna().any(axis=0).sum())
X = X.loc[:, X.notna().all(axis=0)]                      # genes missing in any covered parcel dropped
X.index = [id2label[i] for i in X.index]
X = X.loc[[l for l in X.index if l is not None and l in Y180.index]]
print(f"[{VARIANT}] X: {X.shape[0]} parcels x {X.shape[1]} genes "
      f"({n_partial} genes dropped for missing values)", file=sys.stderr)
Y = Y180.loc[X.index, ["dCT", "CT"]]

# ---- PLS ---------------------------------------------------------------------
sp = pls.spin_test_pls(X, Y180[["dCT", "CT"]], n_perm=N_SPIN, seed=SEED, centroids=pls.HCP_CENTROIDS)
bt = pls.bootstrap_weights(X, Y, n_boot=N_BOOT, seed=SEED)
fit = sp["fit"]; sal = fit.saliences()
assert sal.loc["CT", "PLS1"] > 0.9 and sal.loc["dCT", "PLS2"] > 0.9, sal
sc = fit.scores()
S1, S2 = -sc["PLS1_gene_scores"], -sc["PLS2_gene_scores"]
Z1, Z2 = -bt["Z"]["PLS1"], -bt["Z"]["PLS2"]

c123w = pd.read_csv(REPO / "data" / "weights.csv", index_col=0)
c123s = pd.read_csv(REPO / "data" / "ahba_dme_hcp_top8kgenes_scores.csv")
C = c123s.assign(label=[id2label[i] for i in c123s.id]).set_index("label")[["C1", "C3"]]
# PLS1 must align with C1 (the orientation the slide is built on)
assert stats.spearmanr(S1, C.C1.reindex(S1.index), nan_policy="omit").statistic > 0.5

if VARIANT == "3d_ds5":                                  # reproduce 12_hcp_pls.py exactly
    saved = pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)
    for k in ("PLS1", "PLS2"):
        assert np.abs(fit.weights()[k] - saved[f"hcp_opt2_dCT_CT_{k}_U"].loc[X.columns]).max() < 1e-4, k
        assert np.abs(bt["Z"][k] - saved[f"hcp_opt2_dCT_CT_{k}_Z"].loc[X.columns]).max() < 1e-3, k   # saved at 6 sig. figs
    old = pd.read_csv(RES / "hcp_pls_scores.csv", index_col=0)["thinning_score"]
    assert np.corrcoef(S2.loc[old.index], old)[0, 1] > 0.9999

pd.DataFrame([dict(variant=VARIANT, matrix=VARIANTS[VARIANT], component=k, n_parcels=X.shape[0],
                   n_genes=X.shape[1], cov_explained=fit.cov_explained[j],
                   p_spin_singular=sp["p_singular"][j], boot_reproducibility=bt["reproducibility"][k],
                   sal_dCT=sal.loc["dCT", k], sal_CT=sal.loc["CT", k])
              for j, k in enumerate(("PLS1", "PLS2"))]
             ).to_csv(out("components.tsv"), sep="\t", index=False, float_format="%.5g")

# ---- maps, systems -----------------------------------------------------------
maps = Y180[["CT", "dCT"]].join(pd.DataFrame({"PLS1": S1, "PLS2": S2})).join(C)
syst = pd.read_csv(ROOT / "data" / "reference" / "hcp_cortices" / "hcp_parcel_systems.csv")
syst = syst.assign(key=syst.label.str.lower()).set_index("key")
maps["system"] = syst.system.reindex(maps.index.str.lower()).to_numpy()
maps["cortex"] = syst.cortex.reindex(maps.index.str.lower()).to_numpy()
assert maps.system.notna().all(), maps.index[maps.system.isna()].tolist()
maps.rename_axis("label").to_csv(out("maps.csv"), float_format="%.6g")

rows = []
for comp, s in (("PLS1", S1), ("PLS2", S2)):
    for ref in ("CT", "dCT", "C1", "C3"):
        rho, p, _ = pls.spin_corr(s, maps[ref], n_perm=N_SPIN, seed=SEED, centroids=pls.HCP_CENTROIDS)
        rows.append(dict(x=ref, y=comp, rho=rho, p_spin=p,
                         n=int(s.index.intersection(maps[ref].dropna().index).size)))
rho, p, _ = pls.spin_corr(maps["CT"], maps["dCT"], n_perm=N_SPIN, seed=SEED, centroids=pls.HCP_CENTROIDS)
rows.append(dict(x="CT", y="dCT", rho=rho, p_spin=p, n=int(maps[["CT", "dCT"]].dropna().shape[0])))
MP = pd.DataFrame(rows)
MP.to_csv(out("map_pairs.tsv"), sep="\t", index=False, float_format="%.4g")

# ---- gene vectors ------------------------------------------------------------
Z = pd.DataFrame({"PLS1": Z1, "PLS2": Z2, "C1": c123w["C1"], "C3": c123w["C3"]}).dropna()
Z.rename_axis("gene").to_csv(out("gene_weights.tsv"), sep="\t", float_format="%.6g")
GP = pd.DataFrame([dict(x=r, y=c, rho=stats.spearmanr(Z[c], Z[r]).statistic, n=len(Z))
                   for c in ("PLS1", "PLS2") for r in ("C1", "C3")])
GP.to_csv(out("gene_pairs.tsv"), sep="\t", index=False, float_format="%.4g")

# ---- cell classes and layers (shared universe) --------------------------------
cell = pd.read_csv(AHBA / "seidlitz_cell_genes.csv")
gcols = [c for c in cell.columns if c not in ("Type", "Paper", "Cluster", "Class")]
SETS = {("cell class", k): {x for x in pd.unique(g[gcols].to_numpy().ravel())
                            if isinstance(x, str) and x.strip()}
        for k, g in cell.groupby("Class")}
may = pd.read_csv(AHBA / "maynard_layers.csv", encoding="utf-8-sig")
for lay, col in [("L1", "Layer1"), ("L2", "Layer2"), ("L3", "Layer3"), ("L4", "Layer4"),
                 ("L5", "Layer5"), ("L6", "Layer6"), ("WM", "WM")]:
    SETS[("layer", lay)] = set(may.loc[(may[f"fdr_{col}"] < 0.05) & (may[f"t_stat_{col}"] > 0), "gene"])
res = []
for (kind, name), gset in SETS.items():
    idx = np.flatnonzero(Z.index.isin(gset))
    if len(idx) < 10:
        continue
    rng = np.random.default_rng(SEED)                  # same draws for every vector
    draws = np.stack([rng.choice(len(Z), len(idx), replace=False) for _ in range(N_PERM)])
    for v in Z.columns:
        z = ((Z[v] - Z[v].mean()) / Z[v].std()).to_numpy()
        obs, null = z[idx].mean(), z[draws].mean(1)
        res.append(dict(kind=kind, set=name, vector=v, n_genes=len(idx), mean_z=obs,
                        z=(obs - null.mean()) / null.std(),
                        p_perm=(1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (1 + N_PERM)))
ST = pd.DataFrame(res); ST["q_bh"] = np.nan
for v, g in ST.groupby("vector"):                      # BH within vector, across all sets
    pv = g.p_perm.to_numpy(); o = np.argsort(pv); m = len(pv)
    q = np.minimum.accumulate((pv[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    ST.loc[g.index[o], "q_bh"] = np.minimum(q, 1)
ST.to_csv(out("sets.tsv"), sep="\t", index=False, float_format="%.4g")

# ---- MAGMA -------------------------------------------------------------------
GA = ["SCZ25_META", "MDD", "ASD", "BIP", "ADHD", "ALZ", "EA", "INT", "HEIGHT", "global_slope_1lmm_hcp", "global_slope_1lmm_dk",
      "global_slope_perregion_hcp", "global_slope_perregion_dk"]
NAME = {"MDD": "MDD_div"}                               # MDD.genes.raw = MDD2025 div (multi-ancestry)
VEC = {"ABCD_PLS1_HCP": Z1, "ABCD_PLS2_HCP": Z2, "AHBA_C1": c123w["C1"], "AHBA_C3": c123w["C3"]}
run_dir = RES / "magma_runs" / P; run_dir.mkdir(parents=True, exist_ok=True)
s2e = magma_utils.sym2entrez()
mrows, missing = [], []
for ga in GA:
    raw = magma_utils.INP / f"{ga}.genes.raw"
    if not raw.exists():
        missing.append(ga); continue
    for v, w in VEC.items():
        t = magma_utils.run(raw, magma_utils.write_covar({v: w}, run_dir / f"{v}.covar", s2e), run_dir / f"{ga}__{v}")
        mrows.append(dict(gene_analysis=NAME.get(ga, ga), **t.iloc[0][["VARIABLE", "NGENES", "BETA", "BETA_STD",
                                                                        "SE", "se_std", "P"]].to_dict()))
MG = pd.DataFrame(mrows)
MG.to_csv(out("magma.tsv"), sep="\t", index=False, float_format="%.5g")

pd.set_option("display.width", 200)
print(pd.read_csv(out("components.tsv"), sep="\t").round(4).to_string(index=False))
print(MP.round(3).to_string(index=False)); print(GP.round(3).to_string(index=False))
MG["cell"] = MG.apply(lambda r: f"{r.BETA_STD:+.3f} ({r.P:.2g})", axis=1)
print(MG.pivot(index="VARIABLE", columns="gene_analysis", values="cell").to_string())
print("MAGMA files missing:", missing)
