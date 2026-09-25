"""
27_tier_profile.py -- what the fast-, middle- and slow-thinning thirds of cortex
ARE, transcriptomically: the mechanism panel of the base-AHBA summary slide.

Question it answers: symptom-linked extra thinning is in normally slow-thinning
cortex (26_gradient_scores.py), while the SCZ/MDD gene-property enrichment is
carried by PLS2 / AHBA C3. Is slow-thinning cortex the glial / white-matter pole
of that axis, and fast-thinning cortex the neuronal / L2-3 pole where the risk
genes are expressed?

Inputs (group-level only):
  results/gradient_tiers.csv           parcel -> fast / mid / slow third (from 26)
  results/<P>_maps.csv                 PLS2, C3 scores for the chosen AHBA matrix
  AHBA abagen hcp_<variant>.csv        expression (same matrix the PLS used)
Usage: python code/27_tier_profile.py [3d_ds5|ds5|3d|base]   (as 22_hcp_summary.py)
  Seidlitz 2020 cell classes, Maynard 2021 layers (as in 22_hcp_summary.py)
  genetic_analysis/inputs/magma/{SCZ25_META,MDD}.genes.out  gene-level ZSTAT
Parcel feature scores: every gene z-scored across parcels, then the mean over a
set's genes (marker sets; the top-N genes by MAGMA ZSTAT for the risk sets).
Each feature map is then z-scored across parcels.
Inference: Spearman rho of each feature map with the normative thinning rate,
p from 5,000 spin rotations (HCP centroids). Tier means are descriptive.
Outputs: results/tier_profile_<variant>.tsv, results/tier_profile_<variant>_parcels.csv
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls                                                    # noqa: E402

ROOT = HERE.parent; REPO = ROOT.parent; RES = ROOT / "results"
AHBA = Path.home() / "Git" / "AHBA" / "data"; EXPR = AHBA / "abagen-data" / "expression"
MAG = REPO / "genetic_analysis" / "inputs" / "magma"
N_SPIN, SEED, TOPN = 5000, 0, (250, 500, 1000)
VARIANTS = {"3d_ds5": "hcp_3d_ds5.csv", "ds5": "hcp_ds5.csv", "3d": "hcp_3d.csv", "base": "hcp_base.csv"}
VARIANT = sys.argv[1] if len(sys.argv) > 1 else "base"
assert VARIANT in VARIANTS, f"variant must be one of {list(VARIANTS)}"
P = "hcp_summary" if VARIANT == "3d_ds5" else f"hcp_summary_{VARIANT}"

# ---- expression, exactly as 22_hcp_summary.py loads the base matrix ----------
Y180 = pd.read_csv(RES / "hcp_y_maps_180.csv", index_col=0)
ul = pd.read_csv(ROOT / "data" / "reference" / "hcp_cortices" / "HCP-MMP1_UniqueRegionList.csv", encoding="utf-8-sig")
ylab = {l.lower(): l for l in Y180.index}
id2label = {int(r.regionID): ylab.get(f"lh_{r.region}".lower()) for r in ul[ul.LR == "L"].itertuples()}
X = pd.read_csv(EXPR / VARIANTS[VARIANT], index_col=0)
X = X.loc[[i for i in X.index if i in id2label]].dropna(how="all")
X = X.loc[:, X.notna().all(axis=0)]
X.index = [id2label[i] for i in X.index]
X = X.loc[[l for l in X.index if l is not None and l in Y180.index]]
Z = (X - X.mean()) / X.std()
print(f"[{VARIANT}] X:", X.shape, file=sys.stderr)

# ---- gene sets -----------------------------------------------------------------
cell = pd.read_csv(AHBA / "seidlitz_cell_genes.csv")
gcols = [c for c in cell.columns if c not in ("Type", "Paper", "Cluster", "Class")]
SETS = {("cell class", k): {x for x in pd.unique(g[gcols].to_numpy().ravel()) if isinstance(x, str) and x.strip()}
        for k, g in cell.groupby("Class") if k != "Neuro"}
may = pd.read_csv(AHBA / "maynard_layers.csv", encoding="utf-8-sig")
for lay, col in [("L1", "Layer1"), ("L2", "Layer2"), ("L3", "Layer3"), ("L4", "Layer4"),
                 ("L5", "Layer5"), ("L6", "Layer6"), ("WM", "WM")]:
    SETS[("layer", lay)] = set(may.loc[(may[f"fdr_{col}"] < 0.05) & (may[f"t_stat_{col}"] > 0), "gene"])
loc = pd.read_csv(REPO / "tools" / "bin" / "NCBI37.3.gene.loc", sep="\t", header=None,
                  names=["GENE", "chr", "start", "end", "strand", "symbol"])
loc = loc[~loc.symbol.duplicated(keep=False)]
e2s = loc.set_index("GENE").symbol
for name, f in (("SCZ", "SCZ25_META"), ("MDD", "MDD")):
    go = pd.read_csv(MAG / f"{f}.genes.out", sep=r"\s+")
    go["symbol"] = go.GENE.map(e2s)
    go = go[go.symbol.isin(Z.columns)].sort_values("ZSTAT", ascending=False)
    for n in TOPN:
        SETS[("risk genes", f"{name} top {n}")] = set(go.symbol.head(n))
    print(name, "genes with MAGMA result in X:", len(go), file=sys.stderr)

# ---- parcel maps ---------------------------------------------------------------
F = pd.DataFrame({name: Z[[c for c in Z.columns if c in gs]].mean(axis=1) for (_, name), gs in SETS.items()})
nset = {name: len([c for c in Z.columns if c in gs]) for (_, name), gs in SETS.items()}
kind = {name: k for (k, name) in SETS}
maps = pd.read_csv(RES / f"{P}_maps.csv", index_col=0)
F = F.join(maps[["PLS2", "C3"]]); kind |= {"PLS2": "axis", "C3": "axis"}; nset |= {"PLS2": np.nan, "C3": np.nan}
F = (F - F.mean()) / F.std()
T = pd.read_csv(RES / "gradient_tiers.csv", index_col=0)
lab = {l.lower(): l for l in T.index}
F.index = [lab[l.lower()] for l in F.index]
F = F.join(T)
assert F.tier.notna().all() and len(F) == len(X), (len(F), F.tier.isna().sum())
F.rename_axis("label").to_csv(RES / f"tier_profile_{VARIANT}_parcels.csv", float_format="%.4g")

rows = []
g = F.normative_thinning_um_yr
for c in [c for c in F.columns if c not in ("normative_thinning_um_yr", "tier")]:
    rho, p, _ = pls.spin_corr(g, F[c], n_perm=N_SPIN, seed=SEED, centroids=pls.HCP_CENTROIDS)
    tm = F.groupby("tier")[c].mean()
    rows.append(dict(feature=c, kind=kind[c], n_genes=nset[c], rho_with_thinning=rho, p_spin=p,
                     mean_fast=tm["fast"], mean_mid=tm["mid"], mean_slow=tm["slow"], n_parcels=int(F[c].notna().sum())))
R = pd.DataFrame(rows)
R.insert(0, "variant", VARIANT)
R.to_csv(RES / f"tier_profile_{VARIANT}.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 200)
print(R.round(3).to_string(index=False))
print("tier sizes in the expression parcels:", F.tier.value_counts().to_dict())
