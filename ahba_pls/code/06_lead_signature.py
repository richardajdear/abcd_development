"""
06_lead_signature.py -- characterise the lead component: option 2 (dCT + CT), PLS2.

Orientation: pls.py fixes the dCT salience positive, so positive Lx / Z mean
"associated with LESS thinning".  Here everything is re-oriented to the
THINNING direction (thinning_Z = -Z, thinning_scores = -Lx) so that positive
weight = higher expression where adolescent thinning is FASTER, matching the
sign in which NSPN-PLS2 and AHBA-C3 correlate positively with it.

Outputs
  results/lead_signature_weights.tsv    gene, thinning_Z at ds0/ds25/ds50, U, rank, decile flags
  results/lead_signature_scores.csv     33 regions, thinning-oriented scores + references
  results/lead_top_bottom_genes.tsv     top/bottom 25 genes (ds25)
  results/lead_overlap_with_references.tsv  top/bottom-decile overlaps with C3 and NSPN-PLS2 (hypergeometric)
  results/lead_celltype_enrichment.tsv  Seidlitz 2020 cell-class marker enrichment (permutation)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls

ROOT = HERE.parent
DATA, RES, REF, FIG = ROOT / "data", ROOT / "results", ROOT / "data" / "reference", ROOT / "figures"
FIG.mkdir(exist_ok=True)
LEAD_OPT, LEAD_COMP, LEAD_DS = "opt2_dCT_CT", "PLS2", "ds25"
rng = np.random.default_rng(0)

# ---- weights, thinning-oriented, all DS levels --------------------------------
W = {}
for ds in ("ds0", "ds25", "ds50"):
    t = pd.read_csv(RES / "pls_weights" / f"{LEAD_OPT}_{ds}.tsv", sep="\t", index_col=0)
    W[ds] = -t[f"{LEAD_COMP}_Z"]
lead = pd.DataFrame({f"thinning_Z_{ds}": W[ds] for ds in W})
lead["U_ds25_thinning"] = -pd.read_csv(RES / "pls_weights" / f"{LEAD_OPT}_ds25.tsv", sep="\t", index_col=0)[f"{LEAD_COMP}_U"]
lead["rank_ds25"] = lead["thinning_Z_ds25"].rank(ascending=False)
n25 = lead["thinning_Z_ds25"].notna().sum()
lead["decile_ds25"] = pd.qcut(lead["thinning_Z_ds25"], 10, labels=False) + 1
lead.index.name = "gene"
lead.sort_values("thinning_Z_ds25", ascending=False).to_csv(
    RES / "lead_signature_weights.tsv", sep="\t", float_format="%.5g")

z25 = lead["thinning_Z_ds25"].dropna().sort_values(ascending=False)
tb = pd.concat([z25.head(25).rename("thinning_Z").to_frame().assign(tail="top"),
                z25.tail(25).rename("thinning_Z").to_frame().assign(tail="bottom")])
tb.to_csv(RES / "lead_top_bottom_genes.tsv", sep="\t", float_format="%.3g")

# ---- scores --------------------------------------------------------------------
sc = pd.read_csv(RES / "pls_scores" / f"{LEAD_OPT}_{LEAD_DS}.csv", index_col=0)
Y34 = pd.read_csv(DATA / "y_maps_bilateral_34.csv", index_col=0)
nspn = pd.read_csv(REF / "nspn_dk_maps_bilateral_34.csv", index_col=0)
c3r = pd.read_csv(REF / f"ahba_c123_scores_recomputed_{LEAD_DS}.csv", index_col=0)
S = pd.DataFrame({"thinning_scores_gene_side": -sc[f"{LEAD_COMP}_gene_scores"],
                  "thinning_scores_imaging_side": -sc[f"{LEAD_COMP}_imaging_scores"]})
S["dCT"] = Y34.loc[S.index, "dCT"]; S["CT"] = Y34.loc[S.index, "CT"]
S["NSPN_PLS2"] = nspn.loc[S.index, "PLS2"]; S["NSPN_CT_delta"] = nspn.loc[S.index, "CT_delta"]
S["C3_recomputed"] = c3r.loc[S.index, "C3"]
S.to_csv(RES / "lead_signature_scores.csv", float_format="%.5g")

# ---- overlap of top/bottom deciles with reference top/bottom deciles ------------
nspn_w = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")["PLS2_z"]
c3_w = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)["C3"]
orows = []
for rn, rw in (("NSPN_PLS2_z", nspn_w), ("C3", c3_w)):
    shared = z25.index.intersection(rw.dropna().index)
    a, b = z25.loc[shared], rw.loc[shared]
    n = len(shared); k = n // 10
    for tail in ("top", "bottom"):
        A = set(a.nlargest(k).index if tail == "top" else a.nsmallest(k).index)
        B = set(b.nlargest(k).index if tail == "top" else b.nsmallest(k).index)
        ov = len(A & B)
        p = stats.hypergeom.sf(ov - 1, n, k, k)
        orows.append(dict(reference=rn, tail=tail, n_shared=n, decile_size=k, overlap=ov,
                          expected=k * k / n, fold=ov / (k * k / n), p_hypergeom=p,
                          rho_all_shared=stats.spearmanr(a, b).statistic))
pd.DataFrame(orows).to_csv(RES / "lead_overlap_with_references.tsv", sep="\t", index=False, float_format="%.4g")

# ---- cell-type marker enrichment (Seidlitz et al. 2020 compilation) -----------
cell = pd.read_csv("/Users/richard/Git/AHBA/data/seidlitz_cell_genes.csv")
gene_cols = [c for c in cell.columns if c not in ("Type", "Paper", "Cluster", "Class")]
classes = {}
for cls, sub in cell.groupby("Class"):
    g = set(pd.unique(sub[gene_cols].to_numpy().ravel()))
    g = {x for x in g if isinstance(x, str) and x.strip()}
    classes[cls] = sorted(g & set(z25.index))
univ = z25.index.to_numpy(); zv = z25.to_numpy()
N_PERM = 20000
crows = []
for cls, genes in classes.items():
    if len(genes) < 10:
        continue
    obs = z25.loc[genes].mean()
    null = np.array([zv[rng.choice(len(zv), len(genes), replace=False)].mean() for _ in range(N_PERM)])
    p = (1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (1 + N_PERM)
    crows.append(dict(cell_class=cls, n_markers=len(genes), mean_thinning_Z=obs,
                      null_mean=null.mean(), null_sd=null.std(), z=(obs - null.mean()) / null.std(), p_perm=p))
ct = pd.DataFrame(crows).sort_values("z", ascending=False)
ct.to_csv(RES / "lead_celltype_enrichment.tsv", sep="\t", index=False, float_format="%.4g")

# ---- figure ---------------------------------------------------------------------
# Figure 1 is built in R: code/fig1_lead_signature.R (ggplot2 + patchwork + ggseg).
# This script writes only tables; nothing here renders.

print("top 15:", ", ".join(z25.head(15).index))
print("bottom 15:", ", ".join(z25.tail(15).index[::-1]))
print(pd.DataFrame(orows)[["reference", "tail", "overlap", "expected", "fold", "p_hypergeom"]].round(3).to_string(index=False))
print(ct[["cell_class", "n_markers", "mean_thinning_Z", "z", "p_perm"]].round(3).to_string(index=False))
print("DS stability of thinning_Z (Spearman):", lead[[f"thinning_Z_{d}" for d in W]].corr("spearman").round(3).iloc[0].to_dict())
