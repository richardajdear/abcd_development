"""
13_celltype_compare.py -- cell-class marker enrichment for three gene rankings:
the ABCD thinning signature in DK, the same in HCP-MMP, and AHBA C3.

Same test as 06_lead_signature.py (Seidlitz et al. 2020 marker compilation,
mean weight of a class against 20,000 random gene sets of equal size drawn from
the vector's own universe), run on all three so the profiles are comparable.
Each vector keeps its own universe -- the three matrices do not share a gene
list -- so compare the z pattern across classes, not the absolute z.

Note the z values are descriptive: marker genes are co-expressed, so the
permutation null (independent genes) is anti-conservative.  Sign convention:
positive = the class's markers are expressed more where thinning is faster
(for C3, where the C3 score is higher).

Output: results/celltype_compare.tsv
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES, REF = ROOT / "results", ROOT / "data" / "reference"
CELL = Path.home() / "Git" / "AHBA" / "data" / "seidlitz_cell_genes.csv"
N_PERM, SEED = 20000, 0
rng = np.random.default_rng(SEED)

vectors = {
    "ABCD PLS2, DK": pd.read_csv(RES / "lead_signature_weights.tsv", sep="\t",
                                 index_col=0)["thinning_Z_ds25"].dropna(),
    "ABCD PLS2, HCP-MMP": -pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t",
                                       index_col=0)["hcp_opt2_dCT_CT_PLS2_Z"].dropna(),
    "AHBA C3": pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)["C3"].dropna(),
}

cell = pd.read_csv(CELL)
gene_cols = [c for c in cell.columns if c not in ("Type", "Paper", "Cluster", "Class")]
classes = {}
for cls, sub in cell.groupby("Class"):
    g = {x for x in pd.unique(sub[gene_cols].to_numpy().ravel()) if isinstance(x, str) and x.strip()}
    classes[cls] = g

rows = []
for vname, w in vectors.items():
    zv = w.to_numpy()
    for cls, gset in classes.items():
        genes = sorted(gset & set(w.index))
        if len(genes) < 10:
            continue
        obs = w.loc[genes].mean()
        null = np.array([zv[rng.choice(len(zv), len(genes), replace=False)].mean()
                         for _ in range(N_PERM)])
        p = (1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (1 + N_PERM)
        rows.append(dict(vector=vname, cell_class=cls, n_markers=len(genes),
                         n_universe=len(w), mean_weight=obs, null_mean=null.mean(),
                         null_sd=null.std(), z=(obs - null.mean()) / null.std(), p_perm=p))
    print(vname, len(w), "genes", file=sys.stderr)

T = pd.DataFrame(rows)
T.to_csv(RES / "celltype_compare.tsv", sep="\t", index=False, float_format="%.4g")
piv = T.pivot(index="cell_class", columns="vector", values="z").round(1)
piv["n_markers"] = T.groupby("cell_class").n_markers.min()
print(piv.sort_values("ABCD PLS2, DK", ascending=False).to_string())
print("\nrank agreement of the class profiles (Spearman):")
p2 = T.pivot(index="cell_class", columns="vector", values="z")
print(p2.corr(method="spearman").round(2).to_string())
