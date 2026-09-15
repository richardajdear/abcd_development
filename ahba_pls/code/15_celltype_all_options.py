"""
15_celltype_all_options.py -- cell-class marker enrichment for ALL eight gene
rankings that appear in panel a of the enrichment figure, so the cell-type
profiles can be compared the way the GWAS betas are.

Supersedes 13_celltype_compare.py (three rankings) for figure purposes; that
script is kept because its own-universe values are what the earlier README
section quotes.

Two universes per vector, because z depends on what the null draws from:
  own     -- the vector's own genes (matches 13_celltype_compare.py)
  shared  -- the intersection of all eight vectors' genes, so a difference
             between two columns cannot come from a difference in universe.
The figure uses `shared`; both are written.

Test: Seidlitz et al. 2020 marker compilation; for each class, the mean weight
of its markers against 20,000 random gene sets of the same size drawn from the
universe. Marker genes are co-expressed, so the independent-gene null is
anti-conservative -- the pattern and sign are the result, not the magnitude.

Sign convention: positive = the class's markers are expressed more where
adolescent thinning is faster (published components keep their own sign, so for
C1/C3 and NSPN-PLS2 positive means "higher where the component score is high").

Output: results/celltype_all_options.tsv
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

dk2 = pd.read_csv(RES / "pls_weights" / "opt2_dCT_CT_ds25.tsv", sep="\t", index_col=0)
dk1 = pd.read_csv(RES / "pls_weights" / "opt1_dCT_ds25.tsv", sep="\t", index_col=0)
dk3 = pd.read_csv(RES / "pls_weights" / "opt3_dCT_dT1T2_ds25.tsv", sep="\t", index_col=0)
hcp = pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)
c123 = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)
nspn = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")

# same names and order as panel a of the enrichment figure
# "DK, HCP gene basis" is the DK fit on dk_3d.csv with the SAME 7,973 genes as
# the HCP matrix (12_hcp_pls.py). Without it the HCP-vs-DK contrast in this
# table confounds parcellation with the abagen build and DS gene set, which
# matters because the astrocyte sign difference is read off exactly that pair.
dkm = pd.read_csv(RES / "dk_matched_weights.tsv", sep="\t", index_col=0)
vectors = {
    "ABCD PLS2, HCP-MMP": -hcp["hcp_opt2_dCT_CT_PLS2_Z"],
    "ABCD PLS2, DK (HCP gene basis)": dkm["DK_PLS2_matchedX"],
    "ABCD PLS2, DK": -dk2["PLS2_Z"],
    "ABCD dCT alone, HCP-MMP": -hcp["hcp_opt1_dCT_PLS1_Z"],
    "ABCD dCT alone, DK": -dk1["PLS1_Z"],
    "ABCD dCT+dT1T2 PLS2, DK": -dk3["PLS2_Z"],
    "AHBA C3": c123["C3"],
    "NSPN PLS2": nspn["PLS2_z"],
    "AHBA C1": c123["C1"],
}
vectors = {k: v.dropna() for k, v in vectors.items()}
shared = set.intersection(*(set(v.index) for v in vectors.values()))
print(f"shared universe: {len(shared)} genes; own universes: "
      + ", ".join(f"{k} {len(v)}" for k, v in vectors.items()), file=sys.stderr)
assert len(shared) > 3000, len(shared)

cell = pd.read_csv(CELL)
gene_cols = [c for c in cell.columns if c not in ("Type", "Paper", "Cluster", "Class")]
classes = {cls: {x for x in pd.unique(sub[gene_cols].to_numpy().ravel())
                 if isinstance(x, str) and x.strip()}
           for cls, sub in cell.groupby("Class")}

rows = []
for universe in ("own", "shared"):
    for vname, w in vectors.items():
        wv = w.loc[sorted(shared)] if universe == "shared" else w
        zv = wv.to_numpy()
        rng = np.random.default_rng(SEED)     # same draws for every column
        for cls, gset in classes.items():
            genes = sorted(gset & set(wv.index))
            if len(genes) < 10:
                continue
            obs = wv.loc[genes].mean()
            null = np.array([zv[rng.choice(len(zv), len(genes), replace=False)].mean()
                             for _ in range(N_PERM)])
            p = (1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (1 + N_PERM)
            rows.append(dict(universe=universe, vector=vname, cell_class=cls,
                             n_markers=len(genes), n_universe=len(wv), mean_weight=obs,
                             null_mean=null.mean(), null_sd=null.std(),
                             z=(obs - null.mean()) / null.std(), p_perm=p))
        print(universe, vname, file=sys.stderr)

T = pd.DataFrame(rows)
T.to_csv(RES / "celltype_all_options.tsv", sep="\t", index=False, float_format="%.4g")

S = T[T.universe == "shared"]
piv = S.pivot(index="cell_class", columns="vector", values="z")[list(vectors)]
print(f"\nshared universe (n = {S.n_universe.iloc[0]} genes), marker-enrichment z:")
print(piv.round(1).to_string())
print("\nprofile agreement across rankings (Spearman over the 9 classes):")
print(piv.corr(method="spearman").round(2).to_string())
print("\nAstro / Oligo detail:")
print(S[S.cell_class.isin(["Astro", "Oligo"])]
      .pivot(index="vector", columns="cell_class", values=["z", "p_perm"]).round(3).to_string())
