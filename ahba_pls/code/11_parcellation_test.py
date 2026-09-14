"""
11_parcellation_test.py -- is the ABCD signature's weaker SCZ/MDD enrichment a
consequence of the DK parcellation (33 regions) rather than of the phenotype?

The shipped AHBA C3 weights were derived on the HCP-MMP matrix
(AHBA/data/abagen-data/expression/hcp_3d_ds5.csv: 137 left parcels x 7,973
genes).  This project's PLS ran on the DK matrix at the same DS filter
(dk_3d_ds5.csv: 33 x 7,973) -- 4.2x fewer observations per gene weight.  Fewer
regions means noisier weights, and noise in a regressor attenuates the
gene-property beta, so granularity is a live explanation for C3's ~2x larger
effect.

Controlled comparison (identical genes, identical method, parcellation the only
difference):

  1. PCA on z-scored genes in HCP-137 and in DK-33; take the component whose
     gene loadings best match the shipped C3 weights.  Gene weights = loadings.
  2. An "oracle" PLS: Y = the HCP-derived C3 score map re-expressed in DK space
     (data/reference/ahba_c123_scores_recomputed_ds25.csv is the DK projection
     of the shipped weights).  This asks what a 33-region PLS can recover when
     the imaging map is *exactly* C3 -- the ceiling set by DK resolution and by
     the PLS step itself, with no phenotype noise at all.
  3. All weight vectors, plus the shipped C3 and the ABCD PLS2 signature, go
     into one MAGMA gene-property run on a single shared universe so the betas
     are directly comparable.

Reads:  AHBA/data/abagen-data/expression/{hcp,dk}_3d_ds5.csv, data/weights.csv,
        results/lead_signature_weights.tsv, hpc/work/results/magma/*.genes.raw
Writes: results/parcellation_test_weights.tsv, results/parcellation_test.tsv
"""
from __future__ import annotations
import subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls

ROOT = HERE.parent
REPO = ROOT.parent
RES, REF = ROOT / "results", ROOT / "data" / "reference"
GS = REF / "gene_sets"
RUN = RES / "magma_runs"; RUN.mkdir(exist_ok=True)
EXPR = Path.home() / "Git" / "AHBA" / "data" / "abagen-data" / "expression"
MAGMA = REPO / "tools" / "bin" / ("magma_mac/magma" if (REPO / "tools/bin/magma_mac/magma").exists() else "magma")
RAW = {"SCZ": REPO / "hpc/work/results/magma/SCZ.genes.raw",
       "MDD": REPO / "hpc/work/results/magma/MDD.genes.raw"}

# ---------------------------------------------------------------- inputs -----
# The ds5 (top-8k differential stability) gene sets are parcellation-SPECIFIC --
# hcp_3d_ds5 and dk_3d_ds5 share only ~7.0k of their 7,973 genes -- so the
# unfiltered matrices are subset to the exact gene list the shipped C3 weights
# were fitted on. Parcellation is then the only difference between the two fits.
Xh = pd.read_csv(EXPR / "hcp_3d.csv", index_col=0)
Xd = pd.read_csv(EXPR / "dk_3d.csv", index_col=0)
c123 = pd.read_csv(REPO / "data" / "weights.csv", index_col=0)
genes = [g for g in c123.index if g in Xh.columns and g in Xd.columns]
assert len(genes) == len(c123), f"{len(genes)} of {len(c123)} shipped genes found"
Xh, Xd, c123 = Xh[genes], Xd[genes], c123.loc[genes]
print(f"HCP {Xh.shape}  DK {Xd.shape}  shipped weights {c123.shape}", file=sys.stderr)


def pca_weights(X: pd.DataFrame, k: int = 5) -> tuple[pd.DataFrame, np.ndarray]:
    """PCA (via SVD of the centred, gene-z-scored matrix); gene loadings per PC."""
    Xz = pls.zscore_cols(X.to_numpy(float))
    Xz = Xz - Xz.mean(0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    var = S ** 2 / (S ** 2).sum()
    W = pd.DataFrame(Vt[:k].T, index=X.columns, columns=[f"PC{i+1}" for i in range(k)])
    return W, var[:k]


def match_to(W: pd.DataFrame, target: pd.Series) -> tuple[str, pd.Series, float]:
    """Pick the component whose loadings best match a target weight vector."""
    rho = {c: stats.spearmanr(W[c], target).statistic for c in W.columns}
    best = max(rho, key=lambda c: abs(rho[c]))
    sgn = np.sign(rho[best])
    return best, W[best] * sgn, rho[best] * sgn


Wh, vh = pca_weights(Xh)
Wd, vd = pca_weights(Xd)
rows_match = []
vectors: dict[str, pd.Series] = {"C3_shipped_DME_HCP137": c123["C3"]}
for comp in ("C1", "C2", "C3"):
    for tag, W, var in (("HCP137", Wh, vh), ("DK33", Wd, vd)):
        name, w, rho = match_to(W, c123[comp])
        rows_match.append(dict(target=comp, parcellation=tag, matched=name,
                               rho_to_shipped=rho, var_explained=var[int(name[-1]) - 1]))
        if comp == "C3":
            vectors[f"C3_PCA_{tag}"] = w
print(pd.DataFrame(rows_match).round(3).to_string(index=False), file=sys.stderr)

# ---- oracle PLS: Y = the C3 map itself, at DK resolution --------------------
c3_dk = pd.read_csv(REF / "ahba_c123_scores_recomputed_ds25.csv", index_col=0)["C3"]
X25 = pd.read_parquet(ROOT / "data" / f"X_ds25.parquet")
regions33 = [l.strip() for l in open(ROOT / "data" / "region_order_33.txt") if l.strip()]
X25 = X25.loc[regions33]
Yor = c3_dk.loc[regions33].to_frame("C3_map")
orc = pls.pls_svd(X25, Yor)
w_or = orc.weights()["PLS1"]
# orient to the shipped C3
if stats.spearmanr(w_or, c123["C3"].reindex(w_or.index)).statistic < 0:
    w_or = -w_or
vectors["C3_oraclePLS_DK33"] = w_or

# ---- the actual ABCD signature (thinning-oriented), for reference ----------
lead = pd.read_csv(RES / "lead_signature_weights.tsv", sep="\t", index_col=0)
vectors["ABCD_PLS2_thinning_DK33"] = lead["thinning_Z_ds25"]

W_all = pd.DataFrame(vectors)
W_all.index.name = "gene"
W_all.to_csv(RES / "parcellation_test_weights.tsv", sep="\t", float_format="%.6g")

# ---------------------------------------------------------------- MAGMA ------
sym2ent = (pd.read_csv(GS / "magma_SCZ_genes.tsv", sep="\t").dropna(subset=["symbol"])
           .drop_duplicates("symbol").set_index("symbol")["GENE"].astype(int))
cov = W_all.dropna(how="any")                      # one shared universe
cov = cov.loc[cov.index.intersection(sym2ent.index)]
cov.index = sym2ent.loc[cov.index].values
cov = cov[~cov.index.duplicated()]
cov.index.name = "GENE"
path = RUN / "parcellation_test.covar"
cov.to_csv(path, sep="\t", float_format="%.6g")
print(f"shared universe: {len(cov)} genes", file=sys.stderr)

out = []
for dis, raw in RAW.items():
    prefix = RUN / f"parcellation_{dis}"
    cmd = [str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(path),
           "--out", str(prefix)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"MAGMA failed for {dis}:\n{r.stdout[-1500:]}\n{r.stderr[-500:]}")
    t = pd.read_csv(f"{prefix}.gsa.out", sep=r"\s+", comment="#")
    if "FULL_NAME" in t.columns:
        t["VARIABLE"] = t["FULL_NAME"].fillna(t["VARIABLE"])
    t.insert(0, "disorder", dis)
    out.append(t)
M = pd.concat(out, ignore_index=True)
M["n_regions"] = M.VARIABLE.map(lambda v: 137 if "HCP137" in v else 33)
M.to_csv(RES / "parcellation_test.tsv", sep="\t", index=False, float_format="%.5g")

order = ["C3_shipped_DME_HCP137", "C3_PCA_HCP137", "C3_PCA_DK33",
         "C3_oraclePLS_DK33", "ABCD_PLS2_thinning_DK33"]
tab = (M[M.VARIABLE.isin(order)]
       .pivot(index="VARIABLE", columns="disorder", values=["BETA_STD", "P"])
       .loc[order].round(4))
print("\nMAGMA gene-property, shared universe of", int(M.NGENES.iloc[0]), "genes:")
print(tab.to_string())
print("\ngene-weight agreement with the shipped C3 (Spearman):")
print({k: round(stats.spearmanr(v.reindex(c123.index).dropna(),
                                c123["C3"].reindex(v.reindex(c123.index).dropna().index)).statistic, 3)
       for k, v in vectors.items() if k != "C3_shipped_DME_HCP137"})
