"""03b -- collect the local PRS screen into matched-cell tables (README_HPC.md 8.5).

table_c3axis_prs_local.tsv        lmer betas from tools/prs_assoc.R, rule-4 matched cells only
                                  (EUR GWAS raw score -> EUR stratum; multi-ancestry / pooled GWAS
                                  _zanc score -> full sample), with and without the global covariate.
table_c3axis_prs_delta_local.tsv  SCZ 2025: OLS with family-clustered SE for each phenotype, and paired
                                  delta-beta as the PRS coefficient on the difference score (exact:
                                  same design matrix, same children). Model as prs_assoc.R, which in
                                  this pipeline has no age term (covar_quant carries baseline_age, not age_c).
Run from the repo root after 03_prs_screen_local.sh.
"""
import glob
import os

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, REPO

C3 = REPO / "genetic_analysis/work/results_70tab_hcp/c3axis"
rows = []
for v in ("", "_adjG"):
    for f in glob.glob(str(C3 / f"prs_local{v}/assoc_*.tsv")):
        t = pd.read_csv(f, sep="\t"); t["cell"] = os.path.basename(f)[6:-4]; t["adjG"] = bool(v); rows.append(t)
P = pd.concat(rows)
P["method"] = P.cell.str.split("__").str[0]; P["arm"] = P.cell.str.split("__").str[1]
P["zanc"] = P.cell.str.contains("_zanc"); eur = P.arm.str.contains("EUR|eur")
M = P[(eur & ~P.zanc & (P.stratum == "EUR")) | (~eur & P.zanc & (P.stratum == "full"))].copy()
M["gwas"] = M.arm.replace({"SCZ25_EUR": "SCZ2025", "SCZ25_META": "SCZ2025", "SCZ_eur": "SCZ_PGC3",
                           "SCZ_pooled": "SCZ_PGC3", "MDD_eur": "MDD", "MDD_pooled": "MDD"})
M["arm_type"] = np.where(M.arm.str.contains("EUR|eur"), "EUR", "pooled")
M[["gwas", "arm_type", "method", "phenotype", "adjG", "n", "beta", "se", "p"]].sort_values(
    ["phenotype", "adjG", "gwas", "arm_type", "method"]).to_csv(HERE / "table_c3axis_prs_local.tsv", sep="\t",
                                                               index=False, float_format="%.5g")

# ---- paired delta-beta, SCZ 2025 ----------------------------------------------------------
ph = pd.read_csv(C3 / "pheno/phenotypes_gcta.txt", sep=" ", dtype={"FID": str, "IID": str})
qc = pd.read_csv(C3 / "pheno/covar_quant.txt", sep=" ", dtype={"FID": str, "IID": str})
cc = pd.read_csv(C3 / "pheno/covar_categorical.txt", sep=" ", dtype={"FID": str, "IID": str})
eur_ids = set(pd.read_csv(REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep", sep=r"\s+", header=None,
                          dtype=str).iloc[:, 1])
base = ph.merge(qc.drop(columns=["baseline_age", "n_visits"]), on=["FID", "IID"]).merge(
    cc[["FID", "IID", "sex"]], on=["FID", "IID"])
PH = ["c3axis_rc", "proj_C3", "proj_dCT", "proj_PLS2", "c1axis_rc", "c2axis_rc", "global_slope_c3axis"]
cells = {"SCZ2025 EUR PRSCS": ("scores_scz2025/PRSCS/SCZ25_EUR", True),
         "SCZ2025 EUR SBayesRC": ("scores_scz2025/SBayesRC/SCZ25_EUR", True),
         "SCZ2025 pooled PRSCS": ("scores_scz2025/PRSCS/SCZ25_META/_zanc", False),
         "SCZ2025 pooled SBayesRC": ("scores_scz2025/SBayesRC/SCZ25_META/_zanc", False)}
pcs = " + ".join(f"PC{i}" for i in range(1, 11))
out = []
for name, (rel, is_eur) in cells.items():
    hits = glob.glob(str(REPO / "genetic_analysis/work" / rel / "score_*.profile"))
    if not hits:
        continue
    pr = pd.read_csv(hits[0], sep=r"\s+", dtype={"FID": str, "IID": str})
    sc = [c for c in pr.columns if c.upper().startswith("SCORE")][0]
    d = base.merge(pr[["IID", sc]].rename(columns={sc: "PRS"}), on="IID")
    if is_eur:
        d = d[d.IID.isin(eur_ids)]
    for c in ["PRS"] + PH:
        d[c] = (d[c] - d[c].mean()) / d[c].std()
    g = d.FID.astype("category").cat.codes

    def fit(y):
        m = smf.ols(f"{y} ~ PRS + sex + {pcs}", data=d).fit(cov_type="cluster", cov_kwds={"groups": g})
        return m.params["PRS"], m.bse["PRS"], m.pvalues["PRS"]
    contrasts = [(y, y) for y in PH] + [(f"{a} - {b}", (a, b)) for a, b in
                                         [("c3axis_rc", "c1axis_rc"), ("c3axis_rc", "c2axis_rc"), ("proj_C3", "c1axis_rc")]]
    for lab, y in contrasts:
        if isinstance(y, tuple):
            d["_diff"] = d[y[0]] - d[y[1]]; y = "_diff"
        b, s, p = fit(y)
        out.append(dict(cell=name, contrast=lab, beta=b, se=s, p=p, ci_lo=b - 1.96 * s, ci_hi=b + 1.96 * s, n=len(d)))
pd.DataFrame(out).to_csv(HERE / "table_c3axis_prs_delta_local.tsv", sep="\t", index=False, float_format="%.5g")
print(len(M), "matched rows;", len(out), "delta rows")
