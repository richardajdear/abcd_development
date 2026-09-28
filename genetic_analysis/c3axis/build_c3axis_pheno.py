"""Export the C3-axis single-LMM slopes as a GCTA/GENESIS/prs_assoc phenotype set.

Multi-label version of orderops/build_1lmm_pheno.py.  Reads the fitted
out/<run>_c3axis/fits/blups.parquet, aligns to the analysis set of an existing
export (phenotypes_gcta.txt with FID = family id), z-scores over that set
(rule 3), and writes:

  <out>/phenotypes_gcta.txt          FID = family id (prs_assoc.R, GENESIS)
  <out>_fidiid/phenotypes_gcta.txt   FID = IID (GCTA GREML, rule 2)
  <out>/phenotype_manifest.tsv       name, mpheno, priority, role
  covar_quant.txt / covar_categorical.txt copied from the reference export
  <out>_adjG/                        same phenotypes; covar_quant.txt carries the
                                     extra column gcov_global_slope (the global_mean
                                     single-LMM slope), for the PRS sensitivity model

Per-subject outputs: never commit (rule 16).

Usage:
    python genetic_analysis/c3axis/build_c3axis_pheno.py \
        out/thickness_hcp_70_aa6e91efba82_c3axis \
        genetic_analysis/work/pheno_70tab_hcp \
        genetic_analysis/work/results_70tab_hcp/c3axis/pheno
"""
import sys
from pathlib import Path

import pandas as pd

ROLES = {  # name in the fit -> (exported name, priority, role)
    "c3axis_rc":   ("c3axis_rc",   1, "PRIMARY: LH row-centred component matched to AHBA C3 (rc5)"),
    "proj_C3":     ("proj_C3",     2, "secondary: row-centred projection on AHBA C3"),
    "proj_dCT":    ("proj_dCT",    2, "secondary: amplitude of the normative thinning contrast (the dCT-aligned component)"),
    "proj_PLS2":   ("proj_PLS2",   3, "secondary: row-centred projection on PLS2 (reliability ~0.10)"),
    "c1axis_rc":   ("c1axis_rc",   4, "specificity control: C1 axis (rc1)"),
    "c2axis_rc":   ("c2axis_rc",   4, "specificity control: C2 axis (rc2)"),
    "global_mean": ("global_slope_c3axis", 0, "check: must reproduce global_slope_1lmm (r > 0.99)"),
}

run, ref_dir, out = (Path(p) for p in sys.argv[1:4])
b = pd.read_parquet(run / "fits/blups.parquet")
W = b.pivot(index="subject", columns="label", values="re_slope")
missing = set(ROLES) - set(W.columns)
assert not missing, f"fit lacks {missing}"
W.index = W.index.str.extract(r"([A-Z0-9]{8})$")[0].to_numpy()

ref = pd.read_csv(ref_dir / "phenotypes_gcta.txt", sep=r"\s+", dtype={"FID": str, "IID": str})
ref["_t"] = ref.IID.str.extract(r"([A-Z0-9]{8})$")[0]
m = ref[["FID", "IID", "_t"]].merge(W, left_on="_t", right_index=True, how="inner")
assert len(m) == len(ref), (len(m), len(ref))
cols = list(ROLES)
m[cols] = (m[cols] - m[cols].mean()) / m[cols].std()
exp = m[["FID", "IID"]].copy()
for k, (name, _, _) in ROLES.items():
    exp[name] = m[k]
names = [ROLES[k][0] for k in ROLES]

out.mkdir(parents=True, exist_ok=True)
fidiid = out.parent / f"{out.name}_fidiid"; fidiid.mkdir(exist_ok=True)
adjg = out.parent / f"{out.name}_adjG"; adjg.mkdir(exist_ok=True)
exp.to_csv(out / "phenotypes_gcta.txt", sep=" ", index=False)
exp.assign(FID=exp.IID).to_csv(fidiid / "phenotypes_gcta.txt", sep=" ", index=False)
exp.to_csv(adjg / "phenotypes_gcta.txt", sep=" ", index=False)
man = pd.DataFrame([dict(name=ROLES[k][0], mpheno=i + 1, priority=ROLES[k][1], role=ROLES[k][2],
                         n_nonmissing=int(exp[ROLES[k][0]].notna().sum())) for i, k in enumerate(ROLES)])
for d in (out, fidiid, adjg):
    man.to_csv(d / "phenotype_manifest.tsv", sep="\t", index=False)
    for f in ("covar_quant.txt", "covar_categorical.txt"):
        src = pd.read_csv(ref_dir / f, sep=r"\s+", dtype={"FID": str, "IID": str})
        if d == fidiid:
            src["FID"] = src["IID"]
        src.to_csv(d / f, sep=" ", index=False)
q = pd.read_csv(adjg / "covar_quant.txt", sep=r"\s+", dtype={"FID": str, "IID": str})
q.merge(exp[["IID", "global_slope_c3axis"]].rename(columns={"global_slope_c3axis": "gcov_global_slope"}),
        on="IID", how="left").to_csv(adjg / "covar_quant.txt", sep=" ", index=False)
man[man.name != "global_slope_c3axis"].to_csv(adjg / "phenotype_manifest.tsv", sep="\t", index=False)

# checks printed for the log
if (ref_dir / "phenotypes_gcta.txt").exists() and "global_slope" in ref.columns:
    r = m.merge(ref[["_t", "global_slope"]], on="_t").pipe(lambda d: d.global_mean.corr(d.global_slope))
    print(f"r(global_slope_c3axis, reference per-region global_slope) = {r:.3f} (expect ~0.82 on HCP)")
print(exp[names].corr().round(2).to_string())
print(f"wrote {out}, {fidiid}, {adjg}: n = {len(exp)}, {len(names)} phenotypes")
