"""Gate B: per-SNP z (BETA/SE) of the new slope_thk EUR scan against
scan_1lmm/assoc_eur/global_slope_1lmm, plus lambda_GC from gwas_summary.tsv.

  python gate_b.py <results_70tab_hcp dir> [arm: eur|pooled]
"""
import sys
import pandas as pd

R = sys.argv[1]; arm = sys.argv[2] if len(sys.argv) > 2 else "eur"
suf = "_eur" if arm == "eur" else ""
rd = lambda f: pd.read_csv(f, sep="\t", usecols=["SNP", "A1", "BETA", "SE"])
a = rd(f"{R}/scan_1lmm/assoc{suf}/global_slope_1lmm.sumstats.tsv.gz")
b = rd(f"{R}/scan_geom/assoc{suf}/slope_thk.sumstats.tsv.gz")
m = a.merge(b, on="SNP", suffixes=("_ref", "_new"))
m = m[m.A1_ref == m.A1_new]
r = (m.BETA_ref / m.SE_ref).corr(m.BETA_new / m.SE_new)
s = pd.read_csv(f"{R}/scan_geom/assoc{suf}/gwas_summary.tsv", sep="\t").set_index("phenotype")
ref = pd.read_csv(f"{R}/scan_1lmm/assoc{suf}/gwas_summary.tsv", sep="\t").set_index("phenotype")
lam, lam_ref = s.loc["slope_thk", "lambda_gc"], ref.loc["global_slope_1lmm", "lambda_gc"]
print(f"arm {arm}: SNPs ref {len(a)}, new {len(b)}, matched {len(m)}; r(z) = {r:.5f}")
print(f"lambda_GC slope_thk {lam:.4f} vs global_slope_1lmm {lam_ref:.4f}")
ok = r >= 0.98 and abs(lam - lam_ref) <= 0.01
print(f"GATE B ({arm}): {'PASS' if ok else 'FAIL'}")
