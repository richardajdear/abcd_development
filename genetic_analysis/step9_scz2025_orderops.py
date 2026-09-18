#!/usr/bin/env python
"""Order-of-operations table for the 2025-SCZ scores: per-region-BLUP-mean
global_slope (step 9 main tables) vs the single-LMM slope
(prs_scz2025_1lmm/), matched and secondary cells, both atlases.

Written to results_70tab*/prs_scz2025/table_order_of_operations_scz2025.tsv
(whichever atlas has run; rerun after both to complete).
"""
from __future__ import annotations
import glob, re, sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"dk": REPO / "genetic_analysis/work/results_70tab", "hcp": REPO / "genetic_analysis/work/results_70tab_hcp"}
METHODS = ["CT", "PRSCS", "SBayesR", "SBayesRC", "PRSCSX"]
PHENOS = ["global_slope", "baseline_thickness"]
# (trait_arm, target stratum, score version) -> cell label
CELLS = {("SCZ25_EUR", "EUR", "raw"): "EUR GWAS -> EUR (primary)",
         ("SCZ25_META", "full", "zanc"): "multi GWAS -> pooled zanc (primary)",
         ("SCZ25_MATCHED", "full", "raw"): "ancestry-matched composite -> pooled (primary)",
         ("SCZ25_META", "EUR", "raw"): "multi GWAS -> EUR (secondary)",
         ("SCZ25_EUR", "full", "zanc"): "EUR weights -> pooled zanc (secondary)"}

rows = []
for atlas, root in ROOTS.items():
    base_p = root / "prs_scz2025/table_scz2025_main.tsv"
    if not base_p.exists() or not (root / "prs_scz2025_1lmm").exists():
        continue
    base = pd.read_csv(base_p, sep="\t")
    base = base[base.phenotype.isin(PHENOS)]
    for f in glob.glob(str(root / "prs_scz2025_1lmm/assoc_*.tsv")):
        m = re.match(rf".*/assoc_({'|'.join(METHODS)})_(SCZ25_\w+?)(_zanc)?\.tsv$", f)
        if not m:
            continue
        meth, arm, z = m.group(1), m.group(2), ("zanc" if m.group(3) else "raw")
        t = pd.read_csv(f, sep="\t")
        for (a, st, sv), label in CELLS.items():
          if a != arm or sv != z:
            continue
          for ph in PHENOS:
            one = t[(t.phenotype == f"{ph}_1lmm") & (t.stratum == st)]
            b = base[(base.trait_arm == arm) & (base.method == meth) & (base.score == z) & (base.target_stratum == st) & (base.phenotype == ph)]
            if one.empty or b.empty:
                continue
            one = one.sort_values("p").iloc[0]; b = b.sort_values("p").iloc[0]
            rows.append(dict(atlas=atlas, phenotype=ph, cell=label, trait_arm=arm, stratum=st, method=meth, score=z,
                             beta_perregion=b.beta, se_perregion=b.se, p_adj_perregion=b.p_adj, thr_perregion=b.threshold,
                             beta_1lmm=one.beta, se_1lmm=one.se, p_adj_1lmm=one.p_adj, thr_1lmm=one.threshold, n=int(one.n)))
out = pd.DataFrame(rows)
if out.empty:
    sys.exit("no rows")
out["_c"] = out.cell.map({c: i for i, c in enumerate(CELLS.values())}); out["_m"] = out.method.map({m: i for i, m in enumerate(METHODS)})
out = out.sort_values(["atlas", "phenotype", "_c", "_m"]).drop(columns=["_c", "_m"])
for atlas, root in ROOTS.items():
    sub = out[out.atlas == atlas]
    if not sub.empty:
        sub.to_csv(root / "prs_scz2025/table_order_of_operations_scz2025.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
f = out.copy()
for c in ["beta_perregion", "se_perregion", "beta_1lmm", "se_1lmm"]: f[c] = f[c].round(4)
for c in ["p_adj_perregion", "p_adj_1lmm"]: f[c] = f[c].map(lambda v: f"{v:.3g}")
print(f.drop(columns=["trait_arm", "stratum", "score", "n", "thr_perregion", "thr_1lmm"]).to_string(index=False))
d = out.beta_1lmm - out.beta_perregion
print(f"\nbeta_1lmm - beta_perregion: mean {d.mean():+.4f}, max |diff| {d.abs().max():.4f}; SE ratio 1lmm/perregion mean {(out.se_1lmm/out.se_perregion).mean():.3f}")
