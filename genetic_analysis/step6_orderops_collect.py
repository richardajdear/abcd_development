#!/usr/bin/env python
"""Order-of-operations table, full panel: per-region-BLUP-mean phenotypes (step 6,
prs_final/table_main.tsv) vs the single-LMM construction (prs_final_1lmm/
assoc_*.tsv), for global_slope AND baseline_thickness, every trait arm in its
rule-4 matched stratum (SCZ/MDD pooled -> full zanc; everything else -> EUR raw),
both atlases.  Writes results_70tab*/prs_final_1lmm/table_order_of_operations_all.tsv
(the original table_order_of_operations.tsv, SCZ/MDD x global_slope, is left as
the slide script expects it).
"""
from __future__ import annotations
import glob, re, sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"dk": REPO / "genetic_analysis/work/results_70tab", "hcp": REPO / "genetic_analysis/work/results_70tab_hcp"}
METHODS = ["CT", "PRSCS", "SBayesR", "SBayesRC"]
MATCHED = {"SCZ_pooled": "full", "MDD_pooled": "full", "SCZ_eur": "EUR", "MDD_eur": "EUR", "ASD": "EUR", "ALZ": "EUR",
           "ALZ_noAPOE": "EUR", "ALZ_IGAP": "EUR", "ALZ_IGAP_noAPOE": "EUR", "EA": "EUR"}
PHENOS = ["global_slope", "baseline_thickness"]

rows = []
for atlas, root in ROOTS.items():
    base = pd.read_csv(root / "prs_final/table_main.tsv", sep="\t")
    for f in sorted(glob.glob(str(root / "prs_final_1lmm/assoc_*.tsv"))):
        m = re.match(rf".*/assoc_({'|'.join(METHODS)})_({'|'.join(MATCHED)})(_zanc)?\.tsv$", f)
        if not m:
            continue
        meth, arm, z = m.group(1), m.group(2), ("zanc" if m.group(3) else "raw")
        st = MATCHED[arm]
        if z != ("zanc" if st == "full" else "raw"):
            continue
        t = pd.read_csv(f, sep="\t")
        for ph in PHENOS:
            one = t[(t.phenotype == f"{ph}_1lmm") & (t.stratum == st)]
            b = base[(base.trait_arm == arm) & (base.method == meth) & (base.score == z) & (base.target_stratum == st) & (base.phenotype == ph)]
            if one.empty or b.empty:
                continue
            one = one.sort_values("p").iloc[0]; b = b.sort_values("p").iloc[0]
            rows.append(dict(atlas=atlas, phenotype=ph, trait_arm=arm, stratum=st, method=meth, score=z,
                             beta_perregion=b.beta, se_perregion=b.se, p_adj_perregion=b.p_adj,
                             beta_1lmm=one.beta, se_1lmm=one.se, p_adj_1lmm=one.p_adj, n=int(one.n)))
out = pd.DataFrame(rows)
if out.empty:
    sys.exit("no rows")
order = list(MATCHED)
out["_a"] = out.trait_arm.map({a: i for i, a in enumerate(order)}); out["_m"] = out.method.map({m: i for i, m in enumerate(METHODS)})
out = out.sort_values(["atlas", "phenotype", "_a", "_m"]).drop(columns=["_a", "_m"])
for atlas, root in ROOTS.items():
    out[out.atlas == atlas].to_csv(root / "prs_final_1lmm/table_order_of_operations_all.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
f = out.copy()
for c in ["beta_perregion", "se_perregion", "beta_1lmm", "se_1lmm"]: f[c] = f[c].round(4)
for c in ["p_adj_perregion", "p_adj_1lmm"]: f[c] = f[c].map(lambda v: f"{v:.3g}")
print(f.drop(columns=["score", "n", "stratum"]).to_string(index=False))
