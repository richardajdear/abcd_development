"""
09_summarise_enrichment.py -- one table + one figure for H2 across options,
DS levels and benchmarks.

results/enrichment_summary.tsv : per vector x disorder: MAGMA gene-property
    BETA_STD / P (continuous), top-decile BETA / P, and the length-matched
    permutation z / p for the main prioritised and pool sets.
Figures are built separately in R (code/fig1_*.R, fig2_*.R).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parent
RES, FIG = ROOT / "results", ROOT / "figures"

M = pd.read_csv(RES / "magma_gene_property.tsv", sep="\t")
P = pd.read_csv(RES / "permutation_enrichment.tsv", sep="\t")
C = pd.read_csv(RES / "magma_conditional.tsv", sep="\t")

ORDER = [("lead_opt2_PLS2_ds0", "lead_opt2_PLS2", "ds0", "ABCD PLS2 (dCT+CT), ds0"),
         ("lead_opt2_PLS2_ds25", "lead_opt2_PLS2", "ds25", "ABCD PLS2 (dCT+CT), ds25"),
         ("lead_opt2_PLS2_ds50", "lead_opt2_PLS2", "ds50", "ABCD PLS2 (dCT+CT), ds50"),
         ("opt1_dCT_PLS1_ds25", "opt1_dCT_PLS1", "ds25", "dCT alone (opt 1)"),
         ("opt3_dCT_dT1T2_PLS2_ds25", "opt3_dCT_dT1T2_PLS2", "ds25", "dCT+dT1T2 PLS2 (opt 3)"),
         ("opt4_full4_PLS2_ds25", "opt4_full4_PLS2", "ds25", "4-map PLS2 (opt 4)"),
         ("opt5_slopePCs_PLS2_ds25", "opt5_slopePCs_PLS2", "ds25", "slope-PC PLS2 (opt 5)"),
         ("control_opt2_PLS1_static_ds25", "control_opt2_PLS1_static", "ds25", "static PLS1 (control)"),
         ("AHBA_C1", "AHBA_C1", "ds0", "AHBA C1"), ("AHBA_C2", "AHBA_C2", "ds0", "AHBA C2"),
         ("AHBA_C3", "AHBA_C3", "ds0", "AHBA C3 (Dear 2024)"),
         ("NSPN_PLS1_z", "NSPN_PLS1_z", "ds0", "NSPN PLS1"), ("NSPN_PLS2_z", "NSPN_PLS2_z", "ds0", "NSPN PLS2 (Whitaker 2016)")]
SETS = ["SCZ_prioritised", "SCZ_locus_pool", "MDD_highconf", "MDD_pool"]

rows = []
for mv, pv, ds, label in ORDER:
    for dis in ("SCZ", "MDD"):
        m = M[(M.VARIABLE == mv) & (M.disorder == dis)].iloc[0]
        t = M[(M.VARIABLE == mv + "_top10") & (M.disorder == dis)].iloc[0]
        r = dict(vector=mv, label=label, disorder=dis, n_genes=int(m.NGENES),
                 magma_beta_std=m.BETA_STD, magma_se=m.SE, magma_p=m.P,
                 magma_top10_beta=t.BETA, magma_top10_p=t.P)
        for s in SETS:
            if s.startswith(dis):
                q = P[(P.vector == pv) & (P.ds == ds) & (P.gene_set == s)]
                if len(q):
                    q = q.iloc[0]
                    r[f"perm_{s}_n"] = int(q.n_set); r[f"perm_{s}_z"] = q.z_lengthmatched; r[f"perm_{s}_p2"] = q.p2_lengthmatched
        rows.append(r)
S = pd.DataFrame(rows)
S.to_csv(RES / "enrichment_summary.tsv", sep="\t", index=False, float_format="%.4g")

# Figure 2 is built in R: code/fig2_enrichment.R (ggplot2 + patchwork).
# This script writes only the summary table.

pd.set_option("display.width", 250)
print(S[["label", "disorder", "n_genes", "magma_beta_std", "magma_p", "perm_SCZ_prioritised_z", "perm_SCZ_locus_pool_z", "perm_MDD_highconf_z", "perm_MDD_pool_z"]].round(3).to_string(index=False))
