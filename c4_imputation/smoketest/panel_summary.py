"""Reference-panel C4 structure frequencies, GREx distribution, and the power table.

  python c4_imputation/smoketest/panel_summary.py

Writes results/reference_panel_c4.tsv (structure frequencies and per-person
GREx moments among the 111 CEU panel individuals) and results/power.tsv (the
minimum detectable effect for the primary test at the EUR-arm n, given the
leave-fold-out imputation accuracy of C4A GREx in results/smoketest_accuracy.tsv).
"""
from __future__ import annotations

import gzip
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

C4 = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("c4grex", C4 / "03_c4_grex.py")
c4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c4)  # noqa: E702
PANEL = C4 / "resources/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37.vcf.gz"

with gzip.open(PANEL, "rt") as fh:
    row = next(l for l in fh if not l.startswith("#") and l.split("\t")[2] == "C4").rstrip("\n").split("\t")
alleles = ["REF"] + [x.strip("<>") for x in row[4].split(",")]
people = []
for g in row[9:]:
    h = [alleles[int(x)] for x in re.split(r"[|/]", g)]
    e = {s: c4.seg_counts(h[0])[s] + c4.seg_counts(h[1])[s] for s in c4.SEGMENTS}
    people.append({**{g_: sum(w * e[s] for s, w in wt.items()) for g_, wt in c4.GREX_WEIGHTS.items()},
                   "C4A_copies": e["AL"] + e["AS"], "structs": [c4.structure(x) for x in h]})
P = pd.DataFrame(people)
freq = pd.Series([s for ss in P.structs for s in ss]).value_counts()
rows = [dict(metric="structure_count", key=k, value=v) for k, v in freq.items()]
rows += [dict(metric="structure_freq", key=k, value=v / freq.sum()) for k, v in freq.items()]
rows += [dict(metric="n_haplotypes", key="", value=int(freq.sum())),
         dict(metric="frac_people_common5", key="", value=np.mean([set(s) <= c4.COMMON5 for s in P.structs]))]
for c in ("C4A_GREx", "C4B_GREx", "C4A_copies"):
    rows += [dict(metric=f"{c}_mean", key="", value=P[c].mean()), dict(metric=f"{c}_sd", key="", value=P[c].std())]
rows.append(dict(metric="corr_C4A_C4B_GREx", key="", value=P.C4A_GREx.corr(P.C4B_GREx)))
pd.DataFrame(rows).to_csv(C4 / "results/reference_panel_c4.tsv", sep="\t", index=False, float_format="%.4f")

# ---- power: MDE per SD of TRUE C4A GREx, EUR arm -----------------------------
acc = pd.read_csv(C4 / "results/smoketest_accuracy.tsv", sep="\t").set_index("label")
z = norm.ppf(0.975) + norm.ppf(0.80)
prow = []
for n in (4308, 3730):                       # full EUR anchor set; Hernandez-style common-5 subset size
    for lab in ("leave_fold_out_thin3", "leave_fold_out_thin1"):
        r = acc.loc[lab, "r_C4A_GREx"]
        prow.append(dict(n=n, accuracy_source=lab, r_imputed_true=r,
                         mde80_per_sd_true_grex=z / (r * np.sqrt(n)), mde80_if_perfect=z / np.sqrt(n)))
pd.DataFrame(prow).to_csv(C4 / "results/power.tsv", sep="\t", index=False, float_format="%.4f")
print(freq.head(8).to_string()); print(pd.DataFrame(prow).round(4).to_string(index=False))
