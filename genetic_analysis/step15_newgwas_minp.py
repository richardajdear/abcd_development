#!/usr/bin/env python
"""Step 15 min-p permutation (rule 6) for the new C+T cells, matched strata,
both atlases, reusing step6_minp_permutation.run() unchanged.
    python genetic_analysis/step15_newgwas_minp.py [--n-perm 2000]"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "genetic_analysis"))
import step6_minp_permutation as s6  # noqa: E402
CELLS = [("BIP_eur", "EUR"), ("BIP_pooled", "full"), ("ADHD", "EUR"), ("INT", "EUR")]
s6.MATCHED.update(dict(CELLS))   # so the matched flag is "yes" for these arms
ap = argparse.ArgumentParser(); ap.add_argument("--n-perm", type=int, default=2000); a = ap.parse_args()
eur = REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep"
for suf in ("", "_hcp"):
    pheno = REPO / f"genetic_analysis/work/pheno_70tab{suf}"
    parts = []
    for arm, st in CELLS:
        d = REPO / "genetic_analysis/work/scores_newgwas/CT" / arm
        parts.append(s6.run(pheno, eur, d, s6.token_in(d), [st], a.n_perm, 20260908 + 1, "matched", arm))
    res = pd.concat(parts, ignore_index=True)
    out = REPO / f"genetic_analysis/work/results_70tab{suf}/prs_newgwas/table_minp_permutation.tsv"
    res.to_csv(out, sep="\t", index=False)
    print(out); print(res[res.phenotype == "global_slope"][["trait_arm", "stratum", "n", "min_p", "p_perm", "p_bonferroni"]].to_string(index=False))
