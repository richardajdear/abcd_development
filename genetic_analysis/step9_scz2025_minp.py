#!/usr/bin/env python
"""Family-block min-p permutation for the 2025-SCZ C+T cells (rule 6).

Same statistic and module as step6_minp_permutation.py; only the score
directories differ.  Cells (ancestry-matched, rule 4):
    CT/SCZ25_EUR      -> EUR      (primary, EUR arm)
    CT/SCZ25_META     -> full     (primary, pooled arm)
    CT/SCZ25_META     -> EUR      (secondary: multi-ancestry GWAS on the EUR
                                   target is under-powered, not confounded)
    CT/SCZ25_MATCHED  -> full     (ancestry-matched composite; already z within
                                   cluster, so read against the _zanc cells)
Covariates are the module's primary arm (sex + 10 PCs, no age), as in step 6.

    python genetic_analysis/step9_scz2025_minp.py --parc dsk|hcp [--n-perm 2000]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "genetic_analysis" / "setup"))
import prs_paired_delta as m  # noqa: E402

PHENOS = ["global_slope", "baseline_thickness", "slope_PC1", "slope_PC2", "slope_PC3"]
CELLS = [("SCZ25_EUR", "EUR", "primary"), ("SCZ25_META", "full", "primary"),
         ("SCZ25_META", "EUR", "secondary"), ("SCZ25_MATCHED", "full", "primary_matched")]


def token_in(d: Path) -> str:
    toks = {re.match(r"score_(.+)_[^_]+\.profile$", p.name).group(1) for p in d.glob("score_*.profile")}
    if len(toks) != 1:
        raise SystemExit(f"{d}: expected one score token, found {sorted(toks)}")
    return toks.pop()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parc", choices=["dsk", "hcp"], required=True)
    ap.add_argument("--scores-root", default=str(REPO / "genetic_analysis/work/scores_scz2025"))
    ap.add_argument("--eur-ids", default=str(REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep"))
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260908 + 1)
    a = ap.parse_args()
    suf = "" if a.parc == "dsk" else "_hcp"
    pheno = REPO / f"genetic_analysis/work/pheno_70tab{suf}"
    out = REPO / f"genetic_analysis/work/results_70tab{suf}/prs_scz2025/table_minp_permutation.tsv"
    parts = []
    for arm, stratum, design in CELLS:
        d = Path(a.scores_root) / "CT" / arm
        tok = token_in(d)
        m.DISORDERS = [tok]
        m.NEW = [p for p in PHENOS if p != m.REF]
        dd, cov, scores = m.load(pheno, d, Path(a.eur_ids), with_age=False)
        if not scores:
            raise SystemExit(f"no score_{tok}_<thr>.profile under {d}")
        res = m.minp_permutation([m.Stratum(dd, stratum, cov)], scores, a.n_perm, a.seed)
        res.insert(0, "design", design); res.insert(1, "trait_arm", arm)
        res["matched"] = "yes" if design != "secondary" else "no (multi->EUR, not confounded)"
        print(f"  {design:<16} {arm:<14} -> {stratum:<5} {len(scores)} thresholds -> {len(res)} rows", flush=True)
        parts.append(res)
    res = pd.concat(parts, ignore_index=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out} ({len(res)} rows)")
    show = res[res.phenotype == "global_slope"]
    print(show[["design", "trait_arm", "stratum", "n", "min_p", "p_perm", "p_bonferroni"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
