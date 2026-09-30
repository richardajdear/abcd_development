"""Laptop step L1: export the per-child CBCL outcomes for the D7 symptom arm.

    PYTHONPATH=src:genetic_analysis python directions/d7_partitioned_prs/code/L1_export_symptom_outcomes.py

The outcomes are the Figure-1 definitions, built by executing the outcome block
of ahba_pls/code/23_cbcl_explore.py (via genetic_analysis/fig1_prep_cbcl.py), not
re-implemented:  y_LATE ~ PRS + y_BASE + age_LATE + sex + site, log1p raw CBCL
scores, LATE = mean of waves 5-7.  Sample = the imaged children that script
builds (8,716; ~8,596 genotyped), i.e. the same children as the imaging arm.

Writes directions/d7_partitioned_prs/work/d7_symptom_outcomes.tsv
(INDIVIDUAL-LEVEL, gitignored).  Copy it to the same path on CSD3:

    scp directions/d7_partitioned_prs/work/d7_symptom_outcomes.tsv \
        rajd2@login.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/directions/d7_partitioned_prs/work/
"""
from __future__ import annotations

from pathlib import Path

import fig1_prep_cbcl

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "work" / "d7_symptom_outcomes.tsv"
CBCL = ["pfactor", "internal", "external", "depress", "thought"]


def main() -> None:
    ns = fig1_prep_cbcl.build_outcomes()
    D, FU = ns["D"], ns["FU"]
    cols = ["family", "site", "sex", f"age_{FU}"]
    for o in CBCL:
        cols += [f"{o}_{FU}", f"{o}_BASE"]
    X = D[cols].copy()
    X.columns = [c.replace(f"_{FU}", "_LATE") for c in X.columns]
    X.insert(0, "token", [str(i).replace("sub-", "").replace("NDARINV", "").replace("_", "").upper()
                          for i in D.index])
    assert X.token.str.len().eq(8).all(), X.token.str.len().value_counts()
    assert X.token.is_unique
    OUT.parent.mkdir(parents=True, exist_ok=True)
    X.to_csv(OUT, sep="\t", index=False)
    print(OUT, X.shape, {o: int(X[f"{o}_LATE"].notna().sum()) for o in CBCL})


if __name__ == "__main__":
    main()
