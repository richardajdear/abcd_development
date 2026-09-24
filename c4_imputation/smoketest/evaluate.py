"""Score imputed C4 calls against the held-out truth, pooled over folds.

  python evaluate.py --pairs out/fold1_thin3 out/fold2_thin3 ... --label thin3 --out ../results/smoketest_accuracy.tsv

Each prefix needs <prefix>.truth.tsv (make_test_data.py) and <prefix>.calls.tsv
(03_c4_grex.py).  Appends one row per label to --out.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", nargs="+", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = pd.concat([pd.read_csv(f"{p}.calls.tsv", sep="\t").merge(
        pd.read_csv(f"{p}.truth.tsv", sep="\t"), on="IID", validate="1:1") for p in a.pairs])

    def hap_matches(r):  # matched structures between the two unordered pairs (0, 1 or 2)
        c = Counter([r.hap1, r.hap2]); t = Counter([r.true_hap1, r.true_hap2])
        return sum((c & t).values())

    d["m"] = d.apply(hap_matches, axis=1)
    hi = d.post_mean >= 0.7
    row = dict(label=a.label, n_folds=len(a.pairs), n_individuals=len(d),
               genotype_concordance=(d.m == 2).mean(), haplotype_concordance=d.m.sum() / (2 * len(d)),
               frac_post_ge_0_7=hi.mean(), genotype_concordance_post_ge_0_7=(d.m[hi] == 2).mean(),
               r_C4A_GREx=np.corrcoef(d.C4A_GREx, d.true_C4A_GREx)[0, 1],
               r_C4B_GREx=np.corrcoef(d.C4B_GREx, d.true_C4B_GREx)[0, 1],
               r_C4A_copies=np.corrcoef(d.C4A_copies, d.true_AL + d.true_AS)[0, 1],
               r_HERV_copies=np.corrcoef(d.HERV_copies, d.true_AL + d.true_BL)[0, 1],
               median_post_mean=d.post_mean.median())
    out = Path(a.out)
    prev = pd.read_csv(out, sep="\t") if out.exists() else pd.DataFrame()
    prev = prev[prev.get("label", pd.Series(dtype=str)) != a.label] if len(prev) else prev
    pd.concat([prev, pd.DataFrame([row])]).to_csv(out, sep="\t", index=False, float_format="%.4f")
    print(" ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}" for k, v in row.items()))


if __name__ == "__main__":
    main()
