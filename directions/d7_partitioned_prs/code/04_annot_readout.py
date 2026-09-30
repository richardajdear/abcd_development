"""Step 4c (CSD3, secondary): disorder-side enrichment of the D7 sets.

    $D7PY code/04_annot_readout.py --annot $WORK/annot_d7.txt \
        --fit SCZ25_META:<baseline .snpRes>:<annotated .snpRes> --fit ... \
        --out directions/d7_partitioned_prs/results/table_d7_annot_enrichment.tsv

For every set column of the annotation and each fit, the posterior share of
genetic variance  sum_{j in set} 2 p_j (1 - p_j) b_j^2 / sum_j (same)  over the
share of SNPs, from the SBayesRC posterior means (A1Frq, A1Effect in .snpRes).
Posterior means shrink, so this is a descriptive enrichment; the comparison
that matters is baseline-prior fit vs set-annotated fit of the same GWAS, and
the set vs its matched random sets (step 3, column f_enrich, same quantity in
the target sample).  Also copy any annotation-level output GCTB writes for the
annotated fit (list the files next to the .snpRes) into results/ by hand if it
reports per-annotation heritability enrichment; record its name in the README.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd


def vshare(res: pd.DataFrame, flag: pd.Series) -> tuple[float, float]:
    p = res.A1Frq.values
    v = 2 * p * (1 - p) * res.A1Effect.values ** 2
    f = flag.reindex(res.Name).fillna(0).values.astype(bool)
    return float(v[f].sum() / v.sum()), float(f.mean())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot", required=True)
    ap.add_argument("--fit", action="append", required=True, help="ARM:baseline.snpRes:annot.snpRes")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    an = pd.read_csv(a.annot, sep=r"\s+")
    an = an.set_index(an.columns[0])
    setcols = [c for c in an.columns if c.startswith("d7_")]
    rows = []
    for spec in a.fit:
        arm, base, ann = spec.split(":")
        for label, path in (("baseline_prior", base), ("d7_annotated", ann)):
            res = pd.read_csv(path, sep=r"\s+", usecols=["Name", "A1Frq", "A1Effect"])
            for c in setcols:
                h, s = vshare(res, an[c])
                rows.append(dict(arm=arm, fit=label, set=c[3:], snp_share=s, var_share=h,
                                 enrichment=h / s if s > 0 else np.nan, n_snp=len(res)))
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False, float_format="%.5g")
    print(pd.DataFrame(rows).pivot_table(index=["arm", "set"], columns="fit", values="enrichment").round(3))


if __name__ == "__main__":
    main()
