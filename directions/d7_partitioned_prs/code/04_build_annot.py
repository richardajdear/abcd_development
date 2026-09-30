"""Step 4a (CSD3, secondary): SBayesRC annotation file with the D7 primary sets.

    $D7PY code/04_build_annot.py --annot $ANNOT --bim $GENO.bim --gene-loc $GENELOC \
        --out $WORK/annot_d7.txt

Appends one binary column per primary set (K1, K2 x2, K3 x2) to the baseline
2.2 annotation that every Figure-1 SBayesRC fit used.  A SNP is in a set when it
lies in the MAGMA 35/10 kb window of ANY set gene (window overlap, the usual
annotation convention; the unique assignment of step 1 is only needed for
additive partitions).  SNPs of the annotation file absent from the genotype bim
have no position here and get 0 (counted in the log).

Why refit: step 3 partitions the posterior of a model whose prior knew nothing
about these sets.  Refitting with the sets as annotations lets SBayesRC learn
per-set enrichment from the disorder GWAS (n ~ 1e5), which (i) is the
disorder-side single-cell -> genetics test at full discovery power and (ii)
gives weights that borrow that enrichment; step 3 is then re-run on them.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
UP, DOWN = 35_000, 10_000
PRIMARY = ["K1_brain_expressed", "K2_neuronal", "K2_glial", "K3_syngo", "K3_oligodendrocyte"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot", required=True)
    ap.add_argument("--bim", required=True)
    ap.add_argument("--gene-loc", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    an = pd.read_csv(a.annot, sep=r"\s+")
    snpcol = an.columns[0]
    bim = pd.read_csv(a.bim, sep=r"\s+", header=None, usecols=[0, 1, 3], names=["chr", "snp", "bp"],
                      dtype={"chr": str}).drop_duplicates("snp").set_index("snp")
    pos = bim.reindex(an[snpcol])
    gl = pd.read_csv(a.gene_loc, sep=r"\s+", header=None,
                     names=["entrez", "chr", "start", "stop", "strand", "symbol"], dtype={"chr": str})
    gl["w0"] = np.where(gl.strand == "+", gl.start - UP, gl.start - DOWN)
    gl["w1"] = np.where(gl.strand == "+", gl.stop + DOWN, gl.stop + UP)
    sets = pd.read_csv(HERE / "gene_sets" / "d7_gene_sets.tsv", sep="\t")
    chr_ = pos.chr.values
    bp = pos.bp.values
    print(f"annotation SNPs {len(an):,}; without a bim position {int(np.isnan(bp).sum()):,}")
    for s in PRIMARY:
        g = gl[gl.entrez.isin(sets[sets.set == s].entrez)]
        flag = np.zeros(len(an), dtype=np.int8)
        for ch, gg in g.groupby("chr"):
            m = np.flatnonzero(chr_ == ch)
            o = m[np.argsort(bp[m])]
            b = bp[o]
            for w0, w1 in zip(gg.w0.values, gg.w1.values):
                lo, hi = np.searchsorted(b, w0, "left"), np.searchsorted(b, w1, "right")
                flag[o[lo:hi]] = 1
        an[f"d7_{s}"] = flag
        print(f"  {s}: {len(g)} genes, {flag.sum():,} SNPs ({flag.mean():.3%})")
    an.to_csv(a.out, sep=" ", index=False)


if __name__ == "__main__":
    main()
