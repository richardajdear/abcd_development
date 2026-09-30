"""Step 1 (CSD3, ~10 min): assign every genotyped SNP to at most one gene.

    $D7PY code/01_assign_snps.py --bim $GENO.bim --gene-loc $GENELOC --out $WORK

Windows are MAGMA's window=35,10 (35 kb upstream, 10 kb downstream, strand-
aware), NCBI37.3 / GRCh37, the annotation of steps 7, 10 and 14 and of README_HPC
8.3 C3-D.  MAGMA lets a SNP belong to several genes; a partition must not, or
the set scores would double-count and no longer add up to the full score.  So a
SNP inside several windows goes to ONE gene: the one whose body is nearest
(distance 0 inside the body), ties to the lower Entrez id.  SNPs in no window
are INTERGENIC.  chr6:25,000,000-34,000,000 (the extended MHC) is its own bucket
and is excluded from every gene and set (C3-D), because it would dominate the
SCZ score; the full score is reported with and without it.

Outputs (gitignored, SNP-level):
    snp2gene.tsv.gz   snp, chr, bp, col   (col = gene column; -1 intergenic, -2 MHC)
    genes_index.tsv   col, entrez, symbol, chr, start, stop, n_snp_bim
Bucket columns are appended by step 2 after the gene columns.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

UP, DOWN = 35_000, 10_000
MHC = (6, 25_000_000, 34_000_000)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bim", required=True)
    ap.add_argument("--gene-loc", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    bim = pd.read_csv(a.bim, sep=r"\s+", header=None, usecols=[0, 1, 3],
                      names=["chr", "snp", "bp"], dtype={"chr": str})
    bim = bim[bim.chr.isin([str(i) for i in range(1, 23)])].copy()
    bim["chr"] = bim.chr.astype(int)
    assert bim.snp.is_unique, "duplicate SNP ids in bim"
    gl = pd.read_csv(a.gene_loc, sep=r"\s+", header=None,
                     names=["entrez", "chr", "start", "stop", "strand", "symbol"], dtype={"chr": str})
    gl = gl[gl.chr.isin([str(i) for i in range(1, 23)])].copy()
    gl["chr"] = gl.chr.astype(int)
    gl = gl.sort_values(["chr", "start", "entrez"]).reset_index(drop=True)
    gl["w0"] = np.where(gl.strand == "+", gl.start - UP, gl.start - DOWN)
    gl["w1"] = np.where(gl.strand == "+", gl.stop + DOWN, gl.stop + UP)

    col = np.full(len(bim), -1, dtype=np.int64)
    for ch, gg in gl.groupby("chr"):
        m = (bim.chr == ch).values
        idx = np.flatnonzero(m)
        order = np.argsort(bim.bp.values[idx], kind="stable")
        idx = idx[order]
        bp = bim.bp.values[idx]
        best_d = np.full(len(idx), np.iinfo(np.int64).max)
        best_c = np.full(len(idx), -1, dtype=np.int64)
        best_e = np.full(len(idx), np.iinfo(np.int64).max)
        for gi, r in gg.iterrows():
            lo, hi = np.searchsorted(bp, r.w0, "left"), np.searchsorted(bp, r.w1, "right")
            if hi <= lo:
                continue
            x = bp[lo:hi]
            dist = np.where(x < r.start, r.start - x, np.where(x > r.stop, x - r.stop, 0))
            better = (dist < best_d[lo:hi]) | ((dist == best_d[lo:hi]) & (r.entrez < best_e[lo:hi]))
            best_d[lo:hi] = np.where(better, dist, best_d[lo:hi])
            best_c[lo:hi] = np.where(better, gi, best_c[lo:hi])
            best_e[lo:hi] = np.where(better, r.entrez, best_e[lo:hi])
        col[idx] = best_c
    mhc = ((bim.chr == MHC[0]) & (bim.bp >= MHC[1]) & (bim.bp <= MHC[2])).values
    col[mhc] = -2

    bim["col"] = col
    bim[["snp", "chr", "bp", "col"]].to_csv(out / "snp2gene.tsv.gz", sep="\t", index=False)
    n_snp = pd.Series(col[col >= 0]).value_counts().reindex(range(len(gl)), fill_value=0)
    gl.insert(0, "col", np.arange(len(gl)))
    gl["n_snp_bim"] = n_snp.values
    gl[["col", "entrez", "symbol", "chr", "start", "stop", "n_snp_bim"]].to_csv(
        out / "genes_index.tsv", sep="\t", index=False)
    print(f"SNPs {len(bim):,}: genic {np.sum(col >= 0):,}, intergenic {np.sum(col == -1):,}, "
          f"MHC {np.sum(col == -2):,}; genes {len(gl):,}, with >=1 SNP {int((n_snp > 0).sum()):,}")


if __name__ == "__main__":
    main()
