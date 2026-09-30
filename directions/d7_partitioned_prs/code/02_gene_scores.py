"""Step 2 (CSD3, one pass over the genotypes, ~2-4 h): per-gene partial scores.

    $D7PY code/02_gene_scores.py --bfile $GENO --work $WORK --threads 16 \
        --arm 'SCZ25_META|<weights>|<profile>|pooled' --arm ...

For every score arm, writes an n x (g + 2) float32 matrix whose column j is
    sum over SNPs assigned to gene j of  (SBayesRC posterior weight x A1 dosage)
with the last two columns INTERGENIC and MHC, so every row sums to the full
genome-wide score.  This subsets the existing genome-wide SBayesRC posterior
(fitted once, jointly, with LD) -- the valid partitioned score of README_HPC
8.3 C3-D -- instead of refitting per set.  Genotypes are read once and all arms
are scored in the same pass.

GATE (hard stop): the row sums must reproduce the PLINK 1.9 `--score ... sum`
.profile SCORESUM that Figure 1 used, r > 0.9999, per arm.  PLINK's handling of
missing calls (mean imputation from the loaded sample, or none) is tried both
ways and the matching one is recorded; if neither passes the arm is not written.

Outputs (gitignored; individual-level):
    gene_scores/<ARM>.npy        float32, rows = .fam order
    gene_scores/<ARM>_cols.tsv   col, entrez, n_snp (weighted SNPs), sum_abs_w
    gene_scores/iids.txt         .fam IIDs in row order
    gene_scores/gate.tsv         arm, n_snp_weights, n_snp_matched, r_profile, imputation
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed

CHUNK = 20_000


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bfile", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--arm", action="append", required=True, help="NAME|weights|profile|cell")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    work = Path(a.work); od = work / "gene_scores"; od.mkdir(parents=True, exist_ok=True)

    s2g = pd.read_csv(work / "snp2gene.tsv.gz", sep="\t")
    genes = pd.read_csv(work / "genes_index.tsv", sep="\t")
    g = len(genes); INTER, MHC = g, g + 1
    bim = pd.read_csv(a.bfile + ".bim", sep=r"\s+", header=None,
                      names=["chr", "snp", "cm", "bp", "a1", "a2"], dtype={"chr": str})
    fam = pd.read_csv(a.bfile + ".fam", sep=r"\s+", header=None, usecols=[1], names=["iid"])
    n = len(fam)
    colmap = s2g.set_index("snp").col
    bim["col"] = colmap.reindex(bim.snp).values          # NaN = non-autosomal
    bim["col"] = np.where(bim.col == -1, INTER, np.where(bim.col == -2, MHC, bim.col))

    arms = [x.split("|") for x in a.arm]
    # W is the weight on the count of bim-A1.  A weight on bim-A2 is a "flip":
    # beta*(2 - x) = 2*beta - beta*x, so W = -beta and a constant 2*beta per
    # non-missing call is added to the same gene column.
    W = np.zeros((len(bim), len(arms)), dtype=np.float64)
    FLIP = np.zeros((len(bim), len(arms)), dtype=bool)
    gate = []
    for k, (name, wf, prof, cell) in enumerate(arms):
        w = pd.read_csv(wf, sep=r"\s+", header=None, names=["snp", "a1w", "beta"]).drop_duplicates("snp")
        m = bim[["snp", "a1", "a2"]].merge(w, on="snp", how="left")
        same = (m.a1w == m.a1).values
        flip = (m.a1w == m.a2).values & ~same
        W[:, k] = np.where(same, m.beta, np.where(flip, -m.beta, 0.0))
        FLIP[:, k] = flip
        gate.append(dict(arm=name, cell=cell, n_snp_weights=len(w), n_snp_matched=int((same | flip).sum())))

    use = np.flatnonzero(bim.col.notna().values & (W != 0).any(1))
    print(f"n={n:,}; SNPs with a weight in any arm: {len(use):,}", flush=True)
    acc = {"imp": np.zeros((len(arms), n, g + 2)), "raw": np.zeros((len(arms), n, g + 2))}
    nsnp = np.zeros((len(arms), g + 2), dtype=np.int64)
    sabs = np.zeros((len(arms), g + 2))

    bed = open_bed(a.bfile + ".bed", count_A1=True, num_threads=a.threads)
    for c0 in range(0, len(use), CHUNK):
        ix = use[c0:c0 + CHUNK]
        X = bed.read(index=np.s_[:, ix], dtype="float32")          # count of bim A1, NaN = missing
        miss = np.isnan(X)
        mu = np.nanmean(X, axis=0)
        cols = bim.col.values[ix].astype(np.int64)
        order = np.argsort(cols, kind="stable")
        cs = cols[order]
        starts = np.flatnonzero(np.r_[True, cs[1:] != cs[:-1]])
        ucol = cs[starts]
        X_imp = np.where(miss, mu, X)
        X_raw = np.where(miss, 0.0, X)
        called = (~miss).astype(np.float32)
        for k in range(len(arms)):
            wk = W[ix, k]
            const = np.where(FLIP[ix, k], -2.0 * wk, 0.0)          # = 2*beta on flipped SNPs
            for mode, Xm, cnt in (("imp", X_imp, 1.0), ("raw", X_raw, called)):
                contrib = (Xm * wk + cnt * const)[:, order]
                acc[mode][k][:, ucol] += np.add.reduceat(contrib, starts, axis=1)
            nsnp[k, ucol] += np.add.reduceat((wk[order] != 0).astype(np.int64), starts)
            sabs[k, ucol] += np.add.reduceat(np.abs(wk[order]), starts)
        if (c0 // CHUNK) % 25 == 0:
            print(f"  {c0 + len(ix):,} / {len(use):,}", flush=True)

    fam.iid.to_csv(od / "iids.txt", index=False, header=False)
    for k, (name, wf, prof, cell) in enumerate(arms):
        p = pd.read_csv(prof, sep=r"\s+")
        p = p.set_index("IID").reindex(fam.iid)
        best = None
        for mode in ("imp", "raw"):
            r = np.corrcoef(acc[mode][k].sum(1), p.SCORESUM.values)[0, 1]
            gate[k][f"r_{mode}"] = r
            if best is None or r > best[1]:
                best = (mode, r)
        gate[k]["imputation"], gate[k]["r_profile"] = best
        print(f"{name}: r(profile) imp={gate[k]['r_imp']:.6f} raw={gate[k]['r_raw']:.6f}", flush=True)
        if best[1] < 0.9999:
            print(f"  GATE FAILED for {name}: not written", flush=True)
            continue
        np.save(od / f"{name}.npy", acc[best[0]][k].astype(np.float32))
        pd.DataFrame(dict(col=np.arange(g + 2),
                          entrez=list(genes.entrez) + ["INTERGENIC", "MHC"],
                          n_snp=nsnp[k], sum_abs_w=sabs[k])).to_csv(
            od / f"{name}_cols.tsv", sep="\t", index=False)
    pd.DataFrame(gate).to_csv(od / "gate.tsv", sep="\t", index=False)
    if any(gg["r_profile"] < 0.9999 for gg in gate):
        raise SystemExit("FATAL: at least one arm failed the profile gate (see gate.tsv)")


if __name__ == "__main__":
    main()
