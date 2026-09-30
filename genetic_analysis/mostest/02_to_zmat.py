#!/usr/bin/env python
"""Step 16.3b: k per-trait REGENIE step-2 outputs of one (run, chr) -> one
z-matrix file.

Input:  <prefix>_<trait>.regenie.gz for every trait in --traits (in order).
Output: .npz with z (n_snp x k, float32, z = BETA / SE, sign w.r.t. ALLELE1),
        snp / chrom / pos / a0 / a1 / freq / n arrays.  SNP-level only; no
        per-subject content, but gitignored with the rest of work/.

Columns are resolved BY NAME (repo convention).  The variant list must be
identical across traits (it is: same genotypes, same sample, no missing
phenotypes, same --minMAC), and this is asserted rather than assumed.
MAF < 1 % is dropped here (rule: ASSOC_MAF = 0.01; --minMAC 20 already ran).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

MAF = 0.01
COLS = {"CHROM": "chrom", "GENPOS": "pos", "ID": "snp", "ALLELE0": "a0",
        "ALLELE1": "a1", "A1FREQ": "freq", "N": "n", "BETA": "beta", "SE": "se"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--chr", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    traits = Path(a.traits).read_text().split()
    ref = None
    zs = []
    for t in traits:
        f = Path(f"{a.prefix}_{t}.regenie.gz")
        if not f.exists():
            raise SystemExit(f"FATAL: missing {f}")
        d = pd.read_csv(f, sep=r"\s+", usecols=list(COLS)).rename(columns=COLS)
        if "TEST" in d.columns:
            d = d[d.TEST == "ADD"]
        if ref is None:
            ref = d[["snp", "chrom", "pos", "a0", "a1", "freq", "n"]].reset_index(drop=True)
            if (ref.chrom.astype(int) != a.chr).any():
                raise SystemExit(f"FATAL: {f} contains chromosomes other than {a.chr}")
        elif len(d) != len(ref) or not np.array_equal(d.snp.values, ref.snp.values):
            raise SystemExit(f"FATAL: variant list of {t} differs from {traits[0]}")
        zs.append((d.beta / d.se).to_numpy(np.float32))
    z = np.column_stack(zs)
    keep = (np.minimum(ref.freq, 1 - ref.freq) >= MAF).to_numpy() & np.isfinite(z).all(axis=1)
    ref = ref[keep]
    np.savez(a.out, z=z[keep], traits=np.array(traits),
             **{c: (ref[c].astype(str).to_numpy(dtype="U") if ref[c].dtype == object
                    else ref[c].to_numpy()) for c in ref.columns})   # no pickled objects
    print(f"chr{a.chr}: {keep.sum()} of {keep.size} variants kept (MAF >= {MAF}), k = {len(traits)}")


if __name__ == "__main__":
    main()
