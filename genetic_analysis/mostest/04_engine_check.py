#!/usr/bin/env python
"""Step 16.5 (gate G1): REGENIE reproduces the existing GENESIS pooled scan.

The `global` run puts global_slope_1lmm and baseline_thickness_1lmm through the
same REGENIE pipeline as the regional families.  Its z must track the GENESIS
single-LMM pooled scan (steps 3-4; <GENESIS_1LMM>/<trait>.sumstats.tsv.gz,
columns SNP CHR POS A1 A2 AF1 N BETA SE P) for the same children.  The two
engines differ (REGENIE LOCO ridge vs GENESIS sparse-kinship LMM; RINT vs
z-scored phenotype), so agreement is expected to be high, not exact.

Writes <root>/engine_check.tsv: trait, n_common, r_z (allele-aligned),
r_abs_z, lambda_regenie, lambda_genesis, pass.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R_MIN = 0.90
LAMBDA_TOL = 0.05


def lam(z):
    return float(np.median(z ** 2) / stats.chi2.ppf(0.5, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--genesis-dir", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    parts = [np.load(root / "zmat" / f"global_chr{c}.npz") for c in range(1, 23)]
    traits = list(parts[0]["traits"])
    rg = pd.concat([pd.DataFrame({"SNP": d["snp"], "a1": d["a1"], "a0": d["a0"],
                                  **{t: d["z"][:, j] for j, t in enumerate(traits)}}) for d in parts])
    rows = []
    for t in traits:
        f = Path(a.genesis_dir) / f"{t}.sumstats.tsv.gz"
        if not f.exists():
            raise SystemExit(f"FATAL: missing {f}")
        g = pd.read_csv(f, sep="\t", usecols=["SNP", "A1", "A2", "BETA", "SE"])
        g["zg"] = g.BETA / g.SE
        m = rg[["SNP", "a1", "a0", t]].merge(g, on="SNP")
        same = (m.a1 == m.A1) & (m.a0 == m.A2); flip = (m.a1 == m.A2) & (m.a0 == m.A1)
        m = m[same | flip].copy()
        m["zg_al"] = np.where(same[same | flip].to_numpy(), m.zg, -m.zg)
        r = float(np.corrcoef(m[t], m.zg_al)[0, 1])
        r_abs = float(np.corrcoef(m[t].abs(), m.zg.abs())[0, 1])
        lr, lg = lam(m[t].to_numpy()), lam(m.zg.to_numpy())
        # sign convention of GENESIS A1 is not asserted: |r| is the gate
        ok = abs(r) >= R_MIN and abs(lr - lg) <= LAMBDA_TOL
        rows.append(dict(trait=t, n_common=len(m), r_z=round(r, 4), r_abs_z=round(r_abs, 4),
                         lambda_regenie=round(lr, 4), lambda_genesis=round(lg, 4), pass_=ok))
    out = pd.DataFrame(rows).rename(columns={"pass_": "pass"})
    out.to_csv(root / "engine_check.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    if not out["pass"].all():
        raise SystemExit("GATE G1 FAILED: REGENIE does not reproduce the GENESIS scan -- stop and diagnose")


if __name__ == "__main__":
    main()
