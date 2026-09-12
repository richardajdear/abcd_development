"""Z-score each phenotype column in phenotypes_gcta.txt.

Site glue, but it fixes a real blocker rather than a path mismatch.

`03_gwas` died on the second phenotype with

    Error: the Vp is below 1e-5. Please check: 1. Is there a scaling issue with
    the phenotype? ...

fastGWA refuses any phenotype whose phenotypic variance is under 1e-5, and
`global_slope` — a per-subject rate of cortical thinning in mm/year — has
Var = 1.06e-6.  It is not a degenerate phenotype; it is a well-behaved one
measured in units that make its variance tiny.  The five exported phenotypes
span eight orders of magnitude in variance:

    baseline_thickness      3.77e-03
    global_slope            1.06e-06     <-- below fastGWA's floor
    slope_PC3               2.35e+02
    slope_PC2               2.69e+02
    slope_PC1               2.60e+02

GCTA-REML has no such floor, which is why `02_reml` produced estimates for
`global_slope` while `03_gwas` could not run it at all.

Standardising is safe, not a fudge.  A linear model's p-values are invariant
under an affine transform of the outcome: BETA and SE rescale by the same
factor, so t, chi-square and P are unchanged.  h2 is a variance ratio and is
likewise unchanged.  What does change is the unit of BETA, which becomes "per
SD of the phenotype" — arguably the more useful unit anyway, since it makes
effect sizes comparable across five phenotypes that otherwise are not.

Mean and SD are computed over the non-missing values actually present in the
file, i.e. after the genotype intersection, so they describe the analysis
sample rather than the full export.

    python standardise_pheno.py <pheno-dir>

Writes phenotypes_gcta.txt in place, keeping a .prescale backup, and prints the
scaling factors so a BETA can be converted back to native units.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd

FLOOR = 1e-5   # fastGWA's Vp guard


def main(pheno_dir: str) -> int:
    p = Path(pheno_dir) / "phenotypes_gcta.txt"
    backup = p.with_suffix(".txt.prescale")
    if not backup.exists():
        shutil.copy2(p, backup)
    df = pd.read_csv(backup, sep=" ")

    cols = [c for c in df.columns if c not in ("FID", "IID")]
    print(f"{p}: {len(df)} rows, {len(cols)} phenotypes")
    print(f"{'phenotype':<22}{'mean':>14}{'sd':>14}{'var':>14}  fastGWA")
    for c in cols:
        m, s = df[c].mean(), df[c].std()
        v = s * s
        flag = "would FAIL" if v < FLOOR else "ok"
        print(f"{c:<22}{m:>14.6g}{s:>14.6g}{v:>14.6g}  {flag}")
        if s > 0:
            df[c] = (df[c] - m) / s
        else:
            print(f"  WARNING: {c} has zero variance; left unscaled")

    df.to_csv(p, sep=" ", index=False, na_rep="NA")
    print(f"\nwrote {p} (all phenotypes now mean 0, SD 1)")
    print(f"backup of the native-unit values: {backup}")
    print("To convert a standardised BETA back to native units, multiply by the "
          "SD above.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
