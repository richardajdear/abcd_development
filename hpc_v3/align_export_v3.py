"""Align the hpc_v3 GCTA export to the v1/v2 site convention, in place.

``make_phenotypes_v3.py`` writes the export straight from the analysis code,
which is correct but not what the cluster consumes.  Three differences were
found by comparing its output against ``hpc/work/pheno_allanc_prs/`` (the
$PHENO every v2 step actually ran on), and each one is silent:

1. **IID spelling.**  The export writes ``NDAR_INV005V6D2C``; every ``.fam``,
   ``.grm.id`` and ``.profile`` on CSD3 spells the same subject
   ``sub-005V6D2C``.  GENESIS normalises to the 8-char token and would not
   care, but GCTA (step 05) matches FID+IID literally and would silently
   return a zero-subject REML.

2. **Scale.**  v1 z-scored every phenotype column
   (``hpc/work/standardise_pheno.py``) because ``global_slope`` has
   Var = 1.06e-06 in native units.  Betas are therefore per phenotype SD
   throughout the project, and the v3 columns must use the same unit or the
   paired Delta-beta against ``global_slope`` compares two different scales.

3. **Sample.**  The export carries all 8,192 phenotyped subjects; the v2
   analysis sample is the 8,082 that are also genotyped.  Standardising over
   8,192 would put ``global_slope`` on a slightly different scale than v2's
   published number, so we subset first and standardise after -- which
   reproduces v2's five columns to within rounding.

The ancestry PCs are also missing from the fresh export (the release adapter's
``genetic_pcs`` lookup raises and ``gcta_export.build`` swallows it into
``n_pcs = 0``).  The covariate files are otherwise value-identical to v2's for
all 8,082 shared subjects (verified: baseline_age, n_visits, sex, site all
agree exactly), so we copy v2's covariate files rather than re-deriving them --
same covariates, same PCs, same order as every v2 result.

    python hpc_v3/align_export_v3.py <src-gcta_inputs_v3> <v2-pheno-dir> <out-dir>
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd

SETTLED = ["baseline_thickness", "global_slope", "slope_PC3", "slope_PC2",
           "slope_PC1"]


def token(s: pd.Series) -> pd.Series:
    """The 8-char NDAR token every ABCD id spelling ends with."""
    t = s.astype(str).str.extract(r"([A-Z0-9]{8})$")[0]
    assert t.notna().all(), "some ids do not end in an 8-char token"
    return t


def main(src: str, v2dir: str, out: str) -> int:
    src, v2dir, out = Path(src), Path(v2dir), Path(out)
    out.mkdir(parents=True, exist_ok=True)

    phen = pd.read_csv(src / "phenotypes_gcta.txt", sep=" ")
    ref = pd.read_csv(v2dir / "phenotypes_gcta.txt", sep=" ")
    phen["_t"] = token(phen.IID)
    ref["_t"] = token(ref.IID)

    # 1 + 3: adopt the reference's FID/IID spelling and its subject set, in its
    # row order, so the file is a column-wise superset of v2's $PHENO.
    keep = ref[["FID", "IID", "_t"]].merge(
        phen.drop(columns=["FID", "IID"]), on="_t", how="inner")
    assert len(keep) == len(ref), f"{len(keep)} of {len(ref)} v2 subjects matched"
    assert (keep.FID.astype(str).values == ref.FID.astype(str).values).all()

    cols = [c for c in phen.columns if c not in ("FID", "IID", "_t")]
    keep = keep[["FID", "IID"] + cols]

    # 2: z-score over the analysis sample, and print the factor that converts a
    # standardised beta back to native units.
    print(f"{'phenotype':<22}{'mean':>14}{'sd (native)':>16}")
    for c in cols:
        m, s = keep[c].mean(), keep[c].std()
        assert s > 0, f"{c} has zero variance"
        print(f"{c:<22}{m:>14.6g}{s:>16.6g}")
        keep[c] = (keep[c] - m) / s

    # the five settled columns must now reproduce v2's file
    md = max(float((keep[c].values - ref[c].values).__abs__().max())
             for c in SETTLED)
    print(f"\nmax |v3 - v2| over the five settled z-scored columns: {md:.3g}")
    assert md < 1e-6, "settled phenotypes do not reproduce v2's export"

    keep.to_csv(out / "phenotypes_gcta.txt", sep=" ", index=False, na_rep="NA")
    for f in ("covar_quant.txt", "covar_categorical.txt"):
        shutil.copy2(v2dir / f, out / f)          # identical values + the PCs
    for f in ("phenotype_manifest.tsv", "phenotype_manifest_new_only.tsv"):
        shutil.copy2(src / f, out / f)
    print(f"wrote {out}: {keep.shape[0]} rows x {keep.shape[1]} cols")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(*sys.argv[1:]))
