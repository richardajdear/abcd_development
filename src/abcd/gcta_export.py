"""Export subject-level phenotypes and covariates in GCTA/PLINK format.

Run as ``python -m abcd.gcta_export <run-dir>``.  Writes three files that the
scripts in ``hpc/`` expect:

* ``phenotypes_gcta.txt``   FID IID <phenotype columns>
* ``covar_quant.txt``       FID IID <continuous covariates>
* ``covar_categorical.txt`` FID IID <discrete covariates>

Two things here are easy to get wrong and expensive to discover on the cluster.

**Subject IDs.** ABCD writes ``NDAR_INV...`` in the genetics tables and
``sub-NDARINV...`` in the imaging tables.  A join without normalisation returns
zero rows silently.  The genotype ``.fam`` follows the genetics convention, so
that is the target form.

**Family IDs.** GCTA's ``--grm-cutoff`` and fastGWA use the GRM, not FID, to
handle relatedness, so FID is bookkeeping.  But writing IID into FID (a common
shortcut) makes every subject look like their own family, which breaks any
downstream tool that *does* read FID.  We write the real family ID.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import paths


def _to_genetics_id(s: pd.Series) -> pd.Series:
    """Normalise imaging-form subject ids to the genetics/.fam form."""
    out = s.astype(str).str.replace(r"^sub-", "", regex=True)
    # sub-NDARINVABC123 -> NDAR_INVABC123
    return out.str.replace(r"^NDARINV", "NDAR_INV", regex=True)


def build(run_dir: str | Path) -> dict[str, pd.DataFrame]:
    run_dir = Path(run_dir)
    ph = pd.read_parquet(run_dir / "phenotypes" / "phenotypes.parquet")

    # phenotypes.parquet is long: subject x label x phenotype.  GCTA needs one
    # column per phenotype, so pivot on the whole-cortex mean.  Regional GWAS
    # is a separate (much larger) job; see hpc/README.md.
    if "label" in ph.columns:
        glob = ph[ph.label.isin(["mean", "global", "whole_cortex"])]
        if glob.empty:
            # No explicit global row: average regions per subject.
            glob = (ph.groupby(["subject", "phenotype"], observed=True)
                      .value.mean().reset_index())
        ph = glob
    wide = ph.pivot_table(index="subject", columns="phenotype",
                          values="value", observed=True)

    meta_path = run_dir / "model_table.parquet"
    cols = ["subject", "family_id", "site", "sex", "age_c"]
    meta = pd.read_parquet(meta_path, columns=[c for c in cols if c])
    meta = (meta.sort_values("age_c").groupby("subject", observed=True).first()
                .reset_index())

    df = wide.reset_index().merge(meta, on="subject", how="inner")
    df["IID"] = _to_genetics_id(df["subject"])
    df["FID"] = df["family_id"].astype(str)

    pheno_cols = [c for c in wide.columns if c in ("intercept", "slope")]
    if not pheno_cols:
        raise ValueError(f"no intercept/slope columns in {run_dir}; got {list(wide.columns)}")

    phen = df[["FID", "IID"] + pheno_cols]
    qcov = df[["FID", "IID", "age_c"]].rename(columns={"age_c": "baseline_age"})
    ccov = df[["FID", "IID", "sex", "site"]]
    return {"phenotypes_gcta": phen, "covar_quant": qcov, "covar_categorical": ccov}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_dir")
    ap.add_argument("--out-dir", default=None,
                    help="default: <run-dir>/gcta_inputs")
    a = ap.parse_args(argv)
    out = Path(a.out_dir) if a.out_dir else Path(a.run_dir) / "gcta_inputs"
    out.mkdir(parents=True, exist_ok=True)
    for name, frame in build(a.run_dir).items():
        p = out / f"{name}.txt"
        frame.to_csv(p, sep=" ", index=False, na_rep="NA")
        print(f"{p}  {frame.shape[0]} rows x {frame.shape[1]} cols")
    print("\nNote: PCs of ancestry are NOT included. Add them to "
          "covar_quant.txt before running 02_reml/03_gwas -- ABCD is "
          "multi-ancestry and fastGWA does not absorb structure for you.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
