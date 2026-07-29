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
import yaml

from . import paths


def _to_genetics_id(s: pd.Series) -> pd.Series:
    """Normalise imaging-form subject ids to the genetics/.fam form."""
    out = s.astype(str).str.replace(r"^sub-", "", regex=True)
    # sub-NDARINVABC123 -> NDAR_INVABC123
    return out.str.replace(r"^NDARINV", "NDAR_INV", regex=True)


class FamilyEffectConflict(ValueError):
    """Raised when a run's BLUPs cannot support a genetic analysis.

    Fitting ``(1 | family_id)`` partitions the between-family variance into its
    own random effect.  The subject-level BLUPs that remain are *within-family*
    deviations, and the between-family component -- which is exactly what
    genetic relatedness explains -- has been removed from them.

    Measured on release 7.0, whole-cortex thickness, 8,192 subjects (624
    same-sex twin pairs, 578 sibling pairs):

        with (1 | family_id):   r_sib = -0.117   <-- impossible; forced negative
        without:                r_sib = +0.340

    Cortical thickness is strongly familial, so a *negative* sibling
    correlation is not a weak result, it is a structurally impossible one: the
    family random effect centres each family at zero, so with two members per
    family their deviations must sum to zero and anti-correlate by
    construction.  Falconer's estimator differences two correlations and so
    partly cancels the artefact -- it returned a plausible-looking h2 = 0.485
    for baseline thickness from those invalid correlations -- which is why this
    guard checks the model specification rather than trusting the h2 value to
    look wrong.

    Use a config with ``family_effect: false`` for genetic work (see
    ``configs/ct_70_genetic.yaml``) and handle relatedness where it belongs:
    in the GRM for GCTA, or via a mixed-model GWAS.
    """


def _check_genetic_validity(run_dir: Path) -> None:
    """Refuse to export BLUPs that had the between-family variance removed."""
    cfg_path = run_dir / "config.yaml"
    if not cfg_path.exists():          # older runs predate the config dump
        return
    import yaml
    cfg = yaml.safe_load(cfg_path.read_text()) or {}
    if cfg.get("family_effect"):
        raise FamilyEffectConflict(
            f"{run_dir.name} was fitted with family_effect: true, so its BLUPs are "
            "within-family deviations and carry no between-family genetic variance. "
            "Re-fit with family_effect: false (configs/ct_70_genetic.yaml) before "
            "exporting for GCTA/GWAS. See FamilyEffectConflict.__doc__."
        )


def build(run_dir: str | Path) -> dict[str, pd.DataFrame]:
    run_dir = Path(run_dir)
    _check_genetic_validity(run_dir)
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

    # Ancestry PCs must be quantitative covariates: ABCD is multi-ancestry, and
    # neither GCTA-REML nor fastGWA absorbs population structure on its own.
    # 7.0 tabulates 32 PCs in the static table; 5.1's local copy has none, in
    # which case we emit the file without them and say so loudly.
    n_pcs = 0
    try:
        from . import io as _io
        release = (yaml.safe_load((run_dir / "config.yaml").read_text()) or {}).get("release")
        pcs = _io.get_adapter(release).genetic_pcs(n=10)
        pcs = pcs.copy()
        pcs["IID"] = _to_genetics_id(pcs["subject"] if "subject" in pcs else pcs.index.to_series())
        pc_cols = [c for c in pcs.columns if c.lower().startswith("pc")]
        qcov = qcov.merge(pcs[["IID"] + pc_cols], on="IID", how="left")
        n_pcs = len(pc_cols)
    except Exception as exc:                       # noqa: BLE001 - release-dependent
        qcov.attrs["pc_warning"] = str(exc)

    qcov.attrs["n_pcs"] = n_pcs
    return {"phenotypes_gcta": phen, "covar_quant": qcov, "covar_categorical": ccov}


def main(argv: list[str] | None = None) -> int:
    from .config import active_run_dir

    ap = argparse.ArgumentParser(
        description=__doc__,
        epilog="run_dir defaults to the run for $ABCD_CONFIG.",
    )
    ap.add_argument("run_dir", nargs="?", default=None,
                    help="run directory; omit to use $ABCD_CONFIG")
    ap.add_argument("--config", default=None,
                    help="config name, overriding $ABCD_CONFIG")
    ap.add_argument("--out-dir", default=None,
                    help="default: <run-dir>/gcta_inputs")
    a = ap.parse_args(argv)
    run_dir = Path(a.run_dir) if a.run_dir else active_run_dir(a.config)
    out = Path(a.out_dir) if a.out_dir else run_dir / "gcta_inputs"
    out.mkdir(parents=True, exist_ok=True)
    built = build(run_dir)
    for name, frame in built.items():
        p = out / f"{name}.txt"
        frame.to_csv(p, sep=" ", index=False, na_rep="NA")
        print(f"{p}  {frame.shape[0]} rows x {frame.shape[1]} cols")

    qcov = built["covar_quant"]
    n_pcs = qcov.attrs.get("n_pcs", 0)
    if n_pcs:
        miss = qcov[[c for c in qcov.columns if c.lower().startswith("pc")]].isna().any(axis=1).sum()
        print(f"\ncovar_quant.txt includes {n_pcs} ancestry PCs "
              f"({miss} subjects missing PCs -> NA; GCTA drops them).")
    else:
        print("\nWARNING: no ancestry PCs in covar_quant.txt "
              f"({qcov.attrs.get('pc_warning', 'unavailable for this release')}). "
              "ABCD is multi-ancestry and neither GCTA-REML nor fastGWA absorbs "
              "population structure -- add PCs before 02_reml/03_gwas.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
