"""Export subject-level phenotypes and covariates in GCTA/PLINK format.

Run as ``python -m abcd.gcta_export <run-dir>``.  Writes four files that the
scripts in ``hpc/`` expect:

* ``phenotypes_gcta.txt``   FID IID <phenotype columns>
* ``covar_quant.txt``       FID IID <continuous covariates>
* ``covar_categorical.txt`` FID IID <discrete covariates>
* ``phenotype_manifest.tsv``  name, mpheno column index, priority, role

The manifest exists so the sbatch scripts do not have to rediscover a
phenotype's column position by grepping the header of a space-delimited file --
a step that silently produced the wrong ``--mpheno`` whenever a phenotype name
was a prefix of another.  It is written by the same code that writes the
phenotype columns, so the two cannot disagree.

Three things here are easy to get wrong and expensive to discover on the cluster.

**Subject IDs.** ABCD writes ``NDAR_INV...`` in the genetics tables and
``sub-NDARINV...`` in the imaging tables.  A join without normalisation returns
zero rows silently.  The genotype ``.fam`` follows the genetics convention, so
that is the target form.

**Family IDs.** GCTA's ``--grm-cutoff`` and fastGWA use the GRM, not FID, to
handle relatedness, so FID is bookkeeping.  But writing IID into FID (a common
shortcut) makes every subject look like their own family, which breaks any
downstream tool that *does* read FID.  We write the real family ID.

**Which phenotypes.** The five in :data:`PHENOTYPES`, in the priority order
REPORT_7.0 section 10 derives -- not the ``intercept``/``slope`` pair an earlier
version of this module exported.  ``baseline_thickness`` is carried as a
positive control rather than as a target: it is the phenotype that exposed the
family-effect artefact (section 4.2), and a REML run that fails to recover a
high h2 for it means something is wrong with the GRM or the covariates, not
with the developmental phenotype.
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


#: The phenotypes exported for genetic analysis, in the priority order derived
#: in REPORT_7.0 section 10.  ``priority`` 0 marks the positive control, which
#: is not a target but must be run: it is the only phenotype on which a
#: specification error is *visible* (section 4.2).
#:
#: The ordering is not by heritability.  Held-out h2 is nearly flat across the
#: four slope phenotypes (0.263-0.464, per-split SD 0.10-0.18) while AHBA
#: coupling is not (|rho| 0.51-0.85), so transcriptional association breaks the
#: tie -- which is why PC3 and PC2 outrank PC1 despite PC1 scoring higher on h2.
PHENOTYPES: tuple[dict, ...] = (
    dict(name="baseline_thickness", priority=0, role="positive control",
         display="baseline thickness (control)"),
    dict(name="global_slope", priority=1, role="primary",
         display="global mean slope"),
    dict(name="slope_PC3", priority=2, role="target", display="slope PC3"),
    dict(name="slope_PC2", priority=3, role="target", display="slope PC2"),
    dict(name="slope_PC1", priority=4, role="secondary", display="slope PC1"),
)

PHENOTYPE_NAMES: tuple[str, ...] = tuple(p["name"] for p in PHENOTYPES)


def subject_phenotypes(run_dir: str | Path) -> pd.DataFrame:
    """The five settled subject-level phenotypes, one column each.

    Single source of truth for what these phenotypes *are*: both this module
    and ``tools/regen_h2_tables.py`` build them here, so the vector a GWAS is
    run on is by construction the same vector whose heritability the report
    publishes.  They were previously defined twice, in two files.

    ``global_slope`` and ``baseline_thickness`` are the cross-region means of
    the per-subject slope and intercept BLUPs.  The PC scores project each
    subject onto loadings fitted on the whole sample -- correct here, and *not*
    the same procedure as section 10's held-out evaluation, which refits
    loadings per split precisely so that a phenotype definition never sees the
    twins it is scored on.  For a GWAS there is no such circularity to avoid:
    the definition does not consume the genotypes.
    """
    run_dir = Path(run_dir)
    ph = pd.read_parquet(run_dir / "phenotypes" / "phenotypes.parquet")
    from . import covariance as cov

    out = {
        "global_slope": ph[ph.phenotype == "slope"].groupby("subject")["value"].mean(),
        "baseline_thickness": ph[ph.phenotype == "intercept"].groupby("subject")["value"].mean(),
    }
    W = cov.slope_matrix(run_dir)
    scores = cov.subject_scores(W, cov.slope_pcs(W, 3).loadings)
    for c in ("PC1", "PC2", "PC3"):
        out[f"slope_{c}"] = scores[c]

    frame = pd.DataFrame(out)
    frame.index.name = "subject"
    return frame[list(PHENOTYPE_NAMES)]


def build(run_dir: str | Path) -> dict[str, pd.DataFrame]:
    run_dir = Path(run_dir)
    _check_genetic_validity(run_dir)
    wide = subject_phenotypes(run_dir)

    meta_path = run_dir / "model_table.parquet"
    cols = ["subject", "family_id", "site", "sex", "age_c", "n_visits"]
    meta = pd.read_parquet(meta_path, columns=cols)
    meta = (meta.sort_values("age_c").groupby("subject", observed=True).first()
                .reset_index())

    df = wide.reset_index().merge(meta, on="subject", how="inner")
    if df.empty:
        raise ValueError(
            f"phenotypes and model_table for {run_dir.name} share no subjects"
        )
    df["IID"] = _to_genetics_id(df["subject"])
    df["FID"] = df["family_id"].astype(str)

    pheno_cols = list(PHENOTYPE_NAMES)
    phen = df[["FID", "IID"] + pheno_cols]

    # n_visits is a covariate because the phenotype is a shrunken BLUP: a
    # 2-visit subject's slope is pulled harder toward the fixed effect than a
    # 4-visit subject's, and attrition is not random.  This adjusts the mean,
    # not the variance -- the heteroscedasticity it leaves behind is a known
    # limitation, not something a linear covariate can absorb.
    qcov = (df[["FID", "IID", "age_c", "n_visits"]]
            .rename(columns={"age_c": "baseline_age"}))
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

    # --mpheno is 1-based over the phenotype columns only, i.e. it ignores
    # FID/IID.  Deriving it here from the frame that is about to be written is
    # the whole point: a script that greps the header instead has to reproduce
    # that offset convention, and got it wrong.
    manifest = pd.DataFrame([
        {"name": p["name"],
         "mpheno": phen.columns.get_loc(p["name"]) - 1,
         "priority": p["priority"],
         "role": p["role"],
         "n_nonmissing": int(phen[p["name"]].notna().sum())}
        for p in PHENOTYPES
    ]).sort_values("priority", ignore_index=True)

    return {"phenotypes_gcta": phen, "covar_quant": qcov,
            "covar_categorical": ccov, "phenotype_manifest": manifest}


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
        # The manifest is read by shell (cut/awk), so it is tab-separated with a
        # .tsv suffix; the GCTA inputs are space-delimited as GCTA expects.
        if name == "phenotype_manifest":
            p = out / "phenotype_manifest.tsv"
            frame.to_csv(p, sep="\t", index=False)
        else:
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
