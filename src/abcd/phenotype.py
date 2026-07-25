"""
Subject-level developmental phenotypes from fitted mixed models.

Three things this module insists on, each a response to a specific way the
5.1 pipeline could mislead:

1. **A phenotype is the random slope, not a difference of predictions.**
   The 5.1 code differenced model predictions at ages 9 and 16 while the
   oldest observed scan was 15.75.  Because the model is linear in age, that
   difference is exactly ``7 x (fixed slope + random slope)`` -- the same
   quantity rescaled, dressed up as a change score, and extrapolated.  Here
   the slope is the phenotype and its units are stated.

2. **Reliability travels with the phenotype.**  Every subject-by-region slope
   carries the conditional SD from which reliability is computed.  A GWAS run
   on slopes whose reliability is 0.15 has its heritability attenuated by
   roughly that factor; that is a fact about the design, and it should be
   visible in the phenotype file rather than discovered afterwards.

3. **Reliability is confounded with attrition.**  Subjects with three visits
   have systematically more reliable slopes than subjects with two, and
   whether a subject returns is not random.  ``reliability_weight`` and the
   ``n_visits`` column exist so downstream analyses can condition on or weight
   by design, instead of silently treating all slopes as equivalent.

The output is deliberately long/tidy and metric-agnostic:

    subject | metric | label | phenotype | value | se | reliability | n_visits
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import paths


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_fits(run_dir: str | Path, fits_name: str = "fits") -> dict[str, pd.DataFrame]:
    """Load the four Parquet outputs written by R/fit_lmm.R."""
    d = Path(run_dir) / fits_name
    out = {}
    for k in ("blups", "fixed", "varcomp", "diagnostics"):
        f = d / f"{k}.parquet"
        if not f.exists():
            raise FileNotFoundError(
                f"{f} not found -- run: Rscript R/fit_lmm.R --run-dir {run_dir}"
            )
        out[k] = pd.read_parquet(f)
    return out


# --------------------------------------------------------------------------
# Reliability
# --------------------------------------------------------------------------

def slope_reliability(blups: pd.DataFrame, varcomp: pd.DataFrame) -> pd.DataFrame:
    """Attach per-subject slope reliability to the BLUP table.

    For a random effect ``b`` with prior variance ``tau2`` and posterior
    (conditional) variance ``v``, the reliability of the BLUP is

        rho = 1 - v / tau2

    which is 0 when the data say nothing about that subject (the BLUP is
    shrunk entirely to the fixed effect) and approaches 1 when the subject's
    own data dominate.  This is the same quantity as the ratio of BLUP
    variance to prior variance, computed per subject rather than in aggregate.
    """
    tau2 = (varcomp[(varcomp.grp == "subject") & (varcomp.var1 == "age_c")
                    & (varcomp.var2.isna())][["label", "vcov"]]
            .rename(columns={"vcov": "prior_var_slope"}))
    tau2_i = (varcomp[(varcomp.grp == "subject") & (varcomp.var1 == "(Intercept)")
                      & (varcomp.var2.isna())][["label", "vcov"]]
              .rename(columns={"vcov": "prior_var_intercept"}))
    b = blups.merge(tau2, on="label", how="left").merge(tau2_i, on="label", how="left")
    b["reliability_slope"] = 1.0 - (b.se_re_slope ** 2) / b.prior_var_slope
    b["reliability_intercept"] = 1.0 - (b.se_re_intercept ** 2) / b.prior_var_intercept
    # numerical guard: a singular fit gives prior_var == 0 -> undefined
    for c in ("reliability_slope", "reliability_intercept"):
        b[c] = b[c].where(np.isfinite(b[c])).clip(0.0, 1.0)
    return b


def design_reliability(times: np.ndarray | list[float], tau2: float,
                       sigma2: float) -> float:
    """Analytic reliability of an OLS slope for one subject's visit schedule.

    rho = S_tt * tau2 / (S_tt * tau2 + sigma2),  S_tt = sum (t - tbar)^2

    Useful for projecting what a future release's design will buy before the
    data exist: reliability grows with the *spread* of visit times, so a
    longer span helps more than an extra visit squeezed into the same window.
    """
    t = np.asarray(times, dtype=float)
    if t.size < 2:
        return 0.0
    stt = float(((t - t.mean()) ** 2).sum())
    return stt * tau2 / (stt * tau2 + sigma2)


# --------------------------------------------------------------------------
# Phenotype construction
# --------------------------------------------------------------------------

def build_phenotypes(run_dir: str | Path, fits_name: str = "fits",
                     model_table: pd.DataFrame | None = None) -> pd.DataFrame:
    """Build the tidy subject-level phenotype table for a fitted run.

    Phenotypes produced, per subject and region:

    ``slope``
        Random slope + fixed slope: the subject's own rate of change in
        metric units per year.  This is the primary phenotype.
    ``slope_dev``
        Random slope alone: deviation from the sample-average rate.  Identical
        to ``slope`` up to an additive constant, so it gives the same GWAS;
        provided because it is the natural scale for spatial comparisons.
    ``intercept``
        Random intercept + fixed intercept at the centring age: the subject's
        level, i.e. the phenotype a cross-sectional study would measure.
    """
    fits = load_fits(run_dir, fits_name)
    b = slope_reliability(fits["blups"], fits["varcomp"])

    fx = fits["fixed"]
    fixed_slope = (fx[fx.term == "age_c"].set_index("label").estimate
                   .rename("fixed_slope"))
    fixed_int = (fx[fx.term == "(Intercept)"].set_index("label").estimate
                 .rename("fixed_intercept"))
    b = b.merge(fixed_slope, on="label", how="left").merge(fixed_int, on="label", how="left")

    diag = fits["diagnostics"][["label", "singular", "converged"]]
    b = b.merge(diag, on="label", how="left")

    if model_table is None:
        mt = pd.read_parquet(Path(run_dir) / "model_table.parquet",
                             columns=["subject", "metric", "n_visits", "age_span",
                                      "age_first", "age_last"])
        model_table = mt.drop_duplicates("subject")
    b = b.merge(model_table, on="subject", how="left")

    manifest = json.loads((Path(run_dir) / "manifest.json").read_text())
    age_centre = manifest.get("age_centre")

    long = []
    for name, value, se, rel in (
        ("slope", b.re_slope + b.fixed_slope, b.se_re_slope, b.reliability_slope),
        ("slope_dev", b.re_slope, b.se_re_slope, b.reliability_slope),
        ("intercept", b.re_intercept + b.fixed_intercept, b.se_re_intercept,
         b.reliability_intercept),
    ):
        long.append(pd.DataFrame({
            "subject": b.subject, "metric": b.metric, "label": b.label,
            "phenotype": name, "value": value, "se": se, "reliability": rel,
            "n_visits": b.n_visits, "age_span": b.age_span,
            "age_first": b.age_first, "age_last": b.age_last,
            "singular_fit": b.singular, "converged": b.converged,
        }))
    out = pd.concat(long, ignore_index=True)
    out.attrs["age_centre"] = age_centre
    out.attrs["run_id"] = manifest["run_id"]
    return out


def reliability_summary(pheno: pd.DataFrame) -> pd.DataFrame:
    """Reliability by phenotype and visit count -- the attrition confound."""
    return (pheno.groupby(["phenotype", "n_visits"], observed=True)
            .agg(n=("value", "size"),
                 mean_reliability=("reliability", "mean"),
                 median_reliability=("reliability", "median"),
                 sd_value=("value", "std"))
            .round(4).reset_index())


def write_phenotypes(pheno: pd.DataFrame, run_dir: str | Path,
                     name: str = "phenotypes") -> Path:
    """Write the tidy phenotype table plus its reliability summary."""
    d = Path(run_dir) / "phenotypes"
    d.mkdir(parents=True, exist_ok=True)
    pheno.to_parquet(d / f"{name}.parquet", index=False)
    reliability_summary(pheno).to_csv(d / f"{name}_reliability.csv", index=False)
    return d


# --------------------------------------------------------------------------
# Export for GCTA / PLINK
# --------------------------------------------------------------------------

def to_gcta(pheno: pd.DataFrame, phenotype: str = "slope",
            labels: list[str] | None = None,
            min_reliability: float | None = None,
            id_map: pd.Series | None = None) -> pd.DataFrame:
    """Wide FID/IID table for GCTA ``--pheno`` / ``--mpheno``.

    ``id_map`` maps ``subject`` to genotype IID when the imaging and genetic
    identifiers differ (they do in ABCD).  Without it the subject id is used
    for both columns and the caller is responsible for the join.
    """
    p = pheno[pheno.phenotype == phenotype]
    if labels is not None:
        p = p[p.label.isin(labels)]
    if min_reliability is not None:
        p = p[p.reliability >= min_reliability]
    wide = p.pivot_table(index="subject", columns="label", values="value")
    iid = wide.index.map(id_map) if id_map is not None else wide.index
    out = pd.DataFrame({"FID": iid, "IID": iid}).set_index(wide.index)
    return pd.concat([out, wide], axis=1).reset_index(drop=True)


if __name__ == "__main__":  # pragma: no cover
    import argparse

    ap = argparse.ArgumentParser(description="Build subject-level phenotypes.")
    ap.add_argument("run_dir")
    ap.add_argument("--fits-name", default="fits")
    a = ap.parse_args()
    ph = build_phenotypes(a.run_dir, a.fits_name)
    d = write_phenotypes(ph, a.run_dir)
    print(f"{len(ph):,} rows -> {d}")
    print(reliability_summary(ph).to_string(index=False))
