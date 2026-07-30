"""
Regional maps derived from a fitted run.

The distinction this module exists to enforce
---------------------------------------------
When ``global_covariate`` is set (the default, and the right default), each
regional model is

    value ~ sex + age_c + global_between_c + global_within + (1 + age_c | subject) + ...

so the region's ``age_c`` coefficient is **not** its thinning rate.  It is the
rate *relative to* whatever the cortex as a whole is doing, because
``global_within`` absorbs the shared component.  On 7.0 the adjusted
coefficient is positive in 24 of 68 regions -- which does not mean those
regions thicken.  It means they thin more slowly than average.

Both quantities are scientifically useful and they answer different questions:

* ``adjusted``  -- where is maturation fast or slow *for this cortex*?  This is
  the right map for asking about regional specialisation, and it is what the
  random slopes (and hence the GWAS phenotype) are conditioned on.
* ``total``     -- how fast does this region actually thin, in mm/year?  This
  is the map that is comparable to published developmental maps, including
  the differenced-prediction map from the 5.1 thesis analysis.

Getting this wrong changes conclusions.  The correlation with AHBA C3 is
-0.57 (p_spin = 0.003) for the total rate and -0.32 (p_spin = 0.10) for the
adjusted coefficient: the transcriptional association lives substantially in
the *global* thinning component that the covariate removes.  A pipeline that
only ever exposes the adjusted coefficient would report a null and be wrong
about why.

Reconstruction, not refitting
-----------------------------
``total_age_slope`` recovers the unadjusted rate from the adjusted fit by the
chain rule, rather than refitting without the covariate:

    d(value)/d(age) = beta_age + beta_global_within * d(global_within)/d(age)

The last term is one number for the whole run -- the within-subject regression
of global thickness on centred age -- so this costs nothing.  Validation: on
7.0 the reconstruction correlates r = 0.979 with an independent no-covariate
refit on 5.1, and is negative in 67 of 68 regions as an absolute rate must be.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .phenotype import load_fits


def global_age_slope(model_table: pd.DataFrame) -> float:
    """Within-subject regression slope of global thickness on centred age.

    ``global_within`` is already subject-demeaned by ``assemble``, so an OLS
    fit against ``age_c`` on one region's rows (every region carries the same
    subject-level global columns) is the within-subject slope.

    Returns mm/year, negative during adolescence.
    """
    if "global_within" not in model_table.columns:
        raise KeyError(
            "model_table has no 'global_within' column -- this run was "
            "assembled without a global covariate, so its age coefficients "
            "are already total rates and need no reconstruction."
        )
    one = model_table[model_table.label == model_table.label.iloc[0]]
    g = one[["age_c", "global_within"]].dropna()
    if g.global_within.abs().max() == 0:
        raise ValueError("global_within is identically zero")
    X = np.c_[np.ones(len(g)), g.age_c.to_numpy()]
    return float(np.linalg.lstsq(X, g.global_within.to_numpy(), rcond=None)[0][1])


def total_age_slope(fixed: pd.DataFrame, dglobal_dage: float) -> pd.Series:
    """Total (unadjusted) thinning rate per region, mm/year.

    ``fixed`` is the ``fixed.parquet`` written by ``R/fit_lmm.R``;
    ``dglobal_dage`` comes from :func:`global_age_slope`.
    """
    piv = fixed.pivot(index="label", columns="term", values="estimate")
    missing = {"age_c", "global_within"} - set(piv.columns)
    if missing:
        raise KeyError(f"fixed.parquet has no term(s) {sorted(missing)}")
    return (piv["age_c"] + piv["global_within"] * dglobal_dage).rename("slope_total")


def regional_maps(run_dir: str | Path, fits_name: str = "fits") -> pd.DataFrame:
    """Region-level summary maps for one run, indexed by ``label``.

    Columns
    -------
    slope_adjusted : age coefficient with the global covariate in the model --
                     rate *relative* to cortex-wide thinning.
    slope_total    : absolute rate, mm/year (see module docstring).
    slope_t        : t statistic for the fitted age coefficient.
    tau_slope      : between-subject SD of the random slope, mm/year.

    Both run types are handled, and the distinction is decided by the *fitted
    terms*, not by the model table: ``assemble`` writes ``global_within`` into
    ``model_table.parquet`` regardless of whether the fit used it, so keying on
    the table's columns would silently misclassify a no-global run.

    For a run fitted *without* the global covariate the age coefficient is
    already the absolute rate, so ``slope_total`` and ``slope_adjusted`` are
    the same column.  They are both returned, identical, so that downstream
    code and figures can name the quantity they mean without branching on the
    config -- and ``global_covariate`` in the returned frame's ``attrs``
    records which case produced it.
    """
    run_dir = Path(run_dir)
    f = load_fits(run_dir, fits_name)
    piv = f["fixed"].pivot(index="label", columns="term", values="estimate")

    out = pd.DataFrame({"slope_adjusted": piv["age_c"]})
    if "global_within" in piv.columns:
        mt = pd.read_parquet(run_dir / "model_table.parquet",
                             columns=["label", "age_c", "global_within"])
        out["slope_total"] = total_age_slope(f["fixed"], global_age_slope(mt))
        out.attrs["global_covariate"] = True
    else:
        out["slope_total"] = out["slope_adjusted"]
        out.attrs["global_covariate"] = False

    stat = f["fixed"].pivot(index="label", columns="term", values="statistic")
    if "age_c" in stat.columns:
        out["slope_t"] = stat["age_c"]

    vc = f["varcomp"]
    tau = vc[(vc.grp == "subject") & (vc.var1 == "age_c") & (vc.var2.isna())]
    if len(tau):
        out["tau_slope"] = tau.set_index("label")["sdcor"]
    return out
