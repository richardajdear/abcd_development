"""02 -- reliability and coupling of the single-LMM C3-axis phenotypes (README_HPC.md 8.2).

Reads the fit of out/<run>_c3axis (build_c3axis_model_table.py, then R/fit_lmm.R) and writes
c3axis_1lmm_summary.tsv: slope variance, residual variance, slope reliability (1 - mean PEV / tau^2),
the fixed age slope, and correlations of each child slope with the global-mean slope, sex, age at
first scan, and (rc components) the BLUP-level row-centred score from the per-parcel fit.
Run from the repo root:  python genetic_analysis/c3axis/02_c3axis_1lmm_summary.py
"""
import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, RUN, child_covariates, load_inputs, pca_raw, row_center

FIT = RUN.parent / f"{RUN.name}_c3axis" / "fits"
b3 = pd.read_parquet(FIT / "blups.parquet"); vc = pd.read_parquet(FIT / "varcomp.parquet")
fx = pd.read_parquet(FIT / "fixed.parquet")
tau = vc[(vc.grp == "subject") & (vc.var1 == "age_c") & vc.var2.isna()].set_index("label").vcov
tauI = vc[(vc.grp == "subject") & (vc.var1 == "(Intercept)") & vc.var2.isna()].set_index("label").vcov
sig = vc[vc.grp == "Residual"].set_index("label").vcov
rel = 1 - b3.groupby("label").se_re_slope.apply(lambda s: (s ** 2).mean()) / tau

Wh, WL, REFL, maps, ylab = load_inputs()
cv = child_covariates(WL.index)
Sl = b3.pivot(index="subject", columns="label", values="re_slope").loc[WL.index]
WLc = row_center(WL); Gc, _ = pca_raw(WLc)
sc = pd.DataFrame(WLc.to_numpy() @ Gc.to_numpy(), index=WL.index, columns=Gc.columns)
blup = {"c1axis_rc": sc.c1, "c2axis_rc": sc.c2, "c3axis_rc": sc.c5}

summ = pd.DataFrame({"tau2_slope": tau, "sigma2_resid": sig, "slope_reliability": rel,
                     "fixed_age_slope": fx[fx.term == "age_c"].set_index("label").estimate,
                     "slope_SD_over_baseline_SD": np.sqrt(tau) / np.sqrt(tauI)})
summ["r_with_global_mean_slope"] = Sl.corrwith(Sl.global_mean)
summ["r_male"] = Sl.corrwith(cv.male)
summ["r_age_first"] = Sl.corrwith(cv.age_first)
summ["r_vs_BLUPlevel_LHscore"] = [abs(Sl[c].corr(blup[c])) if c in blup else np.nan for c in summ.index]
summ.to_csv(HERE / "c3axis_1lmm_summary.tsv", sep="\t", float_format="%.4g")
print(summ.round(3).to_string())
