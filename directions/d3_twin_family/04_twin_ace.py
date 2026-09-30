"""D3 part 1: twin correlations, univariate ACE, and bivariate ACE of the thinning rate
with each candidate co-developing trait.

    python directions/d3_twin_family/04_twin_ace.py [n_boot]     (env abcd-spatial; default 500)

Inputs : work/traits.parquet, work/pairs.csv (03_build_traits_pairs.py)
Outputs: results/twin_correlations.tsv  rMZ, rDZ, rSIB (double-entry ICC), Falconer h2
         results/ace_univariate.tsv     ACE / AE / CE / E fits, -2LL, LRTs, bootstrap CIs
         results/ace_bivariate.tsv      dCT x trait, ACE and AE: h2s, rA, rC, rE, rP and the
                                        share of rP due to A, bootstrap CIs. AE is the primary
                                        bivariate model: C for dCT is estimated at 0 (boundary),
                                        and a correlation is flagged not identified when either
                                        trait's variance component is < 0.02 (e.g. rA hits +-1).
                                        dCT_rq (scan-quality adjusted) is run against the
                                        puberty, thickness and CBCL traits as a sensitivity.
All traits are the covariate-residualised z-scores (<trait>_r). CIs are 2.5/97.5
percentiles of a pair bootstrap (resampling within zygosity), warm-started at the optimum.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2  # noqa: E402

import ace  # noqa: E402
from common import RES, WORK  # noqa: E402

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 500
T = pd.read_parquet(WORK / "traits.parquet")
pairs = pd.read_csv(WORK / "pairs.csv")
TRAITS = {
    "dCT_r": "thinning rate (dCT; higher = slower)",
    "dCT_rq": "thinning rate, scan-quality adjusted", "CT0_r": "cortical thickness (CT0)",
    "dT1T2_r": "T1w/T2w slope", "pub_timing_r": "puberty timing", "pub_tempo_r": "puberty tempo",
    "depress_chg_r": "CBCL depressive change", "internal_chg_r": "CBCL internalising change",
    "external_chg_r": "CBCL externalising change", "pfactor_chg_r": "CBCL p-factor change",
    "fluid_gain_r": "fluid cognition gain", "cryst_gain_r": "crystallised cognition gain",
    "total_gain_r": "total cognition gain",
}
q = lambda s: (float(np.nanpercentile(s, 2.5)), float(np.nanpercentile(s, 97.5)))

# ---------------------------------------------------------------- univariate ----------
corr_rows, uni_rows = [], []
for t, lab in TRAITS.items():
    mz = ace.pair_arrays(T, pairs, [t], "MZ"); dz = ace.pair_arrays(T, pairs, [t], "DZ")
    sib = ace.pair_arrays(T, pairs, [t], "SIB")
    corr_rows.append(dict(trait=t, label=lab, **ace.twin_corrs(mz, dz, sib)))
    fits = {m: ace.fit(mz, dz, m) for m in ("ACE", "AE", "CE", "E")}
    bs = ace.bootstrap(mz, dz, "ACE", NB, seed=11, start=fits["ACE"]["theta"], keys=["h2_1", "c2_1", "e2_1"])
    bsae = ace.bootstrap(mz, dz, "AE", NB, seed=12, start=fits["AE"]["theta"], keys=["h2_1"])
    for m, f in fits.items():
        row = dict(trait=t, label=lab, model=m, n_mz=f["n_mz"], n_dz=f["n_dz"], minus2ll=f["minus2ll"],
                   npar=f["npar"], h2=f["h2_1"], c2=f["c2_1"], e2=f["e2_1"], converged=f["converged"])
        if m != "ACE":
            d = f["minus2ll"] - fits["ACE"]["minus2ll"]; df = fits["ACE"]["npar"] - f["npar"]
            # C and A bounded at 0: the one-parameter tests are 50:50 chi2(0)/chi2(1) mixtures
            row.update(lrt_vs_ACE=d, df=df, p_vs_ACE=(0.5 if df == 1 else 1.0) * chi2.sf(max(d, 0), df))
        if m == "ACE":
            row.update(h2_lo=q(bs.h2_1)[0], h2_hi=q(bs.h2_1)[1], c2_lo=q(bs.c2_1)[0], c2_hi=q(bs.c2_1)[1])
        if m == "AE":
            row.update(h2_lo=q(bsae.h2_1)[0], h2_hi=q(bsae.h2_1)[1])
        uni_rows.append(row)
    print(t, {k: round(v, 3) for k, v in corr_rows[-1].items() if isinstance(v, float)}, flush=True)
pd.DataFrame(corr_rows).to_csv(RES / "twin_correlations.tsv", sep="\t", index=False, float_format="%.4g")
pd.DataFrame(uni_rows).to_csv(RES / "ace_univariate.tsv", sep="\t", index=False, float_format="%.4g")

# ---------------------------------------------------------------- bivariate -----------
NBB = max(NB // 2, 100)
KEYS = ["h2_1", "h2_2", "c2_1", "c2_2", "rA", "rC", "rE", "rP", "rP_A"]
bi_rows = []
QSET = {"CT0_r", "pub_timing_r", "pub_tempo_r", "depress_chg_r", "internal_chg_r", "external_chg_r",
        "pfactor_chg_r"}
jobs = [("dCT_r", t) for t in TRAITS if t not in ("dCT_r", "dCT_rq")] + [("dCT_rq", t) for t in TRAITS if t in QSET]
for t1, t in jobs:
    lab = TRAITS[t]
    tr = [t1, t]
    mz = ace.pair_arrays(T, pairs, tr, "MZ"); dz = ace.pair_arrays(T, pairs, tr, "DZ")
    for m in ("ACE", "AE"):
        f = ace.fit(mz, dz, m)
        keys = [k for k in KEYS if not (m == "AE" and k in ("c2_1", "c2_2", "rC"))]
        bs = ace.bootstrap(mz, dz, m, NBB, seed=21, start=f["theta"], keys=keys)
        row = dict(trait1=t1, trait2=t, label2=lab, model=m, n_mz=f["n_mz"], n_dz=f["n_dz"],
                   minus2ll=f["minus2ll"], npar=f["npar"], n_boot=NBB)
        for k in keys:
            row[k] = f[k]; row[f"{k}_lo"], row[f"{k}_hi"] = q(bs[k])
        # bootstrap two-sided p for rA != 0 (share of resamples on the other side of 0)
        row["rA_identified"] = bool(min(f["h2_1"], f["h2_2"]) >= 0.02)
        row["rC_identified"] = bool(m == "ACE" and min(f.get("c2_1", 0), f.get("c2_2", 0)) >= 0.02)
        row["p_rA_boot"] = float(min(1.0, 2 * min((bs.rA <= 0).mean(), (bs.rA >= 0).mean())))
        bi_rows.append(row)
    print("bivariate", t1, t, round(bi_rows[-1]["rA"], 3), bi_rows[-1]["rA_lo"].__round__(3),
          bi_rows[-1]["rA_hi"].__round__(3), flush=True)
pd.DataFrame(bi_rows).to_csv(RES / "ace_bivariate.tsv", sep="\t", index=False, float_format="%.4g")
