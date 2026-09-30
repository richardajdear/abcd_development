"""D3 part 3: co-twin control. Does the faster-thinning member of a pair have more
symptoms at 15-17, holding constant everything the pair shares?

    python directions/d3_twin_family/06_cotwin_control.py        (env abcd-spatial)

Outcome model (the Figure-1 CBCL change model, genetic_analysis/fig1_prep_cbcl.py):
  y_LATE ~ dCT + y_BASE + age_LATE + sex + site      (log1p raw CBCL; beta / SD(y_LATE))
  population : all imaged children, family-clustered SE (reproduces Figure 1g)
  within-pair: same model with a pair fixed effect (= regression of within-pair
               differences), HC1 SE; run separately in MZ, DZ and sibling pairs.
               Sex and site drop out for MZ pairs (constant within pair).
If the within-MZ coefficient matches the population one, the slope-symptom association
is not explained by genes or family environment shared by co-twins; if it falls to
zero it is. DZ and sibling pairs sit in between (they share half the segregating genes).
Sensitivities: + CT0 (baseline thickness); + scan quality (the child's mean log1p surface
topological-defect count and its slope on age, as in dCT_rq); both.
Measurement diagnostics (results/cotwin_measurement.tsv), because co-twins are scanned in
the same session and a classical independent-error reliability argument does not hold:
  r_dCT_raw / r_dCT_resid   co-twin correlation of dCT before / after the covariates
  r_scan_resid     co-twin correlation of scan-level residuals from each child's own
                   line (children with >= 3 scans), matched visit
  r_logtopo        co-twin correlation of per-scan log1p topological defects
  age_gap_days     median |age difference| at a matched visit (0 = same-day scans)
Writes results/cotwin_control.tsv and results/cotwin_measurement.tsv.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402

from common import RES, WORK  # noqa: E402

T = pd.read_parquet(WORK / "traits.parquet")
pairs = pd.read_csv(WORK / "pairs.csv")
CBCL = ["depress", "internal", "external", "pfactor"]

rows = []
for o in CBCL:
    y, yb = f"{o}_LATE", f"{o}_BASE"
    sdy = T[y].std()
    QC = " + qc_ltopo_mean + qc_ltopo_slope"
    for extra in ("", " + CT0", QC, " + CT0" + QC):
        d = T.dropna(subset=[y, yb, "dCT", "age_LATE", "CT0", "qc_ltopo_mean", "qc_ltopo_slope"]).copy()
        d["family"] = d.family.astype(str)
        m = smf.ols(f"{y} ~ dCT + {yb} + age_LATE + C(sex) + C(site){extra}", d).fit(
            cov_type="cluster", cov_kwds={"groups": d.family})
        adj = {"": "none", " + CT0": "CT0", QC: "quality", " + CT0" + QC: "CT0+quality"}[extra]
        rows.append(dict(outcome=o, adjust=adj, design="population",
                         n_children=int(m.nobs), n_pairs=np.nan, beta=m.params["dCT"] / sdy,
                         se=m.bse["dCT"] / sdy, p=m.pvalues["dCT"]))
        for ty in ("MZ", "DZ", "SIB"):
            P = pairs[pairs.type == ty].reset_index(drop=True)
            L = pd.concat([T.reindex(P.id1).reset_index().assign(pair=P.index),
                           T.reindex(P.id2).reset_index().assign(pair=P.index)])
            L = L.dropna(subset=[y, yb, "dCT", "age_LATE", "CT0", "qc_ltopo_mean", "qc_ltopo_slope"])
            L = L[L.groupby("pair").pair.transform("size") == 2]
            cols = ["dCT", yb, "age_LATE", "CT0", "male", "qc_ltopo_mean", "qc_ltopo_slope", y]
            W = L[cols] - L.groupby("pair")[cols].transform("mean")
            W["site_diff"] = (L.groupby("pair").site.transform("nunique") > 1).astype(float)
            xs = ["dCT", yb, "age_LATE"] + (["CT0"] if "CT0" in extra else []) + \
                 (["qc_ltopo_mean", "qc_ltopo_slope"] if "qc_" in extra else []) + ([] if ty == "MZ" else ["male"])
            mw = smf.ols(f"{y} ~ 0 + " + " + ".join(xs), W).fit(cov_type="HC1")
            # HC1 on demeaned data: correct dof for the absorbed pair means
            npair = L.pair.nunique(); k = len(xs)
            se = mw.bse["dCT"] * np.sqrt((len(W) - k) / (len(W) - npair - k))
            from scipy.stats import t as tdist
            tval = mw.params["dCT"] / se
            rows.append(dict(outcome=o, adjust=adj, design=f"within {ty}",
                             n_children=len(W), n_pairs=int(npair), beta=mw.params["dCT"] / sdy,
                             se=se / sdy, p=float(2 * tdist.sf(abs(tval), len(W) - npair - k))))
R_ = pd.DataFrame(rows)
R_["lo"], R_["hi"] = R_.beta - 1.96 * R_.se, R_.beta + 1.96 * R_.se
R_.to_csv(RES / "cotwin_control.tsv", sep="\t", index=False, float_format="%.4g")
print(R_[R_.adjust == "none"].round(4).to_string(index=False))

# ---------------------------------------------------------------- measurement ---------
sm = pd.read_csv(WORK / "scan_means_ct.csv")
qc = pd.read_csv(HERE.parents[1] / "abcd-data-release-7.0/y/mr_y_qc__post__aut.tsv", sep="\t", low_memory=False,
                 usecols=["participant_id", "session_id", "mr_y_qc__post__aut__smri__topodfct_count"])
qc["subject"] = "sub-NDARINV" + qc.participant_id.str.replace("sub-", "", regex=False)
qc["visit"] = qc.session_id.map({"ses-00A": "v0", "ses-02A": "v2", "ses-04A": "v4", "ses-06A": "v6"})
qc["ltopo"] = np.log1p(pd.to_numeric(qc.mr_y_qc__post__aut__smri__topodfct_count, errors="coerce"))
sm = sm.merge(qc[["subject", "visit", "ltopo"]], on=["subject", "visit"], how="left")


def own_resid(g):
    if len(g) < 3:
        return pd.Series(np.nan, index=g.index)
    bb = np.polyfit(g.age, g.mean_val, 1)
    return g.mean_val - np.polyval(bb, g.age)


sm["res"] = sm.groupby("subject", group_keys=False)[["age", "mean_val"]].apply(own_resid)
S = sm.set_index(["subject", "visit"])
de = lambda u, w: float(np.corrcoef(np.r_[u, w], np.r_[w, u])[0, 1])
mm = []
for ty in ("MZ", "DZ", "SIB"):
    P = pairs[pairs.type == ty]
    rec = []
    for v in ("v0", "v2", "v4", "v6"):
        a = S.reindex(list(zip(P.id1, [v] * len(P)))); b = S.reindex(list(zip(P.id2, [v] * len(P))))
        rec.append(pd.DataFrame({"ra": a.res.to_numpy(), "rb": b.res.to_numpy(), "ta": a.ltopo.to_numpy(),
                                 "tb": b.ltopo.to_numpy(), "gap": (a.age.to_numpy() - b.age.to_numpy())}))
    X = pd.concat(rec); xr = X.dropna(subset=["ra", "rb"]); xt = X.dropna(subset=["ta", "tb"])
    A1, A2 = T.reindex(P.id1), T.reindex(P.id2)
    ok = A1.dCT.notna().to_numpy() & A2.dCT.notna().to_numpy()
    mm.append(dict(pair_type=ty, n_pairs=len(P),
                   r_dCT_raw=de(A1.dCT.to_numpy()[ok], A2.dCT.to_numpy()[ok]),
                   r_dCT_resid=de(A1.dCT_r.to_numpy()[ok], A2.dCT_r.to_numpy()[ok]),
                   r_dCT_resid_quality=de(A1.dCT_rq.to_numpy()[ok], A2.dCT_rq.to_numpy()[ok]),
                   n_scan_pairs=len(xr), r_scan_resid=de(xr.ra, xr.rb),
                   r_logtopo=de(xt.ta, xt.tb), age_gap_days=float(np.nanmedian(np.abs(X.gap)) * 365.25),
                   mean_blup_reliability=float(T.dCT_rel.mean())))
pd.DataFrame(mm).to_csv(RES / "cotwin_measurement.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(mm).round(3).to_string(index=False))
