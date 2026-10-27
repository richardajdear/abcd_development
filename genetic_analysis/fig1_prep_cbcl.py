"""Figure 1 panel h: symptoms vs the two Figure-1 phenotypes (HCP-MMP, 7.0).

Usage (repo root, after `Rscript genetic_analysis/fig1_prep_1lmm.R`):
    python genetic_analysis/fig1_prep_cbcl.py

Brain predictors -- the SAME per-child quantities the genetics panels use (the
single LMM on the per-scan cortical mean, orderops/build_1lmm_pheno.py
definition), read from /tmp/hcp70_1lmm_blups.csv (individual-level, never
committed):
    global_slope        random slope (per SD; negative = faster thinning)
    baseline_thickness  random intercept (per SD; thickness at the sample mean
                        age, the positive-control phenotype of the PRS grid)
Outcomes and models are taken verbatim from ahba_pls/code/23_cbcl_explore.py
(its outcome-construction block is executed, not re-implemented), so the
definitions match fig_cbcl_candidate.png:
    CBCL change   y_LATE ~ brain + y_BASE + age_LATE + sex + site   (log1p raw
                  scores; LATE = mean of waves 5-7; beta / SD(y_LATE))
    KSADS         lifetime dx ~ brain + age_LAST + sex + site  (logit)
    OLS / logit, site fixed effects, family-clustered SEs.
Each brain predictor is fitted ALONE (one model per phenotype), and jointly
(both in the model) so the slope's association can be read net of baseline.
Writes fig1_inputs/hcp70_cbcl_assoc.tsv (coefficients only).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SRC_SCRIPT = REPO / "ahba_pls/code/23_cbcl_explore.py"
BLUPS = Path("/tmp/hcp70_1lmm_blups.csv")
OUT = REPO / "genetic_analysis/fig1_inputs/hcp70_cbcl_assoc.tsv"

CBCL = ["pfactor", "internal", "external", "depress", "thought"]
# PQ-BC: number of psychotic-like experiences endorsed (0-21) and distress score
# (sum of distress ratings over endorsed items); youth report, all waves
PPS = {"ple_count": "mh_y_pps_count", "ple_distress": "mh_y_pps__severity_score"}
DXS = ["mdd_youth_DX", "mdd_parent_DX", "psychosis_parent_DX"]


def build_outcomes() -> dict:
    """Execute 23_cbcl_explore.py up to its model fits; return its namespace."""
    src = SRC_SCRIPT.read_text()
    cut = src.index("\ndef fit(")
    ns = {"__file__": str(SRC_SCRIPT), "__name__": "cbcl_prefix"}
    argv, sys.argv = sys.argv, [str(SRC_SCRIPT), "70"]
    try:
        exec(compile(src[:cut], str(SRC_SCRIPT), "exec"), ns)
    finally:
        sys.argv = argv
    return ns


def main() -> None:
    ns = build_outcomes()
    D = ns["D"].copy()
    b = pd.read_csv(BLUPS).set_index("subject")
    z = lambda s: (s - s.mean()) / s.std()
    D["global_slope"] = z(b.re_slope).reindex(D.index)
    D["baseline_thickness"] = z(b.re_intercept).reindex(D.index)
    assert D.global_slope.notna().mean() > 0.99, D.global_slope.notna().mean()
    import statsmodels.formula.api as smf

    def fit(y, xs, extra, kind="ols"):
        f = f"{y} ~ " + " + ".join(xs) + " + C(sex) + C(site)" + "".join(f" + {e}" for e in extra)
        d = D.dropna(subset=[y, *xs, *extra])
        m = (smf.logit if kind == "logit" else smf.ols)(f, d).fit(
            cov_type="cluster", cov_kwds={"groups": d.family},
            **({"disp": 0} if kind == "logit" else {}))
        return m, d

    # PQ-BC (youth self-report psychotic-like experiences, 7.0 mh_y_pps.tsv), built
    # like the CBCL scores: log1p, BASE = ses-00A, LATE = mean of ses-05A..07A
    pps = pd.read_csv(ns["R70"] / "mh_y_pps.tsv", sep="\t",
                      usecols=["participant_id", "session_id", "mh_y_pps_age", *PPS.values()])
    pps["subject"] = "sub-NDARINV" + pps.participant_id.str.replace("sub-", "", regex=False)
    pps["wave"] = pps.session_id.str.slice(4, 6).astype(int)
    pps = pps[pps.subject.isin(D.index)].drop_duplicates(["subject", "wave"])
    for k, c in PPS.items():
        pps[k] = np.log1p(pd.to_numeric(pps[c], errors="coerce"))
    for name, waves in ns["WAVES"].items():
        g = pps[pps.wave.isin(waves)].groupby("subject")[list(PPS)].mean()
        D = D.join(g.add_suffix(f"_{name}"))

    FU = ns["FU"]
    rows = []
    for model_set, xs_list in (("alone", [["global_slope"], ["baseline_thickness"]]),
                               ("joint", [["global_slope", "baseline_thickness"]])):
        for xs in xs_list:
            for o in CBCL + list(PPS):
                y = f"{o}_{FU}"
                m, d = fit(y, xs, [f"age_{FU}", f"{o}_BASE"])
                for x in xs:
                    rows.append(dict(outcome=o, model="change", fit=model_set, brain=x,
                                     beta=m.params[x], se=m.bse[x], p=m.pvalues[x],
                                     n=int(m.nobs), n_cases=np.nan, y_sd=D[y].std()))
            for o in DXS:
                m, d = fit(o, xs, ["age_LAST"], kind="logit")
                for x in xs:
                    rows.append(dict(outcome=o, model="ksads_logit", fit=model_set, brain=x,
                                     beta=m.params[x], se=m.bse[x], p=m.pvalues[x],
                                     n=int(m.nobs), n_cases=int(d[o].sum()), y_sd=np.nan))
    # ---- panel h: depressive symptoms at 15-17 y beyond baseline vs ΔCT ----
    # y = residual of log1p DSM-depressive (LATE) on baseline score + age_LATE +
    # sex + site, in SD of y_LATE (so its OLS slope on z(ΔCT) reproduces panel g's
    # estimate up to the Frisch-Waugh step on the covariates); x = ΔCT in µm/yr
    # from hcp70_child_traits.csv via the same sid key fig1_prep_1lmm.R writes.
    yv = f"depress_{FU}"
    d = D.dropna(subset=[yv, "depress_BASE", f"age_{FU}", "global_slope"]).copy()
    r = smf.ols(f"{yv} ~ depress_BASE + age_{FU} + C(sex) + C(site)", d).fit()
    d["dep_resid_sd"] = r.resid / D[yv].std()
    order = sorted(b.index)                                   # codepoint sort, as sids in R
    sid = pd.Series(range(len(order)), index=order)
    tr = pd.read_csv(OUT.parent / "hcp70_child_traits.csv").set_index("sid")
    d["dct_um_per_yr"] = tr.dct_um_per_yr.reindex(sid.reindex(d.index).values).values
    assert d.dct_um_per_yr.notna().all()
    assert abs(np.corrcoef(d.dct_um_per_yr, d.global_slope)[0, 1] - 1) < 1e-6
    # equal-width 1 µm/yr bins centred on the mean ΔCT (edges mean ± 4.5); the few
    # children beyond the outer edges (~3-4%) are left out of the bin means but kept
    # in the fit, so the fastest-thinning bin still holds ~250+ children
    W_BIN = 1.0
    c0 = np.round(d.dct_um_per_yr.mean())
    edges = c0 + W_BIN * np.arange(-4.5, 4.5 + 1e-9, 1.0)
    q = pd.cut(d.dct_um_per_yr, edges, right=False, labels=False)
    bins = d.groupby(q).agg(x=("dct_um_per_yr", "mean"), y=("dep_resid_sd", "mean"),
                            se=("dep_resid_sd", "sem"), n=("dep_resid_sd", "size"))
    lf = smf.ols("dep_resid_sd ~ dct_um_per_yr", d).fit(cov_type="cluster", cov_kwds={"groups": d.family})
    bins["fit_slope_per_um"] = lf.params["dct_um_per_yr"]; bins["fit_intercept"] = lf.params["Intercept"]
    bins["fit_p"] = lf.pvalues["dct_um_per_yr"]; bins["n_total"] = len(d)
    bins["lo_edge"] = edges[bins.index.astype(int)]; bins["width"] = W_BIN
    bins["n_outside"] = int(q.isna().sum())
    print("bin n:", bins.n.tolist(), "outside:", int(q.isna().sum()))
    bins.to_csv(OUT.parent / "hcp70_dep_bins.csv", index=False, float_format="%.5g")   # aggregate: committed
    print("panel h: n", len(d), "slope per SD ΔCT", round(lf.params["dct_um_per_yr"] * d.dct_um_per_yr.std(), 4),
          "p", round(lf.pvalues["dct_um_per_yr"], 4))

    A = pd.DataFrame(rows)
    A["est"] = np.where(A.model == "change", A.beta / A.y_sd, np.exp(A.beta))
    A["lo"] = np.where(A.model == "change", (A.beta - 1.96 * A.se) / A.y_sd,
                       np.exp(A.beta - 1.96 * A.se))
    A["hi"] = np.where(A.model == "change", (A.beta + 1.96 * A.se) / A.y_sd,
                       np.exp(A.beta + 1.96 * A.se))
    A["source"] = "7.0 p/mh_p_cbcl.tsv + y/mh_{p,y}_ksads; single-LMM HCP phenotypes"
    A.drop(columns="y_sd").to_csv(OUT, sep="\t", index=False, float_format="%.5g")
    show = A[A.fit == "alone"].pivot_table(index="outcome", columns="brain", values=["est", "p"])
    print(show.round(4).to_string())
    print(OUT)


if __name__ == "__main__":
    main()
