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

    FU = ns["FU"]
    rows = []
    for model_set, xs_list in (("alone", [["global_slope"], ["baseline_thickness"]]),
                               ("joint", [["global_slope", "baseline_thickness"]])):
        for xs in xs_list:
            for o in CBCL:
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
