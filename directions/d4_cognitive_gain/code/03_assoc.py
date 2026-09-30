"""D4 (docs/DIRECTIONS.md): does the individual thinning rate predict longitudinal
cognitive gain?  Tests prediction T5 (faster thinning -> smaller gain, given baseline).

Brain predictors -- the Figure 1 phenotypes (single LMM on the per-scan HCP-MMP cortical
mean, 02_fit_1lmm.R), z-scored over the 8,716 imaging children:
    global_slope        random slope; NEGATIVE = FASTER thinning (Figure 1 sign)
    baseline_thickness  random intercept (thickness at the sample-mean age)
So T5 predicts beta > 0 for global_slope: slower thinning, larger gain.

Outcomes -- NIH Toolbox uncorrected standard scores (nc_y_nihtb). The fluid and total
composites exist only at baseline (ses-00A) and year 6 (ses-06A): Card Sort was not
administered at years 2 or 4, List Sorting not at year 2. Primary interval is therefore baseline ->
year 6 for every measure; year 4 is a secondary interval for the tasks given then.

Models (OLS, site fixed effects, family-clustered SE; beta in SD of y_late per SD brain):
  M1 ANCOVA     y_late ~ brain + y_base + age_base + age_late + sex + site
  M2 + SES      M1 + caregiver education (5 lvl) + household income (6 lvl + missing)
                  + ADI national percentile (baseline address)
  M3 + joint    M2 + baseline_thickness + log1p(mean topological defects)  [slope only]
  DIFF          (y_late - y_base) ~ brain + age_base + age_late + sex + site  (Lord check)
  LEVEL_BASE    y_base ~ brain + age_base + sex + site   (is the slope related to level?)
  SEX_INT       M1 + brain x sex
Multiple testing: BH-FDR within each model over the 10 outcomes (3 composites + 7 tasks).

Run from repo root (env abcd-spatial):
    python directions/d4_cognitive_gain/code/03_assoc.py
Writes results/d4_assoc.tsv, results/d4_sample.tsv (aggregates only).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[3]
R70 = REPO / "abcd-data-release-7.0"
D4 = REPO / "directions/d4_cognitive_gain"
W, RES = D4 / "work", D4 / "results"

COMPOSITES = {"fluid": "comp__fluid__uncor", "cryst": "comp__cryst__uncor", "total": "comp__tot__uncor"}
TASKS = {"flanker": "flnkr__uncor", "pattern": "pttcp__uncor", "picseq": "picsq__uncor",
         "picvocab": "picvcb__uncor", "reading": "readr__uncor", "cardsort": "crdst__uncor",
         "listsort": "lswmt__uncor"}
AGECORR = {"fluid": "comp__fluid__agecorr", "cryst": "comp__cryst__agecorr", "total": "comp__tot__agecor"}
OUTCOMES = {**COMPOSITES, **TASKS}
BASE, LATE, MID = "ses-00A", "ses-06A", "ses-04A"


def key(s: pd.Series) -> pd.Series:
    return s.str.extract(r"([A-Z0-9]{8})$")[0]


def load() -> pd.DataFrame:
    z = lambda s: (s - s.mean()) / s.std()
    b = pd.read_csv(W / "hcp70_1lmm_blups.csv")
    b["k"] = key(b.subject)
    D = b.set_index("k")[["re_slope", "re_intercept", "slope_reliability"]]
    D["global_slope"], D["baseline_thickness"] = z(D.re_slope), z(D.re_intercept)

    mt = pd.read_parquet(REPO / "out/thickness_hcp_70_aa6e91efba82/model_table.parquet",
                         columns=["subject", "visit", "sex", "site", "family_id", "n_visits"])
    mt["k"] = key(mt.subject)
    dem = mt.sort_values("visit").groupby("k").agg(sex=("sex", "first"), site=("site", "first"),
                                                    family=("family_id", "first"), n_scans=("n_visits", "first"))
    D = D.join(dem)

    # NIH Toolbox, uncorrected (and age-corrected composites for a sensitivity)
    cols = {f"nc_y_nihtb__{v}_score": n for n, v in OUTCOMES.items()}
    cols |= {f"nc_y_nihtb__{v}_score": f"{n}_ac" for n, v in AGECORR.items()}
    nt = pd.read_csv(R70 / "nc_y_nihtb.tsv", sep="\t", usecols=["participant_id", "session_id", *cols])
    nt = nt.rename(columns=cols)
    nt["k"] = key(nt.participant_id)
    age = pd.read_csv(R70 / "g/dyn/ab_g_dyn.tsv", sep="\t",
                      usecols=["participant_id", "session_id", "ab_g_dyn__visit_age",
                               "ab_g_dyn__cohort_edu__cgs", "ab_g_dyn__cohort_income__hhold__6lvl"])
    age["k"] = key(age.participant_id)
    nt = nt.merge(age[["k", "session_id", "ab_g_dyn__visit_age"]], on=["k", "session_id"], how="left")
    for ses, tag in ((BASE, "BASE"), (LATE, "LATE"), (MID, "MID")):
        w = nt[nt.session_id == ses].set_index("k")
        D = D.join(w[[*OUTCOMES, *[f"{n}_ac" for n in AGECORR], "ab_g_dyn__visit_age"]]
                   .add_suffix(f"_{tag}"), how="left")
        D = D.rename(columns={f"ab_g_dyn__visit_age_{tag}": f"age_{tag}"})

    # SES at baseline
    s0 = age[age.session_id == BASE].set_index("k")
    D["edu"] = s0.ab_g_dyn__cohort_edu__cgs.reindex(D.index)
    inc = s0.ab_g_dyn__cohort_income__hhold__6lvl.reindex(D.index)
    D["income"] = inc.where(inc.between(1, 6)).fillna(0).astype(int).astype(str)  # 0 = missing/declined
    adi = pd.read_csv(R70 / "le_l_adi.tsv", sep="\t",
                      usecols=["participant_id", "session_id", "le_l_adi__addr1__national_prcnt"])
    adi["k"] = key(adi.participant_id)
    D["adi"] = adi[adi.session_id == BASE].set_index("k").le_l_adi__addr1__national_prcnt.reindex(D.index)

    # image quality: mean log1p topological defects over the child's imaging sessions
    q = pd.read_csv(R70 / "y/mr_y_qc__post__aut.tsv", sep="\t",
                    usecols=["participant_id", "session_id", "mr_y_qc__post__aut__smri__topodfct_count"])
    q["k"] = key(q.participant_id)
    sess = mt.assign(session_id=mt.visit.str.replace("v", "ses-0", regex=False) + "A")[["k", "session_id"]]
    q = sess.merge(q, on=["k", "session_id"], how="left")
    D["topo"] = np.log1p(q.groupby("k").mr_y_qc__post__aut__smri__topodfct_count.mean()).reindex(D.index)

    D["site"], D["family"] = D.site.astype(str), D.family.astype(str)
    D["edu"] = D.edu.astype("Int64").astype(str)
    return D


def fit(D, y, xs, extra, cats=("sex", "site")):
    f = f"{y} ~ " + " + ".join(xs) + "".join(f" + C({c})" for c in cats) + "".join(f" + {e}" for e in extra)
    need = [y, *[x.split(":")[0] for x in xs], *[e for e in extra if not e.startswith("C(")], *cats]
    need += [e[2:-1] for e in extra if e.startswith("C(")]
    d = D.dropna(subset=list(dict.fromkeys(need)))
    return smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": d.family}), d


def main() -> None:
    D = load()
    assert len(D) == 8716 and D.global_slope.notna().all()
    rows = []
    ses_terms = ["C(edu)", "C(income)", "adi"]

    def add(model, o, m, d, x, ysd, interval):
        rows.append(dict(model=model, outcome=o, interval=interval, term=x, beta=m.params[x],
                         se=m.bse[x], p=m.pvalues[x], n=int(m.nobs), n_families=d.family.nunique(),
                         est=m.params[x] / ysd, lo=(m.params[x] - 1.96 * m.bse[x]) / ysd,
                         hi=(m.params[x] + 1.96 * m.bse[x]) / ysd))

    for o in OUTCOMES:
        yb, yl = f"{o}_BASE", f"{o}_LATE"
        ysd = D[yl].std()
        for x in ("global_slope", "baseline_thickness"):
            m, d = fit(D, yl, [x], [yb, "age_BASE", "age_LATE"]);                add("M1", o, m, d, x, ysd, "0-6")
            m, d = fit(D, yl, [x], [yb, "age_BASE", "age_LATE", *ses_terms]);   add("M2", o, m, d, x, ysd, "0-6")
            D["_diff"] = D[yl] - D[yb]
            m, d = fit(D, "_diff", [x], ["age_BASE", "age_LATE"]);               add("DIFF", o, m, d, x, ysd, "0-6")
            m, d = fit(D, yb, [x], ["age_BASE"]);                                add("LEVEL_BASE", o, m, d, x, D[yb].std(), "0")
        m, d = fit(D, yl, ["global_slope", "baseline_thickness"],
                   [yb, "age_BASE", "age_LATE", *ses_terms, "topo"])
        add("M3", o, m, d, "global_slope", ysd, "0-6")
        m, d = fit(D, yl, ["global_slope", "global_slope:C(sex)"], [yb, "age_BASE", "age_LATE"])
        t = [k for k in m.params.index if k.startswith("global_slope:")][0]
        add("SEX_INT", o, m, d, t, ysd, "0-6")
        if D[f"{o}_MID"].notna().sum() > 3000:                                  # year-4 interval
            ym = f"{o}_MID"
            m, d = fit(D, ym, ["global_slope"], [yb, "age_BASE", "age_MID"])
            add("M1_Y4", o, m, d, "global_slope", D[ym].std(), "0-4")
    for o in AGECORR:                                                            # age-corrected scores
        yb, yl = f"{o}_ac_BASE", f"{o}_ac_LATE"
        m, d = fit(D, yl, ["global_slope"], [yb, "age_BASE", "age_LATE"])
        add("M1_AGECORR", o, m, d, "global_slope", D[yl].std(), "0-6")
    # restriction: children with >= 3 scans (more reliable slopes)
    D3 = D[D.n_scans >= 3]
    for o in COMPOSITES:
        m, d = fit(D3, f"{o}_LATE", ["global_slope"], [f"{o}_BASE", "age_BASE", "age_LATE"])
        add("M1_3SCANS", o, m, d, "global_slope", D[f"{o}_LATE"].std(), "0-6")

    A = pd.DataFrame(rows)
    A["tier"] = np.where(A.outcome.isin(list(COMPOSITES)), "composite", "task")
    A["p_fdr"] = np.nan
    for (mod, term), g in A.groupby(["model", "term"]):
        A.loc[g.index, "p_fdr"] = multipletests(g.p, method="fdr_bh")[1]
    RES.mkdir(exist_ok=True)
    A.to_csv(RES / "d4_assoc.tsv", sep="\t", index=False, float_format="%.5g")

    # sample description (aggregates only)
    S = pd.DataFrame([dict(
        n_imaging=len(D), n_fluid_both=int(D[["fluid_BASE", "fluid_LATE"]].notna().all(1).sum()),
        n_cryst_both=int(D[["cryst_BASE", "cryst_LATE"]].notna().all(1).sum()),
        age_base_mean=D.age_BASE.mean(), age_late_mean=D.age_LATE.mean(),
        interval_mean=(D.age_LATE - D.age_BASE).mean(),
        fluid_gain_mean=(D.fluid_LATE - D.fluid_BASE).mean(), fluid_gain_sd=(D.fluid_LATE - D.fluid_BASE).std(),
        cryst_gain_mean=(D.cryst_LATE - D.cryst_BASE).mean(), cryst_gain_sd=(D.cryst_LATE - D.cryst_BASE).std(),
        r_fluid_base_late=D.fluid_BASE.corr(D.fluid_LATE), r_cryst_base_late=D.cryst_BASE.corr(D.cryst_LATE),
        slope_reliability_median=D.slope_reliability.median(),
        r_slope_fluid_gain=D.global_slope.corr(D.fluid_LATE - D.fluid_BASE),
        r_slope_topo=D.global_slope.corr(D.topo), r_slope_adi=D.global_slope.corr(D.adi))]).T
    S.columns = ["value"]
    S.to_csv(RES / "d4_sample.tsv", sep="\t", float_format="%.4g")
    show = A[(A.term == "global_slope") & A.model.isin(["M1", "M2", "M3", "DIFF", "LEVEL_BASE"])]
    print(show.pivot_table(index="outcome", columns="model", values="est").round(3).to_string())
    print(show.pivot_table(index="outcome", columns="model", values="p").round(4).to_string())
    print(S.round(3).to_string())


if __name__ == "__main__":
    main()
