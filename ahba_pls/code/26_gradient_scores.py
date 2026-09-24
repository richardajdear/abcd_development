"""
26_gradient_scores.py -- EXPLORATORY: does how closely a child's regional
thinning follows the normative gradient predict later symptoms?

Motivation: the case-control (24) and association (25) maps both showed that
children whose anxiety / general symptoms rise thin relatively more where cortex
normally thins least and less where it thins most -- a FLATTER version of the
group pattern. Those are map-level tests on 179 parcel-wise fits; this is the
direct child-level test.

Per child i (bilateral, 179 parcels), with g_r = normative thinning in parcel r
(the group mean of the per-child thinning, = -dCT, positive = faster):
  gradient slope   b_i : OLS slope of the child's thinning on g across parcels,
                         thin_ir = a_i + b_i g_r + e_ir.
                         b = 1 follows the group gradient, b < 1 is flatter, b > 1 steeper.
                         It is separate from the child's global rate (mean_r thin_ir = a_i + b_i mean(g)).
  gradient fit     r_i : Pearson r across parcels between thin_i and g (pattern
                         similarity, independent of scale)
  PLS2 slope       c_i : in the 137 AHBA parcels, thin_ir = a + b g_r + c z(PLS2_r)
                         -- whether the child's thinning follows PLS2 BEYOND the gradient
  (baseline analogues, CT gradient slope: the same regression on first-scan thickness vs the
   group baseline CT map, as a control for a general "how typical is this child's cortex" trait)

Models (OLS, site FE, sex, family-clustered SE; scores standardised):
  CBCL_15-17 ~ score + CBCL_base + sex + site + age_late + age_first + age_span + n_visits
               [+ global thinning] [+ global baseline CT + CT gradient slope]
  outcomes: log1p raw CBCL -- the 7 scales of 25 plus totprob and external.
Caveat: slope BLUPs are shrunk toward the fixed effects, which ARE the group
gradient, so b_i is pulled toward 1 more for children with fewer / closer scans --
hence n_visits and age_span as covariates and a >= 3-scan sensitivity check.
Group-level outputs only: results/gradient_scores_assoc.tsv, gradient_scores_bins.tsv.
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent; REPO = ROOT.parent; RES = ROOT / "results"
RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
R70 = REPO / "abcd-data-release-7.0"

# ------------------------------------------------------------------ brain -----
bl = pd.read_parquet(RUN / "fits" / "blups.parquet")
fx = pd.read_parquet(RUN / "fits" / "fixed.parquet")
thin = -(bl.pivot(index="subject", columns="label", values="re_slope")
         + fx[fx.term == "age_c"].set_index("label").estimate)
bil = lambda X: (lambda Y: Y.set_axis("lh_" + Y.columns, axis=1))(X.T.groupby(X.columns.str.slice(3)).mean().T)
glob_thin = thin.mean(axis=1)
thin = bil(thin)
tt = pd.read_parquet(RUN / "model_table.parquet", columns=["subject", "label", "value", "age"])
ct0 = bil(tt.sort_values("age").groupby(["subject", "label"]).value.first().unstack()).loc[thin.index, thin.columns]
del tt
mt = pd.read_parquet(RUN / "model_table.parquet",
                     columns=["subject", "sex", "site", "family_id", "age_first", "age_span", "n_visits"]
                     ).drop_duplicates("subject").set_index("subject").loc[thin.index]

def slope_on(M, ref):
    x = ref - ref.mean(); Mc = M - M.mean(axis=1).to_numpy()[:, None]
    return (Mc.to_numpy() @ x.to_numpy()) / (x @ x)

g = thin.mean(axis=0)                                         # normative thinning (mm/yr, >0 = faster)
S = pd.DataFrame(index=thin.index)
S["grad_slope"] = slope_on(thin, g)
S["grad_fit"] = [np.corrcoef(row, g)[0, 1] for row in thin.to_numpy()]
S["ct_grad_slope"] = slope_on(ct0, ct0.mean(axis=0))
S["glob_thin"] = glob_thin; S["glob_ct"] = ct0.mean(axis=1)
maps = pd.read_csv(RES / "hcp_summary_maps.csv", index_col=0)
lab = {l.lower(): l for l in thin.columns}; maps.index = [lab[l.lower()] for l in maps.index]
cov = maps.PLS2.dropna().index
Xg = np.column_stack([np.ones(len(cov)), g.loc[cov], (maps.loc[cov, "PLS2"] - maps.loc[cov, "PLS2"].mean()) / maps.loc[cov, "PLS2"].std()])
B = np.linalg.lstsq(Xg, thin[cov].to_numpy().T, rcond=None)[0]
S["pls2_beyond_grad"] = B[2]
print("score summary:\n", S.describe().round(3).to_string(), file=sys.stderr)
print("score correlations:\n", S.corr().round(2).to_string(), file=sys.stderr)
print("grad_slope vs n_visits:", S.grad_slope.groupby(mt.n_visits).mean().round(3).to_dict(), file=sys.stderr)

# ------------------------------------------------------------------ CBCL ------
RAW = {"totprob": None, "internal": "synd__int", "external": "synd__ext", "depress": "dsm__dep",
       "anxdisord": "dsm__anx", "anxdep": "synd__anxdep", "withdep": "synd__wthdep", "rulebreak": "synd__rule"}
SYN = ["anxdep", "wthdep", "som", "soc", "tho", "attn", "rule", "aggr"]
cols = {k: ("mh_p_cbcl_sum" if v is None else f"mh_p_cbcl__{v}_sum") for k, v in RAW.items()}
cols |= {f"s_{s}": f"mh_p_cbcl__synd__{s}_sum" for s in SYN}
cb = pd.read_csv(R70 / "p" / "mh_p_cbcl.tsv", sep="\t", usecols=["participant_id", "session_id", "mh_p_cbcl_age", *set(cols.values())])
cb["subject"] = "sub-NDARINV" + cb.participant_id.str.replace("sub-", "", regex=False)
cb["wave"] = cb.session_id.str.slice(4, 6).astype(int); cb["age"] = cb.mh_p_cbcl_age
cb = cb[cb.subject.isin(thin.index)]
for k, c in cols.items():
    cb[k] = np.log1p(pd.to_numeric(cb[c], errors="coerce"))
S_ = [f"s_{s}" for s in SYN]
lm = cb[cb.wave.isin([5, 6, 7])].groupby("subject")[S_].mean().dropna(); mu, sd = lm.mean(), lm.std()
_, _, vt = np.linalg.svd(((lm - mu) / sd).to_numpy(), full_matrices=False)
cb["pfactor"] = ((cb[S_] - mu) / sd).to_numpy() @ (vt[0] * np.sign(vt[0].sum()))
OUT = list(RAW) + ["pfactor"]
late = cb[cb.wave.isin([5, 6, 7])].groupby("subject")[OUT + ["age"]].mean().add_suffix("_late")
base = cb[cb.wave == 0].groupby("subject")[OUT].mean().add_suffix("_base")

# image quality: FreeSurfer topological defect count (the 7.0 stand-in for the Euler
# number), log1p, averaged over the child's imaging sessions used in the fit
qc = pd.read_csv(R70 / "y" / "mr_y_qc__post__aut.tsv", sep="\t",
                 usecols=["participant_id", "session_id", "mr_y_qc__post__aut__smri__topodfct_count"])
qc["subject"] = "sub-NDARINV" + qc.participant_id.str.replace("sub-", "", regex=False)
qc["v"] = "v" + qc.session_id.str.slice(5, 6)
used = pd.read_parquet(RUN / "model_table.parquet", columns=["subject", "visit"]).drop_duplicates()
qc = qc.merge(used, left_on=["subject", "v"], right_on=["subject", "visit"])
S["qc_defects"] = np.log1p(pd.to_numeric(qc.mr_y_qc__post__aut__smri__topodfct_count, errors="coerce")).groupby(qc.subject).mean()
print("qc coverage:", S.qc_defects.notna().sum(), "| corr with scores:",
      S[["grad_slope", "grad_fit", "glob_thin", "qc_defects"]].corr().qc_defects.round(3).to_dict(), file=sys.stderr)

D = S.join(mt).join(late).join(base)
SC = ["grad_slope", "grad_fit", "pls2_beyond_grad", "ct_grad_slope", "glob_thin"]
for c in SC + ["glob_ct", "qc_defects"]:
    D[c + "_z"] = (D[c] - D[c].mean()) / D[c].std()
D["site"] = D.site.astype(str); D["family"] = D.family_id.astype(str)

ADJ = {"base": [], "+global": ["glob_thin_z"], "+global+CT": ["glob_thin_z", "glob_ct_z", "ct_grad_slope_z"],
       "+global+CT+QC": ["glob_thin_z", "glob_ct_z", "ct_grad_slope_z", "qc_defects_z"]}
rows = []
for sample, dd in (("all", D), (">=3 scans", D[D.n_visits >= 3])):
    for o in OUT:
        for sc in SC:
            for an, extra in ADJ.items():
                if sc in ("glob_thin",) and an != "base":
                    continue
                if sc == "ct_grad_slope" and an.startswith("+global+CT"):
                    continue
                ex = [e for e in extra if e != f"{sc}_z"]
                f = (f"{o}_late ~ {sc}_z + {o}_base + C(sex) + C(site) + age_late + age_first + age_span + n_visits"
                     + "".join(f" + {e}" for e in ex))
                d = dd.dropna(subset=[f"{o}_late", f"{o}_base", "age_late", f"{sc}_z", *ex])
                m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": d.family})
                rows.append(dict(sample=sample, outcome=o, score=sc, adjust=an, beta=m.params[f"{sc}_z"],
                                 se=m.bse[f"{sc}_z"], p=m.pvalues[f"{sc}_z"], n=int(m.nobs),
                                 y_sd_beta=m.params[f"{sc}_z"] / d[f"{o}_late"].std()))
A = pd.DataFrame(rows)
A.to_csv(RES / "gradient_scores_assoc.tsv", sep="\t", index=False, float_format="%.5g")

# dose-response (group-level): mean residualised later score by decile of each
# gradient score, residualised on everything in the fully adjusted model except the score
brows = []
for o in ["rulebreak", "external", "pfactor", "totprob", "anxdisord"]:
    d = D.dropna(subset=[f"{o}_late", f"{o}_base", "age_late"])
    m = smf.ols(f"{o}_late ~ {o}_base + C(sex) + C(site) + age_late + age_first + age_span + n_visits"
                " + glob_thin_z + glob_ct_z + ct_grad_slope_z + qc_defects_z", d).fit()
    res = m.resid / d[f"{o}_late"].std()
    for sc in ("grad_fit", "grad_slope"):
        dec = pd.qcut(d[sc], 10, labels=False) + 1
        for k, v in res.groupby(dec):
            brows.append(dict(outcome=o, score=sc, decile=int(k), score_mid=float(d[sc][dec == k].median()),
                              mean_resid_sd=v.mean(), se=v.std() / np.sqrt(len(v)), n=len(v)))
pd.DataFrame(brows).to_csv(RES / "gradient_scores_bins.tsv", sep="\t", index=False, float_format="%.5g")

pd.set_option("display.width", 220)
q = A[A["sample"] == "all"].assign(c=lambda x: x.apply(lambda r: f"{r.y_sd_beta:+.3f} ({r.p:.2g})", axis=1),
                                   col=lambda x: x.score + " " + x.adjust)
print(q.pivot(index="outcome", columns="col", values="c").to_string())
s3 = A[(A["sample"] == ">=3 scans") & (A.score == "grad_slope")]
print("\n>=3 scans, grad_slope:\n", s3.pivot(index="outcome", columns="adjust", values="y_sd_beta").round(3).to_string())
f3 = A[(A["sample"] == ">=3 scans") & (A.score == "grad_fit")].assign(c=lambda x: x.apply(lambda r: f"{r.y_sd_beta:+.3f} ({r.p:.2g})", axis=1))
print("\n>=3 scans, grad_fit:\n", f3.pivot(index="outcome", columns="adjust", values="c").to_string())
qq = A[(A["sample"] == "all") & (A.adjust == "+global+CT+QC")].assign(c=lambda x: x.apply(lambda r: f"{r.y_sd_beta:+.3f} ({r.p:.2g})", axis=1))
print("\nall, +global+CT+QC:\n", qq.pivot(index="outcome", columns="score", values="c").to_string())
