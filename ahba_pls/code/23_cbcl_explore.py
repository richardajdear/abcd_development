"""
23_cbcl_explore.py -- EXPLORATORY: do children whose cortical thinning follows
the PLS2 (C3-like) pattern have more psychiatric symptoms?

    python code/23_cbcl_explore.py        # 7.0 CBCL + KSADS (default)
    python code/23_cbcl_explore.py 51     # 5.1 CBCL, baseline..year 4 (sensitivity)

Data (individual-level, read locally, NEVER written to the repo -- ABCD DUC):
  out/thickness_hcp_70_aa6e91efba82/fits/{blups,fixed}.parquet  per-child slope and
      intercept BLUPs, 179 LH HCP-MMP parcels, 8,716 children (7.0 imaging)
  out/.../model_table.parquet                                   sex, site, family
  7.0: abcd-data-release-7.0/p/mh_p_cbcl.tsv (ses-00A..07A, ages ~10-17; raw sums),
       y/mh_{p,y}_ksads__dep.tsv, y/mh_p_ksads__psych.tsv (diagnoses; 555 = not given)
  5.1: ABCD/abcd-data-release-5.1/core/mental-health/mh_p_cbcl.csv (baseline..year 4)
Only coefficient tables are written: results/cbcl_explore[_51]_{assoc,spin}.tsv.

Brain phenotypes, one per child, oriented so HIGHER = MORE THINNING (thin = -slope):
  global_thin       mean thinning over all 179 parcels
  pls2_proj         mean over the 137 AHBA parcels of thin_r * z(PLS2_r): a
                    mean-zero spatial contrast (thinning concentrated where PLS2 is high)
  pls2_top          mean thinning in the top-decile PLS2 parcels
  pls2_topbottom    top-decile minus bottom-decile mean thinning
  pls1_proj         thin projected on z(PLS1)                 (static-axis control)
  pls1_thick_proj   baseline-thickness BLUP projected on z(PLS1), sign so that
                    higher = relatively thinner cortex where PLS1 is high
  c3_proj           thin projected on z(C3)                    (published axis)

CBCL (parent report): raw scale sums, log1p (syndrome T-scores are floored at 50):
  3 broadband, 8 syndrome, 6 DSM-oriented scales; pfactor = PC1 of the 8 z-scored
  log1p syndrome scores (fitted at the primary follow-up wave).
Waves -- 7.0: BASE ses-00A, Y3 ses-03A, LATE = mean of ses-05A..07A (ages ~15-17);
         5.1: BASE, Y3, Y4.
Models (OLS, site fixed effects, sex, family-clustered SE; age of the wave):
  base / y3 / y4 / late      y_wave ~ brain
  change                     y_FU ~ brain + y_BASE   (FU = LATE in 7.0, Y3 in 5.1)
  trajectory (7.0)           per-child OLS slope of log1p score on age over all
                             waves (>= 4 waves) ~ brain + age_BASE
  each spatial phenotype also "+ global_thin" (pattern beyond the global rate)
  clinical (logistic): T >= 65 at FU on total / internalising / externalising / DSM dep
  KSADS (logistic, 7.0): lifetime (any wave, present or past) MDD by youth report,
     MDD by parent report, and parent-reported psychosis spectrum (attenuated
     psychosis, schizophrenia, schizophreniform, other psychotic) ~ brain + age_last
Spin null (pls2_proj, c3_proj): partial correlation with the outcome after the same
  covariates + global_thin, against 1,000 projections onto spin-rotated maps.
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls                                                    # noqa: E402

SRC = sys.argv[1] if len(sys.argv) > 1 else "70"
assert SRC in ("70", "51")
ROOT = HERE.parent; REPO = ROOT.parent; RES = ROOT / "results"
RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
R70 = REPO / "abcd-data-release-7.0"
R51 = Path.home() / "Git" / "ABCD" / "abcd-data-release-5.1" / "core"
TAG = "cbcl_explore" if SRC == "70" else "cbcl_explore_51"
N_SPIN, SEED = 1000, 0

# ---------------------------------------------------------------- brain ------
bl = pd.read_parquet(RUN / "fits" / "blups.parquet")
fx = pd.read_parquet(RUN / "fits" / "fixed.parquet")
fslope = fx[fx.term == "age_c"].set_index("label").estimate
fint = fx[fx.term == "(Intercept)"].set_index("label").estimate
slope = bl.pivot(index="subject", columns="label", values="re_slope") + fslope   # mm/yr
inter = bl.pivot(index="subject", columns="label", values="re_intercept") + fint
thin = -slope

maps = pd.read_csv(RES / "hcp_summary_maps.csv", index_col=0)           # 3d_ds5 fit
lab = {l.lower(): l for l in thin.columns}
maps.index = [lab[l.lower()] for l in maps.index]
cov = maps.dropna(subset=["PLS2"]).index                                 # 137 parcels
z = lambda s: (s - s.mean()) / s.std()
zP2, zP1, zC3 = z(maps.loc[cov, "PLS2"]), z(maps.loc[cov, "PLS1"]), z(maps.loc[cov, "C3"])
top = maps.loc[cov, "PLS2"].nlargest(int(round(0.1 * len(cov)))).index
bot = maps.loc[cov, "PLS2"].nsmallest(int(round(0.1 * len(cov)))).index
B = pd.DataFrame({
    "global_thin": thin.mean(axis=1),
    "pls2_proj": thin[cov] @ zP2 / len(cov),
    "pls2_top": thin[top].mean(axis=1),
    "pls2_topbottom": thin[top].mean(axis=1) - thin[bot].mean(axis=1),
    "pls1_proj": thin[cov] @ zP1 / len(cov),
    "pls1_thick_proj": -(inter[cov] @ zP1) / len(cov),
    "c3_proj": thin[cov] @ zC3 / len(cov),
})
BRAIN = list(B.columns)
SPATIAL = [b for b in BRAIN if b != "global_thin"]
mt = pd.read_parquet(RUN / "model_table.parquet",
                     columns=["subject", "visit", "age", "sex", "site", "family_id"]).drop_duplicates(["subject", "visit"])
demo = mt.sort_values("visit").groupby("subject").agg(sex=("sex", "first"), site=("site", "first"),
                                                      family=("family_id", "first"))
B = B.join(demo)

# ---------------------------------------------------------------- CBCL -------
SYN = ["anxdep", "withdep", "somatic", "social", "thought", "attention", "rulebreak", "aggressive"]
BROAD = ["totprob", "internal", "external"]
DSM = ["depress", "anxdisord", "somaticpr", "adhd", "opposit", "conduct"]
OUT = SYN + BROAD + DSM
CLIN = ["totprob", "internal", "external", "depress"]

if SRC == "70":
    k70 = dict(anxdep="synd__anxdep", withdep="synd__wthdep", somatic="synd__som", social="synd__soc",
               thought="synd__tho", attention="synd__attn", rulebreak="synd__rule", aggressive="synd__aggr",
               internal="synd__int", external="synd__ext", depress="dsm__dep", anxdisord="dsm__anx",
               somaticpr="dsm__somat", adhd="dsm__adhd", opposit="dsm__opp", conduct="dsm__cond")
    raw = {k: f"mh_p_cbcl__{v}_sum" for k, v in k70.items()} | {"totprob": "mh_p_cbcl_sum"}
    tsc = {"totprob": "mh_p_cbcl_tscore", "internal": "mh_p_cbcl__synd__int_tscore",
           "external": "mh_p_cbcl__synd__ext_tscore", "depress": "mh_p_cbcl__dsm__dep_tscore"}
    cb = pd.read_csv(R70 / "p" / "mh_p_cbcl.tsv", sep="\t",
                     usecols=["participant_id", "session_id", "mh_p_cbcl_age", *raw.values(), *tsc.values()])
    cb["subject"] = "sub-NDARINV" + cb.participant_id.str.replace("sub-", "", regex=False)
    cb["wave"] = cb.session_id.str.slice(4, 6).astype(int)            # ses-06A -> 6
    cb["age"] = cb.mh_p_cbcl_age
    src = "7.0 p/mh_p_cbcl.tsv"
    WAVES = {"BASE": [0], "Y3": [3], "LATE": [5, 6, 7]}
    FU = "LATE"
else:
    raw = {s: f"cbcl_scr_syn_{s}_r" for s in SYN + BROAD} | {s: f"cbcl_scr_dsm5_{s}_r" for s in DSM}
    tsc = {s: f"cbcl_scr_syn_{s}_t" for s in BROAD} | {"depress": "cbcl_scr_dsm5_depress_t"}
    cb = pd.read_csv(R51 / "mental-health" / "mh_p_cbcl.csv", low_memory=False,
                     usecols=["src_subject_id", "eventname", *raw.values(), *tsc.values()])
    lt = pd.read_csv(R51 / "abcd-general" / "abcd_y_lt.csv", usecols=["src_subject_id", "eventname", "interview_age"])
    cb = cb.merge(lt, on=["src_subject_id", "eventname"], how="left")
    cb["subject"] = "sub-" + cb.src_subject_id.str.replace("_", "")
    EVN = {"baseline_year_1_arm_1": 0, "1_year_follow_up_y_arm_1": 1, "2_year_follow_up_y_arm_1": 2,
           "3_year_follow_up_y_arm_1": 3, "4_year_follow_up_y_arm_1": 4}
    cb["wave"] = cb.eventname.map(EVN); cb["age"] = cb.interview_age / 12.0
    src = "5.1 core/mental-health/mh_p_cbcl.csv"
    WAVES = {"BASE": [0], "Y3": [3], "Y4": [4]}
    FU = "Y3"

cb = cb[cb.subject.isin(B.index)].drop_duplicates(["subject", "wave"])
for k, c in raw.items():
    cb[k] = np.log1p(pd.to_numeric(cb[c], errors="coerce"))
for k, c in tsc.items():
    t = pd.to_numeric(cb[c], errors="coerce")
    cb[f"{k}_clin"] = (t >= 65).astype(float).where(t.notna())
W = {}
for name, waves in WAVES.items():
    g = cb[cb.wave.isin(waves)].groupby("subject")
    w = g[OUT + ["age"]].mean()                                 # multi-wave: mean of available waves
    w[[f"{k}_clin" for k in CLIN]] = g[[f"{k}_clin" for k in CLIN]].max()   # ... any clinical-range wave
    W[name] = w
fu = W[FU][SYN].dropna(); mu, sd = fu.mean(), fu.std()
_, _, vt = np.linalg.svd(((fu - mu) / sd).to_numpy(), full_matrices=False)
load = pd.Series(vt[0] * np.sign(vt[0].sum()), index=SYN)
for k in W:
    W[k]["pfactor"] = ((W[k][SYN] - mu) / sd) @ load
cb["pfactor"] = ((cb[SYN] - mu) / sd) @ load
OUTP = OUT + ["pfactor"]

D = B.copy()
for k in W:
    D = D.join(W[k].add_suffix(f"_{k}"))
if SRC == "70":                                                 # per-child symptom trajectory
    def traj(g):
        if g.age.notna().sum() < 4:
            return pd.Series(np.nan, index=OUTP)
        a = g.age - g.age.mean()
        return (g[OUTP].mul(a, axis=0).sum() / (a ** 2).sum()).where(g[OUTP].notna().sum() >= 4)
    TR = cb.groupby("subject")[["age", *OUTP]].apply(traj)
    D = D.join(TR.add_suffix("_TRAJ"))
    # KSADS lifetime diagnoses (1 at any wave, present or past; 555 = not administered)
    def ever(path, cols):
        k = pd.read_csv(path, sep="\t", usecols=["participant_id", "session_id", *cols], dtype=str)
        v = k[cols].apply(pd.to_numeric, errors="coerce").where(lambda x: x.isin([0, 1]))
        k["subject"] = "sub-NDARINV" + k.participant_id.str.replace("sub-", "", regex=False)
        k["any"] = v.max(axis=1); k["seen"] = v.notna().any(axis=1)
        s = k[k.seen].groupby("subject")["any"].max()
        return s
    Y = R70 / "y"
    D["mdd_youth_DX"] = ever(Y / "mh_y_ksads__dep.tsv", ["mh_y_ksads__dep__mdd__pres_dx", "mh_y_ksads__dep__mdd__past_dx"])
    D["mdd_parent_DX"] = ever(Y / "mh_p_ksads__dep.tsv", ["mh_p_ksads__dep__mdd__pres_dx", "mh_p_ksads__dep__mdd__past_dx"])
    D["psychosis_parent_DX"] = ever(Y / "mh_p_ksads__psych.tsv",
        ["mh_p_ksads__psych__aps__pres_dx", "mh_p_ksads__psych__aps__past_dx", "mh_p_ksads__psych__schiz__pres_dx",
         "mh_p_ksads__psych__schizform__pres_dx", "mh_p_ksads__psych__schiz__oth__pres_dx",
         "mh_p_ksads__psych__schiz__oth__past_dx", "mh_p_ksads__psych__oth__pres_dx"])
    D["age_LAST"] = cb.groupby("subject").age.max()
    DX = ["mdd_youth_DX", "mdd_parent_DX", "psychosis_parent_DX"]
else:
    DX = []
for b in BRAIN:
    D[b] = (D[b] - D[b].mean()) / D[b].std()
D["site"] = D.site.astype(str); D["family"] = D.family.astype(str)
print(f"CBCL source: {src}; n per wave {({k: int(D[f'pfactor_{k}'].notna().sum()) for k in W})}; "
      f"p-factor loadings {load.round(2).to_dict()}"
      + (f"; KSADS cases {({d: int(D[d].sum()) for d in DX})}" if DX else ""), file=sys.stderr)

def fit(y, x, extra, kind="ols"):
    f = f"{y} ~ {x} + C(sex) + C(site)" + "".join(f" + {e}" for e in extra)
    d = D.dropna(subset=[y, x, *extra])
    fitter = smf.logit if kind == "logit" else smf.ols
    kw = dict(disp=0) if kind == "logit" else {}
    m = fitter(f, d).fit(cov_type="cluster", cov_kwds={"groups": d.family}, **kw)
    return dict(beta=m.params[x], se=m.bse[x], p=m.pvalues[x], n=int(m.nobs),
                n_cases=int(d[y].sum()) if kind == "logit" else np.nan)

MODELS = [(k.lower(), k, [f"age_{k}"]) for k in W] + [("change", FU, [f"age_{FU}", "BASE"])]
if SRC == "70":
    MODELS += [("trajectory", "TRAJ", ["age_BASE"])]
rows = []
for b in BRAIN:
    for g in ([[], ["global_thin"]] if b in SPATIAL else [[]]):
        for o in OUTP:
            for model, wave, extra in MODELS:
                ex = [f"{o}_BASE" if e == "BASE" else e for e in extra]
                rows.append(dict(outcome=o, brain=b, adj_global=bool(g), model=model,
                                 **fit(f"{o}_{wave}", b, ex + g)))
        for o in CLIN:
            rows.append(dict(outcome=f"{o}_clin", brain=b, adj_global=bool(g), model="clinical_logit",
                             **fit(f"{o}_clin_{FU}", b, [f"age_{FU}"] + g, kind="logit")))
        for d_ in DX:
            rows.append(dict(outcome=d_, brain=b, adj_global=bool(g), model="ksads_logit",
                             **fit(d_, b, ["age_LAST"] + g, kind="logit")))
A = pd.DataFrame(rows)
ycol = {m: w for m, w, _ in MODELS}
A["y_sd_beta"] = [r.beta / D[f"{r.outcome}_{ycol[r.model]}"].std() if r.model in ycol else np.nan
                  for r in A.itertuples()]
A["source"] = src
A.to_csv(RES / f"{TAG}_assoc.tsv", sep="\t", index=False, float_format="%.5g")

# ---------------------------------------------------------------- spin null --
perm = pls.cached_spin_permutations(list(cov), n_perm=N_SPIN, seed=SEED, path=pls.HCP_CENTROIDS)
T = thin[cov].loc[D.index].to_numpy()
srows = []
spin_models = [(FU.lower(), FU, [f"age_{FU}"]), ("change", FU, [f"age_{FU}", "BASE"])]
if SRC == "70":
    spin_models.append(("trajectory", "TRAJ", ["age_BASE"]))
for name, zmap in (("pls2_proj", zP2), ("c3_proj", zC3)):
    zm = zmap.to_numpy()
    Pn = T @ np.column_stack([zm] + [zm[perm[i]] for i in range(N_SPIN)]) / len(cov)
    for o in ["totprob", "internal", "external", "pfactor", "depress", "thought", "withdep", "anxdisord", "adhd"]:
        for model, wave, extra in spin_models:
            y = f"{o}_{wave}"; ex = [f"{o}_BASE" if e == "BASE" else e for e in extra]
            ok = D[[y, *ex, "global_thin"]].notna().all(axis=1).to_numpy()
            Xc = pd.get_dummies(D.loc[ok, ["sex", "site"]], drop_first=True).astype(float)
            Xc = np.column_stack([np.ones(ok.sum()), Xc, D.loc[ok, ex + ["global_thin"]].to_numpy()])
            H = lambda v: v - Xc @ np.linalg.lstsq(Xc, v, rcond=None)[0]
            ry = H(D.loc[ok, y].to_numpy()); rP = H(Pn[ok])
            r = (rP * ry[:, None]).sum(0) / np.sqrt((rP ** 2).sum(0) * (ry ** 2).sum())
            srows.append(dict(brain=name, outcome=o, model=model, r_obs=r[0], n=int(ok.sum()),
                              null_mean=r[1:].mean(), null_sd=r[1:].std(),
                              p_spin=(1 + (np.abs(r[1:]) >= abs(r[0])).sum()) / (1 + N_SPIN)))
S = pd.DataFrame(srows); S["source"] = src
S.to_csv(RES / f"{TAG}_spin.tsv", sep="\t", index=False, float_format="%.5g")

from scipy.stats import false_discovery_control
A["q"] = false_discovery_control(A.p)
print(f"{len(A)} fits | p<0.05: {(A.p < 0.05).sum()} | min q {A.q.min():.3f} | spin p<0.05: "
      f"{(S.p_spin < 0.05).sum()}/{len(S)} (min {S.p_spin.min():.3f})")
print(A.nsmallest(15, "p")[["outcome", "brain", "adj_global", "model", "y_sd_beta", "beta", "p", "q", "n", "n_cases"]]
      .round(4).to_string(index=False))
