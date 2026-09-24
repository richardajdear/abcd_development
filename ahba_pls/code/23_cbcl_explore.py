"""
23_cbcl_explore.py -- EXPLORATORY: do children whose cortical thinning follows
the PLS2 (C3-like) pattern have more CBCL symptoms?

Data (individual-level, read locally, NEVER written to the repo -- ABCD DUC):
  out/thickness_hcp_70_aa6e91efba82/fits/blups.parquet   per-child random slope
      and intercept BLUPs, 179 LH HCP-MMP parcels, 8,716 children (7.0 tables)
  out/.../fits/fixed.parquet                             population age slope per parcel
  out/.../model_table.parquet                            sex, site, family, ages
  CBCL: 7.0 mh_p_cbcl if present in abcd-data-release-7.0/, else the 5.1
        core/mental-health/mh_p_cbcl.csv (same child IDs; baseline..year 4)
Only coefficient tables are written (results/cbcl_explore_*.tsv).

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

CBCL outcomes (parent report), raw scores log1p-transformed (syndrome T-scores
are floored at 50, so raw scores keep the variance):
  broadband totprob, internal, external; 8 syndromes; 6 DSM-oriented scales;
  pfactor = PC1 of the 8 log1p syndrome scores (fitted at follow-up)
Timepoints: FU = year-3 follow-up (the most complete adolescent wave in 5.1);
  BASE = baseline.  Models (OLS, site fixed effects, family-clustered SE):
  fu      y_FU   ~ brain + sex + age_FU + site
  change  y_FU   ~ brain + y_BASE + sex + age_FU + site
  base    y_BASE ~ brain + sex + age_BASE + site
  each also "+ global_thin" for the spatial phenotypes (pattern beyond global rate)
  clinical (logistic): T >= 65 at FU on totprob / internal / external / DSM depress
Spin null (pls2_proj, c3_proj): partial correlation with the outcome, after the
  same covariates + global_thin, compared with 1,000 projections onto
  spin-rotated versions of the same map (HCP centroids, 137 covered parcels).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls                                                    # noqa: E402

ROOT = HERE.parent; REPO = ROOT.parent; RES = ROOT / "results"
RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
CB70 = sorted((REPO / "abcd-data-release-7.0").rglob("mh_p_cbcl*.tsv"))
CB51 = Path.home() / "Git" / "ABCD" / "abcd-data-release-5.1" / "core" / "mental-health" / "mh_p_cbcl.csv"
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
    "global_thin": thin.mean(1),
    "pls2_proj": thin[cov] @ zP2 / len(cov),
    "pls2_top": thin[top].mean(1),
    "pls2_topbottom": thin[top].mean(1) - thin[bot].mean(1),
    "pls1_proj": thin[cov] @ zP1 / len(cov),
    "pls1_thick_proj": -(inter[cov] @ zP1) / len(cov),
    "c3_proj": thin[cov] @ zC3 / len(cov),
})
BRAIN = list(B.columns)
SPATIAL = [b for b in BRAIN if b != "global_thin"]

mt = pd.read_parquet(RUN / "model_table.parquet",
                     columns=["subject", "visit", "age", "sex", "site", "family_id"]).drop_duplicates(["subject", "visit"])
demo = mt.sort_values("visit").groupby("subject").agg(sex=("sex", "first"), site=("site", "first"),
                                                      family=("family_id", "first"), age_v0=("age", "min"))
B = B.join(demo)

# ---------------------------------------------------------------- CBCL -------
SYN = ["anxdep", "withdep", "somatic", "social", "thought", "attention", "rulebreak", "aggressive"]
BROAD = ["totprob", "internal", "external"]
DSM = ["depress", "anxdisord", "somaticpr", "adhd", "opposit", "conduct"]
if CB70:
    sys.exit(f"7.0 CBCL found ({CB70[0].name}): add the 7.0 column mapping before running on it")
src = "5.1 core/mental-health/mh_p_cbcl.csv"
cols = ([f"cbcl_scr_syn_{s}_r" for s in SYN + BROAD] + [f"cbcl_scr_dsm5_{s}_r" for s in DSM]
        + [f"cbcl_scr_syn_{s}_t" for s in BROAD] + ["cbcl_scr_dsm5_depress_t"])
cb = pd.read_csv(CB51, usecols=["src_subject_id", "eventname", *cols], low_memory=False)
lt = pd.read_csv(CB51.parents[1] / "abcd-general" / "abcd_y_lt.csv",
                 usecols=["src_subject_id", "eventname", "interview_age"])
cb = cb.merge(lt, on=["src_subject_id", "eventname"], how="left")
cb["subject"] = "sub-" + cb.src_subject_id.str.replace("_", "")
EV = {"BASE": "baseline_year_1_arm_1", "FU": "3_year_follow_up_y_arm_1", "FU4": "4_year_follow_up_y_arm_1"}
ren = {f"cbcl_scr_syn_{s}_r": s for s in SYN + BROAD} | {f"cbcl_scr_dsm5_{s}_r": s for s in DSM}
OUT = SYN + BROAD + DSM
W = {}
for k, ev in EV.items():
    d = cb[cb.eventname == ev].set_index("subject")
    d = d[~d.index.duplicated()]
    w = np.log1p(d[list(ren)].rename(columns=ren).apply(pd.to_numeric, errors="coerce"))
    w["age"] = d.interview_age / 12.0
    for s in BROAD + ["depress"]:
        tcol = f"cbcl_scr_syn_{s}_t" if s in BROAD else "cbcl_scr_dsm5_depress_t"
        w[f"{s}_clin"] = (pd.to_numeric(d[tcol], errors="coerce") >= 65).astype(float).where(d[tcol].notna())
    W[k] = w
# p-factor: PC1 of the z-scored log1p syndrome scores at FU, applied to every wave
fu = W["FU"][SYN].dropna(); mu, sd = fu.mean(), fu.std()
_, _, vt = np.linalg.svd(((fu - mu) / sd).to_numpy(), full_matrices=False)
load = pd.Series(vt[0] * np.sign(vt[0].sum()), index=SYN)
for k in W:
    W[k]["pfactor"] = ((W[k][SYN] - mu) / sd) @ load
OUT = OUT + ["pfactor"]
print(f"CBCL source: {src}; FU n = {W['FU'].pfactor.notna().sum()}, "
      f"p-factor loadings {load.round(2).to_dict()}", file=sys.stderr)

# ---------------------------------------------------------------- models -----
D = B.copy()
for k in W:
    D = D.join(W[k].add_suffix(f"_{k}"))
for b in BRAIN:
    D[b] = (D[b] - D[b].mean()) / D[b].std()
D["site"] = D.site.astype(str); D["family"] = D.family.astype(str)

def fit(y, x, extra, kind="ols"):
    f = f"{y} ~ {x} + C(sex) + C(site)" + "".join(f" + {e}" for e in extra)
    d = D.dropna(subset=[y, x, *extra])
    if kind == "logit":
        m = smf.logit(f, d).fit(disp=0, cov_type="cluster", cov_kwds={"groups": d.family})
    else:
        m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": d.family})
    return dict(beta=m.params[x], se=m.bse[x], p=m.pvalues[x], n=int(m.nobs))

rows = []
for o in OUT:
    for b in BRAIN:
        for glob_ in ([False, True] if b in SPATIAL else [False]):
            g = ["global_thin"] if glob_ else []
            for model, y, extra in [("fu", f"{o}_FU", ["age_FU"]),
                                    ("change", f"{o}_FU", ["age_FU", f"{o}_BASE"]),
                                    ("base", f"{o}_BASE", ["age_BASE"]),
                                    ("fu_year4", f"{o}_FU4", ["age_FU4"])]:
                r = fit(y, b, extra + g)
                rows.append(dict(outcome=o, brain=b, adj_global=glob_, model=model, **r))
for o in BROAD + ["depress"]:
    for b in BRAIN:
        for glob_ in ([False, True] if b in SPATIAL else [False]):
            r = fit(f"{o}_clin_FU", b, ["age_FU"] + (["global_thin"] if glob_ else []), kind="logit")
            rows.append(dict(outcome=f"{o}_clin", brain=b, adj_global=glob_, model="fu_logit", **r))
A = pd.DataFrame(rows)
A["y_sd_beta"] = np.nan                       # OLS betas are per SD brain on the log1p scale;
for (o, m), g in A[A.model != "fu_logit"].groupby(["outcome", "model"]):
    ycol = f"{o}_{'BASE' if m == 'base' else 'FU4' if m == 'fu_year4' else 'FU'}"
    A.loc[g.index, "y_sd_beta"] = g.beta / D[ycol].std()     # ... and per SD outcome
A["source"] = src
A.to_csv(RES / "cbcl_explore_assoc.tsv", sep="\t", index=False, float_format="%.5g")

# ---------------------------------------------------------------- spin null --
perm = pls.cached_spin_permutations(list(cov), n_perm=N_SPIN, seed=SEED, path=pls.HCP_CENTROIDS)
T = thin[cov].loc[D.index].to_numpy()
srows = []
for name, zmap in (("pls2_proj", zP2), ("c3_proj", zC3)):
    zm = zmap.to_numpy()
    Pn = T @ np.column_stack([zm] + [zm[perm[i]] for i in range(N_SPIN)]) / len(cov)
    for o in ["totprob", "internal", "external", "pfactor", "depress", "thought", "withdep", "anxdisord", "adhd"]:
        for model in ("fu", "change"):
            y = f"{o}_FU"; extra = ["age_FU"] + ([f"{o}_BASE"] if model == "change" else [])
            ok = D[[y, *extra, "global_thin"]].notna().all(1).to_numpy()
            Xc = pd.get_dummies(D.loc[ok, ["sex", "site"]], drop_first=True).astype(float)
            Xc = np.column_stack([np.ones(ok.sum()), Xc, D.loc[ok, extra + ["global_thin"]].to_numpy()])
            H = lambda v: v - Xc @ np.linalg.lstsq(Xc, v, rcond=None)[0]
            ry = H(D.loc[ok, y].to_numpy()); rP = H(Pn[ok])
            r = (rP * ry[:, None]).sum(0) / np.sqrt((rP ** 2).sum(0) * (ry ** 2).sum())
            srows.append(dict(brain=name, outcome=o, model=model, r_obs=r[0], n=int(ok.sum()),
                              null_mean=r[1:].mean(), null_sd=r[1:].std(),
                              p_spin=(1 + (np.abs(r[1:]) >= abs(r[0])).sum()) / (1 + N_SPIN)))
S = pd.DataFrame(srows); S["source"] = src
S.to_csv(RES / "cbcl_explore_spin.tsv", sep="\t", index=False, float_format="%.5g")

# ---------------------------------------------------------------- summary ----
pd.set_option("display.width", 220, "display.max_rows", 200)
q = A[(A.model == "fu") & A.outcome.isin(BROAD + ["pfactor", "depress", "thought"])]
q = q.assign(cell=q.apply(lambda r: f"{r.y_sd_beta:+.3f} ({r.p:.2g})", axis=1),
             col=q.brain + np.where(q.adj_global, "|g", ""))
print(q.pivot(index="outcome", columns="col", values="cell").to_string())
print(S.round(4).to_string(index=False))
print("brain phenotype correlations:\n", B[BRAIN].corr().round(2).to_string())
