"""
24_casecontrol_maps.py -- EXPLORATORY: regional case-control maps of adolescent
cortical thinning for children who DEVELOPED symptoms vs those who did not, and
their spatial correlation with the PLS2 map.

Hypothesis (user): the case-control difference maps correlate with PLS2.

Symptom definitions -- the CBCL / KSADS outcomes that 23_cbcl_explore.py found
associated with GLOBAL thinning in the developmental models (change, years 5-7,
symptom slope), 7.0 tables:
  CBCL T-score scales: case = below threshold at baseline (ses-00A) AND at/above it
    at any of ses-05A..07A (ages ~15-17); control = below threshold at EVERY wave.
    Inclusion = p < 0.05 for global thinning in any developmental model of
    results/cbcl_explore_assoc.tsv (so opposit, significant only at baseline, is out).
    Thresholds: Achenbach borderline cut-offs -- T >= 65 for syndrome / DSM scales,
    T >= 60 for the internalising broadband scale.
      depress   DSM depressive problems      mh_p_cbcl__dsm__dep_tscore
      withdep   withdrawn/depressed          mh_p_cbcl__synd__wthdep_tscore
      anxdep    anxious/depressed            mh_p_cbcl__synd__anxdep_tscore
      anxdisord DSM anxiety problems         mh_p_cbcl__dsm__anx_tscore
      internal  internalising                mh_p_cbcl__synd__int_tscore
      rulebreak rule-breaking                mh_p_cbcl__synd__rule_tscore
  pfactor (no T-score): PC1 of the 8 log1p syndrome raw sums (as in 23); case =
    top decile at any of 05A..07A and not top decile at baseline; control = below
    the 75th percentile at every wave.
  mdd_parent: KSADS parent-report MDD, lifetime (present or past dx at any wave)
    vs never (among children with >= 1 administered interview).

Brain: per-child regional thinning thin_ir = -(fixed slope_r + random slope_ir)
(mm/yr) from the HCP-MMP fit (out/thickness_hcp_70_aa6e91efba82), averaged over
hemispheres to 179 bilateral parcels as the group maps are; global = mean of all 358.

Map, per parcel r (OLS, all parcels at once):
  thin_ir ~ case_i + sex + site (FE) + age_first + n_visits + age_span  [+ global_i]
  d_r = beta_case / residual SD  (Cohen's d, positive = cases thin FASTER)
  "absolute" map without global_i; "relative" map with global thinning as a
  covariate (where cases thin more than their own global rate predicts).
  n_visits and age_span are covariates because slope BLUPs are shrunk more for
  children with fewer / closer scans, and scan count may differ by case status.

Map-level tests (Spearman rho):
  vs PLS2 (137 AHBA parcels), AHBA C3, and the normative thinning map dCT
  (179 parcels; a proportional "cases thin more everywhere" effect would make the
  absolute map mirror dCT, which itself correlates with PLS2).
  p_spin  : 5,000 spin rotations (HCP centroids)
  p_label : 1,000 permutations of case/control labels within site -- is the
            observed map more PLS2-like than maps from random groups of the same size?
  partial : rho(map, PLS2 | dCT), both ranked.

Only group-level summaries (per-parcel d, t, p; map correlations) are written.
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls                                                    # noqa: E402

ROOT = HERE.parent; REPO = ROOT.parent; RES = ROOT / "results"
RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
R70 = REPO / "abcd-data-release-7.0"
N_SPIN, N_LABEL, SEED = 5000, 1000, 0

# ------------------------------------------------------------------ brain -----
bl = pd.read_parquet(RUN / "fits" / "blups.parquet")
fx = pd.read_parquet(RUN / "fits" / "fixed.parquet")
thin = -(bl.pivot(index="subject", columns="label", values="re_slope")
         + fx[fx.term == "age_c"].set_index("label").estimate)
glob = thin.mean(axis=1)                                      # both hemispheres, 358 parcels
# bilateral per child (lh/rh mean), as the group Y maps and PLS2 are; labelled lh_*
thin = thin.T.groupby(thin.columns.str.slice(3)).mean().T
thin.columns = "lh_" + thin.columns
assert thin.shape[1] == 179, thin.shape
mt = pd.read_parquet(RUN / "model_table.parquet",
                     columns=["subject", "visit", "sex", "site", "age_first", "age_span", "n_visits"]
                     ).drop_duplicates("subject").set_index("subject")
C = mt.loc[thin.index, ["sex", "site", "age_first", "age_span", "n_visits"]]

maps = pd.read_csv(RES / "hcp_summary_maps.csv", index_col=0)
lab = {l.lower(): l for l in thin.columns}
maps.index = [lab[l.lower()] for l in maps.index]
REF = {"PLS2": maps.PLS2.dropna(), "C3": maps.C3.dropna(), "dCT": maps.dCT.dropna(), "PLS1": maps.PLS1.dropna()}

# ------------------------------------------------------------------ cases -----
T = {"depress": "mh_p_cbcl__dsm__dep_tscore", "withdep": "mh_p_cbcl__synd__wthdep_tscore",
     "anxdep": "mh_p_cbcl__synd__anxdep_tscore",
     "anxdisord": "mh_p_cbcl__dsm__anx_tscore", "internal": "mh_p_cbcl__synd__int_tscore",
     "rulebreak": "mh_p_cbcl__synd__rule_tscore"}
THR = {"internal": 60}
SYN = ["anxdep", "wthdep", "som", "soc", "tho", "attn", "rule", "aggr"]
cb = pd.read_csv(R70 / "p" / "mh_p_cbcl.tsv", sep="\t",
                 usecols=["participant_id", "session_id", *T.values(), *[f"mh_p_cbcl__synd__{s}_sum" for s in SYN]])
cb["subject"] = "sub-NDARINV" + cb.participant_id.str.replace("sub-", "", regex=False)
cb["wave"] = cb.session_id.str.slice(4, 6).astype(int)
cb = cb[cb.subject.isin(thin.index)]
for c in [*T.values(), *[f"mh_p_cbcl__synd__{s}_sum" for s in SYN]]:
    cb[c] = pd.to_numeric(cb[c], errors="coerce")
# p-factor as in 23: PC1 of z-scored log1p syndrome sums, fitted on years 5-7 means
L = np.log1p(cb[[f"mh_p_cbcl__synd__{s}_sum" for s in SYN]])
late_mean = L[cb.wave.isin([5, 6, 7])].groupby(cb.subject[cb.wave.isin([5, 6, 7])]).mean().dropna()
mu, sd = late_mean.mean(), late_mean.std()
_, _, vt = np.linalg.svd(((late_mean - mu) / sd).to_numpy(), full_matrices=False)
w = vt[0] * np.sign(vt[0].sum())
cb["pfactor"] = ((L - mu) / sd).to_numpy() @ w

def define(col, thr, late=(5, 6, 7), ctrl_thr=None):
    d = cb[["subject", "wave", col]].dropna()
    base = d[d.wave == 0].set_index("subject")[col]
    lmax = d[d.wave.isin(late)].groupby("subject")[col].max()
    amax = d.groupby("subject")[col].max()
    ctrl_thr = thr if ctrl_thr is None else ctrl_thr
    case = base.index[(base < thr)].intersection(lmax.index[lmax >= thr])
    ctrl = base.index.intersection(lmax.index).intersection(amax.index[amax < ctrl_thr])
    return pd.Series(1.0, index=case), pd.Series(0.0, index=ctrl)

DEF = {}
for k, c in T.items():
    DEF[k] = define(c, THR.get(k, 65))
pf = cb["pfactor"]
DEF["pfactor"] = define("pfactor", pf.quantile(0.90), ctrl_thr=pf.quantile(0.75))
k = pd.read_csv(R70 / "y" / "mh_p_ksads__dep.tsv", sep="\t", dtype=str,
                usecols=["participant_id", "mh_p_ksads__dep__mdd__pres_dx", "mh_p_ksads__dep__mdd__past_dx"])
k["subject"] = "sub-NDARINV" + k.participant_id.str.replace("sub-", "", regex=False)
v = k[["mh_p_ksads__dep__mdd__pres_dx", "mh_p_ksads__dep__mdd__past_dx"]].apply(pd.to_numeric, errors="coerce")
v = v.where(v.isin([0, 1])); k["dx"] = v.max(axis=1)
e = k.dropna(subset=["dx"]).groupby("subject").dx.max()
e = e[e.index.isin(thin.index)]
DEF["mdd_parent"] = (e[e == 1], e[e == 0])

# ------------------------------------------------------------------ maps ------
def design(idx, with_global):
    X = pd.get_dummies(C.loc[idx, ["sex", "site"]].astype(str), drop_first=True).astype(float)
    X[["age_first", "age_span", "n_visits"]] = C.loc[idx, ["age_first", "age_span", "n_visits"]].to_numpy(float)
    if with_global:
        X["global"] = glob.loc[idx].to_numpy()
    return np.column_stack([np.ones(len(idx)), X.to_numpy(float)])

rng = np.random.default_rng(SEED)
maprows, corrrows = [], []
for name, (cas, ctl) in DEF.items():
    y = pd.concat([cas, ctl]); y = y[~y.index.duplicated()]
    idx = y.index[C.loc[y.index].notna().all(axis=1).to_numpy()]
    y = y.loc[idx]; Yb = thin.loc[idx].to_numpy()
    sites = C.loc[idx, "site"].astype(str).to_numpy()
    for adj in (False, True):
        Z = design(idx, adj)
        H = lambda M: M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]
        R = H(Yb)                                             # residualised thinning (FWL)
        def dmap(lbl):
            rc = H(lbl.astype(float))
            b = (R * rc[:, None]).sum(0) / (rc @ rc)
            res = R - np.outer(rc, b)
            dof = len(lbl) - Z.shape[1] - 1
            s = np.sqrt((res ** 2).sum(0) / dof)
            return b / s, b / (s / np.sqrt(rc @ rc)), dof
        d, t, dof = dmap(y.to_numpy())
        pv = 2 * stats.t.sf(np.abs(t), dof)
        kind = "relative" if adj else "absolute"
        m = pd.Series(d, index=thin.columns)
        maprows.append(pd.DataFrame({"outcome": name, "map": kind, "label": thin.columns, "d": d, "t": t, "p": pv,
                                     "n_case": int(y.sum()), "n_ctrl": int((y == 0).sum())}))
        # label-permutation null (within site)
        null = {r: [] for r in ("PLS2", "C3", "dCT")}
        for _ in range(N_LABEL):
            lp = y.to_numpy().copy()
            for s_ in np.unique(sites):
                ii = np.flatnonzero(sites == s_); lp[ii] = rng.permutation(lp[ii])
            dp = pd.Series(dmap(lp)[0], index=thin.columns)
            for r in null:
                sh = REF[r].index
                null[r].append(stats.spearmanr(dp.loc[sh], REF[r]).statistic)
        for r in ("PLS2", "C3", "dCT", "PLS1"):
            rho, psp, _ = pls.spin_corr(REF[r], m, n_perm=N_SPIN, seed=SEED, centroids=pls.HCP_CENTROIDS)
            pl = ((1 + (np.abs(np.array(null[r])) >= abs(rho)).sum()) / (1 + N_LABEL)) if r in null else np.nan
            corrrows.append(dict(outcome=name, map=kind, reference=r, rho=rho, p_spin=psp, p_label=pl,
                                 n_parcels=len(REF[r]), n_case=int(y.sum()), n_ctrl=int((y == 0).sum())))
        # partial rho with PLS2 given dCT (ranks)
        sh = REF["PLS2"].index
        rk = lambda s: stats.rankdata(s)
        a, b_, c_ = rk(m.loc[sh]), rk(REF["PLS2"]), rk(REF["dCT"].loc[sh])
        res_ = lambda u, v: u - np.polyval(np.polyfit(v, u, 1), v)
        corrrows.append(dict(outcome=name, map=kind, reference="PLS2 | dCT",
                             rho=np.corrcoef(res_(a, c_), res_(b_, c_))[0, 1], p_spin=np.nan, p_label=np.nan,
                             n_parcels=len(sh), n_case=int(y.sum()), n_ctrl=int((y == 0).sum())))
        print(f"{name:10s} {kind:8s} cases {int(y.sum()):5d} ctrl {int((y==0).sum()):5d} | "
              f"mean d {d.mean():+.3f} | parcels p<.05 {(pv<.05).sum():3d}/179 | "
              + " ".join(f"{r['reference']} {r['rho']:+.2f} (spin {r['p_spin']:.3f}, label {r['p_label']:.3f})"
                         for r in corrrows[-5:-1]), file=sys.stderr)

MAPS = pd.concat(maprows, ignore_index=True)
MAPS.to_csv(RES / "casecontrol_maps.tsv", sep="\t", index=False, float_format="%.5g")
CR = pd.DataFrame(corrrows)
CR.to_csv(RES / "casecontrol_map_corr.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 200)
print(CR.pivot_table(index=["outcome", "map"], columns="reference", values="rho").round(2).to_string())
