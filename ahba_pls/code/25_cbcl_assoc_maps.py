"""
25_cbcl_assoc_maps.py -- EXPLORATORY: regional maps of the LINEAR association
between each child's thinning rate and later CBCL symptoms, controlling for
baseline symptoms; do they resemble PLS2?  (Continuous counterpart of the
case-control maps in 24_casecontrol_maps.py -- same children, same outcomes,
same nulls, no cut-offs.)

Outcomes: log1p CBCL raw sums (as in 23), the outcomes associated with global
thinning in a developmental model: depress (DSM depressive), anxdisord (DSM
anxiety), anxdep (anxious/depressed), withdep (withdrawn/depressed), internal
(internalising), rulebreak (rule-breaking), pfactor (PC1 of the 8 syndromes,
fitted at years 5-7).
  y_late = mean over available ses-05A..07A (ages ~15-17);  y_base = ses-00A.

Per parcel r (bilateral, 179), all children with y_late and y_base:
  y_late ~ thin_r + y_base + sex + site (FE) + age_late + age_first + age_span
           + n_visits  [+ global thinning]  [+ CT0_r]
  map value = partial correlation of y_late with thin_r (sign: positive = faster
  thinning in r goes with more symptoms at 15-17 than baseline predicts);
  t and p from the same fit. CT0_r = the child's observed first-scan thickness
  in r (lh/rh mean). Variants: absolute, relative (+ global), absolute_ct, relative_ct.

Map-level tests (Spearman rho) vs PLS2 (137 AHBA parcels), AHBA C3, normative
dCT (mm/yr, all negative: higher = SLOWER thinning), PLS1:
  p_spin  : 5,000 spin rotations of the reference map (HCP centroids)
  p_perm  : 1,000 within-site permutations of the covariate-residualised outcome
            (Freedman-Lane) -- is the map more PLS2-like than maps built from a
            random outcome with the same covariate structure?
Group-level outputs only: results/cbcl_assoc_maps.tsv, cbcl_assoc_map_corr.tsv.
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
N_SPIN, N_PERM, SEED = 5000, 1000, 0

# ------------------------------------------------------------------ brain -----
bl = pd.read_parquet(RUN / "fits" / "blups.parquet")
fx = pd.read_parquet(RUN / "fits" / "fixed.parquet")
thin = -(bl.pivot(index="subject", columns="label", values="re_slope")
         + fx[fx.term == "age_c"].set_index("label").estimate)
glob = thin.mean(axis=1)                                      # 358 parcels
bil = lambda X: (lambda Y: Y.set_axis("lh_" + Y.columns, axis=1))(X.T.groupby(X.columns.str.slice(3)).mean().T)
thin = bil(thin); assert thin.shape[1] == 179
tt = pd.read_parquet(RUN / "model_table.parquet", columns=["subject", "label", "value", "age"])
ct0 = bil(tt.sort_values("age").groupby(["subject", "label"]).value.first().unstack()).loc[thin.index, thin.columns]
del tt
mt = pd.read_parquet(RUN / "model_table.parquet",
                     columns=["subject", "sex", "site", "age_first", "age_span", "n_visits"]
                     ).drop_duplicates("subject").set_index("subject")
C = mt.loc[thin.index]

maps = pd.read_csv(RES / "hcp_summary_maps.csv", index_col=0)
lab = {l.lower(): l for l in thin.columns}
maps.index = [lab[l.lower()] for l in maps.index]
REF = {"PLS2": maps.PLS2.dropna(), "C3": maps.C3.dropna(), "dCT": maps.dCT.dropna(), "PLS1": maps.PLS1.dropna()}

# ------------------------------------------------------------------ CBCL ------
RAW = {"depress": "dsm__dep", "anxdisord": "dsm__anx", "anxdep": "synd__anxdep", "withdep": "synd__wthdep",
       "internal": "synd__int", "rulebreak": "synd__rule"}
SYN = ["anxdep", "wthdep", "som", "soc", "tho", "attn", "rule", "aggr"]
cols = {**{k: f"mh_p_cbcl__{v}_sum" for k, v in RAW.items()}, **{f"s_{s}": f"mh_p_cbcl__synd__{s}_sum" for s in SYN}}
cb = pd.read_csv(R70 / "p" / "mh_p_cbcl.tsv", sep="\t", usecols=["participant_id", "session_id", "mh_p_cbcl_age", *cols.values()])
cb["subject"] = "sub-NDARINV" + cb.participant_id.str.replace("sub-", "", regex=False)
cb["wave"] = cb.session_id.str.slice(4, 6).astype(int); cb["age"] = cb.mh_p_cbcl_age
cb = cb[cb.subject.isin(thin.index)]
for k, c in cols.items():
    cb[k] = np.log1p(pd.to_numeric(cb[c], errors="coerce"))
S_ = [f"s_{s}" for s in SYN]
lm = cb[cb.wave.isin([5, 6, 7])].groupby("subject")[S_].mean().dropna()
mu, sd = lm.mean(), lm.std()
_, _, vt = np.linalg.svd(((lm - mu) / sd).to_numpy(), full_matrices=False)
cb["pfactor"] = ((cb[S_] - mu) / sd).to_numpy() @ (vt[0] * np.sign(vt[0].sum()))
OUT = list(RAW) + ["pfactor"]
late = cb[cb.wave.isin([5, 6, 7])].groupby("subject")[OUT + ["age"]].mean()
base = cb[cb.wave == 0].groupby("subject")[OUT].mean()

# ------------------------------------------------------------------ maps ------
rng = np.random.default_rng(SEED)
maprows, corrrows = [], []
for o in OUT:
    d = pd.DataFrame({"y": late[o], "y0": base[o], "age_late": late["age"]}).join(C).join(glob.rename("glob"))
    d = d.dropna(); d = d.loc[d.index.intersection(thin.index)]
    idx = d.index; sites = d.site.astype(str).to_numpy()
    Yb, CTb = thin.loc[idx].to_numpy(), ct0.loc[idx].to_numpy()
    for adj, adj_ct in ((False, False), (True, False), (False, True), (True, True)):
        X = pd.get_dummies(d[["sex", "site"]].astype(str), drop_first=True).astype(float)
        X[["y0", "age_late", "age_first", "age_span", "n_visits"]] = d[["y0", "age_late", "age_first", "age_span", "n_visits"]]
        if adj:
            X["glob"] = d.glob
        Z = np.column_stack([np.ones(len(d)), X.to_numpy(float)])
        H = lambda M: M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]
        R = H(Yb)                                             # residualised regional thinning
        if adj_ct:
            Q = H(CTb); R = R - Q * ((R * Q).sum(0) / (Q * Q).sum(0))
        ry = H(d.y.to_numpy())
        def rmap(v):                                          # partial r per parcel
            if adj_ct:
                V = v[:, None] - Q * ((v[:, None] * Q).sum(0) / (Q * Q).sum(0))
            else:
                V = v[:, None]
            return (R * V).sum(0) / np.sqrt((R * R).sum(0) * (V * V).sum(0))
        r = rmap(ry)
        dof = len(d) - Z.shape[1] - 1 - int(adj_ct)
        t = r * np.sqrt(dof / (1 - r ** 2)); pv = 2 * stats.t.sf(np.abs(t), dof)
        kind = ("relative" if adj else "absolute") + ("_ct" if adj_ct else "")
        m = pd.Series(r, index=thin.columns)
        maprows.append(pd.DataFrame({"outcome": o, "map": kind, "label": thin.columns, "r": r, "t": t, "p": pv, "n": len(d)}))
        null = {k: [] for k in ("PLS2", "C3", "dCT")}
        for _ in range(N_PERM):                               # Freedman-Lane: permute residualised y within site
            yp = ry.copy()
            for s_ in np.unique(sites):
                ii = np.flatnonzero(sites == s_); yp[ii] = yp[rng.permutation(ii)]
            mp = pd.Series(rmap(yp), index=thin.columns)
            for k in null:
                null[k].append(stats.spearmanr(mp.loc[REF[k].index], REF[k]).statistic)
        for k in ("PLS2", "C3", "dCT", "PLS1"):
            rho, psp, _ = pls.spin_corr(REF[k], m, n_perm=N_SPIN, seed=SEED, centroids=pls.HCP_CENTROIDS)
            pp = (1 + (np.abs(np.array(null[k])) >= abs(rho)).sum()) / (1 + N_PERM) if k in null else np.nan
            corrrows.append(dict(outcome=o, map=kind, reference=k, rho=rho, p_spin=psp, p_perm=pp, n=len(d)))
        print(f"{o:9s} {kind:11s} n {len(d)} | mean r {r.mean():+.4f} | parcels p<.05 {(pv<.05).sum():3d}/179 | "
              + " ".join(f"{c['reference']} {c['rho']:+.2f} (spin {c['p_spin']:.3f}, perm {c['p_perm']:.3f})"
                         for c in corrrows[-4:-1]), file=sys.stderr)

pd.concat(maprows, ignore_index=True).to_csv(RES / "cbcl_assoc_maps.tsv", sep="\t", index=False, float_format="%.5g")
CR = pd.DataFrame(corrrows); CR.to_csv(RES / "cbcl_assoc_map_corr.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 200)
print(CR.pivot_table(index=["outcome", "map"], columns="reference", values="rho").round(2).to_string())
both = CR[(CR.p_spin < .05) & (CR.p_perm < .05)]
print("\npass both nulls:\n", both[["outcome", "map", "reference", "rho", "p_spin", "p_perm"]].round(3).to_string(index=False))
