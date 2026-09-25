"""
30_cbcl_assoc_cortices.py -- EXPLORATORY: the linear thinning-symptom
association maps of 25_cbcl_assoc_maps.py rebuilt on the 22 Glasser 2016
cortices instead of the 179 bilateral HCP-MMP parcels, to average out
parcel-level noise. Outcomes: p-factor and total problems (log1p CBCL total).

Cortices: HCP-MMP1_UniqueRegionList.csv `Cortex_ID` 1-22 (Glasser et al. 2016,
Suppl. Neuroanatomical Results; lookup data/reference/hcp_cortices/
hcp_parcel_systems.csv, vendored for fig5_hcp_summary.R).

Per child, thinning in a cortex = mean of their bilateral-parcel thinning over
the cortex's parcels (likewise baseline CT). Model and variants exactly as in
25: per region, partial r of CBCL_late (mean of ses-05A..07A) with thinning,
given CBCL_base + sex + site + age_late + age_first + age_span + n_visits
[+ global thinning] [+ that region's baseline CT].

Reference maps per cortex: mean of PLS2 / C3 over that cortex's AHBA-covered
parcels (hcp_3d_ds5 fit, 137 parcels) and of normative dCT over all its parcels.
Nulls, both at the level being tested:
  p_spin : 5,000 spin rotations of the PARCEL-level reference map, each
           re-averaged into cortices (so the null keeps the real spatial
           autocorrelation; spinning 22 centroids would be far too coarse)
  p_perm : 1,000 within-site Freedman-Lane permutations of the residualised
           outcome, maps rebuilt each time
The parcel-level maps for the same two outcomes are recomputed here with the
same code and nulls, so the two resolutions are compared like-for-like.

Group-level outputs only: results/cbcl_assoc_cortex_maps.tsv (per-cortex r),
results/cbcl_assoc_cortex_corr.tsv (map-vs-reference, both levels).
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
C = pd.read_parquet(RUN / "model_table.parquet", columns=["subject", "sex", "site", "age_first", "age_span", "n_visits"]
                    ).drop_duplicates("subject").set_index("subject").loc[thin.index]

lab = {l.lower(): l for l in thin.columns}
maps = pd.read_csv(RES / "hcp_summary_maps.csv", index_col=0)
maps.index = [lab[l.lower()] for l in maps.index]
REF = {"PLS2": maps.PLS2.dropna(), "C3": maps.C3.dropna(), "dCT": maps.dCT.dropna()}

sysl = pd.read_csv(ROOT / "data" / "reference" / "hcp_cortices" / "hcp_parcel_systems.csv", index_col=0)
sysl = sysl[sysl.index.str.lower().isin(lab)]                # parcel H is excluded from the fit (medial wall)
sysl.index = [lab[l.lower()] for l in sysl.index]
cid = sysl.cortex_id.loc[thin.columns]
cname = sysl.groupby("cortex_id").cortex.first().str.replace("Temporo-Parieto_Occipital", "Temporo-Parieto-Occipital")
assert cid.notna().all() and cid.nunique() == 22, cid.nunique()
CORT = sorted(cid.unique())

def agg_matrix(parcels) -> tuple[np.ndarray, list]:
    """Row-normalised (n_cortex x n_parcels) averaging matrix over the given parcels."""
    c = cid.loc[parcels].to_numpy()
    keep = [k for k in CORT if (c == k).any()]
    A = np.array([(c == k).astype(float) / (c == k).sum() for k in keep])
    return A, keep

thin_c = pd.DataFrame(thin.to_numpy() @ agg_matrix(thin.columns)[0].T, index=thin.index, columns=CORT)
ct0_c = pd.DataFrame(ct0.to_numpy() @ agg_matrix(thin.columns)[0].T, index=thin.index, columns=CORT)

# reference maps and spin nulls at both levels: (observed ref vector, null matrix n_spin x n_regions, region labels)
def ref_level(k: str, level: str):
    r = REF[k]
    perms = pls.cached_spin_permutations(r.index, n_perm=N_SPIN, seed=SEED, path=pls.HCP_CENTROIDS)
    P = r.to_numpy()[perms]                                   # (n_spin, n_parcels) rotated reference
    if level == "parcel":
        return r.to_numpy(), P, list(r.index)
    A, keep = agg_matrix(r.index)
    return A @ r.to_numpy(), P @ A.T, keep

REFL = {(k, lev): ref_level(k, lev) for k in REF for lev in ("parcel", "cortex")}
cov = pd.Series({k: int((cid.loc[REF["C3"].index] == k).sum()) for k in CORT})
print("AHBA-covered parcels per cortex:", cov.to_dict(), file=sys.stderr)

def spear_rows(v: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Spearman of vector v with every row of M."""
    rv = stats.rankdata(v); rv = (rv - rv.mean()) / rv.std()
    RM = np.apply_along_axis(stats.rankdata, 1, M); RM = (RM - RM.mean(1, keepdims=True)) / RM.std(1, keepdims=True)
    return (RM * rv).mean(1)

# ------------------------------------------------------------------ CBCL ------
RAW = {"totprob": "mh_p_cbcl_sum"}
SYN = ["anxdep", "wthdep", "som", "soc", "tho", "attn", "rule", "aggr"]
cols = RAW | {f"s_{s}": f"mh_p_cbcl__synd__{s}_sum" for s in SYN}
cb = pd.read_csv(R70 / "p" / "mh_p_cbcl.tsv", sep="\t", usecols=["participant_id", "session_id", "mh_p_cbcl_age", *cols.values()])
cb["subject"] = "sub-NDARINV" + cb.participant_id.str.replace("sub-", "", regex=False)
cb["wave"] = cb.session_id.str.slice(4, 6).astype(int); cb["age"] = cb.mh_p_cbcl_age
cb = cb[cb.subject.isin(thin.index)]
for k, c in cols.items():
    cb[k] = np.log1p(pd.to_numeric(cb[c], errors="coerce"))
S_ = [f"s_{s}" for s in SYN]
lm = cb[cb.wave.isin([5, 6, 7])].groupby("subject")[S_].mean().dropna(); mu, sd = lm.mean(), lm.std()
_, _, vt = np.linalg.svd(((lm - mu) / sd).to_numpy(), full_matrices=False)
cb["pfactor"] = ((cb[S_] - mu) / sd).to_numpy() @ (vt[0] * np.sign(vt[0].sum()))
OUT = ["pfactor", "totprob"]
late = cb[cb.wave.isin([5, 6, 7])].groupby("subject")[OUT + ["age"]].mean()
base = cb[cb.wave == 0].groupby("subject")[OUT].mean()

# ------------------------------------------------------------------ maps ------
rng = np.random.default_rng(SEED)
maprows, corrrows = [], []
LEVELS = {"parcel": (thin, ct0), "cortex": (thin_c, ct0_c)}
for o in OUT:
    d = pd.DataFrame({"y": late[o], "y0": base[o], "age_late": late["age"]}).join(C).join(glob.rename("glob")).dropna()
    idx = d.index; sites = d.site.astype(str).to_numpy()
    for adj, adj_ct in ((False, False), (True, False), (False, True), (True, True)):
        kind = ("relative" if adj else "absolute") + ("_ct" if adj_ct else "")
        X = pd.get_dummies(d[["sex", "site"]].astype(str), drop_first=True).astype(float)
        X[["y0", "age_late", "age_first", "age_span", "n_visits"]] = d[["y0", "age_late", "age_first", "age_span", "n_visits"]]
        if adj:
            X["glob"] = d.glob
        Z = np.column_stack([np.ones(len(d)), X.to_numpy(float)])
        H = lambda M: M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]
        ry = H(d.y.to_numpy())
        # per level: residualised thinning (and CT) and the partial-r map function
        fns = {}
        for lev, (TH, CT0) in LEVELS.items():
            R = H(TH.loc[idx].to_numpy())
            Q = H(CT0.loc[idx].to_numpy()) if adj_ct else None
            if adj_ct:
                R = R - Q * ((R * Q).sum(0) / (Q * Q).sum(0))
            def rmap(v, R=R, Q=Q):
                V = v[:, None] - Q * ((v[:, None] * Q).sum(0) / (Q * Q).sum(0)) if Q is not None else v[:, None]
                return (R * V).sum(0) / np.sqrt((R * R).sum(0) * (V * V).sum(0))
            fns[lev] = (rmap, list(TH.columns))
        obs = {lev: pd.Series(f(ry), index=cols_) for lev, (f, cols_) in fns.items()}
        dof = len(d) - Z.shape[1] - 1 - int(adj_ct)
        rc = obs["cortex"]; tc = rc * np.sqrt(dof / (1 - rc ** 2))
        maprows.append(pd.DataFrame({"outcome": o, "map": kind, "cortex_id": rc.index, "cortex": cname.loc[rc.index].values,
                                     "r": rc.values, "t": tc.values, "p": 2 * stats.t.sf(np.abs(tc.values), dof), "n": len(d)}))
        null = {(k, lev): [] for k in REF for lev in LEVELS}
        for _ in range(N_PERM):                               # Freedman-Lane: permute residualised y within site
            yp = ry.copy()
            for s_ in np.unique(sites):
                ii = np.flatnonzero(sites == s_); yp[ii] = yp[rng.permutation(ii)]
            for lev, (f, cols_) in fns.items():
                mp = pd.Series(f(yp), index=cols_)
                for k in REF:
                    rv, _, regs = REFL[(k, lev)]
                    null[(k, lev)].append(stats.spearmanr(mp.loc[regs], rv).statistic)
        for lev in LEVELS:
            for k in REF:
                rv, P, regs = REFL[(k, lev)]
                m = obs[lev].loc[regs].to_numpy()
                rho = stats.spearmanr(m, rv).statistic
                ns = spear_rows(m, P)
                npm = np.array(null[(k, lev)])
                corrrows.append(dict(outcome=o, map=kind, level=lev, reference=k, rho=rho, n_regions=len(regs),
                                     p_spin=(1 + (np.abs(ns) >= abs(rho)).sum()) / (1 + N_SPIN),
                                     p_perm=(1 + (np.abs(npm) >= abs(rho)).sum()) / (1 + N_PERM), n=len(d)))
        last = [c for c in corrrows[-6:]]
        print(f"{o:8s} {kind:11s} n {len(d)} | " + " ".join(f"{c['level'][:3]} {c['reference']} {c['rho']:+.2f} "
              f"({c['p_spin']:.3f}/{c['p_perm']:.3f})" for c in last), file=sys.stderr)

pd.concat(maprows, ignore_index=True).to_csv(RES / "cbcl_assoc_cortex_maps.tsv", sep="\t", index=False, float_format="%.5g")
CR = pd.DataFrame(corrrows); CR.to_csv(RES / "cbcl_assoc_cortex_corr.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 220)
CR["c"] = CR.apply(lambda r: f"{r.rho:+.2f} ({r.p_spin:.3f}, {r.p_perm:.3f})", axis=1)
print(CR.pivot_table(index=["outcome", "map"], columns=["level", "reference"], values="c", aggfunc="first").to_string())
