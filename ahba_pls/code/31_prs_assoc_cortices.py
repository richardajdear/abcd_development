"""
31_prs_assoc_cortices.py -- EXPLORATORY: regional maps of the association
between SCZ / MDD polygenic scores and adolescent thinning rate, on the 22
Glasser cortices and on the 179 bilateral HCP-MMP parcels, tested against
normative dCT, PLS2 and AHBA C3.

Scores (PRS-CS and SBayesRC, the Figure 1f cells; per-child .profile files copied from CSD3
into genetic_analysis/work/, gitignored, never committed):
  SCZ25_META  pooled arm,  within-ancestry standardised score (_zanc)
  SCZ25_EUR   EUR arm,     raw score, EUR subset (needs inputs/ancestry/eur_anchor.keep)
  MDD_pooled  pooled arm,  _zanc      (legacy/hpc_v2/work/results_v2/prs_final/PRSCS/)
  MDD_eur     EUR arm,     raw, EUR subset
Arms whose files are absent are skipped and listed.

Model per region (fixed part of tools/prs_assoc.R, the Figure 1 model):
  thinning_r ~ PRS + sex + baseline_age + PC1..PC10 [+ global thinning] [+ CT0_r]
Map value = partial r of thinning_r with the score. OLS, so no family random
effect: that changes per-region SEs (anti-conservative per-region p), not the
map or the map-level tests. Thinning is oriented positive = faster, so a
positive value = higher score, faster thinning.

Whole-cortex check: the same model on global_slope_1lmm (phenotypes_gcta.txt)
is printed next to the Figure 1 PRS-CS beta, to confirm IDs and covariates.

Nulls for map vs reference (as 30_cbcl_assoc_cortices.py):
  p_spin : 5,000 spins of the PARCEL-level reference, re-averaged into cortices
  p_perm : 1,000 permutations of the residualised score across children
Group-level outputs only: results/prs_assoc_cortex_maps.tsv,
results/prs_assoc_cortex_corr.tsv, results/prs_assoc_global_check.tsv.
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
GW = REPO / "genetic_analysis" / "work"
PH = GW / "results_70tab_hcp" / "prs_final_1lmm" / "pheno"
# the EUR list and the MDD scores live in the legacy tree the Figure 1 cluster jobs read
# (orderops/prs_1lmm.sbatch, step6_prs_1lmm_controls.sbatch); mirrored locally at the same relative paths
EURF = next((p for p in (REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep",
                         GW / "inputs/ancestry/eur_anchor.keep") if p.exists()),
            REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep")
LEG = REPO / "legacy/hpc_v2/work/results_v2/prs_final"
N_SPIN, N_PERM, SEED = 5000, 1000, 0
SC25 = GW / "scores_scz2025"
ARMS = {  # arm: (profile, EUR subset?, trait_arm in the cluster tables, score type)
    "PRSCS_SCZ25_META":    (SC25 / "PRSCS/SCZ25_META/_zanc/score_SCZ25META_prscs.profile", False, "SCZ25_META", "zanc"),
    "PRSCS_SCZ25_EUR":     (SC25 / "PRSCS/SCZ25_EUR/score_SCZ25EUR_prscs.profile", True, "SCZ25_EUR", "raw"),
    "PRSCS_MDD_pooled":    (LEG / "PRSCS/MDD_pooled/_zanc/score_MDDpooled_prscs.profile", False, "MDD_pooled", "zanc"),
    "PRSCS_MDD_eur":       (LEG / "PRSCS/MDD_eur/score_MDDeur_prscs.profile", True, "MDD_eur", "raw"),
    "SBayesRC_SCZ25_META": (SC25 / "SBayesRC/SCZ25_META/_zanc/score_SCZ25META_sbrc.profile", False, "SCZ25_META", "zanc"),
    "SBayesRC_SCZ25_EUR":  (SC25 / "SBayesRC/SCZ25_EUR/score_SCZ25EUR_sbrc.profile", True, "SCZ25_EUR", "raw"),
    # PGC3 (Trubetskoy 2022) SCZ, used only while the SCZ 2025 SBayesRC scores are absent locally
    "SBayesRC_SCZpgc3_pooled": (LEG / "SBayesRC/SCZ_pooled/_zanc/score_SCZpooled_sbrc.profile", False, "SCZ_pooled", "zanc"),
    "SBayesRC_SCZpgc3_eur":    (LEG / "SBayesRC/SCZ_eur/score_SCZeur_sbrc.profile", True, "SCZ_eur", "raw"),
    "SBayesRC_MDD_pooled": (LEG / "SBayesRC/MDD_pooled/_zanc/score_MDDpooled_sbrc.profile", False, "MDD_pooled", "zanc"),
    "SBayesRC_MDD_eur":    (LEG / "SBayesRC/MDD_eur/score_MDDeur_sbrc.profile", True, "MDD_eur", "raw"),
}

def find_profile(p: Path) -> Path | None:
    if p.exists():
        return p
    c = sorted(p.parent.glob("score_*" + p.name.rsplit("_", 1)[-1])) if p.parent.exists() else []
    return c[0] if c else None

present = {a: (find_profile(v[0]),) + v[1:] for a, v in ARMS.items()}
skip = [a for a, v in present.items() if v[0] is None or (v[1] and not EURF.exists())]
present = {a: v for a, v in present.items() if a not in skip}
if "SBayesRC_SCZ25_META" in present:
    present = {a: v for a, v in present.items() if "pgc3" not in a}
print("arms run:", list(present), "| skipped (files absent):", skip, file=sys.stderr)
assert present, "no PRS profile files found"

# ------------------------------------------------------------------ brain -----
bl = pd.read_parquet(RUN / "fits" / "blups.parquet")
fx = pd.read_parquet(RUN / "fits" / "fixed.parquet")
thin = -(bl.pivot(index="subject", columns="label", values="re_slope")
         + fx[fx.term == "age_c"].set_index("label").estimate)
glob = thin.mean(axis=1)
bil = lambda X: (lambda Y: Y.set_axis("lh_" + Y.columns, axis=1))(X.T.groupby(X.columns.str.slice(3)).mean().T)
thin = bil(thin); assert thin.shape[1] == 179
tt = pd.read_parquet(RUN / "model_table.parquet", columns=["subject", "label", "value", "age"])
ct0 = bil(tt.sort_values("age").groupby(["subject", "label"]).value.first().unstack()).loc[thin.index, thin.columns]
del tt

lab = {l.lower(): l for l in thin.columns}
maps = pd.read_csv(RES / "hcp_summary_maps.csv", index_col=0)
maps.index = [lab[l.lower()] for l in maps.index]
REF = {"PLS2": maps.PLS2.dropna(), "C3": maps.C3.dropna(), "dCT": maps.dCT.dropna()}
sysl = pd.read_csv(ROOT / "data" / "reference" / "hcp_cortices" / "hcp_parcel_systems.csv", index_col=0)
sysl = sysl[sysl.index.str.lower().isin(lab)]
sysl.index = [lab[l.lower()] for l in sysl.index]
cid = sysl.cortex_id.loc[thin.columns]
cname = sysl.groupby("cortex_id").cortex.first()
CORT = sorted(cid.unique()); assert len(CORT) == 22

def agg_matrix(parcels):
    c = cid.loc[parcels].to_numpy(); keep = [k for k in CORT if (c == k).any()]
    return np.array([(c == k).astype(float) / (c == k).sum() for k in keep]), keep

A_all = agg_matrix(thin.columns)[0]
thin_c = pd.DataFrame(thin.to_numpy() @ A_all.T, index=thin.index, columns=CORT)
ct0_c = pd.DataFrame(ct0.to_numpy() @ A_all.T, index=thin.index, columns=CORT)

def ref_level(k, level):
    r = REF[k]
    P = r.to_numpy()[pls.cached_spin_permutations(r.index, n_perm=N_SPIN, seed=SEED, path=pls.HCP_CENTROIDS)]
    if level == "parcel":
        return r.to_numpy(), P, list(r.index)
    A, keep = agg_matrix(r.index)
    return A @ r.to_numpy(), P @ A.T, keep
REFL = {(k, lev): ref_level(k, lev) for k in REF for lev in ("parcel", "cortex")}

def spear_rows(v, M):
    rv = stats.rankdata(v); rv = (rv - rv.mean()) / rv.std()
    RM = np.apply_along_axis(stats.rankdata, 1, M); RM = (RM - RM.mean(1, keepdims=True)) / RM.std(1, keepdims=True)
    return (RM * rv).mean(1)

# ------------------------------------------------------------ covariates -----
qc = pd.read_csv(PH / "covar_quant.txt", sep=r"\s+"); cc = pd.read_csv(PH / "covar_categorical.txt", sep=r"\s+")
ph = pd.read_csv(PH / "phenotypes_gcta.txt", sep=r"\s+")
cov = qc.merge(cc, on=["FID", "IID"]).merge(ph, on=["FID", "IID"])
cov["subject"] = "sub-NDARINV" + cov.IID.str.replace("sub-", "", regex=False)
cov = cov.set_index("subject")
PCS = [f"PC{i}" for i in range(1, 11)]
eur = set(pd.read_csv(EURF, sep=r"\s+", header=None).iloc[:, -1]) if EURF.exists() else set()
fig1 = pd.concat([pd.read_csv(REPO / "genetic_analysis" / "fig1_inputs" / "hcp70_prs_key_arms.tsv", sep="\t"),
                  pd.read_csv(GW / "results_70tab_hcp" / "prs_final_1lmm" / "table_order_of_operations_all.tsv", sep="\t")
                  .assign(phenotype=lambda x: x.phenotype.str.replace("_1lmm", ""))], ignore_index=True)
fig1 = fig1[fig1.atlas == "hcp"].drop_duplicates(["phenotype", "trait_arm", "method", "score"])

rng = np.random.default_rng(SEED)
maprows, corrrows, grows = [], [], []
for arm, (prof, is_eur, f1arm, stype) in present.items():
    sc = pd.read_csv(prof, sep=r"\s+")
    sc["subject"] = "sub-NDARINV" + sc.IID.str.replace("sub-", "", regex=False)
    d = cov.join(sc.set_index("subject").SCORESUM.rename("prs"), how="inner")
    d = d[d.index.isin(thin.index)]
    if is_eur:
        d = d[d.IID.isin(eur)]
    d = d.join(glob.rename("glob")).dropna(subset=["prs", "sex", "baseline_age", "glob", *PCS])
    idx = d.index
    X0 = np.column_stack([np.ones(len(d)), (d.sex.astype(str) == "M").astype(float), d[["baseline_age", *PCS]].to_numpy(float)])
    # whole-cortex check against Figure 1 (global_slope_1lmm, standardised, as prs_assoc.R)
    y = d.global_slope_1lmm; ok = y.notna()
    Xg = np.column_stack([((d.prs - d.prs.mean()) / d.prs.std())[ok], X0[ok.to_numpy()]])
    bg = np.linalg.lstsq(Xg, ((y - y[ok].mean()) / y[ok].std())[ok], rcond=None)[0][0]
    f1 = fig1[(fig1.phenotype == "global_slope") & (fig1.trait_arm == f1arm) & (fig1.method == arm.split("_")[0]) & (fig1.score == stype)]
    grows.append(dict(arm=arm, n=len(d), beta_global_slope_ols=bg,
                      beta_fig1_lmm=f1.beta_1lmm.iloc[0] if len(f1) else np.nan))
    for adj, adj_ct in ((False, False), (True, False), (False, True), (True, True)):
        kind = ("relative" if adj else "absolute") + ("_ct" if adj_ct else "")
        Z = np.column_stack([X0, d.glob.to_numpy()]) if adj else X0
        H = lambda M: M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]
        rs = H(d.prs.to_numpy(float))
        fns = {}
        for lev, (TH, CT0) in {"parcel": (thin, ct0), "cortex": (thin_c, ct0_c)}.items():
            R = H(TH.loc[idx].to_numpy()); Q = H(CT0.loc[idx].to_numpy()) if adj_ct else None
            if adj_ct:
                R = R - Q * ((R * Q).sum(0) / (Q * Q).sum(0))
            def rmap(v, R=R, Q=Q):
                V = v[:, None] - Q * ((v[:, None] * Q).sum(0) / (Q * Q).sum(0)) if Q is not None else v[:, None]
                return (R * V).sum(0) / np.sqrt((R * R).sum(0) * (V * V).sum(0))
            fns[lev] = (rmap, list(TH.columns))
        obs = {lev: pd.Series(f(rs), index=c) for lev, (f, c) in fns.items()}
        dof = len(d) - Z.shape[1] - 1 - int(adj_ct)
        for lev in ("cortex", "parcel"):
            r = obs[lev]; t = r * np.sqrt(dof / (1 - r ** 2))
            maprows.append(pd.DataFrame({"arm": arm, "map": kind, "level": lev, "region": r.index.astype(str),
                                         "cortex": cname.reindex(r.index).values if lev == "cortex" else None,
                                         "r": r.values, "t": t.values, "p": 2 * stats.t.sf(np.abs(t.values), dof), "n": len(d)}))
        null = {(k, lev): [] for k in REF for lev in fns}
        for _ in range(N_PERM):
            vp = rng.permutation(rs)
            for lev, (f, c) in fns.items():
                mp = pd.Series(f(vp), index=c)
                for k in REF:
                    rv, _, regs = REFL[(k, lev)]
                    null[(k, lev)].append(stats.spearmanr(mp.loc[regs], rv).statistic)
        for lev in fns:
            for k in REF:
                rv, P, regs = REFL[(k, lev)]
                m = obs[lev].loc[regs].to_numpy(); rho = stats.spearmanr(m, rv).statistic
                ns = spear_rows(m, P); npm = np.array(null[(k, lev)])
                corrrows.append(dict(arm=arm, map=kind, level=lev, reference=k, rho=rho, n_regions=len(regs),
                                     mean_r=obs[lev].mean(), p_spin=(1 + (np.abs(ns) >= abs(rho)).sum()) / (1 + N_SPIN),
                                     p_perm=(1 + (np.abs(npm) >= abs(rho)).sum()) / (1 + N_PERM), n=len(d)))
        print(f"{arm:10s} {kind:11s} n {len(d)} mean r(cortex) {obs['cortex'].mean():+.4f} | " + " ".join(
              f"{c['level'][:3]} {c['reference']} {c['rho']:+.2f} ({c['p_spin']:.3f}/{c['p_perm']:.3f})" for c in corrrows[-6:]),
              file=sys.stderr)

pd.concat(maprows, ignore_index=True).to_csv(RES / "prs_assoc_cortex_maps.tsv", sep="\t", index=False, float_format="%.5g")
CR = pd.DataFrame(corrrows); CR.to_csv(RES / "prs_assoc_cortex_corr.tsv", sep="\t", index=False, float_format="%.4g")
G = pd.DataFrame(grows); G.to_csv(RES / "prs_assoc_global_check.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 220)
print(G.round(4).to_string(index=False))
CR["c"] = CR.apply(lambda r: f"{r.rho:+.2f} ({r.p_spin:.3f}, {r.p_perm:.3f})", axis=1)
print(CR.pivot_table(index=["arm", "map"], columns=["level", "reference"], values="c", aggfunc="first").to_string())
