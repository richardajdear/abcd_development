"""D1 step 6: is there one disease-relevant process behind thinning and diffusion PC4?

    python directions/d1_multimodal/08_single_factor.py [n_boot]      (env abcd-spatial; default 1000 / 300 twin)

If one latent process P drives both measures (x_i = a_i P + e_i), four things follow:
  (1) the child slopes correlate (r = a1 a2 / sqrt(var1 var2), beyond shared artefact);
  (2) a composite has higher reliability than either part;
  (3) every polygenic score acts through P, so the PRS effects are proportional across scores:
      b(x2, SCZ) / b(x1, SCZ) = b(x2, MDD) / b(x1, MDD)  <=>  tetrad t = b2S*b1M - b2M*b1S = 0;
  (4) their genetic correlation (twin rA) is high.
(1) and (2) were answered in 03/05 (r(thinning, dPC4) = -0.07). This script tests (3) and (4) and shows what
combining does to the PRS signal. A composite can raise the PRS z without a shared process, simply by adding two
independent signals; that is aggregation, not evidence for one process.

PRS model: OLS, y_z ~ prs_z + sex + baseline_age + PC1-10 + site, family-clustered SE (D3's population model; it
matches the Figure-1 LMM to ~0.002). PRS-CS scores. Arms: EUR (EUR weights, z within arm), pooled (_zanc), and
non-EUR (pooled scores restricted to eur == 0) -- the only subsample not used when dPC4 was picked out.
Twin: D3 ace.py AE models on site/covariate-residualised traits (as 03_build_traits_pairs.py), pair bootstrap.
Outputs: tables/single_factor_prs.tsv, single_factor_tetrad.tsv, single_factor_twin.tsv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
D3 = HERE.parent / "d3_twin_family"
sys.path.insert(0, str(D3))
import ace  # noqa: E402
from common import residualise  # noqa: E402

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
W, TAB = HERE / "work", HERE / "tables"
rng = np.random.default_rng(6)

B = pd.read_csv(W / "blups_dmri.csv.gz")
thin = (-B[B.metric == "smri__thk_hcp"].set_index("subject").slope).rename("thin")
pcs = pd.read_csv(W / "pls_child.csv.gz").set_index("subject")[["PC2", "PC4"]].rename(columns={"PC2": "dPC2", "PC4": "dPC4"})
ch = pd.read_csv(REPO / "genetic_analysis/work/results_70tab_hcp/d9_sex_moderation/children.csv").set_index("subject")
D = ch.join(thin, how="inner").join(pcs, how="inner")
D["male"] = (D.sex == "M").astype(float)
PCS = [f"PC{i}" for i in range(1, 11)]
SITE = pd.get_dummies(D.site.astype(str), prefix="s", drop_first=True, dtype=float)
COV = pd.concat([D[["male", "baseline_age", *PCS]], SITE], axis=1)
z = lambda s: (s - s.mean()) / s.std()
print(f"children with thinning + diffusion PCs + covariates: {len(D)}")


def arm_index(arm, score):
    ok = D[score].notna() & D.thin.notna() & D.dPC4.notna() & COV.notna().all(axis=1)
    if arm == "EUR":
        ok &= D.eur == 1
    if arm == "nonEUR":
        ok &= D.eur == 0
    return D.index[ok]


def ols(y, x, idx):
    X = sm.add_constant(pd.concat([z(x.loc[idx]).rename("prs"), COV.loc[idx]], axis=1))
    m = sm.OLS(z(y.loc[idx]), X).fit(cov_type="cluster", cov_kwds={"groups": D.loc[idx, "family"].astype(str)})
    return m.params["prs"], m.bse["prs"], m.pvalues["prs"]


SC = {("SCZ", "EUR"): "scz_prscs_eur", ("SCZ", "pooled"): "scz_prscs_pooled", ("SCZ", "nonEUR"): "scz_prscs_pooled",
      ("MDD", "EUR"): "mdd_prscs_eur", ("MDD", "pooled"): "mdd_prscs_pooled", ("MDD", "nonEUR"): "mdd_prscs_pooled"}

# ---- (a) composites -------------------------------------------------------------------------------
rows = []
for (dis, arm), sc in SC.items():
    idx = arm_index(arm, sc)
    zt, z4 = z(D.thin.loc[idx]), z(D.dPC4.loc[idx])
    r = np.corrcoef(zt, z4)[0, 1]
    for w in (0, 0.25, 0.5, 1, 2, np.inf):
        y = z4 if np.isinf(w) else z(zt + w * z4)
        b, se, p = ols(y.reindex(D.index), D[sc], idx)
        rows.append(dict(disorder=dis, arm=arm, weight_dPC4=w, est=b, se=se, z=b / se, p=p, n=len(idx), r_thin_dPC4=r))
P = pd.DataFrame(rows)
P.to_csv(TAB / "single_factor_prs.tsv", sep="\t", index=False)
P["cell"] = P.apply(lambda q: f"{q.est:+.3f} (z {q.z:+.1f})", axis=1)
print(P.pivot_table(index=["disorder", "arm"], columns="weight_dPC4", values="cell", aggfunc="first").to_string())

# ---- (b) proportionality across scores (tetrad), family bootstrap -----------------------------------
TT = []
for arm in ("EUR", "pooled"):
    sS, sM = SC[("SCZ", arm)], SC[("MDD", arm)]
    idx = arm_index(arm, sS).intersection(arm_index(arm, sM))
    fams = D.loc[idx, "family"].astype(str); uf = fams.unique(); members = {f: fams.index[fams == f] for f in uf}

    def four(ii):
        b1S = ols(D.thin, D[sS], ii)[0]; b1M = ols(D.thin, D[sM], ii)[0]
        b2S = ols(D.dPC4, D[sS], ii)[0]; b2M = ols(D.dPC4, D[sM], ii)[0]
        return b1S, b1M, b2S, b2M, b2S * b1M - b2M * b1S

    est = four(idx); bs = []
    for _ in range(NB):
        pick = rng.choice(uf, len(uf), replace=True)
        ii = pd.Index(np.concatenate([members[f] for f in pick]))
        bs.append(four(ii))
    bs = np.array(bs)
    lo, hi = np.percentile(bs[:, 4], [2.5, 97.5]); p_t = 2 * min((bs[:, 4] <= 0).mean(), (bs[:, 4] >= 0).mean())
    TT.append(dict(arm=arm, n=len(idx), b_thin_SCZ=est[0], b_thin_MDD=est[1], b_dPC4_SCZ=est[2], b_dPC4_MDD=est[3],
                   ratio_SCZ=est[2] / est[0], ratio_MDD=est[3] / est[1], tetrad=est[4], tetrad_lo=lo, tetrad_hi=hi,
                   p_tetrad_boot=min(1.0, p_t), n_boot=NB))
    print(f"{arm}: dPC4/thin ratio SCZ {est[2]/est[0]:+.2f}, MDD {est[3]/est[1]:+.2f}; tetrad {est[4]:+.5f} "
          f"[{lo:+.5f}, {hi:+.5f}] p {min(1.0, p_t):.3f}", flush=True)
pd.DataFrame(TT).to_csv(TAB / "single_factor_tetrad.tsv", sep="\t", index=False)

# ---- (c) twins: h2 of dPC4 / dPC2 and rA with the thinning rate --------------------------------------
T = pd.read_parquet(D3 / "work" / "traits.parquet")
pairs = pd.read_csv(D3 / "work" / "pairs.csv")
BR = ["male", "age_first", "age_span", "n_visits"]
for t in ("dPC4", "dPC2"):
    T[t] = pcs[t].reindex(T.index)
    T[f"{t}_r"] = residualise(T, t, BR, ["site"])
T["thin_r"] = -T["dCT_r"]                                   # dCT_r: higher = slower; flip so higher = faster
KEYS = ["h2_1", "h2_2", "rA", "rE", "rP", "rP_A"]
NBT = max(100, NB // 3); TW = []
for t in ("dPC4_r", "dPC2_r"):
    mz = ace.pair_arrays(T, pairs, [t], "MZ"); dz = ace.pair_arrays(T, pairs, [t], "DZ")
    u = ace.fit(mz, dz, "AE"); ub = ace.bootstrap(mz, dz, "AE", NBT, seed=21, start=u["theta"], keys=["h2_1"])
    corr = ace.twin_corrs(mz, dz)
    mz2 = ace.pair_arrays(T, pairs, ["thin_r", t], "MZ"); dz2 = ace.pair_arrays(T, pairs, ["thin_r", t], "DZ")
    f = ace.fit(mz2, dz2, "AE"); bs = ace.bootstrap(mz2, dz2, "AE", NBT, seed=22, start=f["theta"], keys=KEYS)
    q = lambda s: np.nanpercentile(s, [2.5, 97.5])
    row = dict(trait=t, **{k: v for k, v in corr.items()}, h2_AE=u["h2_1"], h2_lo=q(ub.h2_1)[0], h2_hi=q(ub.h2_1)[1],
               rA_with_thin=f["rA"], rA_lo=q(bs.rA)[0], rA_hi=q(bs.rA)[1], rE_with_thin=f["rE"], rP_with_thin=f["rP"],
               h2_thin_biv=f["h2_1"], n_mz=f["n_mz"], n_dz=f["n_dz"], n_boot=NBT)
    TW.append(row)
    print(f"twin {t}: h2 {u['h2_1']:.2f} [{row['h2_lo']:.2f}, {row['h2_hi']:.2f}]; rA with faster thinning "
          f"{f['rA']:+.2f} [{row['rA_lo']:+.2f}, {row['rA_hi']:+.2f}]; rE {f['rE']:+.2f}; rP {f['rP']:+.2f}", flush=True)
pd.DataFrame(TW).to_csv(TAB / "single_factor_twin.tsv", sep="\t", index=False)
