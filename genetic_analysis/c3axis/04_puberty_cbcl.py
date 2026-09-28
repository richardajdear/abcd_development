"""04 -- slope components vs puberty (PDS) and CBCL (README_HPC.md 8.6).

Components (child scores, standardised):
  PC1-5  plain PCA of standardised HCP-MMP slope BLUPs, both hemispheres (358 parcels)
  rc1-7  PCA of LEFT-hemisphere row-centred slopes (identical to DME; see 01)
Signs: PC1 high = faster global thinning; PC2/rc1 +C1, PC3/rc2 +C2, rc5 +C3; rc4/rc7 high = amplified
normative dCT contrast (positive with the child's projection on dCT); the rest by the A-P axis.

Puberty: Pubertal Development Scale, parent (ph_p_pds) and youth (ph_y_pds); sex-specific mean.
  timing = within-sex, within-wave residual on cubic age, standardised; baseline and mean over waves.
  tempo  = per-child OLS slope of PDS on age (>= 3 waves, age SD > 0.5), adjusted for baseline timing.
CBCL (p/mh_p_cbcl.tsv): baseline log1p raw sums + p-factor (PC1 of 8 syndromes, fitted at baseline);
  onset = below threshold at ses-00A and at/above at any of 05A-07A vs below at every wave
  (T >= 60 broadband, >= 65 syndrome/DSM; p-factor: top decile late and not at baseline vs < 75th pct always),
  as in ahba_pls/code/24_casecontrol_maps.py.
Model: outcome ~ component + sex + site + age_first + age_span + n_visits; OLS (standardised outcome)
  or logistic; family-clustered SE; BH-FDR within each family (puberty, cbcl_baseline, cbcl_onset).

Writes table_components_puberty_cbcl.tsv (coefficients only). Individual-level data are read
from the gitignored release copy and never written. The PDS tables are looked for at the
release root and under p/pds, y/pds.
Run from the repo root:  python genetic_analysis/c3axis/04_puberty_cbcl.py
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import HERE, REPO, bilat, child_covariates, cov, load_inputs, orient, pca_raw, row_center, sp, zs

warnings.filterwarnings("ignore")
REL = REPO / "abcd-data-release-7.0"
Wh, WL, REFL, maps, ylab = load_inputs()
cv = child_covariates(WL.index)
fam = pd.read_parquet(REPO / "out/thickness_hcp_70_aa6e91efba82/model_table.parquet", columns=["subject", "family_id"]
                      ).drop_duplicates("subject").set_index("subject").family_id.reindex(WL.index)

# ---- components ----------------------------------------------------------------------------
Zh = zs(Wh); U, S_, Vt = np.linalg.svd(Zh.to_numpy() - Zh.to_numpy().mean(0), full_matrices=False)
SCp = pd.DataFrame(U[:, :5] * S_[:5], index=Wh.index, columns=[f"PC{i}" for i in range(1, 6)])
Lp = pd.DataFrame(Vt[:5].T, index=Wh.columns, columns=SCp.columns)
bm = lambda s: bilat(s.to_frame("v"), ylab).v
sgn = {"PC1": -np.sign(np.corrcoef(SCp.PC1, Wh.mean(1))[0, 1]),
       "PC2": np.sign(sp(bm(Lp.PC2).reindex(maps.index), maps.C1)),
       "PC3": np.sign(sp(bm(Lp.PC3).reindex(maps.index), maps.C2))}
ori = cov._orient(Lp.to_numpy().copy(), Lp.index)
for i, c in enumerate(SCp):
    SCp[c] *= sgn.get(c, np.sign(np.dot(ori[:, i], Lp[c])))

WLc = row_center(WL); Lr, _ = pca_raw(WLc)
SCr = pd.DataFrame(WLc.to_numpy() @ Lr.iloc[:, :7].to_numpy(), index=WL.index, columns=[f"rc{i}" for i in range(1, 8)])
Lr = Lr.iloc[:, :7]; Lr.columns = SCr.columns
m = REFL.dCT_lh.dropna(); X = WL[m.index].to_numpy(); v = (m - m.mean()).to_numpy()
proj_dct = (X - X.mean(1, keepdims=True)) @ v / (v @ v)
tgt = {"rc1": "C1", "rc2": "C2", "rc5": "C3"}
for c in SCr:
    if c in tgt:
        s = np.sign(sp(Lr[c].reindex(REFL.index), REFL[tgt[c]]))
    elif c in ("rc4", "rc7"):
        s = np.sign(np.corrcoef(SCr[c], proj_dct)[0, 1])
    else:
        s = np.sign(np.dot(cov._orient(Lr[[c]].to_numpy().copy(), Lr.index)[:, 0], Lr[c]))
    SCr[c] *= s
COMP = SCp.join(SCr); COMP = (COMP - COMP.mean()) / COMP.std()
CNAME = {"PC1": "PC1 (global, high = faster)", "PC2": "PC2 (C1)", "PC3": "PC3 (C2)", "PC4": "PC4", "PC5": "PC5",
         "rc1": "rc1 (C1)", "rc2": "rc2 (C2)", "rc3": "rc3 (CT)", "rc4": "rc4 (dCT contrast)", "rc5": "rc5 (C3)",
         "rc6": "rc6", "rc7": "rc7 (dCT contrast)"}


# ---- puberty ---------------------------------------------------------------------------------
def find(name: str):
    for p in (REL / name, REL / "p" / "pds" / name, REL / "y" / "pds" / name):
        if p.exists():
            return p
    raise FileNotFoundError(name)


def puberty(fname: str, pre: str, label: str) -> pd.DataFrame:
    d = pd.read_csv(find(fname), sep="\t", low_memory=False,
                    usecols=["participant_id", "session_id", f"{pre}_age", f"{pre}__f_mean", f"{pre}__m_mean"])
    d["subject"] = "sub-NDARINV" + d.participant_id.str.replace("sub-", "", regex=False)
    for c in (f"{pre}_age", f"{pre}__f_mean", f"{pre}__m_mean"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d.subject.isin(WL.index)].copy()
    d["sex"] = cv.sex.reindex(d.subject).to_numpy()
    d["pds"] = np.where(d.sex == "F", d[f"{pre}__f_mean"], d[f"{pre}__m_mean"]); d["age"] = d[f"{pre}_age"]
    d = d.dropna(subset=["pds", "age", "sex"])

    def res(g):
        if len(g) < 50:
            return pd.Series(np.nan, index=g.index)
        Xv = np.vander(g.age - g.age.mean(), 4); b, *_ = np.linalg.lstsq(Xv, g.pds, rcond=None)
        r = g.pds - Xv @ b; return r / r.std()
    d["timing"] = d.groupby(["sex", "session_id"], group_keys=False)[["age", "pds"]].apply(res)
    out = pd.DataFrame(index=WL.index)
    out[f"{label}_timing_base"] = d[d.session_id == "ses-00A"].set_index("subject").timing
    out[f"{label}_timing_mean"] = d.groupby("subject").timing.mean()
    out[f"{label}_tempo"] = d.groupby("subject")[["age", "pds"]].apply(
        lambda g: np.polyfit(g.age, g.pds, 1)[0] if len(g) >= 3 and g.age.std() > 0.5 else np.nan)
    return out


PUB = puberty("ph_p_pds.tsv", "ph_p_pds", "pds_p").join(puberty("ph_y_pds.tsv", "ph_y_pds", "pds_y"))

# ---- CBCL ------------------------------------------------------------------------------------
k70 = dict(anxdep="synd__anxdep", withdep="synd__wthdep", somatic="synd__som", social="synd__soc",
           thought="synd__tho", attention="synd__attn", rulebreak="synd__rule", aggressive="synd__aggr",
           internal="synd__int", external="synd__ext", depress="dsm__dep", anxdisord="dsm__anx", adhd="dsm__adhd")
raw = {k: f"mh_p_cbcl__{v}_sum" for k, v in k70.items()} | {"totprob": "mh_p_cbcl_sum"}
tsc = {k: f"mh_p_cbcl__{v}_tscore" for k, v in k70.items()} | {"totprob": "mh_p_cbcl_tscore"}
cb = pd.read_csv(REL / "p" / "mh_p_cbcl.tsv", sep="\t", low_memory=False,
                 usecols=lambda c: c in {"participant_id", "session_id", *raw.values(), *tsc.values()})
cb["subject"] = "sub-NDARINV" + cb.participant_id.str.replace("sub-", "", regex=False)
cb["wave"] = cb.session_id.str.slice(4, 6).astype(int)
cb = cb[cb.subject.isin(WL.index)].drop_duplicates(["subject", "wave"])
for k, c in raw.items():
    cb[k] = np.log1p(pd.to_numeric(cb[c], errors="coerce"))
for k, c in tsc.items():
    cb[f"{k}_T"] = pd.to_numeric(cb[c], errors="coerce")
SYN = ["anxdep", "withdep", "somatic", "social", "thought", "attention", "rulebreak", "aggressive"]
base = cb[cb.wave == 0].set_index("subject")
mu, sdv = base[SYN].mean(), base[SYN].std()
ld = np.linalg.svd(((base[SYN].dropna() - mu) / sdv).to_numpy(), full_matrices=False)[2][0]
cb["pfactor"] = ((cb[SYN] - mu) / sdv) @ (ld * np.sign(ld.sum()))
base = cb[cb.wave == 0].set_index("subject")
BASEO = base[["totprob", "internal", "external", "depress", "adhd", "pfactor"]].add_prefix("cbcl_base_")
THR = {"totprob": 60, "internal": 60, "external": 60, "depress": 65, "anxdisord": 65, "withdep": 65,
       "anxdep": 65, "attention": 65, "rulebreak": 65}
late = cb[cb.wave.isin([5, 6, 7])]
ONS = pd.DataFrame(index=WL.index)


def onset(v0, vlate, vall, is_case_hi, is_ctrl):
    case = v0.index[is_case_hi[0](v0)].intersection(vlate.index[is_case_hi[1](vlate)])
    s = pd.Series(np.nan, index=WL.index); s[s.index.isin(vall.index[is_ctrl(vall)])] = 0.0
    s[s.index.isin(case)] = 1.0; return s


for k, t in THR.items():
    ONS[f"onset_{k}"] = onset(base[f"{k}_T"], late.groupby("subject")[f"{k}_T"].max(),
                              cb.groupby("subject")[f"{k}_T"].max(), (lambda x: x < t, lambda x: x >= t), lambda x: x < t)
p90, p75 = cb.pfactor.quantile(.9), cb.pfactor.quantile(.75)
ONS["onset_pfactor"] = onset(base.pfactor, late.groupby("subject").pfactor.max(), cb.groupby("subject").pfactor.max(),
                             (lambda x: x < p90, lambda x: x >= p90), lambda x: x < p75)

# ---- models ----------------------------------------------------------------------------------
D = COMP.join(cv[["male", "site", "age_first", "age_span", "n_visits"]]).join(PUB).join(BASEO).join(ONS)
D["family"] = fam.astype(str); D["site"] = D.site.astype(str)
OUTC = {"puberty": ["pds_p_timing_mean", "pds_p_timing_base", "pds_p_tempo", "pds_y_timing_mean", "pds_y_tempo"],
        "cbcl_baseline": list(BASEO.columns), "cbcl_onset": list(ONS.columns)}
covs = "male + C(site) + age_first + age_span + n_visits"
rows = []
for fam_name, outs in OUTC.items():
    for y in outs:
        tb = ["pds_p_timing_base"] if "p_tempo" in y else (["pds_y_timing_base"] if "y_tempo" in y else [])
        for c in COMP.columns:
            d = D[[y, c, "male", "site", "age_first", "age_span", "n_visits", "family"] + tb].dropna()
            f = f"{y} ~ {c} + {covs}" + "".join(f" + {t}" for t in tb)
            g = d.family.astype("category").cat.codes
            if fam_name == "cbcl_onset":
                mfit = smf.logit(f, d).fit(disp=0, cov_type="cluster", cov_kwds={"groups": g})
                eff, metric = np.exp(mfit.params[c]), "OR per SD"
            else:
                d = d.assign(**{y: (d[y] - d[y].mean()) / d[y].std()})
                mfit = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": g})
                eff, metric = mfit.params[c], "beta (SD/SD)"
            rows.append(dict(family=fam_name, outcome=y, component=c, label=CNAME[c], effect=eff, metric=metric,
                             est=mfit.params[c], se=mfit.bse[c], p=mfit.pvalues[c], n=len(d),
                             n_case=int(d[y].sum()) if fam_name == "cbcl_onset" else np.nan))
AS = pd.DataFrame(rows)
AS["q_fdr_family"] = AS.groupby("family").p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
AS.to_csv(HERE / "table_components_puberty_cbcl.tsv", sep="\t", index=False, float_format="%.5g")
print(AS.groupby("family").agg(tests=("p", "size"), p05=("p", lambda p: int((p < .05).sum())),
                               q05=("q_fdr_family", lambda q: int((q < .05).sum()))))
