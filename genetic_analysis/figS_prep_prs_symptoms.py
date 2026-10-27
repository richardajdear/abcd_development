"""SI figure inputs: polygenic scores -> symptoms, and mediation through ΔCT.

Usage (repo root, after `Rscript genetic_analysis/fig1_prep_1lmm.R`, which writes
/tmp/hcp70_1lmm_blups.csv):
    python genetic_analysis/figS_prep_prs_symptoms.py

Outcomes and waves exactly as fig1_prep_cbcl.py (CBCL parent report and PQ-BC youth
report; log1p raw scores; BASE = ses-00A, LATE = mean of ses-05A..07A, ages ~15-17).
Scores as directions/d2_prs_age_puberty/01_build_tables.py: EUR-weight scores
z-scored within the EUR arm (4,308 children); multi-ancestry-weight scores z within
ancestry cluster (_zanc) in the pooled arm. Only SCZ (2025) and MDD scores are on
this machine; PRS-CS primary, SBayesRC as a method check.

1. Direct associations (OLS, family-clustered SEs; β in SD of y per SD of score):
     baseline   y_BASE ~ PRS + age_BASE + sex + site + PC1-10
     late       y_LATE ~ PRS + age_LATE + sex + site + PC1-10
     change     y_LATE ~ PRS + y_BASE + age_LATE + sex + site + PC1-10   (Figure 1g model)
   BH-FDR over the 84 PRS-CS fits.
2. Mediation through ΔCT (product of coefficients; same covariates C as 'change'):
     a      ΔCT_z ~ PRS + C
     b, c'  y_LATE ~ PRS + ΔCT_z + C
     c      y_LATE ~ PRS + C                       (total = the 'change' fit)
     indirect = a·b, proportion = a·b / c.  95% CIs: 1,000 bootstrap replicates
   resampling families (seed 1). ΔCT is the shrunk single-LMM slope BLUP (Figure 1),
   so b, and hence a·b, is attenuated by its unreliability (mean reliability ~0.2-0.4).
Writes coefficient tables only:
   fig1_inputs/hcp70_prs_symptom_assoc.tsv, fig1_inputs/hcp70_prs_mediation.tsv
"""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
IN = REPO / "genetic_analysis/fig1_inputs"
BLUPS = Path("/tmp/hcp70_1lmm_blups.csv")
N_BOOT, SEED = 1000, 1
OUTC = ["pfactor", "internal", "external", "depress", "thought", "ple_count", "ple_distress"]
SC = {("SCZ", "PRS-CS", "EUR"): "genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_EUR/score_SCZ25EUR_prscs.profile",
      ("SCZ", "PRS-CS", "pooled"): "genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_META/_zanc/score_SCZ25META_prscs.profile",
      ("MDD", "PRS-CS", "EUR"): "legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_eur/score_MDDeur_prscs.profile",
      ("MDD", "PRS-CS", "pooled"): "legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_pooled/_zanc/score_MDDpooled_prscs.profile",
      ("SCZ", "SBayesRC", "EUR"): "genetic_analysis/work/scores_scz2025/SBayesRC/SCZ25_EUR/score_SCZ25EUR_sbrc.profile",
      ("SCZ", "SBayesRC", "pooled"): "genetic_analysis/work/scores_scz2025/SBayesRC/SCZ25_META/_zanc/score_SCZ25META_sbrc.profile",
      ("MDD", "SBayesRC", "EUR"): "legacy/hpc_v2/work/results_v2/prs_final/SBayesRC/MDD_eur/score_MDDeur_sbrc.profile",
      ("MDD", "SBayesRC", "pooled"): "legacy/hpc_v2/work/results_v2/prs_final/SBayesRC/MDD_pooled/_zanc/score_MDDpooled_sbrc.profile"}
PCS = [f"PC{i}" for i in range(1, 11)]


def load() -> tuple[pd.DataFrame, str]:
    spec = importlib.util.spec_from_file_location("p", REPO / "genetic_analysis/fig1_prep_cbcl.py")
    p = importlib.util.module_from_spec(spec); spec.loader.exec_module(p)
    ns = p.build_outcomes(); D = ns["D"].copy(); FU = ns["FU"]
    pps = pd.read_csv(ns["R70"] / "mh_y_pps.tsv", sep="\t", usecols=["participant_id", "session_id", *p.PPS.values()])
    pps["subject"] = "sub-NDARINV" + pps.participant_id.str.replace("sub-", "", regex=False)
    pps["wave"] = pps.session_id.str.slice(4, 6).astype(int)
    pps = pps[pps.subject.isin(D.index)].drop_duplicates(["subject", "wave"])
    for k, c in p.PPS.items():
        pps[k] = np.log1p(pd.to_numeric(pps[c], errors="coerce"))
    for name, waves in ns["WAVES"].items():
        D = D.join(pps[pps.wave.isin(waves)].groupby("subject")[list(p.PPS)].mean().add_suffix(f"_{name}"))
    b = pd.read_csv(BLUPS).set_index("subject")
    z = lambda s: (s - s.mean()) / s.std()
    D["dct"] = z(b.re_slope).reindex(D.index)
    D["gid"] = D.index.str.replace("sub-NDARINV", "sub-", regex=False)
    pcs = pd.read_csv(REPO / "genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno/covar_quant.txt",
                      sep=r"\s+").set_index("IID")
    D = D.join(pcs[PCS], on="gid")
    return D, FU


def scored(D, path, arm, eur):
    s = pd.read_csv(REPO / path, sep=r"\s+").set_index("IID").SCORESUM
    d = D.copy(); d["prs"] = s.reindex(d.gid).values
    if arm == "EUR":
        d = d[d.gid.isin(eur)]
    d = d.dropna(subset=["prs", *PCS])
    d["prs"] = (d.prs - d.prs.mean()) / d.prs.std()
    return d


def design(d, cols):
    X = pd.get_dummies(d[cols], columns=[c for c in ("sex", "site") if c in cols], drop_first=True, dtype=float)
    X.insert(0, "const", 1.0)
    return X.to_numpy(float)


def main():
    D, FU = load()
    eur = set(pd.read_csv(REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep", sep=r"\s+", header=None)[1])
    rows, med = [], []
    rng = np.random.default_rng(SEED)
    for (trait, meth, arm), path in SC.items():
        d0 = scored(D, path, arm, eur)
        for o in OUTC:
            yl, yb = f"{o}_{FU}", f"{o}_BASE"
            for model, y, extra in (("baseline", yb, ["age_BASE"]), ("late", yl, [f"age_{FU}"]),
                                    ("change", yl, [f"age_{FU}", yb])):
                d = d0.dropna(subset=[y, *extra])
                m = smf.ols(f"{y} ~ prs + C(sex) + C(site) + " + " + ".join(extra + PCS), d).fit(
                    cov_type="cluster", cov_kwds={"groups": d.family})
                sd = d[y].std()
                rows.append(dict(trait=trait, method=meth, arm=arm, outcome=o, model=model, n=int(m.nobs),
                                 beta=m.params["prs"] / sd, se=m.bse["prs"] / sd, p=m.pvalues["prs"]))
            if meth != "PRS-CS":
                continue
            # ---- mediation through ΔCT
            d = d0.dropna(subset=[yl, yb, f"age_{FU}", "dct"]).copy()
            sd = d[yl].std(); d["y"] = d[yl] / sd
            cov = f" + {yb} + age_{FU} + C(sex) + C(site) + " + " + ".join(PCS)
            fa = smf.ols("dct ~ prs" + cov, d).fit(cov_type="cluster", cov_kwds={"groups": d.family})
            fb = smf.ols("y ~ prs + dct" + cov, d).fit(cov_type="cluster", cov_kwds={"groups": d.family})
            fc = smf.ols("y ~ prs" + cov, d).fit(cov_type="cluster", cov_kwds={"groups": d.family})
            Cm = design(d, [yb, f"age_{FU}", "sex", "site", *PCS])
            P, M, Y = d.prs.to_numpy(), d.dct.to_numpy(), d.y.to_numpy()
            fam = d.family.astype(str).to_numpy()
            ufam, inv = np.unique(fam, return_inverse=True)
            members = np.split(np.argsort(inv, kind="stable"), np.cumsum(np.bincount(inv))[:-1])
            boot = np.empty((N_BOOT, 3))
            for r in range(N_BOOT):
                idx = np.concatenate([members[i] for i in rng.integers(0, len(ufam), len(ufam))])
                Xa = np.column_stack([P[idx], Cm[idx]])
                a = np.linalg.lstsq(Xa, M[idx], rcond=None)[0][0]
                bb = np.linalg.lstsq(np.column_stack([P[idx], M[idx], Cm[idx]]), Y[idx], rcond=None)[0]
                cc = np.linalg.lstsq(Xa, Y[idx], rcond=None)[0][0]
                boot[r] = (a * bb[1], cc, a * bb[1] / cc)
            ab = fa.params["prs"] * fb.params["dct"]
            q = lambda v: np.percentile(v, [2.5, 97.5])
            med.append(dict(trait=trait, arm=arm, outcome=o, n=int(fb.nobs),
                            a=fa.params["prs"], a_se=fa.bse["prs"], a_p=fa.pvalues["prs"],
                            b=fb.params["dct"], b_se=fb.bse["dct"], b_p=fb.pvalues["dct"],
                            c_total=fc.params["prs"], c_se=fc.bse["prs"], c_p=fc.pvalues["prs"],
                            c_direct=fb.params["prs"], c_direct_se=fb.bse["prs"],
                            indirect=ab, ind_lo=q(boot[:, 0])[0], ind_hi=q(boot[:, 0])[1],
                            ind_p_boot=2 * min((boot[:, 0] <= 0).mean(), (boot[:, 0] >= 0).mean()),
                            prop=ab / fc.params["prs"], prop_lo=q(boot[:, 2])[0], prop_hi=q(boot[:, 2])[1]))
            print(trait, arm, o, f"a {fa.params['prs']:+.3f} b {fb.params['dct']:+.3f} ab {ab:+.4f} c {fc.params['prs']:+.3f}",
                  flush=True)
    R = pd.DataFrame(rows)
    prim = R.method == "PRS-CS"
    R["q_fdr"] = np.nan
    R.loc[prim, "q_fdr"] = multipletests(R.loc[prim, "p"], method="fdr_bh")[1]
    R["source"] = "7.0 CBCL + PQ-BC; PRS as directions/d2"
    R.to_csv(IN / "hcp70_prs_symptom_assoc.tsv", sep="\t", index=False, float_format="%.5g")
    Mt = pd.DataFrame(med); Mt["n_boot"] = N_BOOT
    Mt.to_csv(IN / "hcp70_prs_mediation.tsv", sep="\t", index=False, float_format="%.5g")
    print(len(R), "association fits;", len(Mt), "mediation fits")


if __name__ == "__main__":
    main()
