"""D3 part 2: within-family (sibling-difference) polygenic score effects on the thinning rate.

    python directions/d3_twin_family/05_within_family_prs.py   (env abcd-spatial)

Between/within (Mundlak) decomposition in families with >= 2 genotyped, imaged children:
  y ~ PRS_between + PRS_within + sex + baseline_age + n_visits + PC1-10     (OLS, family-clustered SE)
  PRS_between = family mean of the score; PRS_within = child's deviation from it.
The within-family coefficient is free of population stratification, assortative mating
and indirect (genetic-nurture) effects; the between coefficient carries all of them.
MZ co-twins share a score, so one member of each MZ pair is dropped before forming
families (they add nothing within and would double-count the family mean).
Reported per score and outcome:
  pop_full     population beta in the full pooled arm (n = 8,596), same covariates. This is
               the fig1 model by OLS rather than LMM, so it matches Figure 1f to ~0.002.
  pop_fam      population beta in the family subsample (the like-for-like comparator)
  within, between, p(within = between) (cluster Wald)
  ratio        within / pop_full, CI = within-family 95% CI / pop_full (pop_full is
               estimated ~6x more precisely, so its uncertainty is ignored). The quantity
               DIRECTIONS.md D3 names: ~1 = direct genetic effect; ~0.5 = half of the
               population association is stratification / nurture. (A ratio to pop_fam
               was tried first; pop_fam is itself near zero in resamples, so its bootstrap
               CI spans +-20 and says nothing.)
  mde80        within-family effect detectable at 80% power, two-sided alpha .05
               (2.80 x within SE), and mde80 / |pop_full|.
Scores: pooled-arm (_zanc) SCZ 2025 and MDD, PRS-CS and SBayesRC, each z-scored in the
full arm. EA scores are cluster-only; see README for the rsync that adds them.
Outcomes: dCT (primary), CT0 (baseline thickness; the scores are null there by design).
Writes results/within_family_prs.tsv and results/within_family_sample.tsv.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402

from common import RES, SCORES, WORK, zscore  # noqa: E402

T = pd.read_parquet(WORK / "traits.parquet")
pairs = pd.read_csv(WORK / "pairs.csv")
G = T[T.gen_family.notna()].copy()                       # genotyped analysis set (pooled arm)
PCS = [f"PC{i}" for i in range(1, 11)]
for name, path in SCORES.items():
    s = pd.read_csv(path, sep=r"\s+")
    s["subject"] = "sub-NDARINV" + s.IID.str.replace("sub-", "", regex=False)
    G[name] = zscore(s.set_index("subject").SCORESUM.reindex(G.index))
# drop one member of every MZ pair (all MZ pairs in the release, not only the analysed ones)
mz_all = pairs[pairs.type == "MZ"]
G = G[~G.index.isin(mz_all.id2)]
G["fam_n"] = G.groupby("gen_family").gen_family.transform("size")
F = G[G.fam_n >= 2].copy()
sample = pd.DataFrame([dict(set="pooled arm, genotyped + imaged", n_children=len(T[T.gen_family.notna()])),
                       dict(set="after dropping MZ co-twins", n_children=len(G)),
                       dict(set="families with >= 2 children", n_children=len(F),
                            n_families=int(F.gen_family.nunique()),
                            size_2=int((F.groupby("gen_family").size() == 2).sum()),
                            size_3plus=int((F.groupby("gen_family").size() >= 3).sum()))])
sample.to_csv(RES / "within_family_sample.tsv", sep="\t", index=False)
print(sample.to_string(index=False))

COV = ["male", "age_first", "n_visits"] + PCS


def ols(d: pd.DataFrame, y: str, xs: list[str]):
    d = d.dropna(subset=[y, *xs])
    X = sm.add_constant(d[xs])
    return sm.OLS(d[y], X).fit(cov_type="cluster", cov_kwds={"groups": d.gen_family.astype(str)}), d


def decompose(d: pd.DataFrame, sc: str) -> pd.DataFrame:
    d = d.copy()
    d["b"] = d.groupby("gen_family")[sc].transform("mean")
    d["w"] = d[sc] - d.b
    return d


rows = []
for sc in SCORES:
    Fd = decompose(F.dropna(subset=[sc]), sc)
    for y in ("dCT", "CT0"):
        mp, dp = ols(G, y, [sc] + COV)
        mf, df_ = ols(Fd, y, [sc] + COV)
        mw, dw = ols(Fd, y, ["w", "b"] + COV)
        wald = mw.t_test("w - b = 0")
        rows.append(dict(score=sc, outcome=y, n_full=int(mp.nobs), n_fam=int(mw.nobs),
                         n_families=int(dw.gen_family.nunique()),
                         pop_full=mp.params[sc], pop_full_se=mp.bse[sc], pop_full_p=mp.pvalues[sc],
                         pop_fam=mf.params[sc], pop_fam_se=mf.bse[sc], pop_fam_p=mf.pvalues[sc],
                         within=mw.params["w"], within_se=mw.bse["w"], within_p=mw.pvalues["w"],
                         between=mw.params["b"], between_se=mw.bse["b"], between_p=mw.pvalues["b"],
                         p_within_eq_between=float(wald.pvalue),
                         p_within_eq_pop=float(2 * __import__("scipy.stats").stats.norm.sf(
                             abs(mw.params["w"] - mp.params[sc]) / mw.bse["w"])),
                         ratio=mw.params["w"] / mp.params[sc],
                         ratio_lo=min((mw.params["w"] - 1.96 * mw.bse["w"]) / mp.params[sc],
                                      (mw.params["w"] + 1.96 * mw.bse["w"]) / mp.params[sc]),
                         ratio_hi=max((mw.params["w"] - 1.96 * mw.bse["w"]) / mp.params[sc],
                                      (mw.params["w"] + 1.96 * mw.bse["w"]) / mp.params[sc]),
                         mde80=2.80 * mw.bse["w"], mde80_over_pop=2.80 * mw.bse["w"] / abs(mp.params[sc]),
                         sd_within_score=float(Fd.w.std())))
        r = rows[-1]
        print(f"{sc:12s} {y:4s} pop {r['pop_full']:+.3f} (p {r['pop_full_p']:.3g}) | fam {r['pop_fam']:+.3f} | "
              f"within {r['within']:+.3f}±{r['within_se']:.3f} between {r['between']:+.3f} | "
              f"ratio {r['ratio']:.2f} [{r['ratio_lo']:.2f}, {r['ratio_hi']:.2f}] MDE {r['mde80']:.3f}", flush=True)
pd.DataFrame(rows).to_csv(RES / "within_family_prs.tsv", sep="\t", index=False, float_format="%.4g")
