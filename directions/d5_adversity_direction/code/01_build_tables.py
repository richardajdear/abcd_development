"""01_build_tables.py -- D5 inputs: environment, polygenic scores, symptoms, and
per-wave whole-cortex thickness, one table per child and one per imaging wave.

    PY=~/.claude-science/conda/envs/abcd-spatial/bin/python
    $PY directions/d5_adversity_direction/code/01_build_tables.py      # repo root

Writes INDIVIDUAL-LEVEL tables to out/d5_adversity_direction/ (gitignored; ABCD
DUC -- never commit them):
    child.parquet   one row per child in the 8,596-child genetics analysis set
    waves.parquet   one row per child x imaging wave (v0/v2/v4/v6)
and a coverage summary (counts only) to results/table_d5_coverage.tsv.

Sources
  phenotype + covariates  genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno/
      global_slope_1lmm, baseline_thickness_1lmm (single LMM on the per-scan
      cortical mean, the Figure-1 traits), sex, site, baseline_age, n_visits, PC1-10,
      FID = family id.  Same children, same covariates as tools/prs_assoc.R.
  polygenic scores (PRS-CS, the Figure-1f method; matched-arm rule):
      pooled arm (all ancestries, within-ancestry z score '_zanc'):
          SCZ 2025 META    genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_META/_zanc/
          MDD pooled       legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_pooled/_zanc/
      EUR arm (raw score, eur_anchor.keep children only):
          SCZ 2025 EUR     .../PRSCS/SCZ25_EUR/
          MDD EUR          legacy/.../PRSCS/MDD_eur/
      EUR membership: legacy/hpc/work/results/ancestry/eur_anchor.keep (the list
      the cluster association jobs read).
  environment, measured at study baseline (ses-00A) unless stated:
      income        ab_p_demo__income__hhold_001, 10 ordered bands (777/999 -> NA)
      parent_edu    max(respondent, partner) highest grade/degree, codes 0-21
      adi           le_l_adi__addr1__national_prcnt (higher = more deprived)
      conflict_y    fc_y_fes__confl_mean   (Family Environment Scale, youth)
      conflict_p    fc_p_fes__confl_mean   (parent)
      bad_events_y  mh_y_ple__exp__bad_count at ses-01A (first PLE wave; the
      bad_events_p  mh_p_ple__exp__bad_count   questionnaire version changes at
                    year 4, so later waves are not on the same scale)
    composites (mean of available z scores, >= 2 components required):
      ses        = z(income), z(parent_edu), -z(adi)          higher = advantaged
      adversity  = z(conflict_y), z(conflict_p), z(bad_events_y), z(bad_events_p)
  symptoms: parent CBCL, log1p raw scale sums, built by executing the outcome
      block of ahba_pls/code/23_cbcl_explore.py (as genetic_analysis/fig1_prep_cbcl.py
      does), so BASE (ses-00A) and LATE (mean of ses-05A..07A) match Figure 1g.
  per-wave table: whole-cortex mean thickness per scan (mean over the 358 parcels
      of out/thickness_hcp_70_aa6e91efba82/model_table.parquet, as in
      genetic_analysis/fig1_prep_1lmm.R) and CBCL at the same session.
  scan quality: log1p(mr_y_qc__post__aut__smri__topodfct_count), FreeSurfer surface
      topological defects, per scan (waves) and averaged over the child's scans
      (child.qc_defects) -- the 7.0 stand-in for the Euler number.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
D5 = HERE.parent
REPO = D5.parents[1]
R70 = REPO / "abcd-data-release-7.0"
RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
PHENO = REPO / "genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno"
OUT = REPO / "out" / "d5_adversity_direction"
RES = D5 / "results"
OUT.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

SCORES = {
    "SCZ_pooled": REPO / "genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_META/_zanc/score_SCZ25META_prscs.profile",
    "SCZ_eur": REPO / "genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_EUR/score_SCZ25EUR_prscs.profile",
    "MDD_pooled": REPO / "legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_pooled/_zanc/score_MDDpooled_prscs.profile",
    "MDD_eur": REPO / "legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_eur/score_MDDeur_prscs.profile",
}
# Educational attainment (Okbay 2016, EUR discovery -> EUR arm only, as in the cluster
# control panel) -- scored on CSD3 only; picked up once copied (README: "EA step").
EA_DIR = REPO / "legacy/hpc_v2/work/results_v2/prs_final/PRSCS/EA"
_ea = sorted(EA_DIR.glob("score_*.profile")) if EA_DIR.exists() else []
if _ea:
    SCORES["EA_eur"] = _ea[0]
EUR_KEEP = REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep"
CBCL_WAVE = {"depress": "mh_p_cbcl__dsm__dep_sum", "internal": "mh_p_cbcl__synd__int_sum",
             "external": "mh_p_cbcl__synd__ext_sum"}
CHILD_OUTCOMES = ["pfactor", "internal", "external", "depress"]


def tail8(s: pd.Series) -> pd.Series:
    """Join key: the 8-character GUID suffix shared by sub-XXXX and sub-NDARINVXXXX."""
    return s.astype(str).str.extract(r"([A-Z0-9]{8})$")[0]


def read_release(name: str, cols: list[str], session: str | None = None) -> pd.DataFrame:
    p = R70 / f"{name}.tsv"
    if not p.exists():
        p = next(R70.rglob(f"{name}.tsv"))
    d = pd.read_csv(p, sep="\t", usecols=["participant_id", "session_id", *cols], dtype=str)
    if session:
        d = d[d.session_id == session]
    d["_t"] = tail8(d.participant_id)
    for c in cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.drop_duplicates("_t").set_index("_t")[cols] if session else d


def z(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std()


def composite(df: pd.DataFrame, parts: dict[str, int]) -> pd.Series:
    Z = pd.concat({k: sign * z(df[k]) for k, sign in parts.items()}, axis=1)
    return Z.mean(axis=1).where(Z.notna().sum(axis=1) >= 2)


def build_outcomes() -> pd.DataFrame:
    """Execute 23_cbcl_explore.py up to its model fits; return its per-child table D."""
    src_path = REPO / "ahba_pls/code/23_cbcl_explore.py"
    src = src_path.read_text()
    ns = {"__file__": str(src_path), "__name__": "cbcl_prefix"}
    argv, sys.argv = sys.argv, [str(src_path), "70"]
    try:
        exec(compile(src[: src.index("\ndef fit(")], str(src_path), "exec"), ns)
    finally:
        sys.argv = argv
    D, FU = ns["D"].copy(), ns["FU"]
    keep = {f"{o}_BASE": f"{o}_base" for o in CHILD_OUTCOMES}
    keep |= {f"{o}_{FU}": f"{o}_late" for o in CHILD_OUTCOMES}
    keep |= {f"age_{FU}": "age_late"}
    D = D[list(keep)].rename(columns=keep)
    D["_t"] = tail8(pd.Series(D.index, index=D.index))
    return D.set_index("_t")


def main() -> None:
    # ---- genetics analysis set: phenotype + covariates ----------------------------
    rd = lambda f: pd.read_csv(PHENO / f, sep=r"\s+", dtype={"FID": str, "IID": str})
    C = rd("phenotypes_gcta.txt").merge(rd("covar_quant.txt"), on=["FID", "IID"]).merge(
        rd("covar_categorical.txt"), on=["FID", "IID"])
    C["_t"] = tail8(C.IID)
    C = C.set_index("_t").rename(columns={"FID": "family_id"})
    assert len(C) == 8596 and C.index.is_unique, C.shape

    # ---- polygenic scores ----------------------------------------------------------
    eur = set(tail8(pd.read_csv(EUR_KEEP, sep=r"\s+", header=None, dtype=str)[1]))
    for k, p in SCORES.items():
        s = pd.read_csv(p, sep=r"\s+", dtype={"IID": str})
        s["_t"] = tail8(s.IID)
        C[k] = s.set_index("_t").SCORESUM.reindex(C.index)
        if k.endswith("_eur"):
            C.loc[~C.index.isin(eur), k] = np.nan
    C["eur"] = C.index.isin(eur)

    # ---- environment ---------------------------------------------------------------
    dm = read_release("ab_p_demo", ["ab_p_demo__income__hhold_001", "ab_p_demo__edu__slf_001",
                                    "ab_p_demo__edu__prtnr_001"], "ses-00A")
    dm = dm.where(dm < 700)                                   # 777 decline, 999 don't know
    E = pd.DataFrame(index=C.index)
    E["income"] = dm.ab_p_demo__income__hhold_001.reindex(C.index)
    E["parent_edu"] = dm[["ab_p_demo__edu__slf_001", "ab_p_demo__edu__prtnr_001"]].max(axis=1).reindex(C.index)
    E["adi"] = read_release("le_l_adi", ["le_l_adi__addr1__national_prcnt"], "ses-00A").iloc[:, 0].reindex(C.index)
    E["conflict_y"] = read_release("fc_y_fes", ["fc_y_fes__confl_mean"], "ses-00A").iloc[:, 0].reindex(C.index)
    E["conflict_p"] = read_release("fc_p_fes", ["fc_p_fes__confl_mean"], "ses-00A").iloc[:, 0].reindex(C.index)
    E["bad_events_y"] = read_release("mh_y_ple", ["mh_y_ple__exp__bad_count"], "ses-01A").iloc[:, 0].reindex(C.index)
    E["bad_events_p"] = read_release("mh_p_ple", ["mh_p_ple__exp__bad_count"], "ses-01A").iloc[:, 0].reindex(C.index)
    E["ses"] = composite(E, {"income": 1, "parent_edu": 1, "adi": -1})
    E["adversity"] = composite(E, {"conflict_y": 1, "conflict_p": 1, "bad_events_y": 1, "bad_events_p": 1})
    C = C.join(E)

    # ---- symptoms (Figure-1g definitions) -----------------------------------------
    C = C.join(build_outcomes())

    # ---- per-wave table: cortical mean thickness + CBCL at the same session --------
    mt = pd.read_parquet(RUN / "model_table.parquet", columns=["subject", "visit", "value", "age"])
    W = mt.groupby(["subject", "visit"], observed=True).agg(mean_ct=("value", "mean"), age=("age", "first")).reset_index()
    W["_t"] = tail8(W.subject)
    W["session_id"] = "ses-" + W.visit.astype(str).str.slice(1).str.zfill(2) + "A"
    cb = read_release("mh_p_cbcl", list(CBCL_WAVE.values()))
    cb = cb.rename(columns={v: k for k, v in CBCL_WAVE.items()})
    for k in CBCL_WAVE:
        cb[k] = np.log1p(cb[k])
    W = W.merge(cb[["_t", "session_id", *CBCL_WAVE]], on=["_t", "session_id"], how="left")
    # scan quality: FreeSurfer surface topological defects (the 7.0 stand-in for the
    # Euler number), log1p; per scan, and averaged over the child's scans
    qc = read_release("mr_y_qc__post__aut", ["mr_y_qc__post__aut__smri__topodfct_count"])
    qc["qc_defects"] = np.log1p(qc.mr_y_qc__post__aut__smri__topodfct_count)
    W = W.merge(qc[["_t", "session_id", "qc_defects"]].drop_duplicates(["_t", "session_id"]),
                on=["_t", "session_id"], how="left")
    W = W[W._t.isin(C.index)].rename(columns={"_t": "guid8"}).drop(columns="subject")
    W.to_parquet(OUT / "waves.parquet", index=False)
    C["qc_defects"] = W.groupby("guid8").qc_defects.mean().reindex(C.index)
    C.reset_index().rename(columns={"_t": "guid8"}).to_parquet(OUT / "child.parquet", index=False)

    # ---- coverage (counts only) ----------------------------------------------------
    cov = [dict(table="child", variable=c, n_nonmissing=int(C[c].notna().sum()), n_total=len(C))
           for c in [*SCORES, *E.columns, *[f"{o}_{w}" for o in CHILD_OUTCOMES for w in ("base", "late")],
                     "qc_defects"]]
    cov += [dict(table="waves", variable=f"{v}@{s}", n_nonmissing=int(W.loc[W.visit == s, v].notna().sum()),
                 n_total=int((W.visit == s).sum())) for s in ["v0", "v2", "v4", "v6"]
            for v in ["mean_ct", *CBCL_WAVE, "qc_defects"]]
    pd.DataFrame(cov).to_csv(RES / "table_d5_coverage.tsv", sep="\t", index=False)
    print(pd.DataFrame(cov).to_string(index=False))


if __name__ == "__main__":
    main()
