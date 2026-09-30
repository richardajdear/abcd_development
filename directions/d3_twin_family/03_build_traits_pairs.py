"""Per-child trait table and relative-pair tables for the twin and family designs.

    python directions/d3_twin_family/03_build_traits_pairs.py      (env abcd-spatial)

TRAITS (one row per imaged child, 8,716; work/traits.parquet, gitignored)
  dCT      single-LMM random slope of HCP-MMP mean thickness, z  (Figure-1 trait;
           HIGHER = SLOWER thinning, as in Figure 1)
  CT0      single-LMM random intercept (thickness at the sample mean age), z
  dT1T2    single-LMM random slope of DK mean T1w/T2w ratio, z  (higher = faster rise)
  pub_tempo    parent-report PDS slope on age (>= 3 waves, age SD > 0.5); timing =
  pub_timing   within-sex, within-wave residual on cubic age, mean over waves
               (definitions copied from genetic_analysis/c3axis/04_puberty_cbcl.py)
  <cbcl>_chg   CBCL change, <cbcl> in depress/internal/external/pfactor: outcome
               construction executed from ahba_pls/code/23_cbcl_explore.py (log1p raw
               sums; LATE = mean of ses-05A..07A), as genetic_analysis/fig1_prep_cbcl.py does
  <cog>_gain   NIH Toolbox uncorrected composite at ses-06A given ses-00A, <cog> in
               fluid/cryst/total (the only waves with the fluid composite)
Each trait also gets an ACE-ready residual  <trait>_r : OLS residual on the covariates
below, z-scored. Twins share site, age and visit schedule, so these must leave the mean
model or they are absorbed as shared environment (C).
  brain:  sex + age_first + age_span + n_visits + site
  dCT_rq: brain covariates + scan quality (QUALITY-ADJUSTED sensitivity): the child's mean
          log1p surface topological-defect count over the scans in the fit, and its OLS slope
          on age (y/mr_y_qc__post__aut.tsv, smri__topodfct_count; the 7.0 stand-in for the
          Euler number). Co-twins' defect counts correlate 0.31 (MZ) vs 0.16 (DZ), so image
          quality is itself twin-similar and could inflate A or confound slope-symptom links.
  puberty: sex + site (+ timing at baseline for tempo)
  CBCL:   y_BASE + age_LATE + sex + site      (so _r is residualised change)
  cognition: y_00A + age_00A + age_06A + sex + site

PAIRS (work/pairs.csv): genetically inferred relationships from y/gn_y_genrel.tsv
  (zygosity 1 = MZ, 2 = DZ, 3 = non-twin sibling; release data dictionary), both members
  imaged. One pair per birth event (triplets contribute one pair) and, for siblings, one
  pair per family; member order randomised with a fixed seed.
Committed: results/pair_counts.tsv, results/trait_availability.tsv.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import PHENO_DIR, R70, REPO, RES, WORK, residualise, to_img, zscore  # noqa: E402

RNG = np.random.default_rng(3)

# ---------------------------------------------------------------- brain ---------------
sm = pd.read_csv(WORK / "scan_means_ct.csv")
demo = (sm.sort_values("visit").groupby("subject")
          .agg(sex=("sex", "first"), site=("site", "first"), family=("family_id", "first"),
               age_first=("age_first", "first"), n_visits=("n_visits", "first"), age_last=("age", "max")))
demo["age_span"] = demo.age_last - demo.age_first
T = demo.copy()
b = pd.read_csv(WORK / "blups_ct.csv").set_index("subject")
T["dCT"] = zscore(b.re_slope); T["CT0"] = zscore(b.re_intercept); T["dCT_rel"] = b.slope_rel
bt = pd.read_csv(WORK / "blups_t1t2.csv").set_index("subject")
T["dT1T2"] = zscore(bt.re_slope).reindex(T.index)
T["male"] = (T.sex == "M").astype(float)
qc = pd.read_csv(R70 / "y/mr_y_qc__post__aut.tsv", sep="\t", low_memory=False,
                 usecols=["participant_id", "session_id", "mr_y_qc__post__aut__smri__topodfct_count"])
qc["subject"] = to_img(qc.participant_id)
qc["visit"] = qc.session_id.map({"ses-00A": "v0", "ses-02A": "v2", "ses-04A": "v4", "ses-06A": "v6"})
qc["ltopo"] = np.log1p(pd.to_numeric(qc.mr_y_qc__post__aut__smri__topodfct_count, errors="coerce"))
sq = sm[["subject", "visit", "age"]].merge(qc[["subject", "visit", "ltopo"]], on=["subject", "visit"], how="left")
T["qc_ltopo_mean"] = sq.groupby("subject").ltopo.mean()
T["qc_ltopo_slope"] = sq.dropna().groupby("subject")[["age", "ltopo"]].apply(
    lambda g: np.polyfit(g.age, g.ltopo, 1)[0] if len(g) >= 2 else np.nan)
T["site"] = T.site.astype(str)

# ---------------------------------------------------------------- puberty -------------
def puberty(fname: str, pre: str) -> pd.DataFrame:
    p = R70 / fname
    d = pd.read_csv(p, sep="\t", low_memory=False,
                    usecols=["participant_id", "session_id", f"{pre}_age", f"{pre}__f_mean", f"{pre}__m_mean"])
    d["subject"] = to_img(d.participant_id)
    for c in (f"{pre}_age", f"{pre}__f_mean", f"{pre}__m_mean"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d.subject.isin(T.index)].copy()
    d["sex"] = T.sex.reindex(d.subject).to_numpy()
    d["pds"] = np.where(d.sex == "F", d[f"{pre}__f_mean"], d[f"{pre}__m_mean"]); d["age"] = d[f"{pre}_age"]
    d = d.dropna(subset=["pds", "age", "sex"])

    def res(g):
        if len(g) < 50:
            return pd.Series(np.nan, index=g.index)
        Xv = np.vander(g.age - g.age.mean(), 4); bb, *_ = np.linalg.lstsq(Xv, g.pds, rcond=None)
        r = g.pds - Xv @ bb; return r / r.std()
    d["timing"] = d.groupby(["sex", "session_id"], group_keys=False)[["age", "pds"]].apply(res)
    out = pd.DataFrame(index=T.index)
    out["pub_timing_base"] = d[d.session_id == "ses-00A"].set_index("subject").timing
    out["pub_timing"] = d.groupby("subject").timing.mean()
    out["pub_tempo"] = d.groupby("subject")[["age", "pds"]].apply(
        lambda g: np.polyfit(g.age, g.pds, 1)[0] if len(g) >= 3 and g.age.std() > 0.5 else np.nan)
    return out


T = T.join(puberty("ph_p_pds.tsv", "ph_p_pds"))

# ---------------------------------------------------------------- CBCL ----------------
SRC = REPO / "ahba_pls/code/23_cbcl_explore.py"
src = SRC.read_text(); cut = src.index("\ndef fit(")
ns = {"__file__": str(SRC), "__name__": "cbcl_prefix"}
argv, sys.argv = sys.argv, [str(SRC), "70"]
try:
    exec(compile(src[:cut], str(SRC), "exec"), ns)
finally:
    sys.argv = argv
D = ns["D"]
CBCL = ["depress", "internal", "external", "pfactor"]
for o in CBCL:
    T[f"{o}_BASE"] = D[f"{o}_BASE"].reindex(T.index)
    T[f"{o}_LATE"] = D[f"{o}_LATE"].reindex(T.index)
T["age_LATE"] = D["age_LATE"].reindex(T.index)
for dx in ("mdd_youth_DX", "mdd_parent_DX", "psychosis_parent_DX"):
    T[dx] = D[dx].reindex(T.index)
T["age_LAST"] = D["age_LAST"].reindex(T.index)

# ---------------------------------------------------------------- NIH Toolbox ---------
nt = pd.read_csv(R70 / "nc_y_nihtb.tsv", sep="\t", low_memory=False)
nt["subject"] = to_img(nt.participant_id)
agec = [c for c in nt.columns if c.startswith("nc_y_nihtb__") and c.endswith("_age")]
nt["age_nt"] = nt[agec].apply(pd.to_numeric, errors="coerce").mean(axis=1)
COG = {"fluid": "nc_y_nihtb__comp__fluid__uncor_score", "cryst": "nc_y_nihtb__comp__cryst__uncor_score",
       "total": "nc_y_nihtb__comp__tot__uncor_score"}
for k, c in COG.items():
    nt[k] = pd.to_numeric(nt[c], errors="coerce")
for ses, tag in (("ses-00A", "00"), ("ses-06A", "06")):
    w = nt[nt.session_id == ses].drop_duplicates("subject").set_index("subject")
    for k in COG:
        T[f"{k}_{tag}"] = w[k].reindex(T.index)
    T[f"age_nt_{tag}"] = w.age_nt.reindex(T.index)

# ---------------------------------------------------------------- residuals -----------
BR = ["male", "age_first", "age_span", "n_visits"]
for t in ("dCT", "CT0", "dT1T2"):
    T[f"{t}_r"] = residualise(T, t, BR, ["site"])
T["dCT_rq"] = residualise(T, "dCT", BR + ["qc_ltopo_mean", "qc_ltopo_slope"], ["site"])
T["pub_timing_r"] = residualise(T, "pub_timing", ["male"], ["site"])
T["pub_tempo_r"] = residualise(T, "pub_tempo", ["male", "pub_timing_base"], ["site"])
for o in CBCL:
    T[f"{o}_chg_r"] = residualise(T, f"{o}_LATE", [f"{o}_BASE", "age_LATE", "male"], ["site"])
for k in COG:
    T[f"{k}_gain_r"] = residualise(T, f"{k}_06", [f"{k}_00", "age_nt_00", "age_nt_06", "male"], ["site"])

# ---------------------------------------------------------------- PRS covariates ------
ph = pd.read_csv(PHENO_DIR / "phenotypes_gcta.txt", sep=r"\s+", dtype={"FID": str})
cq = pd.read_csv(PHENO_DIR / "covar_quant.txt", sep=r"\s+", dtype={"FID": str})
g = ph.merge(cq, on=["FID", "IID"]); g["subject"] = to_img(g.IID)
g = g.set_index("subject")
T["gen_family"] = g.FID.reindex(T.index)
for c in ["global_slope_1lmm"] + [f"PC{i}" for i in range(1, 11)]:
    T[c] = g[c].reindex(T.index)
chk = T[["dCT", "global_slope_1lmm"]].dropna()
r_chk = chk.dCT.corr(chk.global_slope_1lmm)
assert r_chk > 0.999, f"dCT does not reproduce the genetics trait: r = {r_chk:.5f}"
T.to_parquet(WORK / "traits.parquet")

# ---------------------------------------------------------------- pairs ---------------
gr = pd.read_csv(R70 / "y/gn_y_genrel.tsv", sep="\t")
rows = []
for k in ("01", "02", "03", "04"):
    s = gr[["participant_id", f"gn_y_genrel_id__paired__{k}", f"gn_y_genrel_pihat__{k}",
            f"gn_y_genrel_zyg__{k}", "gn_y_genrel_id__fam", "gn_y_genrel_id__birth"]].dropna(
        subset=[f"gn_y_genrel_id__paired__{k}"])
    s.columns = ["a", "b", "pihat", "zyg", "fam", "birth"]; rows.append(s)
P = pd.concat(rows)
P = P[P.a < P.b].drop_duplicates(["a", "b"])
P["id1"], P["id2"] = to_img(P.a), to_img(P.b)
P["type"] = P.zyg.map({1: "MZ", 2: "DZ", 3: "SIB"})
P = P[P.type.notna()]
counts = [dict(stage="release (all genotyped pairs)", **P.type.value_counts().to_dict())]
P = P[P.id1.isin(T.index) & P.id2.isin(T.index)]
counts.append(dict(stage="both imaged (>= 2 scans)", **P.type.value_counts().to_dict()))
P = P.sample(frac=1, random_state=3)          # random pair choice within birth / family
tw = P[P.type != "SIB"].drop_duplicates("birth")
sib = P[P.type == "SIB"].drop_duplicates("fam")
P = pd.concat([tw, sib])
counts.append(dict(stage="one pair per birth event / family", **P.type.value_counts().to_dict()))
flip = RNG.random(len(P)) < 0.5
P.loc[flip, ["id1", "id2"]] = P.loc[flip, ["id2", "id1"]].to_numpy()
P[["id1", "id2", "type", "pihat", "fam", "birth"]].to_csv(WORK / "pairs.csv", index=False)
P1 = T.reindex(P.id1); P2 = T.reindex(P.id2)
P["opp_sex"] = (P1.sex.to_numpy() != P2.sex.to_numpy())
counts.append(dict(stage="of which opposite-sex", **P[P.opp_sex].type.value_counts().to_dict()))
pd.DataFrame(counts).fillna(0).astype({c: int for c in ("MZ", "DZ", "SIB")}).to_csv(
    RES / "pair_counts.tsv", sep="\t", index=False)

# trait availability per pair type (both members non-missing)
RT = ["dCT_r", "dCT_rq", "CT0_r", "dT1T2_r", "pub_timing_r", "pub_tempo_r"] + [f"{o}_chg_r" for o in CBCL] + \
     [f"{k}_gain_r" for k in COG]
av = []
for t in RT:
    both = P1[t].notna().to_numpy() & P2[t].notna().to_numpy()
    av.append(dict(trait=t, n_children=int(T[t].notna().sum()),
                   **{f"{ty}_pairs": int((both & (P.type == ty).to_numpy()).sum()) for ty in ("MZ", "DZ", "SIB")}))
pd.DataFrame(av).to_csv(RES / "trait_availability.tsv", sep="\t", index=False)
print(pd.DataFrame(counts).to_string(index=False))
print(pd.DataFrame(av).to_string(index=False))
print(f"dCT vs genetics trait r = {r_chk:.5f} on n = {len(chk)}")
