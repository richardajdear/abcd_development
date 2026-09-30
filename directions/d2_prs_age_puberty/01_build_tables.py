"""D2 step 1 -- build the individual-level tables for the one-stage PRS x age / PRS x puberty models.

Everything this writes is SUBJECT-LEVEL and goes to a gitignored work directory
(genetic_analysis/work/results_70tab_hcp/d2_prs_age_puberty/). Only 02_fit_models.R's
aggregate tables are committed.

Inputs (all already on disk; nothing fetched):
  out/thickness_hcp_70_aa6e91efba82/model_table.parquet   settled spec on HCP-MMP; the Figure 1 trait
  genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno/{covar_quant,phenotypes_gcta}.txt
                                                           PCs, and the Figure 1 single-LMM phenotypes
  genetic_analysis/work/scores_scz2025/SBayesRC|PRSCS/...  SCZ 2025 scores (EUR weights raw; META weights
                                                           z-scored within ancestry cluster = _zanc)
  legacy/hpc_v2/work/results_v2/prs_final/...              MDD scores, same construction (genotype-only
                                                           products, independent of imaging vintage)
  legacy/hpc/work/results/ancestry/eur_anchor.keep         the EUR arm (4,308 children in the export)
  abcd-data-release-7.0/ph_y_pds.tsv, ph_p_pds.tsv         Pubertal Development Scale, yearly

Outputs (work dir):
  scans.csv      one row per QC-passing scan: mean_ct (unweighted mean over the 358 HCP parcels, as
                 fig1_prep_1lmm.R), age, age_c (centre 12.797), sex, site, family_id, PC1-10, the PRS
                 columns, eur flag, PDS at the same session (youth report; parent as sensitivity)
  pds_long.csv   all yearly PDS waves for the children in scans.csv (for the PDS-tempo LMM)
  children.csv   one row per child: Figure 1 phenotypes, PRS, covariates

Run from the repo root:  python directions/d2_prs_age_puberty/01_build_tables.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / "out/thickness_hcp_70_aa6e91efba82"
EXP = REPO / "genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno"
WORK = REPO / "genetic_analysis/work/results_70tab_hcp/d2_prs_age_puberty"
REL = REPO / "abcd-data-release-7.0"
AGE_CENTRE = 12.797   # fig1_prep_1lmm.R; equals the sample-mean age the run was centred on

SCORES = {
    # column name            : (profile path, arm)   arm 'pooled' = META weights, z within cluster (all 8,596)
    #                                                arm 'eur'    = EUR weights, raw score (EUR children only)
    "scz_sbrc_pooled": ("genetic_analysis/work/scores_scz2025/SBayesRC/SCZ25_META/_zanc/score_SCZ25META_sbrc.profile", "pooled"),
    "scz_sbrc_eur":    ("genetic_analysis/work/scores_scz2025/SBayesRC/SCZ25_EUR/score_SCZ25EUR_sbrc.profile", "eur"),
    "scz_prscs_pooled": ("genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_META/_zanc/score_SCZ25META_prscs.profile", "pooled"),
    "scz_prscs_eur":   ("genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_EUR/score_SCZ25EUR_prscs.profile", "eur"),
    "mdd_sbrc_pooled": ("legacy/hpc_v2/work/results_v2/prs_final/SBayesRC/MDD_pooled/_zanc/score_MDDpooled_sbrc.profile", "pooled"),
    "mdd_sbrc_eur":    ("legacy/hpc_v2/work/results_v2/prs_final/SBayesRC/MDD_eur/score_MDDeur_sbrc.profile", "eur"),
    "mdd_prscs_pooled": ("legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_pooled/_zanc/score_MDDpooled_prscs.profile", "pooled"),
    "mdd_prscs_eur":   ("legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_eur/score_MDDeur_prscs.profile", "eur"),
}
VISIT2SES = {"v0": "ses-00A", "v2": "ses-02A", "v4": "ses-04A", "v6": "ses-06A"}


def scan_means() -> pd.DataFrame:
    cols = ["subject", "visit", "value", "age", "age_c", "sex", "site", "family_id", "n_visits"]
    mt = pd.read_parquet(RUN / "model_table.parquet", columns=cols)
    g = (mt.groupby(["subject", "visit"], observed=True)
           .agg(mean_ct=("value", "mean"), age=("age", "first"), age_c=("age_c", "first"),
                sex=("sex", "first"), site=("site", "first"), family_id=("family_id", "first"),
                n_visits=("n_visits", "first"))
           .reset_index())
    assert abs((g.age - g.age_c).mean() - AGE_CENTRE) < 0.01, (g.age - g.age_c).mean()
    # imaging tables carry sub-NDARINVxxxx; the genetics export and 7.0 tabulated data carry sub-xxxx
    g["subject"] = g.subject.str.replace("sub-NDARINV", "sub-", regex=False)
    g["session_id"] = g.visit.map(VISIT2SES)
    return g


def prs_table(ids: set) -> pd.DataFrame:
    out = None
    for col, (path, arm) in SCORES.items():
        p = pd.read_csv(REPO / path, sep=r"\s+")[["IID", "SCORESUM"]].rename(columns={"IID": "subject", "SCORESUM": col})
        p = p[p.subject.isin(ids)]
        out = p if out is None else out.merge(p, on="subject", how="outer")
    eur = pd.read_csv(REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep", sep=r"\s+", header=None)[1]
    out["eur"] = out.subject.isin(set(eur)).astype(int)
    # raw EUR-weight scores: standardise within the EUR arm (the arm they are used in);
    # _zanc scores are already z within ancestry cluster over all 11,670 genotyped children
    for col, (_, arm) in SCORES.items():
        if arm == "eur":
            m = out.eur == 1
            out.loc[~m, col] = np.nan
            out.loc[m, col] = (out.loc[m, col] - out.loc[m, col].mean()) / out.loc[m, col].std()
    return out


def pds_long(ids: set) -> pd.DataFrame:
    frames = []
    for who in ["y", "p"]:
        d = pd.read_csv(REL / f"ph_{who}_pds.tsv", sep="\t", low_memory=False)
        d = d[d.participant_id.isin(ids)]
        f, m = f"ph_{who}_pds__f_mean", f"ph_{who}_pds__m_mean"
        fc, mc = f"ph_{who}_pds__f_categ", f"ph_{who}_pds__m_categ"
        x = pd.DataFrame({"subject": d.participant_id, "session_id": d.session_id,
                          "pds_age": d[f"ph_{who}_pds_age"],
                          f"pds_mean_{who}": d[f].fillna(d[m]),
                          f"pds_categ_{who}": d[fc].fillna(d[mc])})
        frames.append(x)
    y, p = frames
    long = y.merge(p.drop(columns="pds_age"), on=["subject", "session_id"], how="outer")
    long["pds_age"] = long.pds_age.fillna(long.groupby("subject").pds_age.transform("mean"))
    long = long.dropna(subset=["pds_mean_y", "pds_mean_p"], how="all")
    long["age_c"] = long.pds_age - AGE_CENTRE
    return long.sort_values(["subject", "session_id"])


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    scans = scan_means()
    ids = set(scans.subject)
    q = pd.read_csv(EXP / "covar_quant.txt", sep=r"\s+").rename(columns={"IID": "subject"})
    pcs = q[["subject"] + [f"PC{i}" for i in range(1, 11)]]
    prs = prs_table(ids)
    scans = scans.merge(pcs, on="subject", how="inner").merge(prs, on="subject", how="inner")
    pl = pds_long(ids)
    scans = scans.merge(pl[["subject", "session_id", "pds_mean_y", "pds_categ_y", "pds_mean_p", "pds_categ_p"]],
                        on=["subject", "session_id"], how="left")
    ph = pd.read_csv(EXP / "phenotypes_gcta.txt", sep=r"\s+").rename(columns={"IID": "subject"}).drop(columns="FID")
    children = (scans.groupby("subject").agg(sex=("sex", "first"), site=("site", "first"), family_id=("family_id", "first"),
                                             n_visits=("n_visits", "first"), age_first=("age", "min"))
                .reset_index().merge(ph, on="subject").merge(pcs, on="subject").merge(prs, on="subject"))
    scans.to_csv(WORK / "scans.csv", index=False)
    pl[pl.subject.isin(set(scans.subject))].to_csv(WORK / "pds_long.csv", index=False)
    children.to_csv(WORK / "children.csv", index=False)
    print(f"scans {len(scans)} rows, {scans.subject.nunique()} children ({int(children.eur.sum())} EUR); "
          f"PDS matched at scan: youth {scans.pds_mean_y.notna().mean():.2f}, parent {scans.pds_mean_p.notna().mean():.2f}; "
          f"pds_long {len(pl)} waves, median {pl.groupby('subject').size().median():.0f} per child")


if __name__ == "__main__":
    main()
