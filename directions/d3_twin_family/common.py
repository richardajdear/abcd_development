"""Shared paths and helpers for directions/d3_twin_family/.

Individual-level files are written ONLY to WORK (gitignored). Committed outputs
in results/ are coefficient and count tables.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"            # individual-level, gitignored
RES = HERE / "results"          # aggregate tables, committed
FIG = HERE / "figures"
R70 = REPO / "abcd-data-release-7.0"
RUN_CT = REPO / "out" / "thickness_hcp_70_aa6e91efba82"     # settled HCP-MMP thickness fit
RUN_T1T2 = REPO / "out" / "t1t2_ratio_dsk_70_9a62dde44370"  # matched T1w/T2w fit (DK)
PHENO_DIR = REPO / "genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno"
SCORES = {
    # name: path -- pooled arm = ancestry-standardised (_zanc) scores, n = 8,596
    "SCZ25_prscs": REPO / "genetic_analysis/work/scores_scz2025/PRSCS/SCZ25_META/_zanc/score_SCZ25META_prscs.profile",
    "SCZ25_sbrc": REPO / "genetic_analysis/work/scores_scz2025/SBayesRC/SCZ25_META/_zanc/score_SCZ25META_sbrc.profile",
    "MDD_prscs": REPO / "legacy/hpc_v2/work/results_v2/prs_final/PRSCS/MDD_pooled/_zanc/score_MDDpooled_prscs.profile",
    "MDD_sbrc": REPO / "legacy/hpc_v2/work/results_v2/prs_final/SBayesRC/MDD_pooled/_zanc/score_MDDpooled_sbrc.profile",
}
for d in (WORK, RES, FIG):
    d.mkdir(exist_ok=True)


def to_img(pid: pd.Series) -> pd.Series:
    """Release participant_id 'sub-XXXXXXXX' -> imaging subject 'sub-NDARINVXXXXXXXX'."""
    return "sub-NDARINV" + pid.str.replace("sub-", "", regex=False)


def zscore(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std()


def residualise(df: pd.DataFrame, y: str, covs: list[str], cat: list[str]) -> pd.Series:
    """OLS residual of y on covariates (categoricals dummy-coded), z-scored; NaN where missing."""
    d = df[[y, *covs, *cat]].dropna()
    X = pd.get_dummies(d[covs + cat], columns=cat, drop_first=True, dtype=float)
    X.insert(0, "const", 1.0)
    b, *_ = np.linalg.lstsq(X.to_numpy(), d[y].to_numpy(), rcond=None)
    r = pd.Series(d[y].to_numpy() - X.to_numpy() @ b, index=d.index)
    return zscore(r).reindex(df.index)
