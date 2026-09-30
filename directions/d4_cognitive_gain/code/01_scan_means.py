"""Per-scan HCP-MMP cortical mean thickness with covariates (individual-level -> work/).

Identical construction to genetic_analysis/fig1_prep_1lmm.R's input (Figure 1 phenotypes):
area-weighting is already in the run's parcel values; the per-scan mean is the mean over
parcels, one row per child x visit.
Run from repo root:  python directions/d4_cognitive_gain/code/01_scan_means.py
"""
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RUN = REPO / "out/thickness_hcp_70_aa6e91efba82"
OUT = REPO / "directions/d4_cognitive_gain/work/hcp70_scan_means_cov.csv"

mt = pd.read_parquet(RUN / "model_table.parquet",
                     columns=["subject", "visit", "value", "age", "age_c", "sex", "site", "n_visits"])
d = (mt.groupby(["subject", "visit"], observed=True)
       .agg(mean_ct=("value", "mean"), age=("age", "first"), age_c=("age_c", "first"),
            sex=("sex", "first"), site=("site", "first"), n_visits=("n_visits", "first"))
       .reset_index())
assert d.subject.nunique() == 8716, d.subject.nunique()
d.to_csv(OUT, index=False)
print(len(d), d.subject.nunique(), OUT)
