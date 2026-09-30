"""Per-scan whole-cortex means for the two single-LMM slope phenotypes.

    python directions/d3_twin_family/01_scan_means.py

Thickness: mean over the 358 HCP-MMP parcels of the settled fit
(out/thickness_hcp_70_aa6e91efba82), exactly as genetic_analysis/fig1_prep_1lmm.R
builds its input, so the slope is the Figure-1 / genetics primary trait.
T1w/T2w: mean over the 68 DK parcels of the matched ratio fit
(out/t1t2_ratio_dsk_70_9a62dde44370; config t1t2_70_noglobal_mv2_genetic).
Writes work/scan_means_{ct,t1t2}.csv (individual-level, gitignored).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402

from common import RUN_CT, RUN_T1T2, WORK  # noqa: E402

COLS = ["subject", "visit", "value", "age", "age_c", "sex", "site", "family_id", "age_first", "n_visits"]
for tag, run in (("ct", RUN_CT), ("t1t2", RUN_T1T2)):
    mt = pd.read_parquet(run / "model_table.parquet", columns=COLS)
    agg = (mt.groupby(["subject", "visit"], observed=True)
             .agg(mean_val=("value", "mean"), n_parcels=("value", "size"), age=("age", "first"),
                  age_c=("age_c", "first"), sex=("sex", "first"), site=("site", "first"),
                  family_id=("family_id", "first"), age_first=("age_first", "first"),
                  n_visits=("n_visits", "first"))
             .reset_index())
    assert agg.n_parcels.nunique() == 1, agg.n_parcels.value_counts()
    agg.to_csv(WORK / f"scan_means_{tag}.csv", index=False)
    print(tag, agg.subject.nunique(), "children", len(agg), "scans", int(agg.n_parcels.iloc[0]), "parcels")
