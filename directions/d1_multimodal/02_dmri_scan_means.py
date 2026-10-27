"""D1 step 2: per-scan whole-cortex means of every DTI and RSI release table (DK), Figure-1 sample rules.

    PYTHONPATH=src python directions/d1_multimodal/02_dmri_scan_means.py      (env abcd-spatial)

Reads every  abcd-data-release-7.0/y/mr_y_{dti,rsi}__*__{gm,wm,gwc}__dsk.tsv  and takes the release's
own whole-cortex column  <stem>_mean  (the release convention, as for thickness). Metric = the table stem
without the mr_y_ prefix and the __dsk suffix, e.g. dti__is__fa__gm, rsi__rnd__wm.

Sample rules (mirroring out/thickness_hcp_70_aa6e91efba82/qc_ledger.csv, with the dMRI release flag in
place of the T1 one; surface-defect QC is inherited because the cortical sampling uses the same surfaces):
  visits ses-00A/02A/04A/06A; non-missing whole-cortex value; mr_y_qc__incl__dmri_indicator == 1;
  scanner manufacturer != Philips (adapter decoding); children with >= 2 remaining scans.
Covariates from the 7.0 adapter (longitudinal: age, site, family_id; demographics: sex). age_c is
centred at 12.797 y, the thickness fit's centre, so intercepts are comparable across modalities.

Writes work/scan_means_dmri.parquet (long: subject, visit, metric, mean_val, age, age_c, sex, site,
family_id, n_visits)  -- individual-level, gitignored -- and tables/dmri_sample.tsv (counts only).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
from abcd import paths  # noqa: E402
from abcd.io import get_adapter  # noqa: E402

# The repo-vendored release is first in the search order. Drop the ~/Git/ABCD fallback, whose
# existence check raises PermissionError inside the analysis sandbox.
paths._ABCD_FALLBACKS = ()

WORK = HERE / "work"; WORK.mkdir(exist_ok=True)
TAB = HERE / "tables"; TAB.mkdir(exist_ok=True)
Y = REPO / "abcd-data-release-7.0" / "y"
CEN = 12.797
VIS = {"ses-00A": "v0", "ses-02A": "v2", "ses-04A": "v4", "ses-06A": "v6"}

A = get_adapter("7.0")
lon = A.longitudinal()                       # subject, visit, age, site, family_id
dem = A.demographics()
sc = A.scanner()
print("adapter tables:", lon.shape, dem.shape, sc.shape, "| demo cols", list(dem.columns)[:6], "| scanner cols", list(sc.columns))
sex = dem[["subject", "sex"]].drop_duplicates("subject")
qc = pd.read_csv(Y / "mr_y_qc__incl.tsv", sep="\t", usecols=["participant_id", "session_id", "mr_y_qc__incl__dmri_indicator"])
qc["subject"] = "sub-NDARINV" + qc.participant_id.str.replace("sub-", "", regex=False)
qc["visit"] = qc.session_id.map(VIS)
qc = qc.dropna(subset=["visit"])[["subject", "visit", "mr_y_qc__incl__dmri_indicator"]]

files = sorted(f for f in Y.glob("mr_y_*__dsk.tsv") if re.match(r"mr_y_(dti|rsi)__.*__(gm|wm|gwc)__dsk\.tsv$", f.name))
print(len(files), "tables")
parts, counts = [], []
for f in files:
    stem = f.name[:-4]
    metric = stem.replace("mr_y_", "").replace("__dsk", "")
    d = pd.read_csv(f, sep="\t", usecols=["participant_id", "session_id", f"{stem}_mean"]).rename(columns={f"{stem}_mean": "mean_val"})
    d["subject"] = "sub-NDARINV" + d.participant_id.str.replace("sub-", "", regex=False)
    d["visit"] = d.session_id.map(VIS)
    c = dict(metric=metric, rows=len(d))
    d = d.dropna(subset=["visit", "mean_val"]); c["imaging_visits_nonmissing"] = len(d)
    d = d.merge(qc, on=["subject", "visit"], how="left")
    d = d[d.mr_y_qc__incl__dmri_indicator == 1]; c["release_dmri_include"] = len(d)
    d = d.merge(sc, on=["subject", "visit"], how="left")
    d = d[~d.manufacturer.astype(str).str.contains("philips", case=False)]; c["not_philips"] = len(d)
    d = d.merge(lon, on=["subject", "visit"], how="inner").merge(sex, on="subject", how="inner")
    d = d.dropna(subset=["age", "site", "sex"])
    d["n_visits"] = d.groupby("subject").visit.transform("size")
    d = d[d.n_visits >= 2]; c["scans_in_model"] = len(d); c["children_in_model"] = d.subject.nunique()
    d["age_c"] = d.age - CEN
    d["metric"] = metric
    parts.append(d[["subject", "visit", "metric", "mean_val", "age", "age_c", "sex", "site", "family_id", "n_visits"]])
    counts.append(c)
    print(f"{metric:22s} " + " ".join(f"{k}={v}" for k, v in c.items() if k != "metric"), flush=True)

out = pd.concat(parts, ignore_index=True)
out.to_csv(WORK / "scan_means_dmri.csv.gz", index=False)
pd.DataFrame(counts).to_csv(TAB / "dmri_sample.tsv", sep="\t", index=False)
print("wrote", len(out), "rows;", out.metric.nunique(), "metrics")
