"""D1 step 7a: per-scan whole-cortex values for the HCP-MMP structural measures (CSD3 tables), plus the scanner
unit and software version of every scan, for the scanner-term LMMs (10_scanner_terms.R).

    python directions/d1_multimodal/09_smri_scan_means.py      (env abcd-spatial)

Structural measures (abcd-data-release-7.0/processed/hcp/, same DAIRC surfaces and HCP annotation as thickness):
  area  mr_y_smri__area__hcp_sum      total cortical surface area (mm2)
  vol   mr_y_smri__vol__hcp_sum       total cortical grey-matter volume (mm3)
  gwc   mr_y_smri__gwc__hcp_mean      grey/white contrast (FreeSurfer pctsurfcon, %)
  t1gm  mr_y_smri__t1__gm50__hcp_mean T1w (nu.mgz) at GM mid-depth
  t1wm  mr_y_smri__t1__wm1__hcp_mean  T1w (nu.mgz) 1 mm below the white surface
Sample: exactly the Figure-1 thickness scans (directions/d3_twin_family/work/scan_means_ct.csv: release T1 include,
surface defects <= 46, no Philips, >= 2 visits), so every structural measure shares sample and surfaces with the
thinning slope. Thickness itself is that file's Figure-1 value (unweighted mean of the 358 parcels).
Also carried, on their own samples: T1w/T2w (DK, release; d3 scan_means_t1t2.csv) and 8 representative diffusion
measures (work/scan_means_dmri.csv.gz).
Scanner terms from ab_g_dyn: serial (scanner unit) and serial x software version ("epoch").
Writes work/scan_means_smri.csv.gz (long; gitignored) and tables/scanner_epochs.tsv (counts only).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
W, TAB = HERE / "work", HERE / "tables"
D3W = HERE.parent / "d3_twin_family" / "work"
R70 = REPO / "abcd-data-release-7.0"
VIS = {"ses-00A": "v0", "ses-02A": "v2", "ses-04A": "v4", "ses-06A": "v6"}
KEEP = ["subject", "visit", "metric", "mean_val", "age", "age_c", "sex", "site", "family_id"]
DMRI = ["dti__is__md__gwc", "dti__is__fa__wm", "dti__is__td__wm", "dti__is__md__gm",
        "rsi__rnd__wm", "rsi__fni__gm", "rsi__hnt__wm", "rsi__rni__gm"]


def ids(d: pd.DataFrame) -> pd.DataFrame:
    d["subject"] = "sub-NDARINV" + d.participant_id.str.replace("sub-", "", regex=False)
    d["visit"] = d.session_id.map(VIS)
    return d.dropna(subset=["visit"])


ct = pd.read_csv(D3W / "scan_means_ct.csv")
base = ct[["subject", "visit", "age", "age_c", "sex", "site", "family_id"]]
parts = [ct.assign(metric="smri__thk_hcp")[KEEP]]
STEMS = {"smri__area_hcp": ("mr_y_smri__area__hcp", "_sum"), "smri__vol_hcp": ("mr_y_smri__vol__hcp", "_sum"),
         "smri__gwc_hcp": ("mr_y_smri__gwc__hcp", "_mean"), "smri__t1gm50_hcp": ("mr_y_smri__t1__gm50__hcp", "_mean"),
         "smri__t1wm1_hcp": ("mr_y_smri__t1__wm1__hcp", "_mean")}
for metric, (stem, agg) in STEMS.items():
    col = stem + agg
    d = ids(pd.read_csv(R70 / "processed/hcp" / f"{stem}.tsv", sep="\t", usecols=["participant_id", "session_id", col]))
    d = base.merge(d.rename(columns={col: "mean_val"})[["subject", "visit", "mean_val"]], on=["subject", "visit"], how="inner")
    d = d.dropna(subset=["mean_val"]); d["metric"] = metric
    print(f"{metric:18s} scans {len(d):6d} of {len(base)} Figure-1 scans; children {d.subject.nunique()}")
    parts.append(d[KEEP])
t12 = pd.read_csv(D3W / "scan_means_t1t2.csv").assign(metric="smri__t1t2_dsk")
parts.append(t12[KEEP])
dm = pd.read_csv(W / "scan_means_dmri.csv.gz")
parts.append(dm[dm.metric.isin(DMRI)][KEEP])
S = pd.concat(parts, ignore_index=True)

dy = ids(pd.read_csv(R70 / "g/dyn/ab_g_dyn.tsv", sep="\t", usecols=[
    "participant_id", "session_id", "ab_g_dyn__design_mr__serial", "ab_g_dyn__design_mr__software"]))
dy = dy.rename(columns={"ab_g_dyn__design_mr__serial": "serial", "ab_g_dyn__design_mr__software": "software"})
S = S.merge(dy[["subject", "visit", "serial", "software"]], on=["subject", "visit"], how="left")
S["epoch"] = S.serial.astype(str) + "|" + S.software.astype(str)
S["n_visits"] = S.groupby(["metric", "subject"]).visit.transform("size")
S = S[S.n_visits >= 2]
print("missing serial", S.serial.isna().mean().round(4), "| epochs", S.epoch.nunique(), "| serials", S.serial.nunique())
ep = (S[S.metric == "smri__thk_hcp"].groupby("epoch").agg(scans=("visit", "size"), children=("subject", "nunique"),
      age_mean=("age", "mean"), age_sd=("age", "std")).reset_index())
ep["serial"] = ep.epoch.str.split("|").str[0]
ep.to_csv(TAB / "scanner_epochs.tsv", sep="\t", index=False)
print(f"thickness sample: {len(ep)} serial x software epochs; median scans/epoch {ep.scans.median():.0f}; "
      f"median within-epoch age SD {ep.age_sd.median():.2f} y (age SD overall {S[S.metric=='smri__thk_hcp'].age.std():.2f} y)")
S.to_csv(W / "scan_means_smri.csv.gz", index=False)
print(S.groupby("metric").agg(scans=("visit", "size"), children=("subject", "nunique")).to_string())
