"""
Compare the faster fMRIPrep-T1w -> DAIRC registrations with the pilot reference, and check
abcd.surf_metrics' LTA-based sampling of fMRIPrep volumes against FreeSurfer mri_vol2surf.

1. For each variant in surf_metrics_pilot/t2reg_speed.tsv: time, bbregister cost, and the
   displacement (mm) of the lh white-surface vertices between the variant's and the reference
   transform (both register.dat, the same tkregister convention).
2. For sessions with surf_metrics_pilot/<s>.lh.fmriprep_T2w.mid.mgh (mri_vol2surf --reg bbr.lta
   on the white-pial midpoint surface): per-vertex r against our t2__gm50 sample.
"""

from __future__ import annotations

import sys

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.freesurfer import read_geometry

from abcd import paths
from abcd.hcp_stats import DEFAULT_FS_ROOT
from abcd.surf_metrics import (DEFAULT_FMRIPREP_ROOT, fmriprep_image, load_registered,
                               sample_hemi)

WORK = paths.release_dir("7.0") / "processed" / "hcp" / "surf_metrics_pilot"


def read_dat(path) -> np.ndarray:
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    return np.array([[float(x) for x in l.split()] for l in lines[4:8]])


def main() -> int:
    t = pd.read_csv(WORK / "t2reg_speed.tsv", sep="\t")
    disp = []
    for r in t.itertuples():
        s = f"{r.participant_id}_{r.session_id}"
        ref, var = WORK / "reg" / s / "bbr.dat", WORK / f"reg_{r.variant}" / s / "bbr.dat"
        if r.status != "ok" or not ref.exists() or not var.exists():
            disp.append((np.nan, np.nan))
            continue
        w, _ = read_geometry(str(DEFAULT_FS_ROOT / r.participant_id / r.session_id / "surf" / "lh.white"))
        x = np.c_[w, np.ones(len(w))]
        d = np.linalg.norm((x @ (read_dat(var) - read_dat(ref)).T)[:, :3], axis=1)
        disp.append((d.mean(), d.max()))
    t[["disp_mean_mm", "disp_max_mm"]] = disp
    t["sec_total"] = t.sec_coreg.fillna(0) + t.sec_bbregister
    print(t.to_string(index=False))
    print("\nby variant (median; max for displacement):")
    print(t.groupby("variant").agg(n_ok=("status", lambda s: (s == "ok").sum()),
                                   sec_coreg=("sec_coreg", "median"),
                                   sec_bbregister=("sec_bbregister", "median"),
                                   sec_total=("sec_total", "median"),
                                   mincost_median=("bbr_mincost", "median"),
                                   mincost_max=("bbr_mincost", "max"),
                                   disp_mean_mm=("disp_mean_mm", "median"),
                                   disp_max_mm=("disp_max_mm", "max")).round(3).to_string())

    print("\nLTA sampling vs mri_vol2surf --reg bbr.lta (fMRIPrep T2w, midpoint surface, lh):")
    for f in sorted(WORK.glob("*.lh.fmriprep_T2w.mid.mgh")):
        s = f.name.split(".")[0]
        sub, ses = s.rsplit("_", 1)
        sdir = DEFAULT_FS_ROOT / sub / ses
        vol = load_registered(fmriprep_image(DEFAULT_FMRIPREP_ROOT, sub, ses, "T2w"),
                              WORK / "reg" / s / "bbr.lta", sdir / "mri" / "orig.mgz")
        vals, _ = sample_hemi(sdir, "lh", ("t2__gm50",), {"fmriprep:T2w": vol})
        fsv = np.asarray(nib.load(str(f)).dataobj).reshape(-1)
        ok = np.isfinite(vals["t2__gm50"]) & (fsv > 0)
        print(f"{s}: n={ok.sum()}, r={np.corrcoef(vals['t2__gm50'][ok], fsv[ok])[0, 1]:.5f}, "
              f"median|d|/median={np.median(np.abs(vals['t2__gm50'][ok] - fsv[ok])) / np.median(fsv[ok]):.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
