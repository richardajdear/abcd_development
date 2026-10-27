"""
Is fMRIPrep's desc-preproc_T2w bias-corrected like its desc-preproc_T1w (N4)?

For a few pilot sessions, fit a 2nd-order polynomial in scanner coordinates to log intensity in
white matter and report the fitted field's 5-95 % spread (in %). An N4-corrected image should be
nearly flat; an uncorrected multichannel-coil image typically shows a spread of tens of percent.
Images: fMRIPrep T1w / T2w (WM = fMRIPrep label-WM_probseg > 0.95), DAIRC rawavg / nu
(WM = aseg cerebral WM labels 2/41, resampled by nearest neighbour through scanner RAS).

    python directions/d1_multimodal/csd3_bias_check.py [REG_ROOT]   (compute node; ~1 min/session)

With REG_ROOT (tools/t2reg.sbatch output), also reports the T2w after dividing by its N4 field
(``fmriprep_T2w_n4``), which should be as flat as the fMRIPrep T1w.
"""

from __future__ import annotations

import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates

from abcd import paths
from abcd.hcp_stats import DEFAULT_FS_ROOT

FMRIPREP = Path("/rds/project/rds-CeXlNYOYMxw/derivatives/abcc-fmriprep_v25.1.4")


def field_spread(img: nib.Nifti1Image | nib.MGHImage, mask_xyz: np.ndarray,
                 field: nib.Nifti1Image | None = None) -> float:
    """5-95 % spread (%) of a quadratic fit to log intensity at scanner-RAS points."""
    pts = np.c_[mask_xyz, np.ones(len(mask_xyz))].T
    ijk = (np.linalg.inv(img.affine) @ pts)[:3]
    v = map_coordinates(np.asarray(img.dataobj, dtype=np.float32), ijk, order=1)
    if field is not None:
        fijk = (np.linalg.inv(field.affine) @ pts)[:3]
        v = v / map_coordinates(np.asarray(field.dataobj, dtype=np.float32), fijk, order=1, mode="nearest")
    ok = v > 0
    x = mask_xyz[ok] - mask_xyz[ok].mean(0)
    x /= x.std(0)
    X = np.c_[np.ones(len(x)), x, x ** 2, x[:, [0]] * x[:, [1]], x[:, [0]] * x[:, [2]], x[:, [1]] * x[:, [2]]]
    beta, *_ = np.linalg.lstsq(X, np.log(v[ok]), rcond=None)
    f = X @ beta
    lo, hi = np.percentile(f, [5, 95])
    return float(100 * (np.exp(hi - lo) - 1))


def wm_points(seg: nib.spatialimages.SpatialImage, sel, step: int = 2) -> np.ndarray:
    d = np.asarray(seg.dataobj)
    ijk = np.argwhere(sel(d))[::step ** 3]
    return (seg.affine @ np.c_[ijk, np.ones(len(ijk))].T)[:3].T


def main() -> int:
    pilot = pd.read_csv(paths.release_dir("7.0") / "processed" / "hcp" / "surf_metrics_pilot_sessions.tsv", sep="\t")
    reg_root = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    rows = []
    for sub, ses, manu in zip(pilot.participant_id, pilot.session_id, pilot.ab_g_dyn__design_mr__manufact):
        anat = FMRIPREP / sub / ses / "anat"
        t1 = sorted(anat.glob("*_desc-preproc_T1w.nii.gz"))
        t2 = sorted(anat.glob("*_desc-preproc_T2w.nii.gz"))
        wm = sorted(anat.glob("*_label-WM_probseg.nii.gz"))
        if not (t1 and t2 and wm):
            continue
        fs = DEFAULT_FS_ROOT / sub / ses / "mri"
        pf = wm_points(nib.load(str(wm[0])), lambda d: d > 0.95)
        pa = wm_points(nib.load(str(fs / "aseg.mgz")), lambda d: np.isin(d, (2, 41)))
        r = {"participant_id": sub, "session_id": ses, "manufacturer": manu,
             "fmriprep_T1w": field_spread(nib.load(str(t1[0])), pf),
             "fmriprep_T2w": field_spread(nib.load(str(t2[0])), pf),
             "dairc_rawavg": field_spread(nib.load(str(fs / "rawavg.mgz")), pa),
             "dairc_nu": field_spread(nib.load(str(fs / "nu.mgz")), pa)}
        fp = reg_root / sub / ses / "t2_biasfield_8mm.nii.gz" if reg_root else None
        if fp is not None and fp.exists():
            r["fmriprep_T2w_n4"] = field_spread(nib.load(str(t2[0])), pf, nib.load(str(fp)))
        print(r, flush=True)
        rows.append(r)
        if len(rows) >= 9:
            break
    df = pd.DataFrame(rows)
    print("\nWM bias-field 5-95 % spread (%), quadratic fit:")
    print(df.round(1).to_string(index=False))
    print("\nmedian by manufacturer:")
    print(df.groupby("manufacturer").median(numeric_only=True).round(1).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
