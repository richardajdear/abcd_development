"""
surface.py -- fsaverage annot helpers shared by 18_nspn_to_hcp.py,
19_resample_check.py and 20_basis_vs_parcellation.py.

Everything here is parcels <-> vertices on the fsaverage (164k) mesh, reading
the annot files in ~/Git/AHBA/data/parcellations/.

One trap this module exists to contain: the annot files do NOT name parcels
consistently. `lh.aparc.annot` calls a parcel `lh_bankssts` while
`rh.aparc.annot` calls the same parcel `bankssts`, and HCPMMP1 uses
`L_V1_ROI` / `R_V1_ROI`. Matching or grouping on the raw name silently splits
the hemispheres and lets `lh_unknown` through a background filter, so every
name goes through `strip_hemi` first.
"""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np, pandas as pd
import nibabel as nib

PARC = Path.home() / "Git" / "AHBA" / "data" / "parcellations"
BACKGROUND = ("unknown", "???", "corpuscallosum", "Background")


def strip_hemi(name: str) -> str:
    """Drop a hemisphere prefix and the HCP `_ROI` suffix from an annot label."""
    return re.sub(r"_ROI$", "", re.sub(r"^(lh_|rh_|L_|R_)", "", name))


def is_bg(name: str) -> bool:
    return strip_hemi(name).startswith(BACKGROUND)


def annot(name: str) -> tuple[np.ndarray, list[str]]:
    """(per-vertex label ids, label names) for an annot file in PARC."""
    lab, _, names = nib.freesurfer.read_annot(str(PARC / name))
    return lab, [n.decode() if isinstance(n, bytes) else n for n in names]


def paint(values: pd.Series, lab: np.ndarray, names: list[str], key=strip_hemi) -> np.ndarray:
    """Parcel values -> per-vertex values (NaN where the parcel has no value)."""
    out = np.full(lab.shape, np.nan)
    for i, n in enumerate(names):
        if is_bg(n):
            continue
        k = key(n)
        if k in values.index and np.isfinite(values[k]):
            out[lab == i] = values[k]
    return out


def collect(vert: np.ndarray, lab: np.ndarray, names: list[str], key=strip_hemi) -> pd.Series:
    """Per-vertex values -> parcel means (parcels with no finite vertex dropped)."""
    out = {}
    for i, n in enumerate(names):
        if is_bg(n):
            continue
        v = vert[lab == i]
        if np.isfinite(v).any():
            out[key(n)] = float(np.nanmean(v))
    return pd.Series(out)


def resample(values: pd.Series, src: str, dst: str,
             src_key=strip_hemi, dst_key=strip_hemi) -> pd.Series:
    """Parcel values from one annot to another, via the shared vertex mesh."""
    ls, ns = annot(src)
    ld, nd = annot(dst)
    assert ls.shape == ld.shape, "annots must share the fsaverage mesh"
    return collect(paint(values, ls, ns, src_key), ld, nd, dst_key)
