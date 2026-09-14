"""Tests for the Desikan parse of FreeSurfer ``aparc.stats`` (abcd.dk_stats)."""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import dk_stats  # noqa: E402

APARC = """# Table of FreeSurfer cortical parcellation anatomical statistics
# Measure Cortex, NumVert, Number of Vertices, 132766, unitless
# Measure Cortex, WhiteSurfArea, White Surface Total Area, 89850.6, mm^2
# Measure Cortex, MeanThickness, Mean Thickness, 2.71093, mm
# ColHeaders StructName NumVert SurfArea GrayVol ThickAvg ThickStd MeanCurv GausCurv FoldInd CurvInd
bankssts                                 1546   1043   3051  2.734 0.556     0.097     0.019       11     1.2
caudalanteriorcingulate                  1282    836   2706  2.510 0.852     0.122     0.022       19     1.0
"""


def test_parse_aparc_stats(tmp_path):
    p = tmp_path / "lh.aparc.stats"
    p.write_text(APARC)
    df, meas = dk_stats.parse_aparc_stats(p)
    assert list(df.region) == ["bankssts", "caudalanteriorcingulate"]
    assert df.iloc[0].thickness_mm == 2.734 and df.iloc[0].nverts == 1546
    assert meas == {"NumVert": 132766.0, "WhiteSurfArea": 89850.6, "MeanThickness": 2.71093}


def test_region_to_code_covers_all_34():
    codes = dk_stats._region_to_code()
    assert len(codes) == 34
    assert codes["bankssts"] == "bstmps" and codes["superiorfrontal"] == "sfrt"


def test_wide_table_uses_release_names_and_area_weighted_means():
    long = pd.DataFrame({
        "participant_id": ["sub-A"] * 4, "session_id": ["ses-00A"] * 4,
        "hemi": ["lh", "lh", "rh", "rh"],
        "region": ["bankssts", "superiorfrontal"] * 2,
        "nverts": [100, 200, 100, 200], "area_mm2": [1.0, 2.0, 3.0, 4.0],
        "gmv_mm3": [10.0, 20.0, 30.0, 40.0], "thickness_mm": [2.0, 3.0, 2.5, 3.5],
        "thickness_sd": 0.1, "mean_curv": 0.1, "gauss_curv": 0.1,
        "folding_index": 1.0, "intrinsic_curv": 0.1,
    })
    qc = pd.DataFrame([{"participant_id": "sub-A", "session_id": "ses-00A",
                        "cortex_nverts_lh": 1000, "cortex_nverts_rh": 3000,
                        "cortex_mean_thickness_lh": 2.6, "cortex_mean_thickness_rh": 3.0}])
    codes = {"bankssts": "bstmps", "superiorfrontal": "sfrt"}
    w = dk_stats.wide_table(long, "thickness_mm", qc, codes).iloc[0]
    assert w["mr_y_smri__thk__dsk__bstmps__lh_mean"] == 2.0
    assert w["mr_y_smri__thk__dsk__sfrt__rh_mean"] == 3.5
    # release convention: surface-area-weighted mean over regions
    assert w["mr_y_smri__thk__dsk__lh_mean"] == pytest.approx((1.0 * 2.0 + 2.0 * 3.0) / 3.0)
    assert w["mr_y_smri__thk__dsk_mean"] == pytest.approx((2.0 + 6.0 + 7.5 + 14.0) / 10.0)
    a = dk_stats.wide_table(long, "area_mm2", qc, codes).iloc[0]
    assert a["mr_y_smri__area__dsk__lh_sum"] == 3.0 and a["mr_y_smri__area__dsk_sum"] == 10.0
