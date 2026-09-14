"""Tests for the HCP-MMP extraction (abcd.hcp_stats) and its 7.0 adapter path.

The parser tests run on a synthetic ``mris_anatomical_stats`` table and need
no data.  The adapter test needs ``abcd-data-release-7.0/processed/hcp/`` and skips
otherwise.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import hcp_stats, io, paths  # noqa: E402

HEADER = """Using TH3 vertex volume calc
Total face volume 280344
Total vertex volume 280169 (mask=0)

table columns are:
    number of vertices
    total surface area (mm^2)
    total gray matter volume (mm^3)
    average cortical thickness +- standard deviation (mm)
    integrated rectified mean curvature
    integrated rectified Gaussian curvature
    folding index
    intrinsic curvature index
    structure name

"""

ROWS_LH = (
    " 9288   5697   3632  0.621 1.081     0.090     0.060     1901    35.8  ???\n"
    " 5852   3949   7030  1.773 0.467     0.132     0.033       74     7.8  L_V1_ROI\n"
    "  386    240    453  2.400 0.406     0.081     0.012        1     0.3  L_MST_ROI\n"
)
ROWS_RH = (
    " 9000   5500   3500  0.600 1.000     0.090     0.060     1900    35.0  ???\n"
    " 6000   4000   7200  1.800 0.470     0.130     0.030       70     7.0  R_V1_ROI\n"
    "  400    250    460  2.500 0.400     0.080     0.010        1     0.3  R_MST_ROI\n"
)


def _write_session(root: Path, sub="sub-TEST0001", ses="ses-00A"):
    d = root / sub / ses / hcp_stats.ATLAS
    d.mkdir(parents=True)
    (d / f"lh.{hcp_stats.ATLAS}.log").write_text(HEADER + ROWS_LH)
    (d / f"rh.{hcp_stats.ATLAS}.log").write_text(HEADER + ROWS_RH)
    return hcp_stats.Session(sub, ses)


def test_parse_anatomical_stats(tmp_path):
    p = tmp_path / "lh.log"
    p.write_text(HEADER + ROWS_LH)
    df = hcp_stats.parse_anatomical_stats(p)
    assert list(df.structure) == ["???", "L_V1_ROI", "L_MST_ROI"]
    v1 = df.set_index("structure").loc["L_V1_ROI"]
    assert v1.nverts == 5852 and v1.thickness_mm == 1.773 and v1.thickness_sd == 0.467
    assert v1.folding_index == 74 and v1.intrinsic_curv == 7.8


def test_parse_rejects_malformed_numeric_line(tmp_path):
    p = tmp_path / "lh.log"
    p.write_text(HEADER + " 5852   3949   7030  1.773  L_V1_ROI\n")
    with pytest.raises(ValueError):
        hcp_stats.parse_anatomical_stats(p)


def test_region_names():
    assert hcp_stats.region_from_structure("L_V1_ROI") == "V1"
    assert hcp_stats.region_from_structure("R_a9-46v_ROI") == "a9-46v"
    assert hcp_stats.region_from_structure("???") is None


def test_parse_aseg_holes(tmp_path):
    p = tmp_path / "aseg.stats"
    p.write_text(
        "# Measure lhSurfaceHoles, lhSurfaceHoles, Number of defect holes in lh surfaces prior to fixing, 13, unitless\n"
        "# Measure rhSurfaceHoles, rhSurfaceHoles, Number of defect holes in rh surfaces prior to fixing, 9, unitless\n"
        "# Measure SurfaceHoles, SurfaceHoles, Total number of defect holes in surfaces prior to fixing, 22, unitless\n"
        "# ColHeaders  Index SegId NVoxels\n  1 4 100\n"
    )
    h = hcp_stats.parse_aseg_holes(p)
    assert h == {"lh_holes": 13.0, "rh_holes": 9.0, "total_holes": 22.0}
    assert all(pd.isna(v) for v in hcp_stats.parse_aseg_holes(tmp_path / "nope").values())


def test_extract_session_enforces_180_parcels(tmp_path):
    """A synthetic session has 2 parcels per hemisphere, so it must be rejected."""
    sess = _write_session(tmp_path / "parc")
    long, qc = hcp_stats.extract_session(sess, tmp_path / "parc", tmp_path / "fs")
    assert long is None
    assert qc["status"].startswith("parcel_count_2_2")
    assert qc["medial_wall_verts"] == 9288 + 9000
    assert qc["recon_done"] is False


def test_extract_session_ok_when_parcel_count_relaxed(tmp_path, monkeypatch):
    monkeypatch.setattr(hcp_stats, "N_PARCELS_PER_HEMI", 2)
    sess = _write_session(tmp_path / "parc")
    long, qc = hcp_stats.extract_session(sess, tmp_path / "parc", tmp_path / "fs")
    assert qc["status"] == "ok"
    assert len(long) == 4 and set(long.label) == {"lh_V1", "lh_MST", "rh_V1", "rh_MST"}
    assert "???" not in set(long.region)


def test_wide_table_is_area_weighted(tmp_path, monkeypatch):
    monkeypatch.setattr(hcp_stats, "N_PARCELS_PER_HEMI", 2)
    sess = _write_session(tmp_path / "parc")
    long, _ = hcp_stats.extract_session(sess, tmp_path / "parc", tmp_path / "fs")
    w = hcp_stats.wide_table(long, "thickness_mm").iloc[0]
    assert w["mr_y_smri__thk__hcp__V1__lh_mean"] == 1.773
    assert w["mr_y_smri__thk__hcp__MST__rh_mean"] == 2.5
    lh = (3949 * 1.773 + 240 * 2.400) / (3949 + 240)          # area-weighted
    assert w["mr_y_smri__thk__hcp__lh_mean"] == pytest.approx(lh)
    tot = (3949 * 1.773 + 240 * 2.4 + 4000 * 1.8 + 250 * 2.5) / (3949 + 240 + 4000 + 250)
    assert w["mr_y_smri__thk__hcp_mean"] == pytest.approx(tot)
    a = hcp_stats.wide_table(long, "area_mm2").iloc[0]
    assert a["mr_y_smri__area__hcp__lh_sum"] == 3949 + 240
    assert a["mr_y_smri__area__hcp_sum"] == 3949 + 240 + 4000 + 250


# --- adapter path (needs the derived file on disk) -------------------------

def _hcp_available() -> bool:
    try:
        return (io.Release70Adapter().hcp_dir / "mr_y_smri__thk__hcp.tsv").exists()
    except paths.DataRootError:
        return False


@pytest.mark.skipif(not _hcp_available(), reason="abcd-data-release-7.0/processed/hcp not present")
def test_release70_hcp_thickness_long_format():
    df = io.Release70Adapter().imaging("thickness", "hcp")
    reg = df[~df.is_global]
    assert reg.region.nunique() == 180 and reg.label.nunique() == 360
    assert set(reg.hemi) == {"lh", "rh"}
    assert df.subject.str.startswith("sub-NDARINV").all()
    assert set(df.visit) <= {"v0", "v2", "v4", "v6"}
    g = df[df.is_global]
    assert (g.label == "global_mean").all() and len(g) == reg.groupby(["subject", "visit"]).ngroups
    # HCP labels must line up with the committed spin-test centroids
    cent = pd.read_csv(paths.DATA_DIR / "hcp_centroids.csv")
    assert set(reg.label) == set(cent.label)
