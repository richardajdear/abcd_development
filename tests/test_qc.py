"""Tests for the QC predicate framework."""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import io, paths, qc  # noqa: E402

try:
    paths.abcd_root()
    HAVE_DATA = True
except paths.DataRootError:
    HAVE_DATA = False

needs_data = pytest.mark.skipif(not HAVE_DATA, reason="ABCD release tree not available")


@pytest.fixture
def toy_scans():
    return pd.DataFrame({
        "subject": ["s1", "s1", "s1", "s2", "s2", "s3"],
        "visit": ["v0", "v2", "v4", "v0", "v2", "v0"],
        "n_regions_present": [68, 68, 67, 68, 68, 68],
        "n_regions_total": [68] * 6,
    })


def test_complete_regions_drops_incomplete(toy_scans):
    r = qc.apply_qc(toy_scans, [qc.complete_regions_predicate(68)], "test")
    assert r.n_scans == 5
    assert r.ledger.iloc[-1].n_dropped == 1


def test_unavailable_predicate_passes_all_rows_but_is_logged(toy_scans):
    def boom(df):
        raise qc.SourceUnavailable("no such table in this release")

    p = qc.Predicate("phantom", "a criterion with no data", boom)
    r = qc.apply_qc(toy_scans, [p], "test")
    assert r.n_scans == len(toy_scans), "unavailable predicate must not drop rows"
    assert r.unavailable == ["phantom"]
    assert r.ledger.iloc[-1].status == "unavailable"
    assert "no such table" in r.ledger.iloc[-1].detail


def test_require_available_raises_for_inert_criterion(toy_scans):
    def boom(df):
        raise qc.SourceUnavailable("missing")

    r = qc.apply_qc(toy_scans, [qc.Predicate("phantom", "", boom)], "test")
    with pytest.raises(qc.SourceUnavailable, match="refusing to proceed"):
        r.require_available("phantom")
    r.require_available()  # no names -> no error


def test_ledger_is_a_funnel(toy_scans):
    preds = [qc.complete_regions_predicate(68), qc.complete_regions_predicate(68)]
    led = qc.apply_qc(toy_scans, preds, "test").ledger
    # each step's input equals the previous step's output
    assert (led.n_scans_in.iloc[1:].values == led.n_scans_out.iloc[:-1].values).all()
    assert (led.n_dropped == led.n_scans_in - led.n_scans_out).all()


def test_visit_pattern_is_temporally_ordered(toy_scans):
    pat = qc.visit_pattern(toy_scans, ("v0", "v2", "v4"))
    assert pat["s1"] == "v0+v2+v4"
    assert pat["s2"] == "v0+v2"
    assert pat["s3"] == "v0"


def test_visit_pattern_order_independent_of_row_order(toy_scans):
    shuffled = toy_scans.iloc[[2, 0, 1, 4, 3, 5]].reset_index(drop=True)
    a = qc.visit_pattern(toy_scans, ("v0", "v2", "v4"))
    b = qc.visit_pattern(shuffled, ("v0", "v2", "v4"))
    pd.testing.assert_series_equal(a.sort_index(), b.sort_index())


@pytest.mark.parametrize("rule,expected", [
    ("all", {"s1", "s2", "s3"}),
    ("min_visits", {"s1", "s2"}),
    ("has_last", {"s1"}),
])
def test_sample_rules(toy_scans, rule, expected):
    subs, info = qc.select_subjects(toy_scans, rule, 2, ("v0", "v2", "v4"))
    assert set(subs) == expected
    assert info["rule"] == rule


def test_unknown_policy_and_rule_raise(toy_scans):
    with pytest.raises(ValueError, match="unknown QC policy"):
        qc.build_policy("nonsense", None)
    with pytest.raises(ValueError, match="unknown sample rule"):
        qc.select_subjects(toy_scans, "nonsense", 2, ("v0",))


def test_none_policy_is_empty():
    assert qc.build_policy("none", None) == []


@needs_data
def test_legacy_policy_drops_expected_scans():
    """Regression lock on the published 5.1 QC numbers."""
    a = io.Release51Adapter()
    scans = qc.scan_table(a.imaging("thickness"))
    r = qc.apply_qc(scans, qc.build_policy("legacy_euler", a, 68), "legacy_euler")
    assert len(scans) == 22854, f"input scan count drifted: {len(scans)}"
    assert r.n_scans == 20162, f"legacy QC output drifted: {r.n_scans}"
    assert r.n_subjects == 11162
    assert r.unavailable == []


@needs_data
def test_coded_policy_is_inert_on_51_and_says_so():
    """The local 5.1 copy has no release-native QC sources.

    This is a documentation test: if a future release drops these files in,
    it will fail and prompt an update rather than silently changing the sample.
    """
    a = io.Release51Adapter()
    scans = qc.scan_table(a.imaging("thickness"))
    r = qc.apply_qc(scans, qc.build_policy("coded", a, 68), "coded")
    assert set(r.unavailable) == {
        "release_qc_include", "surface_defects", "scanner_manufacturer"
    }
