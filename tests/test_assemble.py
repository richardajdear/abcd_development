"""Tests for the assembly step."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import assemble, io, paths  # noqa: E402
from abcd.config import RunConfig  # noqa: E402

def _have_release_tables() -> bool:
    """True only if a release *table* is readable, not merely a release directory.

    ``abcd-data-release-7.0/`` may hold only derived ``processed/hcp/`` output (built on
    CSD3 before the tabulated release was in place); that is not release data.
    """
    try:
        paths.abcd_root()
    except paths.DataRootError:
        return False
    for rel, stem in (("7.0", "ab_g_dyn"), ("5.1", "abcd_y_lt")):
        try:
            paths.find_table(paths.release_dir(rel), stem)
            return True
        except paths.DataRootError:
            continue
    return False


HAVE_DATA = _have_release_tables()

needs_data = pytest.mark.skipif(not HAVE_DATA, reason="ABCD release tree not available")

CFG = Path(__file__).resolve().parents[1] / "configs"


def test_derived_metrics_registered():
    assert "t1t2_ratio" in assemble.DERIVED_METRICS
    requires, _ = assemble.DERIVED_METRICS["t1t2_ratio"]
    assert requires == ("t1_gray", "t2_gray")


def test_register_derived_adds_to_registry():
    @assemble.register_derived("_test_metric", ("thickness",))
    def _fn(frames):
        return frames["thickness"]

    try:
        assert "_test_metric" in assemble.DERIVED_METRICS
    finally:
        del assemble.DERIVED_METRICS["_test_metric"]


def test_t1t2_ratio_arithmetic():
    key = ["subject", "visit", "hemi", "region", "label", "is_global"]
    base = pd.DataFrame({
        "subject": ["s1", "s1"], "visit": ["v0", "v0"],
        "hemi": ["lh", "lh"], "region": ["insula", "cuneus"],
        "label": ["lh_insula", "lh_cuneus"], "is_global": [False, False],
    })
    t1 = base.assign(value=[2.0, 4.0])
    t2 = base.assign(value=[1.0, 0.0])  # a zero denominator must become NaN
    out = assemble.DERIVED_METRICS["t1t2_ratio"][1]({"t1_gray": t1, "t2_gray": t2})
    assert out.loc[out.region == "insula", "value"].iloc[0] == 2.0
    assert np.isnan(out.loc[out.region == "cuneus", "value"].iloc[0])
    assert (out.metric == "t1t2_ratio").all()


def test_unsupported_global_covariate_rejected():
    """The 5.1 'predicted mean' option is deliberately not offered."""
    with pytest.raises(Exception):
        RunConfig(global_covariate="predicted_mean")


# --- integration against the real release ---------------------------------

@pytest.fixture(scope="module")
def baseline_table():
    cfg = RunConfig.from_yaml(CFG / "ct_baseline.yaml")
    table, manifest = assemble.assemble(cfg, verbose=False)
    return table, manifest


@needs_data
def test_schema_and_no_nulls(baseline_table):
    table, _ = baseline_table
    expected = {
        "subject", "visit", "metric", "release", "hemi", "region", "label",
        "value", "age", "age_c", "sex", "site", "family_id",
        "global_value", "global_between", "global_between_c", "global_within",
        "age_first", "age_last", "age_span", "n_visits",
    }
    assert set(table.columns) == expected
    assert not table.isna().any().any(), "assembled table must have no nulls"


@needs_data
def test_baseline_sample_size_locked(baseline_table):
    """Regression lock: changes to the sample must be deliberate."""
    _, m = baseline_table
    f = m["final"]
    assert f["n_subjects"] == 6937, f"sample size drifted: {f['n_subjects']}"
    assert f["n_regions"] == 68
    assert f["n_sites"] == 22
    assert f["n_families"] == 5989


@needs_data
def test_min_visits_rule_is_respected(baseline_table):
    table, _ = baseline_table
    per_subject = table.groupby("subject").visit.nunique()
    assert per_subject.min() >= 2


@needs_data
def test_global_within_sums_to_zero_per_subject(baseline_table):
    """The within/between split must be an exact decomposition."""
    table, _ = baseline_table
    scans = table.drop_duplicates(["subject", "visit"])
    per_subject = scans.groupby("subject").global_within.mean()
    assert per_subject.abs().max() < 1e-9
    # and the parts must reconstruct the whole
    recon = table.global_between + table.global_within
    assert np.allclose(recon, table.global_value)


@needs_data
def test_age_centring_is_recorded_and_applied(baseline_table):
    table, m = baseline_table
    assert np.isclose(table.age.mean() - table.age_c.mean(), m["age_centre"], atol=1e-3)


@needs_data
def test_thickness_declines_across_adolescence(baseline_table):
    """Sanity check on the science: cortex thins between ages 9 and 15."""
    table, _ = baseline_table
    means = table.groupby("visit", observed=True).value.mean()
    assert means["v0"] > means["v2"] > means["v4"], means.to_dict()


@needs_data
def test_manifest_records_qc_and_sample(baseline_table):
    _, m = baseline_table
    assert m["qc_ledger"][0]["criterion"] == "initial"
    assert m["sample"]["rule"] == "min_visits"
    assert m["config_hash"] == m["config"]["hash"] if "hash" in m["config"] else True


@needs_data
def test_legacy_replication_matches_thesis_sample():
    """The legacy config must reproduce the 5.1 analysis sample.

    The thesis used left hemisphere only, the static Euler list, and required
    a year-4 scan.
    """
    cfg = RunConfig.from_yaml(CFG / "ct_legacy_replication.yaml")
    table, m = assemble.assemble(cfg, verbose=False)
    assert set(table.hemi.unique()) == {"lh"}
    assert m["final"]["n_regions"] == 34
    assert (table.visit == "v4").groupby(table.subject).any().all(), \
        "has_last rule must guarantee every subject has a year-4 scan"


@needs_data
def test_derived_metric_assembles_end_to_end():
    cfg = RunConfig.from_yaml(CFG / "t1t2_baseline.yaml")
    table, m = assemble.assemble(cfg, verbose=False)
    assert (table.metric == "t1t2_ratio").all()
    assert m["final"]["n_regions"] == 68
    # T1w/T2w ratio should be positive and O(1)
    assert table.value.min() > 0
    assert 0.5 < table.value.median() < 5
