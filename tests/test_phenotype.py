"""Tests for phenotype construction and reliability."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import phenotype  # noqa: E402

RUN = Path(__file__).resolve().parents[1] / "out" / "thickness_dsk_51_f22d9bdb4145"
needs_fits = pytest.mark.skipif(
    not (RUN / "fits" / "blups.parquet").exists(),
    reason="fitted run not present; run assemble + R/fit_lmm.R first",
)


def _blups_varcomp(corr: float):
    """One region, three subjects, with intercept-slope correlation ``corr``.

    ``tau2`` is set small and the conditional SD smaller still, which is the
    numeric signature of a boundary fit: 1 - v/tau2 comes out near 1 even
    though tau2 is not identified.
    """
    blups = pd.DataFrame({
        "label": ["rA"] * 3, "subject": ["s1", "s2", "s3"],
        "re_slope": [0.01, -0.01, 0.0], "re_intercept": [0.1, -0.1, 0.0],
        "se_re_slope": [1e-4] * 3, "se_re_intercept": [1e-3] * 3,
    })
    varcomp = pd.DataFrame([
        {"label": "rA", "grp": "subject", "var1": "age_c",
         "var2": None, "vcov": 1e-6, "sdcor": 1e-3},
        {"label": "rA", "grp": "subject", "var1": "(Intercept)",
         "var2": None, "vcov": 1e-2, "sdcor": 0.1},
        {"label": "rA", "grp": "subject", "var1": "(Intercept)",
         "var2": "age_c", "vcov": -1e-4, "sdcor": corr},
    ])
    return blups, varcomp


def test_boundary_fit_reliability_is_nan_not_high():
    """A degenerate RE covariance must not report near-perfect reliability.

    Regression test.  Boundary fits gave a finite reliability of ~0.87 against
    ~0.10 in well-identified regions, which inflated every mean over regions and
    made release 5.1 appear MORE reliable than 7.0 despite shorter follow-up.
    """
    b_ok = phenotype.slope_reliability(*_blups_varcomp(-0.45))
    assert b_ok.reliability_slope.notna().all()
    assert b_ok.reliability_slope.iloc[0] > 0.9  # the spurious-looking value...

    b_bad = phenotype.slope_reliability(*_blups_varcomp(-1.0))
    assert b_bad.reliability_slope.isna().all()  # ...is suppressed on the boundary
    assert b_bad.reliability_intercept.isna().all()


def test_boundary_detection_is_not_triggered_by_strong_but_valid_correlation():
    """r = -0.99 is strong but identified; only the boundary itself is masked."""
    b = phenotype.slope_reliability(*_blups_varcomp(-0.99))
    assert b.reliability_slope.notna().all()


def test_design_reliability_monotone_in_span():
    """Reliability must increase with span at fixed visit count.

    This is the formal statement of the project's central design finding:
    a third visit inside the same window buys almost nothing, while widening
    the window buys a lot.  If this ever fails, the formula is wrong.
    """
    tau2, s2 = 0.0045552 ** 2, 0.0257443 ** 2
    r_2v_2y = phenotype.design_reliability([0, 2], tau2, s2)
    r_2v_4y = phenotype.design_reliability([0, 4], tau2, s2)
    r_3v_4y = phenotype.design_reliability([0, 2, 4], tau2, s2)
    r_3v_6y = phenotype.design_reliability([0, 3, 6], tau2, s2)
    assert r_2v_2y < r_2v_4y < r_3v_6y
    # the load-bearing comparison: a mid-window visit is nearly free of gain
    assert abs(r_3v_4y - r_2v_4y) < 0.02
    # widening the window is not
    assert r_3v_6y - r_3v_4y > 0.15


def test_design_reliability_degenerate():
    assert phenotype.design_reliability([1.0], 1.0, 1.0) == 0.0
    assert phenotype.design_reliability([], 1.0, 1.0) == 0.0


def test_design_reliability_bounds():
    tau2, s2 = 0.01, 0.01
    for t in ([0, 1], [0, 2, 4], [0, 1, 2, 3, 4, 5]):
        r = phenotype.design_reliability(t, tau2, s2)
        assert 0.0 <= r <= 1.0


@needs_fits
def test_phenotypes_are_tidy_and_complete():
    ph = phenotype.build_phenotypes(RUN)
    assert set(ph.phenotype) == {"slope", "slope_dev", "intercept"}
    assert ph.value.notna().all()
    # one row per subject x region x phenotype
    assert not ph.duplicated(["subject", "label", "phenotype"]).any()
    n_sub, n_reg = ph.subject.nunique(), ph.label.nunique()
    assert len(ph) == n_sub * n_reg * 3


@needs_fits
def test_slope_and_slope_dev_differ_by_a_constant_per_region():
    """slope = slope_dev + fixed slope, so within a region they must be
    perfectly correlated.  This documents that they give the SAME GWAS --
    the reason the 5.1 pipeline's several 'different' phenotypes were not
    actually independent."""
    ph = phenotype.build_phenotypes(RUN)
    for lab in list(ph.label.unique())[:5]:
        s = ph[(ph.label == lab) & (ph.phenotype == "slope")].set_index("subject").value
        d = ph[(ph.label == lab) & (ph.phenotype == "slope_dev")].set_index("subject").value
        diff = (s - d.reindex(s.index)).round(12)
        assert diff.nunique() == 1, f"{lab}: offset is not constant"


@needs_fits
def test_reliability_in_unit_interval_and_ordered_by_visits():
    ph = phenotype.build_phenotypes(RUN)
    r = ph.reliability.dropna()
    assert r.between(0.0, 1.0).all()
    summ = phenotype.reliability_summary(ph)
    sl = summ[summ.phenotype == "slope"].set_index("n_visits").mean_reliability
    assert sl.loc[3] > sl.loc[2], "3-visit slopes must be more reliable than 2-visit"
    # regression lock on the headline numbers
    assert 0.14 < sl.loc[2] < 0.17
    assert 0.18 < sl.loc[3] < 0.21


@needs_fits
def test_intercept_far_more_reliable_than_slope():
    """The project's framing depends on this gap being large."""
    summ = phenotype.reliability_summary(phenotype.build_phenotypes(RUN))
    mi = summ[summ.phenotype == "intercept"].mean_reliability.mean()
    ms = summ[summ.phenotype == "slope"].mean_reliability.mean()
    assert mi > 3 * ms


@needs_fits
def test_to_gcta_shape_and_id_map():
    ph = phenotype.build_phenotypes(RUN)
    labs = sorted(ph.label.unique())[:4]
    g = phenotype.to_gcta(ph, "slope", labels=labs)
    assert list(g.columns[:2]) == ["FID", "IID"]
    assert len(g.columns) == 2 + len(labs)
    assert (g.FID == g.IID).all()
    # an id map must be applied to BOTH id columns
    idmap = pd.Series({s: s.replace("sub-", "G") for s in ph.subject.unique()})
    g2 = phenotype.to_gcta(ph, "slope", labels=labs, id_map=idmap)
    assert g2.IID.str.startswith("G").all()
    assert (g2.FID == g2.IID).all()


@needs_fits
def test_to_gcta_reliability_filter_reduces_rows():
    ph = phenotype.build_phenotypes(RUN)
    labs = sorted(ph.label.unique())[:2]
    lo = phenotype.to_gcta(ph, "slope", labels=labs)
    hi = phenotype.to_gcta(ph, "slope", labels=labs, min_reliability=0.5)
    assert len(hi) < len(lo)


def test_relatedness_ids_are_normalised():
    """The genetics tables use NDAR_INV...; imaging uses sub-NDARINV....

    Joining unnormalised gives zero rows and looks like a null result, so
    this is guarded explicitly.
    """
    from abcd import io as abcd_io
    from abcd.paths import abcd_root

    if not (abcd_root() / "abcd-data-release-5.1").exists():
        pytest.skip("release 5.1 not available")
    rel = abcd_io.Release51Adapter().relatedness()
    assert rel.subject.str.startswith("sub-NDAR").all()
    assert not rel.subject.str.contains("_").any()
    assert set(rel.pair_type) <= {"MZ", "DZ_or_sib", "other", "unknown"}
    assert (rel.pair_type == "MZ").sum() > 500
