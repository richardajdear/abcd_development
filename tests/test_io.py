"""Tests for the release adapters.

Tests needing the raw release are skipped when ABCD_ROOT is unavailable, so
the suite still runs on a machine with only the repo checked out.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import io, paths  # noqa: E402

try:
    paths.abcd_root()
    HAVE_DATA = True
except paths.DataRootError:
    HAVE_DATA = False

needs_data = pytest.mark.skipif(not HAVE_DATA, reason="ABCD release tree not available")


# --- label table (committed to the repo, always testable) -----------------

def test_region_labels_shape():
    lab = io.region_labels("dsk")
    assert (~lab.is_global).sum() == 68, "expected 34 DK regions x 2 hemispheres"
    assert lab[~lab.is_global].region.nunique() == 34
    assert set(lab[~lab.is_global].hemi) == {"lh", "rh"}


def test_region_labels_have_both_naming_families():
    """ABCD spells 19/34 DK regions differently in T2 vs thickness tables."""
    lab = io.region_labels("dsk")
    differing = lab[lab.stem_a != lab.stem_b]
    assert len(differing) == 38, f"expected 19 regions x 2 hemispheres, got {len(differing)}"
    row = lab[lab.label == "lh_caudalanteriorcingulate"].iloc[0]
    assert row.stem_a == "cdacatelh" and row.stem_b == "cdatcgatelh"


def test_region_labels_unique_keys():
    lab = io.region_labels("dsk")
    assert not lab.duplicated(["hemi", "region"]).any()
    assert not lab.duplicated(["stem_a"]).any()


def test_unknown_parcellation_raises():
    with pytest.raises(ValueError, match="no region labels"):
        io.region_labels("nonexistent")


# --- adapter behaviour ----------------------------------------------------

def test_bids_subject_conversion():
    import pandas as pd
    s = pd.Series(["NDAR_INV00CY2MDM"])
    assert io.ReleaseAdapter._to_bids(s).iloc[0] == "sub-NDARINV00CY2MDM"


def test_derived_metric_rejected_by_adapter():
    a = io.Release51Adapter()
    with pytest.raises(KeyError, match="not a raw 5.1 table"):
        a.imaging("t1t2_ratio")


def test_70_adapter_warns_as_unverified():
    with pytest.warns(UserWarning, match="unverified stub"):
        io.get_adapter("7.0")


def test_51_adapter_does_not_warn():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        io.get_adapter("5.1")


@needs_data
def test_imaging_long_schema():
    df = io.Release51Adapter().imaging("thickness")
    assert set(df.columns) == {
        "subject", "visit", "metric", "hemi", "region", "label", "is_global", "value"
    }
    assert set(df.visit.unique()) <= {"v0", "v2", "v4"}
    assert df.subject.str.startswith("sub-NDAR").all()


@needs_data
def test_t2_alias_resolves_to_same_anatomy():
    """The stem_a/stem_b split must map to the *same* region, not shuffle it.

    Left/right thickness and T2 differ in absolute value, but the spatial
    pattern across regions should correlate strongly if the mapping is right,
    and the region sets must be identical.
    """
    a = io.Release51Adapter()
    thk = a.imaging("thickness")
    t2 = a.imaging("t2_gray")
    assert set(thk.label) == set(t2.label)
    # both should recover insula, which is spelled identically in both families
    assert (thk.label == "lh_insula").any() and (t2.label == "lh_insula").any()


@needs_data
def test_t2_alias_is_anatomically_correct():
    """Positional pairing of stem_a/stem_b must map to the same anatomy.

    The alias table was built by column position, so it needs an independent
    anatomical check. Two hold if the mapping is right and break if the
    regions are shuffled:

    * T1w and T2w grey-matter intensity are inversely related, so their
      regional profiles should be strongly *anti*-correlated (observed
      Spearman -0.76; a shuffled alias gives |rho| < 0.25 at the 95th
      percentile over 200 permutations).
    * Left and right hemisphere profiles should be near-identical
      (observed >= 0.98 for every metric).
    """
    from scipy import stats

    a = io.Release51Adapter()

    def profile(metric):
        df = a.imaging(metric)
        return df[~df.is_global].groupby("label").value.mean()

    t1, t2 = profile("t1_gray"), profile("t2_gray")
    common = t1.index.intersection(t2.index)
    assert len(common) == 68
    rho = stats.spearmanr(t1[common], t2[common]).statistic
    assert rho < -0.5, f"T1/T2 profiles should be strongly anticorrelated, got {rho:.3f}"

    for name, prof in [("t1_gray", t1), ("t2_gray", t2)]:
        lh = prof[[i for i in prof.index if i.startswith("lh_")]]
        rh = prof[[i for i in prof.index if i.startswith("rh_")]]
        rh.index = [i.replace("rh_", "lh_") for i in rh.index]
        sym = stats.spearmanr(lh, rh[lh.index]).statistic
        assert sym > 0.9, f"{name} lh/rh asymmetry suggests bad mapping: {sym:.3f}"


@needs_data
def test_family_id_is_subject_constant():
    lt = io.Release51Adapter().longitudinal()
    per_subject = lt.dropna(subset=["family_id"]).groupby("subject").family_id.nunique()
    assert (per_subject <= 1).all(), "family_id must not vary within a subject"


@needs_data
def test_demographics_one_row_per_subject():
    demo = io.Release51Adapter().demographics()
    assert not demo.subject.duplicated().any()
    assert set(demo.sex.dropna().unique()) <= {"M", "F"}


@needs_data
def test_verify_adapter_all_pass():
    result = io.verify_adapter("5.1")
    failed = result[~result.ok]
    assert failed.empty, f"failed checks:\n{failed.to_string(index=False)}"
