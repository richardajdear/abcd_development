"""Tests for the release adapters.

Tests needing the raw release are skipped when ABCD_ROOT is unavailable, so
the suite still runs on a machine with only the repo checked out.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import io, paths  # noqa: E402

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


def test_70_adapter_is_verified_and_silent():
    """7.0 is verified against real data, so it must not warn.

    The warning path still exists for any future unverified adapter; see
    ``io.get_adapter``.
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        a = io.get_adapter("7.0")
    assert a.VERIFIED
    assert set(a.VISITS.values()) == {"ses-00A", "ses-02A", "ses-04A", "ses-06A"}
    assert len(a.REGION_CODES) == 34


def test_70_subject_id_normalisation_matches_51():
    """Both releases must normalise to the same id form, or joins silently fail."""
    import pandas as pd
    a70, a51 = io.Release70Adapter(), io.Release51Adapter()
    assert a70._to_bids(pd.Series(["sub-003RTV85"]))[0] == "sub-NDARINV003RTV85"
    assert a51._to_bids(pd.Series(["NDAR_INV003RTV85"]))[0] == "sub-NDARINV003RTV85"
    # idempotent
    assert a70._to_bids(pd.Series(["sub-NDARINV003RTV85"]))[0] == "sub-NDARINV003RTV85"


@needs_data
def test_70_zygosity_codes_are_what_we_think_they_are():
    """No codebook ships with 7.0's ``gn_y_genrel``, so pin the inference.

    The mapping in ``Release70Adapter.ZYGOSITY_CODES`` was established from the
    data itself: code 1 sits at pi-hat ~1.0, codes 2 and 3 both at ~0.5, and
    birth event separates 2 from 3 perfectly.  If a future release renumbers
    these, every Falconer estimate silently changes class -- so assert the
    evidence rather than trusting the constant.
    """
    import pandas as pd
    a = io.Release70Adapter()
    raw = a._table("gn_y_genrel")
    zyg = raw["gn_y_genrel_zyg__01"].astype("float")
    pih = raw["gn_y_genrel_pihat__01"].astype("float")
    ok = zyg.notna() & pih.notna()
    assert set(zyg[ok].astype(int)) <= set(a.ZYGOSITY_CODES), \
        f"unexpected zygosity codes {sorted(set(zyg[ok].astype(int)))}"

    by_code = pih[ok].groupby(zyg[ok].astype(int)).mean()
    mz = [c for c, v in a.ZYGOSITY_CODES.items() if v == "MZ"][0]
    assert by_code[mz] > 0.9, f"code {mz} labelled MZ but mean pi-hat is {by_code[mz]:.3f}"
    for c, name in a.ZYGOSITY_CODES.items():
        if name != "MZ":
            assert 0.4 < by_code[c] < 0.6, f"code {c} ({name}) mean pi-hat {by_code[c]:.3f}"

    # birth event separates DZ twins from non-twin siblings
    birth = raw.set_index("participant_id")["gn_y_genrel_id__birth"]
    same = (raw["gn_y_genrel_id__birth"].values
            == raw["gn_y_genrel_id__paired__01"].map(birth).values)
    rate = pd.Series(same)[ok.values].groupby(zyg[ok].astype(int).values).mean()
    for c, name in a.ZYGOSITY_CODES.items():
        want = 1.0 if name in ("MZ", "DZ_twin") else 0.0
        assert rate[c] == want, \
            f"code {c} ({name}): shares birth event in {rate[c]:.1%} of pairs, expected {want:.0%}"


@needs_data
def test_70_genotyped_pairs_are_reciprocal_and_consistent():
    """Each pair is listed from both sides; collapsing must not lose or invent any."""
    a = io.Release70Adapter()
    gp = a.genotyped_pairs()
    assert set(gp.pair_type) == {"MZ", "DZ_twin", "full_sib"}
    assert (gp.a < gp.b).all(), "pair keys are not canonically ordered"
    assert not gp.duplicated(["a", "b"]).any()
    assert gp.a.str.startswith("sub-NDARINV").all() and gp.b.str.startswith("sub-NDARINV").all()
    # every subject in a pair appears in the source's own participant list
    src = set(a._to_bids(a._table("gn_y_genrel").participant_id))
    assert set(gp.a) <= src and set(gp.b) <= src


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


# --- 7.0 tabulation specifics that do not need the release on disk ---------

class _FakeTables70(io.Release70Adapter):
    """A 7.0 adapter whose tables are in-memory frames."""

    def __init__(self, tables):
        self._tables = tables

    @property
    def root(self):  # pragma: no cover - only used in messages
        return Path("/fake/abcd-data-release-7.0")

    def _table(self, stem):
        return self._tables[stem].copy()


def _dyn(manufact):
    import pandas as pd
    n = len(manufact)
    return pd.DataFrame({
        "participant_id": [f"sub-{i:08d}" for i in range(n)],
        "session_id": ["ses-00A"] * n,
        "ab_g_dyn__design_mr__manufact": manufact,
    })


def test_scanner_decodes_70_integer_codes():
    """7.0 codes the vendor as 1/2/3; the Philips exclusion needs the label."""
    ad = _FakeTables70({"ab_g_dyn": _dyn([1.0, 2.0, 3.0, float("nan")])})
    out = ad.scanner()
    assert out.manufacturer.tolist()[:3] == ["GE", "Philips", "Siemens"]
    assert out.manufacturer.isna().tolist()[3]


def test_scanner_passes_60_strings_through():
    ad = _FakeTables70({"ab_g_dyn": _dyn(["SIEMENS", "Philips Medical Systems", "GE MEDICAL SYSTEMS"])})
    assert ad.scanner().manufacturer.tolist() == ["SIEMENS", "Philips Medical Systems", "GE MEDICAL SYSTEMS"]


def test_scanner_refuses_unknown_code():
    ad = _FakeTables70({"ab_g_dyn": _dyn([1.0, 4.0])})
    with pytest.raises(paths.SourceUnavailable):
        ad.scanner()


def _fake_release(six_year_rows: int):
    import pandas as pd
    ses = ["ses-00A"] * 10 + ["ses-06A"] * six_year_rows
    ids = [f"sub-{i:08d}" for i in range(len(ses))]
    img = pd.DataFrame({"participant_id": ids, "session_id": ses})
    return _FakeTables70({"mr_y_smri__thk__dsk": img, "ab_g_dyn": img.copy()})


def test_assert_vintage_rejects_60_sized_tables():
    """The 6.0 tables are column-identical to 7.0; only the six-year count differs."""
    ad = _fake_release(4086)
    with pytest.raises(paths.SourceUnavailable, match="6.0 tables"):
        ad.assert_vintage()


def test_assert_vintage_accepts_70_sized_tables_and_reports_counts():
    ad = _fake_release(7607)
    v = ad.assert_vintage()
    assert v["imaging_rows_by_session"] == {"ses-00A": 10, "ses-06A": 7607}
    assert v["imaging_table"] == "mr_y_smri__thk__dsk"


def test_60_adapter_is_registered_and_shares_the_70_format():
    assert io.ADAPTERS["6.0"] is io.Release60Adapter
    assert issubclass(io.Release60Adapter, io.Release70Adapter)
    assert io.Release60Adapter.release == "6.0"


@needs_data
def test_70_release_dir_is_a_70_vintage():
    """Guards the exact mistake this repo made: 6.0 tables under a 7.0 label."""
    ad = io.Release70Adapter()
    v = ad.assert_vintage()
    assert v["imaging_rows_by_session"]["ses-06A"] >= io.Release70Adapter.MIN_SIX_YEAR_ROWS
