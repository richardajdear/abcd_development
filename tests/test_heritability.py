"""Guards on the twin-design heritability estimator.

The estimator is arithmetically trivial; every bug we have actually hit was in
the *inputs* (ID spellings that silently join to nothing, phenotypes from a
model that partitions out the variance being estimated).  These tests target
that, not the formula.
"""

import numpy as np
import pandas as pd
import pytest

from abcd.heritability import falconer, _normalise_ids, _SHARING_GAP


def _synthetic(h2_true, n_mz=5000, n_dz=5000, seed=0, dz_label="DZ_twin",
               r_sib=None, n_sib=5000):
    """Build pair data with a known h2 in the explicit-pair contract.

    Returns ``(y, pairs)`` where ``pairs`` has one row per unordered pair with
    ``a``/``b`` subject columns -- the shape :func:`abcd.heritability.pair_table`
    produces.  Passing ``r_sib`` adds a ``full_sib`` class at that correlation,
    which is how the pooling-bias test is built.
    """
    rng = np.random.default_rng(seed)
    rows, ys = [], {}
    spec = [("MZ", h2_true, n_mz), (dz_label, h2_true / 2, n_dz)]
    if r_sib is not None:
        spec.append(("full_sib", r_sib, n_sib))
    for pt, r, n in spec:
        # bivariate normal with correlation r
        a = rng.normal(size=n)
        b = r * a + np.sqrt(max(1 - r**2, 0)) * rng.normal(size=n)
        for i in range(n):
            fam = f"{pt}_{i}"
            sa, sb = f"{fam}_0", f"{fam}_1"
            ys[sa], ys[sb] = a[i], b[i]
            rows.append({"a": sa, "b": sb, "family_id": fam,
                         "pihat": np.nan, "pair_type": pt})
    return pd.Series(ys), pd.DataFrame(rows)


@pytest.mark.parametrize("h2_true", [0.2, 0.5, 0.8])
def test_recovers_known_h2(h2_true):
    """Unbiasedness check.

    n=5000 pairs per class is far more than ABCD has; it is chosen so the test
    is not flaky.  Falconer's sampling SE is large: at n=400 pairs sd(h2) ~ 0.13
    and at ABCD's 260 MZ pairs it is wider still, which is why the real
    estimates carry bootstrap intervals spanning ~0.5.
    """
    y, pairs = _synthetic(h2_true)
    got = falconer(y, pairs)["h2"]
    assert abs(got - h2_true) < 0.12, f"expected ~{h2_true}, got {got:.3f}"


def test_reports_both_correlations():
    """h2 alone cannot reveal a bad input -- the pair correlations must surface."""
    y, pairs = _synthetic(0.5)
    res = falconer(y, pairs)
    assert {"r_MZ", "r_DZ", "n_MZ", "n_DZ"} <= res.keys()
    assert res["r_MZ"] > res["r_DZ"]


def test_negative_dz_correlation_is_not_silently_swallowed():
    """The family-random-effect artefact: siblings anti-correlate.

    Falconer differences the correlations, so this returns a *plausible* h2.
    The guard is that r_DZ is exposed so a caller can see it is impossible.
    """
    y, pairs = _synthetic(0.5)
    # force DZ pairs to anti-correlate, as within-family deviations do
    for _, row in pairs[pairs.pair_type == "DZ_twin"].iterrows():
        y[row.b] = -y[row.a]
    res = falconer(y, pairs)
    assert res["r_DZ"] < 0
    assert res["h2"] > 0.9  # plausible-looking, from invalid input


def test_refuses_when_ids_do_not_join():
    y, pairs = _synthetic(0.5)
    y.index = ["nonsense_" + i for i in y.index]
    with pytest.raises(ValueError, match="complete"):
        falconer(y, pairs)


def test_pooling_siblings_into_dz_biases_h2_upward():
    """The 5.1-era pooled DZ class is not a harmless convenience.

    Non-twin full siblings correlate less than DZ twins at the same nominal 0.5
    sharing.  Pooling them lowers r_DZ, which Falconer *subtracts*, so h2 comes
    out too high.  Built here with a true h2 of 0.5 (r_MZ 0.50, r_DZtwin 0.25)
    and siblings at 0.10, the pooled estimate should overshoot while the
    twins-only estimate recovers the truth.
    """
    y, pairs = _synthetic(0.5, r_sib=0.10)
    twins = falconer(y, pairs, dz_class="twins_only")
    pooled = falconer(y, pairs, dz_class="pooled")
    assert abs(twins["h2"] - 0.5) < 0.05, twins["h2"]
    assert pooled["h2"] > twins["h2"] + 0.1, (pooled["h2"], twins["h2"])
    assert twins["dz_class"] == "DZ_twin" and pooled["dz_class"] == "DZ_or_sib"
    # the excluded class must be reported so the choice is auditable
    assert twins["r_sib_excluded"] < twins["r_DZ"]


def test_twins_only_is_the_default():
    """Default must be the unbiased class, not the historically-reported one."""
    y, pairs = _synthetic(0.5, r_sib=0.10)
    assert falconer(y, pairs)["dz_class"] == "DZ_twin"


def test_pooled_source_ignores_dz_class_choice():
    """On a 5.1-style table the classes are already pooled; don't pretend otherwise."""
    y, pairs = _synthetic(0.5, dz_label="DZ_or_sib")
    for dz in ("twins_only", "pooled"):
        assert falconer(y, pairs, dz_class=dz)["dz_class"] == "DZ_or_sib"


def test_rejects_unknown_dz_class():
    y, pairs = _synthetic(0.5)
    with pytest.raises(ValueError, match="dz_class"):
        falconer(y, pairs, dz_class="siblings_only")


def test_id_normalisation_reconciles_release_spellings():
    s = pd.Series(["sub-0A4P0LWM", "NDAR_INV0A4P0LWM", "NDARINV0A4P0LWM"])
    out = _normalise_ids(s)
    assert out.iloc[1] == out.iloc[2] == "NDARINV0A4P0LWM"
    assert out.iloc[0] == "0A4P0LWM"


def test_sharing_gap_is_one_half():
    """The factor of 2 is 1/0.5; if this changes the estimator is not Falconer."""
    assert _SHARING_GAP == 0.5
