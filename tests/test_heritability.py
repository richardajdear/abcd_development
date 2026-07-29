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


def _synthetic(h2_true, n_mz=5000, n_dz=5000, seed=0):
    """Build pair data with a known h2, then flatten to (y, pairs)."""
    rng = np.random.default_rng(seed)
    rows, ys = [], {}
    for pt, r in [("MZ", h2_true), ("DZ_or_sib", h2_true / 2)]:
        n = n_mz if pt == "MZ" else n_dz
        # bivariate normal with correlation r
        a = rng.normal(size=n)
        b = r * a + np.sqrt(max(1 - r**2, 0)) * rng.normal(size=n)
        for i in range(n):
            fam = f"{pt}_{i}"
            for j, v in enumerate((a[i], b[i])):
                sid = f"{fam}_{j}"
                ys[sid] = v
                rows.append({"subject": sid, "family_id": fam, "pair_type": pt})
    pairs = pd.DataFrame(rows).set_index("subject")
    return pd.Series(ys), pairs


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
    dz = pairs[pairs.pair_type == "DZ_or_sib"]
    for fam, grp in dz.groupby("family_id"):
        ids = list(grp.index)
        if len(ids) == 2:
            y[ids[1]] = -y[ids[0]]
    res = falconer(y, pairs)
    assert res["r_DZ"] < 0
    assert res["h2"] > 0.9  # plausible-looking, from invalid input


def test_refuses_when_ids_do_not_join():
    y, pairs = _synthetic(0.5)
    y.index = ["nonsense_" + i for i in y.index]
    with pytest.raises(ValueError, match="complete"):
        falconer(y, pairs)


def test_id_normalisation_reconciles_release_spellings():
    s = pd.Series(["sub-0A4P0LWM", "NDAR_INV0A4P0LWM", "NDARINV0A4P0LWM"])
    out = _normalise_ids(s)
    assert out.iloc[1] == out.iloc[2] == "NDARINV0A4P0LWM"
    assert out.iloc[0] == "0A4P0LWM"


def test_sharing_gap_is_one_half():
    """The factor of 2 is 1/0.5; if this changes the estimator is not Falconer."""
    assert _SHARING_GAP == 0.5
