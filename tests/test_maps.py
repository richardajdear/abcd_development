"""Guards on the adjusted-vs-total slope distinction.

The bug this pins is interpretive, not mechanical, which is why it needs a
test: with a global covariate in the model, a region's ``age_c`` coefficient
is a *relative* rate, and reading it as a thinning rate inverts the sign for
roughly a third of the cortex.  A map that silently mixes the two, or a
reconstruction that drifts from the identity below, will produce a real-looking
map with the wrong meaning -- and downstream that changes whether the AHBA C3
association appears at all (-0.57 vs -0.32).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from abcd.maps import global_age_slope, total_age_slope


def _fixed(terms: dict[str, dict[str, float]]) -> pd.DataFrame:
    rows = [{"label": lab, "term": t, "estimate": v}
            for t, per_label in terms.items() for lab, v in per_label.items()]
    return pd.DataFrame(rows)


def test_global_age_slope_recovers_known_slope():
    """global_within is subject-demeaned, so OLS on age_c is the within slope."""
    rng = np.random.default_rng(0)
    n = 4000
    age = rng.normal(0, 2, n)
    truth = -0.017
    mt = pd.DataFrame({
        "label": "lh_bankssts",
        "age_c": age,
        "global_within": truth * age + rng.normal(0, 0.02, n),
    })
    assert global_age_slope(mt) == pytest.approx(truth, abs=1e-3)


def test_global_age_slope_uses_one_region_not_all():
    """Every region repeats the same subject-level global columns.

    Pooling all regions would multiply the row count without changing the
    estimate -- but would silently succeed on a table where regions carry
    *different* global values, which should never happen.  Pinning the
    one-region behaviour keeps that assumption explicit.
    """
    rng = np.random.default_rng(1)
    age = rng.normal(0, 2, 500)
    gw = -0.02 * age
    mt = pd.concat([
        pd.DataFrame({"label": lab, "age_c": age, "global_within": gw})
        for lab in ("lh_bankssts", "lh_cuneus", "rh_insula")
    ], ignore_index=True)
    assert global_age_slope(mt) == pytest.approx(-0.02, abs=1e-6)


def test_global_age_slope_rejects_run_without_covariate():
    mt = pd.DataFrame({"label": "lh_bankssts", "age_c": [0.0, 1.0]})
    with pytest.raises(KeyError, match="no 'global_within'"):
        global_age_slope(mt)


def test_total_slope_is_chain_rule():
    fixed = _fixed({
        "age_c": {"lh_a": 0.010, "lh_b": -0.005},
        "global_within": {"lh_a": 1.2, "lh_b": 0.4},
    })
    got = total_age_slope(fixed, -0.0166)
    assert got["lh_a"] == pytest.approx(0.010 + 1.2 * -0.0166)
    assert got["lh_b"] == pytest.approx(-0.005 + 0.4 * -0.0166)


def test_positive_adjusted_coefficient_can_be_negative_total():
    """The whole point: 'positive age_c' does not mean 'thickens'.

    A region thinning more slowly than the cortex-wide average gets a positive
    adjusted coefficient while its absolute rate stays negative.  This is the
    case that makes reading the adjusted map as a thinning map wrong.
    """
    fixed = _fixed({"age_c": {"lh_precentral": 0.0171},
                    "global_within": {"lh_precentral": 1.45}})
    total = total_age_slope(fixed, -0.0166)["lh_precentral"]
    assert total < 0, "adjusted-positive region must still thin in absolute terms"


def test_total_slope_requires_both_terms():
    with pytest.raises(KeyError, match="global_within"):
        total_age_slope(_fixed({"age_c": {"lh_a": 0.01}}), -0.0166)


def _write_run(tmp_path, *, with_global: bool):
    """Minimal on-disk run directory: the four fit tables plus a model table.

    ``global_within`` is written into the model table in BOTH cases, because
    ``assemble`` does exactly that regardless of the fit spec -- which is the
    trap ``regional_maps`` has to avoid keying on.
    """
    labels = ["lh_a", "lh_b"]
    rng = np.random.default_rng(1)
    n = 500
    age = rng.normal(0, 2, n)
    terms = {"age_c": {"lh_a": -0.020, "lh_b": -0.010}}
    if with_global:
        terms["global_within"] = {"lh_a": 0.8, "lh_b": 1.2}
    fixed = _fixed(terms)
    fixed["statistic"] = -5.0

    fits = tmp_path / "fits"
    fits.mkdir(parents=True)
    fixed.to_parquet(fits / "fixed.parquet")
    pd.DataFrame({"label": labels, "grp": "subject", "var1": "age_c",
                  "var2": None, "sdcor": [0.006, 0.007]}
                 ).to_parquet(fits / "varcomp.parquet")
    for name in ("blups", "diagnostics"):
        pd.DataFrame({"label": labels}).to_parquet(fits / f"{name}.parquet")

    pd.DataFrame({
        "label": np.repeat(labels, n),
        "age_c": np.tile(age, len(labels)),
        "global_within": np.tile(-0.0166 * age, len(labels)),
    }).to_parquet(tmp_path / "model_table.parquet")
    return tmp_path


def test_regional_maps_noglobal_run_is_already_total(tmp_path):
    """Without the covariate the age coefficient IS the absolute rate.

    The two slope columns must be identical rather than one being a
    reconstruction, and ``attrs`` must say so -- otherwise a figure captioned
    'total' silently plots something else.
    """
    from abcd.maps import regional_maps

    out = regional_maps(_write_run(tmp_path, with_global=False))
    assert out.attrs["global_covariate"] is False
    pd.testing.assert_series_equal(
        out["slope_total"], out["slope_adjusted"], check_names=False)
    assert (out["slope_total"] < 0).all()


def test_regional_maps_global_run_reconstructs_and_differs(tmp_path):
    """With the covariate, total != adjusted, and the difference has the sign
    the chain rule dictates: a positive global loading plus cortex-wide
    thinning makes the absolute rate MORE negative than the coefficient."""
    from abcd.maps import regional_maps

    out = regional_maps(_write_run(tmp_path, with_global=True))
    assert out.attrs["global_covariate"] is True
    assert (out["slope_total"] < out["slope_adjusted"]).all()
    assert out.loc["lh_a", "slope_total"] == pytest.approx(
        -0.020 + 0.8 * -0.0166, abs=1e-4)
