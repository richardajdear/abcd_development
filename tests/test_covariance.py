"""Guards on the structural-covariance PCs.

The claim these tests pin is that "PCs of the structural covariance matrix" and
"region loadings of a PCA on the slope matrix" are the same vectors.  The 5.1
thesis presented them as two analyses; if a future change makes them genuinely
diverge, that is a real finding and this test should fail loudly rather than
have the report quietly keep asserting the identity.

They also pin the variance-explained convention.  PCA on a correlation matrix
reports a share of *squared* eigenvalues, which inflated the thesis's headline
from 7% to 21%.  That inflation is exactly 3x here, so it is not a rounding
detail.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from abcd.covariance import (
    PCResult,
    sc_matrix,
    sc_pcs,
    slope_pcs,
    subject_scores,
)


@pytest.fixture
def slopes() -> pd.DataFrame:
    """Synthetic subject x region slopes with two planted spatial factors."""
    rng = np.random.default_rng(7)
    n_sub = 800
    regions = ([f"lh_{r}" for r in ("superiorfrontal", "rostralmiddlefrontal",
                                    "frontalpole", "medialorbitofrontal")]
               + [f"lh_{r}" for r in ("lateraloccipital", "lingual", "cuneus",
                                      "pericalcarine")]
               + [f"lh_{r}" for r in ("insula", "superiortemporal",
                                      "middletemporal", "bankssts")])
    ant = np.array([1.0]*4 + [0.0]*4 + [0.0]*4)
    post = np.array([0.0]*4 + [1.0]*4 + [0.0]*4)
    f1 = rng.normal(0, 1, n_sub)
    f2 = rng.normal(0, 1, n_sub)
    X = (np.outer(f1, ant - post) * 0.9
         + np.outer(f2, np.r_[np.zeros(8), np.ones(4)]) * 0.7
         + rng.normal(0, 0.35, (n_sub, len(regions))))
    return pd.DataFrame(X, columns=regions,
                        index=[f"s{i:04d}" for i in range(n_sub)])


def test_sc_pcs_equal_slope_pcs(slopes):
    """The identity the report rests on: same eigenvectors, |r| ~ 1."""
    sc = sc_pcs(sc_matrix(slopes), n_components=3)
    dr = slope_pcs(slopes, n_components=3)
    for k in range(3):
        r = np.corrcoef(sc.loadings.iloc[:, k], dr.loadings.iloc[:, k])[0, 1]
        assert abs(r) > 0.99, f"component {k+1} diverged: r={r:.4f}"


def test_variance_explained_is_not_the_squared_convention(slopes):
    """L/sum(L), not L^2/sum(L^2) -- the thesis reported the latter."""
    sc = sc_pcs(sc_matrix(slopes), n_components=3)
    dr = slope_pcs(slopes, n_components=3)
    np.testing.assert_allclose(sc.var_explained, dr.var_explained, atol=1e-8)
    assert sc.var_explained_sq_convention[0] > sc.var_explained[0], (
        "the squared convention must inflate PC1, otherwise this test is "
        "not checking what it claims"
    )


def test_variance_explained_sums_below_one_and_descends(slopes):
    sc = sc_pcs(sc_matrix(slopes), n_components=3)
    assert sc.var_explained.sum() <= 1.0 + 1e-9
    assert np.all(np.diff(sc.var_explained) <= 1e-12)


def test_orientation_is_deterministic_under_sign_flip(slopes):
    """A flipped input must not flip the reported map.

    Eigenvector signs are arbitrary, so without an orientation rule a refit on
    new data can invert a component and silently invert every correlation
    reported against it.
    """
    a = sc_pcs(sc_matrix(slopes), n_components=3).loadings
    b = sc_pcs(sc_matrix(-slopes), n_components=3).loadings
    for k in range(3):
        r = np.corrcoef(a.iloc[:, k], b.iloc[:, k])[0, 1]
        assert r > 0.99, f"component {k+1} flipped under input negation (r={r:.3f})"


def test_orientation_puts_frontal_positive(slopes):
    """PC1 here is planted as frontal-vs-occipital; frontal must come out +."""
    load = sc_pcs(sc_matrix(slopes), n_components=3).loadings
    front = load.loc[["lh_superiorfrontal", "lh_frontalpole"], "SC1"].mean()
    occ = load.loc[["lh_lateraloccipital", "lh_cuneus"], "SC1"].mean()
    assert front > occ


def test_loadings_are_z_scored(slopes):
    load = sc_pcs(sc_matrix(slopes), n_components=3).loadings
    np.testing.assert_allclose(load.mean().to_numpy(), 0, atol=1e-10)
    np.testing.assert_allclose(load.std(ddof=0).to_numpy(), 1, atol=1e-10)


def test_sc_matrix_rejects_thin_overlap():
    """Better to raise than to correlate two regions on a handful of subjects."""
    W = pd.DataFrame({"lh_a": [1.0, 2.0, 3.0, 4.0],
                      "lh_b": [1.0, np.nan, np.nan, 2.0]})
    with pytest.raises(ValueError, match="shared"):
        sc_matrix(W, min_overlap=100)


def test_subject_scores_shape_and_centring(slopes):
    load = sc_pcs(sc_matrix(slopes), n_components=3).loadings
    S = subject_scores(slopes, load)
    assert S.shape == (len(slopes), 3)
    assert list(S.index) == list(slopes.index)
    np.testing.assert_allclose(S.mean().to_numpy(), 0, atol=1e-8)


def test_subject_scores_recover_planted_factor(slopes):
    """SC1 scores must track the planted frontal-occipital factor."""
    load = sc_pcs(sc_matrix(slopes), n_components=3).loadings
    S = subject_scores(slopes, load)
    ant = ["lh_superiorfrontal", "lh_rostralmiddlefrontal",
           "lh_frontalpole", "lh_medialorbitofrontal"]
    post = ["lh_lateraloccipital", "lh_lingual", "lh_cuneus", "lh_pericalcarine"]
    contrast = slopes[ant].mean(axis=1) - slopes[post].mean(axis=1)
    assert abs(np.corrcoef(S["SC1"], contrast)[0, 1]) > 0.9
