"""Tests for spatial null models.

The assertions here are about *statistical validity*, not just plumbing.  The
whole point of this module is that a naive permutation null gives wrong
answers on cortical maps; these tests lock in the properties that make the
spin null different.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd import spatial as sp  # noqa: E402


@pytest.fixture(scope="module")
def geom():
    return sp.load_dk_geometry()


def test_geometry_shape(geom):
    assert len(geom.labels) == 68
    assert geom.coords.shape == (68, 3)
    assert set(geom.hemi) == {"lh", "rh"}


def test_centroids_are_on_a_sphere(geom):
    """Radii must be equal -- rotation is only valid on spherical coords."""
    r = np.linalg.norm(geom.coords, axis=1)
    assert r.std() / r.mean() < 1e-6


def test_distance_matrix_properties(geom):
    d = geom.distance_matrix()
    assert d.shape == (68, 68)
    assert np.allclose(np.diag(d), 0)
    assert np.allclose(d, d.T)
    assert (d >= 0).all()


def test_spin_preserves_permutation_structure(geom):
    """Every rotation must be a permutation: each region used exactly once.

    A nearest-centroid spin can in principle map two regions to the same
    target; the implementation must resolve that, or the null is not a
    permutation of the data and its variance is wrong.
    """
    S = sp.spin_indices(geom, n_perm=50, seed=0)
    assert S.shape == (50, 68)
    for row in S:
        assert len(np.unique(row)) == 68, "spin is not a permutation"


def test_spin_keeps_hemispheres_separate(geom):
    """Left regions must rotate to left targets, right to right."""
    S = sp.spin_indices(geom, n_perm=20, seed=1)
    hemi = np.asarray(geom.hemi)
    for row in S:
        assert (hemi[row] == hemi).all()


def test_spin_inflates_null_for_smooth_maps_but_not_noise(geom):
    """The discriminating property of a valid spatial null.

    On a spatially smooth map the spin null must be WIDER than a naive
    permutation null (that width is the correction).  On white noise the two
    must agree, because there is no spatial structure to preserve.
    """
    d = geom.distance_matrix()
    med = np.median(d[np.triu_indices(len(d), 1)])
    K = np.exp(-(d ** 2) / (2 * (0.35 * med) ** 2))
    np.fill_diagonal(K, 0)
    rng = np.random.default_rng(0)
    smooth = pd.Series(K @ rng.normal(size=68), index=list(geom.labels))
    noise = pd.Series(rng.normal(size=68), index=list(geom.labels))

    for x, expect_inflation in [(smooth, True), (noise, False)]:
        t_spin = sp.spatial_corr(x, x, null="spin", geom=geom, n_perm=300)
        t_naive = sp.spatial_corr(x, x, null="naive", n_perm=300, seed=0)
        ratio = t_spin.null_sd / t_naive.null_sd
        if expect_inflation:
            assert ratio > 1.3, f"spin should inflate smooth-map null, got {ratio:.2f}"
        else:
            assert ratio < 1.2, f"spin should not inflate noise null, got {ratio:.2f}"


def test_p_value_is_never_zero(geom):
    """A permutation test cannot justify p = 0; the +1 offset enforces that."""
    x = pd.Series(np.arange(68.0), index=list(geom.labels))
    t = sp.spatial_corr(x, x, null="spin", geom=geom, n_perm=100)
    assert t.p > 0
    assert t.p == pytest.approx(1 / 101)


def test_spin_requires_geometry():
    """Absent geometry must raise, never silently fall back to a bad null."""
    x = pd.Series(np.arange(10.0))
    with pytest.raises((ValueError, sp.SpatialResourceMissing, TypeError)):
        sp.spatial_corr(x, x, null="spin", geom=None, n_perm=10)


def test_hemisphere_symmetry_detects_real_asymmetry(geom):
    """Symmetry measurement drives automatic pairing choice, so test it."""
    lab = list(geom.labels)
    rng = np.random.default_rng(2)
    lh_vals = rng.normal(size=34)
    sym = pd.Series(
        [lh_vals[[l[3:] for l in lab if l.startswith("lh_")].index(l[3:])] for l in lab],
        index=lab,
    )
    assert sp.hemisphere_symmetry(sym, geom) > 0.99
    asym = pd.Series(rng.normal(size=68), index=lab)
    assert sp.hemisphere_symmetry(asym, geom) < 0.6


def test_unpaired_spin_is_available_and_differs(geom):
    S_p = sp.spin_indices(geom, n_perm=20, seed=0, hemi_paired=True)
    S_u = sp.spin_indices(geom, n_perm=20, seed=0, hemi_paired=False)
    assert not np.array_equal(S_p, S_u)


def test_variogram_surrogates_preserve_variogram(geom):
    """Surrogates must match the map's spatial autocorrelation, not just its
    marginal distribution -- otherwise they are just a shuffle."""
    d = geom.distance_matrix()
    med = np.median(d[np.triu_indices(len(d), 1)])
    K = np.exp(-(d ** 2) / (2 * 0.35 * med) ** 2)
    np.fill_diagonal(K, 0)
    rng = np.random.default_rng(0)
    x = K @ rng.normal(size=68)
    surr = sp.variogram_surrogates(x, d, n_perm=30, seed=0)
    assert surr.shape == (30, 68)
    # marginal distribution preserved exactly (rank-remapping)
    assert np.allclose(np.sort(surr[0]), np.sort(x))
    # neighbour correlation of surrogates should be closer to the map's than
    # a plain shuffle is
    iu = np.triu_indices(68, 1)
    near = d[iu] < np.percentile(d[iu], 10)

    def near_sim(v):
        return np.corrcoef(v[iu[0]][near], v[iu[1]][near])[0, 1]

    obs = near_sim(x)
    surr_sim = np.mean([near_sim(s) for s in surr])
    shuf = np.mean([near_sim(rng.permutation(x)) for _ in range(30)])
    assert abs(surr_sim - obs) < abs(shuf - obs)
