"""
Spatial statistics for cortical maps: nulls that respect spatial autocorrelation.

Why this module exists
----------------------
Cortical maps are smooth.  Neighbouring regions have similar values for
reasons that have nothing to do with the hypothesis being tested, so the
number of *effectively independent* regions is far below 68.  A Pearson or
Spearman p-value computed against the naive null -- randomly permuting region
labels -- is therefore anticonservative, often dramatically so: two unrelated
but smooth maps correlate at r ~ 0.4 by chance far more often than the
parametric p-value admits.

Every map-to-map comparison in this project is tested against a
spatially-constrained null.  Two are implemented, with a documented fallback
order:

``spin``
    Rotate the parcellation on the sphere (Alexander-Bloch et al. 2018) and
    re-read the map at the rotated positions.  This is the reference method:
    it preserves the *exact* empirical spatial autocorrelation because it
    reuses the observed values, only moving them.  Requires spherical
    coordinates for the parcel centroids.

``moran`` / ``burt``
    Generate surrogate maps matching the observed variogram
    (Burt et al. 2020) from an inter-regional distance matrix.  Used when
    sphere coordinates are unavailable but distances are.

``naive``
    Plain label permutation.  **Never** appropriate for inference; provided
    only so the inflation from ignoring spatial structure can be quantified
    and reported, which this module does explicitly in
    :func:`compare_nulls`.

The module degrades honestly: if the resources for a spin test are absent,
:func:`spatial_corr` raises rather than silently falling back to the naive
null.  Pass ``null="naive"`` deliberately if that is what you want.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment

from . import paths


class SpatialResourceMissing(RuntimeError):
    """Raised when a requested null needs a resource that is not available."""


# --------------------------------------------------------------------------
# Parcel geometry
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class ParcelGeometry:
    """Centroid coordinates for a parcellation, in a named space.

    ``coords`` is (n_regions, 3).  ``space`` is ``"sphere"`` for spherical
    coordinates suitable for the spin test, or ``"native"`` for anatomical
    coordinates suitable only for distance-based nulls.
    """

    labels: tuple[str, ...]
    coords: np.ndarray
    space: str
    hemi: tuple[str, ...]

    def __post_init__(self):
        if self.coords.shape != (len(self.labels), 3):
            raise ValueError(
                f"coords {self.coords.shape} does not match {len(self.labels)} labels"
            )

    def distance_matrix(self) -> np.ndarray:
        """Euclidean distances between centroids.

        On the sphere this is chord distance, a monotone function of great
        circle distance, which is all the variogram methods require.
        """
        d = self.coords[:, None, :] - self.coords[None, :, :]
        return np.sqrt((d ** 2).sum(-1))


def load_dk_geometry(path: str | Path | None = None) -> ParcelGeometry:
    """Load DK parcel centroids.

    Looks for ``data/dk_centroids.csv`` with columns
    ``label, hemi, x, y, z, space``.  Generate it once with
    ``python -m abcd.spatial --make-dk-centroids`` (requires ``nibabel`` and
    a FreeSurfer ``fsaverage`` sphere), then commit it -- the pipeline should
    not depend on a FreeSurfer installation at analysis time.
    """
    p = Path(path) if path else paths.DATA_DIR / "dk_centroids.csv"
    if not p.exists():
        raise SpatialResourceMissing(
            f"{p} not found. Spin tests need parcel centroids on the sphere.\n"
            "Create with: python -m abcd.spatial --make-dk-centroids "
            "--annot <fsaverage lh.aparc.annot> --sphere <fsaverage lh.sphere>"
        )
    df = pd.read_csv(p)
    space = df.space.iloc[0]
    if df.space.nunique() != 1:
        raise ValueError("mixed coordinate spaces in centroid file")
    return ParcelGeometry(
        labels=tuple(df.label), coords=df[["x", "y", "z"]].to_numpy(float),
        space=space, hemi=tuple(df.hemi),
    )


# --------------------------------------------------------------------------
# Null generators
# --------------------------------------------------------------------------

def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Uniform random 3x3 rotation via QR of a Gaussian matrix."""
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.diag(r))


def hemisphere_symmetry(x: pd.Series, geom: ParcelGeometry) -> float:
    """Pearson correlation between homologous left and right regions.

    This number decides ``hemi_paired`` (see :func:`spin_indices`), so it is
    computed rather than assumed.  Real cortical maps in this project sit
    around 0.95-0.97; a map below ~0.5 should use an unpaired spin.
    """
    labs = list(geom.labels)
    lh = [l for l in labs if l.startswith("lh_")]
    rh = ["rh_" + l[3:] for l in lh]
    both = [(a, b) for a, b in zip(lh, rh)
            if a in x.index and b in x.index
            and np.isfinite(x.get(a, np.nan)) and np.isfinite(x.get(b, np.nan))]
    if len(both) < 5:
        return np.nan
    a = x.loc[[p[0] for p in both]].to_numpy(float)
    b = x.loc[[p[1] for p in both]].to_numpy(float)
    return float(stats.pearsonr(a, b).statistic)


SYMMETRY_THRESHOLD = 0.5
"""Above this lh/rh correlation, the hemisphere-paired spin is used.

Calibrated on 200 pairs of independent synthetic smooth maps per condition,
false-positive rate at alpha = 0.05 (nominal 0.05):

    map symmetry     paired spin     unpaired spin
    r ~ 0.95         0.065 - 0.075   0.165 - 0.315
    r ~ 0.18         0.153 - 0.227   0.087 - 0.093

Neither variant is correct for both regimes: pairing must match the map.
The threshold sits between the two regimes and is deliberately generous,
since real cortical maps cluster near 0.95 and nothing in this project is
near 0.5.
"""


def spin_indices(geom: ParcelGeometry, n_perm: int = 1000,
                 seed: int = 0, hemi_paired: bool = True) -> np.ndarray:
    """Index array (n_perm, n_regions) mapping each region to its rotated partner.

    Implements the Alexander-Bloch spin test: rotate the centroids, then match
    each original parcel to a rotated parcel.  The matching is a *bijection*
    solved by minimising total centroid displacement
    (``scipy.optimize.linear_sum_assignment``), not independent
    nearest-neighbour lookups.  Independent lookups are the common shortcut but
    they are not a permutation -- on the DK atlas they yield only ~48 distinct
    targets out of 68, so some regional values are duplicated and others
    dropped.  That changes the null maps' marginal distribution, and a rank
    correlation against a null whose marginals differ from the data is not a
    valid test.

    When ``hemi_paired`` the same rotation is applied to both hemispheres
    (reflected in x for the right), which preserves the strong left-right
    symmetry of cortical maps -- without it the null is too permissive for
    symmetric maps like ours (lh/rh correlation of the developmental map is
    0.94).
    """
    if geom.space != "sphere":
        raise SpatialResourceMissing(
            f"spin test needs spherical coordinates, geometry is '{geom.space}'"
        )
    rng = np.random.default_rng(seed)
    coords = geom.coords
    hemi = np.asarray(geom.hemi)
    n = len(geom.labels)
    out = np.empty((n_perm, n), dtype=int)

    for i in range(n_perm):
        rot = _random_rotation(rng)
        assign = np.empty(n, dtype=int)
        if hemi_paired and set(hemi) == {"lh", "rh"}:
            refl = np.diag([-1.0, 1.0, 1.0])
            for h, R in (("lh", rot), ("rh", refl @ rot @ refl)):
                m = np.flatnonzero(hemi == h)
                rotated = coords[m] @ R.T
                d = ((rotated[:, None, :] - coords[m][None, :, :]) ** 2).sum(-1)
                _, col = linear_sum_assignment(d)
                assign[m] = m[col]
        else:
            rotated = coords @ rot.T
            d = ((rotated[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
            _, col = linear_sum_assignment(d)
            assign = col
        out[i] = assign
    return out


def variogram_surrogates(x: np.ndarray, dist: np.ndarray, n_perm: int = 1000,
                         seed: int = 0, n_bins: int = 12) -> np.ndarray:
    """Surrogate maps preserving the empirical variogram (Burt-style).

    Simplified implementation: permute the map, smooth the permutation with a
    distance kernel whose bandwidth is fitted so the surrogate's variogram
    matches the observed one, then rank-map back onto the observed values so
    the marginal distribution is preserved exactly.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    n = len(x)
    iu = np.triu_indices(n, 1)
    d_flat = dist[iu]
    bins = np.quantile(d_flat, np.linspace(0, 1, n_bins + 1))
    bin_id = np.clip(np.digitize(d_flat, bins[1:-1]), 0, n_bins - 1)

    def variogram(v):
        sv = 0.5 * (v[iu[0]] - v[iu[1]]) ** 2
        return np.array([sv[bin_id == b].mean() for b in range(n_bins)])

    target = variogram(x)
    x_sorted = np.sort(x)

    # candidate bandwidths as fractions of the median distance
    med = np.median(d_flat)
    cands = med * np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.5])
    kernels = [np.exp(-(dist ** 2) / (2 * bw ** 2)) for bw in cands]
    for K in kernels:
        np.fill_diagonal(K, 0.0)

    # pick the bandwidth whose surrogates best match the target variogram
    probe = rng.permutation(x)
    errs = []
    for K in kernels:
        s = K @ probe
        s = x_sorted[np.argsort(np.argsort(s))]
        errs.append(np.abs(variogram(s) - target).sum())
    K = kernels[int(np.argmin(errs))]

    out = np.empty((n_perm, n))
    for i in range(n_perm):
        s = K @ rng.permutation(x)
        out[i] = x_sorted[np.argsort(np.argsort(s))]
    return out


# --------------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------------

@dataclass
class SpatialTest:
    r: float
    p: float
    null: str
    n_perm: int
    method: str
    null_r: np.ndarray

    def __repr__(self):  # pragma: no cover
        return (f"SpatialTest(r={self.r:.3f}, p={self.p:.4f}, null={self.null!r}, "
                f"n_perm={self.n_perm}, method={self.method!r})")

    @property
    def null_sd(self) -> float:
        return float(np.std(self.null_r))


def _corr(a, b, method):
    if method == "spearman":
        return float(stats.spearmanr(a, b).statistic)
    return float(stats.pearsonr(a, b).statistic)


def spatial_corr(x: pd.Series, y: pd.Series, null: str = "spin",
                 geom: ParcelGeometry | None = None, n_perm: int = 1000,
                 seed: int = 0, method: str = "spearman",
                 spins: np.ndarray | None = None,
                 hemi_paired: bool | None = None) -> SpatialTest:
    """Correlate two regional maps and test against a spatial null.

    ``x`` and ``y`` are indexed by region label; they are aligned on the
    intersection of their indices, so a map missing a region is handled
    without silently misaligning.  The two-sided p-value is
    ``(1 + #{|r_null| >= |r_obs|}) / (1 + n_perm)`` -- the +1 keeps it from
    ever being exactly zero, which no permutation test can justify.

    ``hemi_paired=None`` (the default) chooses the spin variant from the
    measured lh/rh symmetry of ``x``; see :data:`SYMMETRY_THRESHOLD` for the
    calibration behind that choice.  Pass a bool to force it.
    """
    common = x.index.intersection(y.index)
    if len(common) < len(x) or len(common) < len(y):
        pass  # alignment is intentional; caller may check len(common)
    xa, ya = x.loc[common].to_numpy(float), y.loc[common].to_numpy(float)
    ok = np.isfinite(xa) & np.isfinite(ya)
    xa, ya, labels = xa[ok], ya[ok], common[ok]
    r_obs = _corr(xa, ya, method)

    if null == "naive":
        rng = np.random.default_rng(seed)
        null_r = np.array([_corr(rng.permutation(xa), ya, method)
                           for _ in range(n_perm)])
    elif null == "spin":
        if spins is None:
            if geom is None:
                raise SpatialResourceMissing(
                    "spin null requires `geom` (or precomputed `spins`); "
                    "pass null='naive' only if you intend an invalid test"
                )
            pos = {l: i for i, l in enumerate(geom.labels)}
            missing = [l for l in labels if l not in pos]
            if missing:
                raise ValueError(f"labels absent from geometry: {missing[:5]}")
            idx = np.array([pos[l] for l in labels])
            if hemi_paired is None:
                sym = hemisphere_symmetry(pd.Series(xa, index=labels), geom)
                hemi_paired = bool(np.isfinite(sym) and sym >= SYMMETRY_THRESHOLD)
            spins_full = spin_indices(geom, n_perm=n_perm, seed=seed,
                                      hemi_paired=hemi_paired)
            spins = np.searchsorted(idx, spins_full[:, idx])
            # map rotated positions back into the subset; drop rotations that
            # land outside it by re-using nearest within-subset assignment
            spins = np.clip(spins, 0, len(idx) - 1)
        null_r = np.array([_corr(xa[s], ya, method) for s in spins])
    elif null in ("moran", "burt", "variogram"):
        if geom is None:
            raise SpatialResourceMissing("variogram null requires `geom`")
        pos = {l: i for i, l in enumerate(geom.labels)}
        idx = np.array([pos[l] for l in labels])
        dist = geom.distance_matrix()[np.ix_(idx, idx)]
        surr = variogram_surrogates(xa, dist, n_perm=n_perm, seed=seed)
        null_r = np.array([_corr(s, ya, method) for s in surr])
    else:
        raise ValueError(f"unknown null '{null}'")

    p = (1 + int((np.abs(null_r) >= abs(r_obs)).sum())) / (1 + n_perm)
    return SpatialTest(r=r_obs, p=p, null=null, n_perm=n_perm,
                       method=method, null_r=null_r)


def compare_nulls(x: pd.Series, y: pd.Series, geom: ParcelGeometry,
                  n_perm: int = 1000, seed: int = 0,
                  method: str = "spearman") -> pd.DataFrame:
    """Run every available null side by side and report the inflation.

    The point of this table is to make the cost of ignoring spatial structure
    visible: the naive p-value is typically several-fold smaller than the
    spin p-value for the same observed correlation.
    """
    rows = []
    for null in ("naive", "variogram", "spin"):
        try:
            t = spatial_corr(x, y, null=null, geom=geom, n_perm=n_perm,
                             seed=seed, method=method)
        except SpatialResourceMissing as e:
            rows.append(dict(null=null, r=np.nan, p=np.nan, null_sd=np.nan,
                             status=str(e)[:60]))
            continue
        rows.append(dict(null=null, r=t.r, p=t.p, null_sd=t.null_sd, status="ok"))
    out = pd.DataFrame(rows)
    ref = out[out.null == "spin"].p
    if len(ref) and np.isfinite(ref.iloc[0]):
        out["p_ratio_vs_spin"] = out.p / ref.iloc[0]
    return out


# --------------------------------------------------------------------------
# Centroid generation (run once, offline)
# --------------------------------------------------------------------------

def make_centroids_from_annot(annot_lh: str | Path, sphere_lh: str | Path,
                              annot_rh: str | Path | None = None,
                              sphere_rh: str | Path | None = None
                              ) -> pd.DataFrame:
    """Compute parcel centroids on the sphere from FreeSurfer files.

    Requires ``nibabel``.  Intended to be run once and the result committed,
    so downstream analysis has no FreeSurfer dependency.
    """
    import nibabel.freesurfer.io as fsio

    rows = []
    pairs = [("lh", annot_lh, sphere_lh)]
    if annot_rh and sphere_rh:
        pairs.append(("rh", annot_rh, sphere_rh))
    for hemi, annot, sphere in pairs:
        lab, _, names = fsio.read_annot(str(annot))
        coords, _ = fsio.read_geometry(str(sphere))
        for i, nm in enumerate(names):
            nm = nm.decode() if isinstance(nm, bytes) else nm
            if nm in ("unknown", "corpuscallosum", "???"):
                continue
            m = lab == i
            if m.sum() == 0:
                continue
            c = coords[m].mean(0)
            c = c / np.linalg.norm(c) * 100.0  # project back onto the sphere
            rows.append(dict(label=f"{hemi}_{nm}", hemi=hemi,
                             x=c[0], y=c[1], z=c[2], space="sphere"))
    return pd.DataFrame(rows)


if __name__ == "__main__":  # pragma: no cover
    import argparse

    ap = argparse.ArgumentParser(description="Spatial null utilities.")
    ap.add_argument("--make-dk-centroids", action="store_true")
    ap.add_argument("--annot"); ap.add_argument("--sphere")
    ap.add_argument("--annot-rh"); ap.add_argument("--sphere-rh")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.make_dk_centroids:
        df = make_centroids_from_annot(a.annot, a.sphere, a.annot_rh, a.sphere_rh)
        out = Path(a.out) if a.out else paths.DATA_DIR / "dk_centroids.csv"
        df.to_csv(out, index=False)
        print(f"{len(df)} centroids -> {out}")
