"""MOSTest-style multivariate omnibus test on per-measure GWAS z-scores.

Reimplements the statistic of van der Meer et al. 2020 (Nat Commun 11:3512,
doi 10.1038/s41467-020-17368-1) on summary statistics, so the univariate scans
can come from a mixed-model engine that keeps relatives instead of MOSTest's
own OLS on unrelated individuals.  Step 16 uses REGENIE for this, NOT the
GENESIS scan of steps 3-4: GENESIS fits one null model per phenotype and one
task per phenotype x chromosome (7 runs x 68 measures x 22 = ~10k tasks),
whereas REGENIE scans all measures of a run in one pass.  04_engine_check.py
verifies REGENIE against the existing GENESIS scan on global_slope_1lmm.

The recipe, as in MOSTest:
  1. One univariate GWAS per measure on rank-inverse-normalised phenotypes,
     run twice: on the real data, and once with the genotype-to-child link
     permuted (here: phenotype+covariate rows shuffled jointly, which is
     equivalent and preserves the between-measure correlation).
  2. R = correlation of the PERMUTED z-scores across measures (the null
     correlation of z; for fully overlapping samples it equals the residual
     phenotypic correlation).
  3. MOSTest statistic per SNP: T = z' R_reg^-1 z, with R regularised by
     flooring its eigenvalues (R is near-singular when measures are many and
     correlated).
  4. Null of T: a gamma distribution fitted to T on the permuted scan, which
     calibrates the tail empirically; p = gamma survival function.
  5. MinP (the comparator MOSTest reports): min over measures of the
     univariate two-sided p, calibrated on the permuted minP through an
     effective number of tests (Sidak form, tail-matched; see fit_minp_null).

Every function here is pure numpy/scipy and unit-tested in
tests/test_mostest_core.py against synthetic z ~ MVN(0, R).
"""
from __future__ import annotations

import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Null correlation, streamed over chromosome chunks
# ---------------------------------------------------------------------------
class CorrAccumulator:
    """Accumulate the correlation matrix of z-score columns chunk by chunk.

    Memory is O(k^2) regardless of the number of SNPs, so the permuted scan of
    ~7 M SNPs never has to sit in memory at once.
    """

    def __init__(self, k: int):
        self.k = k
        self.n = 0
        self.s = np.zeros(k)
        self.ss = np.zeros((k, k))

    def add(self, z: np.ndarray) -> None:
        z = np.asarray(z, dtype=np.float64)
        ok = np.isfinite(z).all(axis=1)
        z = z[ok]
        self.n += z.shape[0]
        self.s += z.sum(axis=0)
        self.ss += z.T @ z

    def corr(self) -> np.ndarray:
        if self.n < 2:
            raise ValueError("need at least two complete rows")
        mu = self.s / self.n
        cov = self.ss / self.n - np.outer(mu, mu)
        d = np.sqrt(np.diag(cov))
        return cov / np.outer(d, d)


# ---------------------------------------------------------------------------
# Regularised inverse
# ---------------------------------------------------------------------------
def regularised_inverse(R: np.ndarray, floor_frac: float = 1e-4):
    """Inverse of R with eigenvalues floored at floor_frac * max eigenvalue.

    Returns (Rinv, info) where info records the eigenvalue spectrum summary so
    the run log shows whether the floor bit (n_floored > 0) and how
    ill-conditioned R was.
    """
    R = 0.5 * (R + R.T)
    w, U = np.linalg.eigh(R)
    floor = floor_frac * w.max()
    n_floored = int((w < floor).sum())
    w_reg = np.maximum(w, floor)
    Rinv = (U / w_reg) @ U.T
    # participation-ratio effective dimensionality of R
    eff_dim = float(w.sum() ** 2 / (w ** 2).sum())
    info = dict(k=int(R.shape[0]), eig_max=float(w.max()), eig_min=float(w.min()),
                cond=float(w.max() / max(w.min(), 1e-300)), n_floored=n_floored,
                floor_frac=floor_frac, eff_dim_pr=eff_dim)
    return Rinv, info


def most_stat(z: np.ndarray, Rinv: np.ndarray, chunk: int = 500_000) -> np.ndarray:
    """T_i = z_i' Rinv z_i for every row i (chunked to bound memory)."""
    z = np.asarray(z, dtype=np.float64)
    out = np.empty(z.shape[0])
    for a in range(0, z.shape[0], chunk):
        zz = z[a:a + chunk]
        out[a:a + chunk] = np.einsum("ij,jk,ik->i", zz, Rinv, zz)
    return out


def minp_stat(z: np.ndarray) -> np.ndarray:
    """min over measures of the two-sided univariate p (computed on |z| for
    precision: p = 2 * norm.sf(max |z|))."""
    zmax = np.nanmax(np.abs(z), axis=1)
    return 2.0 * stats.norm.sf(zmax)


# ---------------------------------------------------------------------------
# Null fits on the permuted scan
# ---------------------------------------------------------------------------
def _subsample(x: np.ndarray, n: int, seed: int) -> np.ndarray:
    x = x[np.isfinite(x)]
    if x.size <= n:
        return x
    rng = np.random.default_rng(seed)
    return rng.choice(x, size=n, replace=False)


def fit_gamma_null(t_perm: np.ndarray, n_fit: int = 2_000_000, seed: int = 1):
    """MLE gamma (loc fixed at 0) fitted to the permuted MOSTest statistic."""
    x = _subsample(t_perm, n_fit, seed)
    a, loc, scale = stats.gamma.fit(x, floc=0)
    return dict(shape=float(a), scale=float(scale), n_fit=int(x.size),
                mean=float(x.mean()), var=float(x.var()))


def fit_minp_null(minp_perm: np.ndarray, tail_q: float = 1e-3, n_keep: int = 2_000_000,
                  seed: int = 2):
    """Calibrate minP against the permuted scan.

    Bulk (minP above the permuted tail_q quantile): the empirical CDF of the
    permuted minP.  Tail (below it): Sidak extrapolation
    p = 1 - (1 - minP)^M_eff, with M_eff solved so the two pieces meet at the
    tail_q quantile.  A single Sidak or beta form cannot do both: with
    correlated measures the Sidak form tail-matched at 1e-3 gives lambda 0.57
    in the bulk, and an MLE beta fitted to the bulk is ~15 % conservative at
    p = 0.01 (synthetic tests, tests/test_mostest_core.py).
    """
    x = np.sort(_subsample(minp_perm, n_keep, seed))
    q = float(np.quantile(x, tail_q))
    m_eff = float(np.log1p(-tail_q) / np.log1p(-q))
    return dict(m_eff=m_eff, tail_q=tail_q, q_at_tail=q, n_fit=int(x.size), sorted_perm=x)


def p_most(t: np.ndarray, gamma_fit: dict) -> np.ndarray:
    return stats.gamma.sf(t, gamma_fit["shape"], scale=gamma_fit["scale"])


def p_minp(minp: np.ndarray, minp_fit: dict) -> np.ndarray:
    minp = np.asarray(minp, dtype=np.float64)
    x = minp_fit["sorted_perm"]
    emp = np.searchsorted(x, minp, side="right") / x.size
    # -expm1(M * log1p(-p)) = 1 - (1 - p)^M, accurate for tiny p
    sidak = -np.expm1(minp_fit["m_eff"] * np.log1p(-minp))
    return np.where(minp < minp_fit["q_at_tail"], sidak, emp)


def p_chi2(t: np.ndarray, k: int) -> np.ndarray:
    """Analytic chi-square(k) reference, exact for full-rank R with no floor.
    Reported alongside the gamma p as a consistency check."""
    return stats.chi2.sf(t, k)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
def lambda_gc(p: np.ndarray) -> float:
    p = p[np.isfinite(p) & (p > 0)]
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))


def tail_counts(p: np.ndarray, thresholds=(1e-4, 1e-6, 5e-8)) -> dict:
    p = p[np.isfinite(p)]
    n = p.size
    out = {}
    for t in thresholds:
        out[f"n_p_lt_{t:g}"] = int((p < t).sum())
        out[f"exp_p_lt_{t:g}"] = float(n * t)
    return out


def qq_table(p: np.ndarray, n_bins: int = 400) -> np.ndarray:
    """Binned QQ points (expected, observed -log10 p) on a log-spaced grid of
    ranks, so a genome-wide QQ plot can be drawn from a few hundred rows of a
    committed summary table instead of millions of per-SNP p-values."""
    p = np.sort(p[np.isfinite(p) & (p > 0)])
    n = p.size
    ranks = np.unique(np.round(np.logspace(0, np.log10(n), n_bins)).astype(int)) - 1
    exp = -np.log10((ranks + 0.5) / n)
    obs = -np.log10(p[ranks])
    return np.column_stack([exp, obs])
