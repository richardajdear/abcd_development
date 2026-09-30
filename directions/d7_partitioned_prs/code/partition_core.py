"""Statistics for the D7 partitioned polygenic score (no I/O; unit-tested).

Notation.  G is an n x g matrix of per-gene PARTIAL scores: column j is the sum,
over the SNPs uniquely assigned to gene j, of SBayesRC posterior weight x allele
dosage.  The (non-MHC) full score is S = G.sum(1) + intergenic.  A set T has the
partial score S_T = G[:, T].sum(1).  y is the outcome.  A tilde means
residualised on the covariates (Frisch-Waugh-Lovell), after scaling y and S to
unit marginal SD in the analysis sample, so beta matches tools/prs_assoc.R's
scale(y) ~ scale(PRS) + covariates (OLS here, lmer there; gated in 03).

Exact additive decomposition (no independence assumption):
    s_T = cov(y~, S_T~) / cov(y~, S~)      share of the score's association with y
    f_T = cov(S_T~, S~) / var(S~)          share of the score's variance
    ER_T = s_T / f_T                        enrichment ratio; 1 if the set carries
                                            signal in proportion to its variance
Both s and f sum to 1 over any partition of the columns.  Under "signal is
spread in proportion to disorder-score variance", ER = 1 whatever f is; a set
that carries all the signal in a quarter of the variance has ER = 4 and
beta_T = beta_full * s / sqrt(f) = 2 x beta_full (DIRECTIONS.md D7 table).
Because cov is linear, c_j = cov(y~, G_j~) and d_j = cov(G_j~, S~) are computed
once per gene and any set's s, f are sums of them: random-set nulls are free.
"""
from __future__ import annotations

import numpy as np


def residualise(M: np.ndarray, C: np.ndarray) -> np.ndarray:
    """M minus its least-squares projection on the columns of C (C includes 1)."""
    coef, *_ = np.linalg.lstsq(C, M, rcond=None)
    return M - C @ coef


def within_cluster_transform(G: np.ndarray, full: np.ndarray, cluster: np.ndarray):
    """The `_zanc` transform of setup/standardise_within_ancestry.py, applied
    linearly so the partial scores still add up to the transformed full score:
    centre every column within cluster and divide by the cluster SD of the FULL
    score (ddof=1, as that script).  Returns (G', full')."""
    G = G.astype(np.float64, copy=True)
    full = full.astype(np.float64, copy=True)
    for c in np.unique(cluster):
        m = cluster == c
        sd = full[m].std(ddof=1)
        if not np.isfinite(sd) or sd <= 0:
            G[m] = np.nan
            full[m] = np.nan
            continue
        G[m] = (G[m] - G[m].mean(0)) / sd
        full[m] = (full[m] - full[m].mean()) / sd
    return G, full


def gene_moments(y_t: np.ndarray, G_t: np.ndarray, S_t: np.ndarray):
    """Per-gene c_j, d_j and the full-score cov(y,S), var(S) (all residualised)."""
    n = len(y_t)
    c = G_t.T @ y_t / (n - 1)
    d = G_t.T @ S_t / (n - 1)
    return c, d, float(y_t @ S_t / (n - 1)), float(S_t @ S_t / (n - 1))


def share(idx, c, d, cyS, vS):
    s = c[idx].sum() / cyS
    f = d[idx].sum() / vS
    return s, f, s / f


def ols_cluster(y_t: np.ndarray, x_t: np.ndarray, groups: np.ndarray):
    """Slope of y~ on x~ (both already residualised) with a family-clustered
    sandwich SE (CR1).  Returns beta, se, p (normal)."""
    from scipy.stats import norm
    sxx = x_t @ x_t
    b = (x_t @ y_t) / sxx
    e = y_t - b * x_t
    _, inv = np.unique(groups, return_inverse=True)
    u = np.bincount(inv, weights=x_t * e)
    G = len(u)
    n = len(y_t)
    v = (u @ u) / sxx ** 2 * G / (G - 1) * (n - 1) / (n - 2)
    se = np.sqrt(v)
    return b, se, 2 * norm.sf(abs(b / se))


def bin_genes(universe: np.ndarray, *features: np.ndarray, q=(10,)) -> np.ndarray:
    """Joint quantile bins of per-gene features, computed within the universe.
    Returns an int bin label per gene in `universe` order."""
    lab = np.zeros(len(universe), dtype=np.int64)
    mult = 1
    for feat, k in zip(features, q):
        x = feat[universe]
        r = np.argsort(np.argsort(x, kind="stable"), kind="stable")
        b = np.minimum((r * k) // len(x), k - 1)
        lab += b * mult
        mult *= k
    return lab


def matched_draws(target: np.ndarray, universe: np.ndarray, bins: np.ndarray,
                  n_draw: int, rng: np.random.Generator) -> np.ndarray:
    """n_draw random gene sets from `universe` with the target's bin composition.

    target, universe: column indices (target must be a subset of universe);
    bins: labels aligned to `universe`.  A bin with fewer universe genes than the
    target needs is merged with its neighbour label (rare; reported by caller).
    Returns an (n_draw, len(target)) array of column indices."""
    pos = {g: i for i, g in enumerate(universe)}
    tb = bins[[pos[t] for t in target]]
    labels, need = np.unique(tb, return_counts=True)
    pools = {l: universe[bins == l] for l in np.unique(bins)}
    out = np.empty((n_draw, len(target)), dtype=np.int64)
    k = 0
    merged = 0
    for l, m in zip(labels, need):
        pool = pools[l]
        if len(pool) < m:          # merge with the nearest labels until large enough
            merged += 1
            extra = [pools[x] for x in sorted(pools, key=lambda z: abs(z - l)) if x != l]
            while len(pool) < m and extra:
                pool = np.concatenate([pool, extra.pop(0)])
        # without replacement within a draw: take the m smallest of random keys
        keys = rng.random((n_draw, len(pool)))
        pick = np.argpartition(keys, m - 1, axis=1)[:, :m] if m < len(pool) else \
            np.tile(np.arange(len(pool)), (n_draw, 1))
        out[:, k:k + m] = pool[pick]
        k += m
    matched_draws.last_merged = merged
    return out


def null_distribution(draws: np.ndarray, c, d, cyS, vS):
    s = c[draws].sum(1) / cyS
    f = d[draws].sum(1) / vS
    return s, f, s / f


def p_two_sided(obs: float, null: np.ndarray) -> float:
    n = len(null)
    hi = (1 + np.sum(null >= obs)) / (n + 1)
    lo = (1 + np.sum(null <= obs)) / (n + 1)
    return float(min(1.0, 2 * min(hi, lo)))


def p_upper(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))


def family_bootstrap(y_t, S_t, parts: dict, groups, B: int, rng):
    """Family-resampled bootstrap of s, f, ER for each part (a residualised
    n-vector).  Covariates are not re-partialled within replicates (a stated
    approximation; the residualisation step is a fixed linear map).
    Returns {name: (B, 3) array of s, f, ER}."""
    _, inv = np.unique(groups, return_inverse=True)
    nfam = inv.max() + 1
    members = np.split(np.argsort(inv, kind="stable"), np.cumsum(np.bincount(inv))[:-1])
    out = {k: np.empty((B, 3)) for k in parts}
    for b in range(B):
        pick = rng.integers(0, nfam, nfam)
        ix = np.concatenate([members[i] for i in pick])
        yb, Sb = y_t[ix], S_t[ix]
        cyS = yb @ Sb
        vS = Sb @ Sb
        for k, P in parts.items():
            pb = P[ix]
            s, f = (yb @ pb) / cyS, (pb @ Sb) / vS
            out[k][b] = (s, f, s / f)
    return out
