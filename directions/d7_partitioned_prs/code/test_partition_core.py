"""Synthetic checks of partition_core (laptop or CSD3; ~10 s).

    python -m pytest directions/d7_partitioned_prs/code/test_partition_core.py -q

Not collected by tests/ (which the README test count guards); run it directly.
Planted truth: genes in `hot` carry all of the score's association with y.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import partition_core as pc  # noqa: E402


def _sim(n=6000, g=800, hot_frac=0.25, seed=1, beta=0.25):
    rng = np.random.default_rng(seed)
    # genes with heterogeneous variance (like per-gene partial scores)
    scale = rng.gamma(1.5, 1.0, g)
    G = rng.standard_normal((n, g)) * np.sqrt(scale)
    hot = np.sort(rng.choice(g, int(hot_frac * g), replace=False))
    C = np.column_stack([np.ones(n), rng.standard_normal((n, 3))])
    G += C[:, 1:2] * 0.3                       # a covariate loads on every gene
    y = beta * G[:, hot].sum(1) / np.sqrt(scale[hot].sum()) + rng.standard_normal(n)
    fam = np.repeat(np.arange(n // 2), 2)
    return G, y, C, hot, fam, rng


def _prep(G, y, C):
    S = G.sum(1)
    y = (y - y.mean()) / y.std()
    sdS = S.std()
    Gt = pc.residualise(G / sdS, C)
    St = Gt.sum(1)
    yt = pc.residualise(y[:, None], C)[:, 0]
    return yt, Gt, St


def test_additivity():
    G, y, C, hot, fam, rng = _sim()
    yt, Gt, St = _prep(G, y, C)
    c, d, cyS, vS = pc.gene_moments(yt, Gt, St)
    a = np.arange(G.shape[1] // 3)
    b = np.setdiff1d(np.arange(G.shape[1]), a)
    sa, fa, _ = pc.share(a, c, d, cyS, vS)
    sb, fb, _ = pc.share(b, c, d, cyS, vS)
    assert abs(sa + sb - 1) < 1e-10 and abs(fa + fb - 1) < 1e-10


def test_planted_enrichment_recovered():
    G, y, C, hot, fam, rng = _sim()
    yt, Gt, St = _prep(G, y, C)
    c, d, cyS, vS = pc.gene_moments(yt, Gt, St)
    s, f, er = pc.share(hot, c, d, cyS, vS)
    assert 0.8 < s < 1.2, s                    # all signal planted in `hot`
    assert 2.0 < er < 6.5, (er, f)             # ~1/f with f ~ 0.25
    # a random set of the same size has ER near 1
    cold = rng.choice(np.setdiff1d(np.arange(G.shape[1]), hot), len(hot), replace=False)
    _, _, er0 = pc.share(cold, c, d, cyS, vS)
    assert er0 < 1.0


def test_matched_null_and_pvalue():
    G, y, C, hot, fam, rng = _sim()
    yt, Gt, St = _prep(G, y, C)
    c, d, cyS, vS = pc.gene_moments(yt, Gt, St)
    uni = np.arange(G.shape[1])
    bins = pc.bin_genes(uni, d, q=(5,))
    draws = pc.matched_draws(hot, uni, bins, 500, rng)
    assert draws.shape == (500, len(hot))
    ns, nf, ner = pc.null_distribution(draws, c, d, cyS, vS)
    _, f, er = pc.share(hot, c, d, cyS, vS)
    assert abs(np.median(nf) - f) / f < 0.25    # variance matched through the d bins
    assert pc.p_two_sided(er, ner) < 0.01


def test_within_cluster_transform_is_additive():
    rng = np.random.default_rng(3)
    G = rng.standard_normal((400, 20)) + 2.0
    full = G.sum(1)
    cl = np.repeat(np.array(["a", "b"]), 200)
    G2, f2 = pc.within_cluster_transform(G, full, cl)
    assert np.allclose(G2.sum(1), f2)
    for k in ("a", "b"):
        assert abs(f2[cl == k].std(ddof=1) - 1) < 1e-10


def test_cluster_se_and_bootstrap():
    G, y, C, hot, fam, rng = _sim(n=3000, g=300)
    yt, Gt, St = _prep(G, y, C)
    b, se, p = pc.ols_cluster(yt, St / St.std(), fam)
    assert se > 0 and 0 <= p <= 1
    bs = pc.family_bootstrap(yt, St, {"hot": Gt[:, hot].sum(1)}, fam, 200, rng)["hot"]
    assert np.isfinite(bs).all() and bs[:, 2].std() > 0
