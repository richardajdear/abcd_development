"""Unit tests for genetic_analysis/mostest/mostest_core.py (step 16).

Synthetic only: z ~ MVN(0, R) with an R shaped like the DK regional-slope
correlation (one global factor plus weak local structure, 68 measures).
Checks: (1) the gamma-calibrated MOSTest p and the tail-matched minP are
uniform under the null; (2) with full-rank R the gamma p agrees with the
analytic chi-square(k); (3) a distributed, mixed-sign effect that cancels in
the mean is found by MOSTest and not by minP -- the property the step relies on.
"""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "genetic_analysis" / "mostest"))
import mostest_core as mc  # noqa: E402

K = 68


def _toy_R(k=K, seed=0):
    rng = np.random.default_rng(seed)
    g = 0.42 * np.ones(k)                       # global factor, r ~ 0.18
    L = rng.normal(0, 0.25, size=(k, 4))        # weak local factors
    C = np.outer(g, g) + L @ L.T
    C += np.diag(1 - np.diag(C).clip(max=0.9))
    d = np.sqrt(np.diag(C))
    return C / np.outer(d, d)


@pytest.fixture(scope="module")
def null_draws():
    R = _toy_R()
    rng = np.random.default_rng(1)
    chol = np.linalg.cholesky(R)
    z_perm = rng.standard_normal((300_000, K)) @ chol.T
    z_real = rng.standard_normal((300_000, K)) @ chol.T
    return R, z_perm, z_real


def test_corr_accumulator_matches_numpy(null_draws):
    _, z, _ = null_draws
    acc = mc.CorrAccumulator(K)
    for a in range(0, z.shape[0], 70_000):
        acc.add(z[a:a + 70_000])
    np.testing.assert_allclose(acc.corr(), np.corrcoef(z, rowvar=False), atol=1e-10)


def test_null_calibration(null_draws):
    R, z_perm, z_real = null_draws
    acc = mc.CorrAccumulator(K); acc.add(z_perm)
    Rinv, info = mc.regularised_inverse(acc.corr())
    assert info["n_floored"] == 0
    g = mc.fit_gamma_null(mc.most_stat(z_perm, Rinv))
    p = mc.p_most(mc.most_stat(z_real, Rinv), g)
    assert 0.97 < mc.lambda_gc(p) < 1.03
    n = p.size
    for t in (1e-2, 1e-3):
        obs = (p < t).sum(); exp = n * t
        assert abs(obs - exp) < 5 * np.sqrt(exp), (t, obs, exp)
    # gamma vs analytic chi2(k) for full-rank R
    pc = mc.p_chi2(mc.most_stat(z_real, Rinv), K)
    assert np.corrcoef(-np.log10(p), -np.log10(pc))[0, 1] > 0.999
    # minP calibration (M_eff, tail-matched)
    b = mc.fit_minp_null(mc.minp_stat(z_perm))
    pm = mc.p_minp(mc.minp_stat(z_real), b)
    assert 0.95 < mc.lambda_gc(pm) < 1.05
    for t in (1e-2, 1e-3):
        assert abs((pm < t).sum() - n * t) < 5 * np.sqrt(n * t), t
    assert 1 < b["m_eff"] < K


def test_distributed_mixed_sign_effect(null_draws):
    R, z_perm, _ = null_draws
    acc = mc.CorrAccumulator(K); acc.add(z_perm)
    Rinv, _ = mc.regularised_inverse(acc.corr())
    g = mc.fit_gamma_null(mc.most_stat(z_perm, Rinv))
    b = mc.fit_minp_null(mc.minp_stat(z_perm))
    rng = np.random.default_rng(5)
    # 200 causal SNPs, each shifting all 68 measures by +/-1.2 SD-of-z with
    # random signs: mean effect ~0, no single measure near genome-wide.
    signs = rng.choice([-1.0, 1.0], size=(200, K))
    z = rng.standard_normal((200, K)) @ np.linalg.cholesky(R).T + 1.2 * signs
    p_most = mc.p_most(mc.most_stat(z, Rinv), g)
    p_minp = mc.p_minp(mc.minp_stat(z), b)
    z_mean = z.mean(axis=1) / np.sqrt(R.mean())       # univariate scan of the mean
    p_mean = 2 * stats.norm.sf(np.abs(z_mean))
    assert np.median(p_most) < 1e-8
    assert np.median(p_minp) > 1e-4
    assert np.median(p_mean) > 1e-3


def test_floor_engages_on_singular_R():
    R = _toy_R()
    R2 = np.block([[R, R], [R, R]])                     # exactly duplicated measures
    _, info = mc.regularised_inverse(R2)
    assert info["n_floored"] >= K


def test_qq_table_shape():
    p = np.random.default_rng(3).uniform(size=100_000)
    q = mc.qq_table(p)
    assert q.shape[1] == 2 and q.shape[0] > 100
    assert np.abs(q[:, 0] - q[:, 1]).max() < 1.5
