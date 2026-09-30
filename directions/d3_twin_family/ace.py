"""Classical twin (ACE) models by full-information maximum likelihood.

Why not OpenMx: there is no conda build of r-openmx for this platform, so the models
are written out here (~150 lines, normal-theory FIML, the same likelihood OpenMx
maximises). 07_openmx_check.R cross-checks the univariate fits against OpenMx
installed from CRAN when that is available.

Data: one row per pair, columns for twin 1 and twin 2 of each trait; traits are
covariate-residualised z-scores (03_build_traits_pairs.py), so only a mean per trait is
estimated. Missing values are handled by FIML (per missingness pattern).

Univariate:  var = a^2 + c^2 + e^2;  cov_MZ = a^2 + c^2;  cov_DZ = a^2/2 + c^2.
Bivariate (correlated factors): A = D_a R_A D_a, C = D_c R_C D_c, E = D_e R_E D_e;
  within-person S = A + C + E, cross-twin X = k A + C with k = 1 (MZ), 1/2 (DZ).
  Reported: h2_k, c2_k, e2_k per trait; rA, rC, rE; rP; and the share of rP due to A.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

LOG2PI = np.log(2 * np.pi)


def _loglik(Y: np.ndarray, mu: np.ndarray, S: np.ndarray) -> float:
    """Sum of multivariate-normal log densities with FIML over missingness patterns."""
    miss = np.isnan(Y)
    keys = np.packbits(~miss, axis=1, bitorder="little")[:, 0]
    ll = 0.0
    for k in np.unique(keys):
        rows = keys == k
        obs = ~miss[rows][0]
        if not obs.any():
            continue
        y = Y[rows][:, obs] - mu[obs]
        Si = S[np.ix_(obs, obs)]
        try:
            L = np.linalg.cholesky(Si)
        except np.linalg.LinAlgError:
            return -1e12
        z = np.linalg.solve(L, y.T)
        ll += -0.5 * (y.shape[0] * (obs.sum() * LOG2PI + 2 * np.log(np.diag(L)).sum()) + (z ** 2).sum())
    return ll


def _cov(p: dict, k: float, nv: int) -> np.ndarray:
    def comp(sd, r):
        R = np.eye(nv)
        if nv == 2:
            R[0, 1] = R[1, 0] = r
        return np.outer(sd, sd) * R
    A, C, E = comp(p["a"], p.get("rA", 0)), comp(p["c"], p.get("rC", 0)), comp(p["e"], p.get("rE", 0))
    S = A + C + E
    X = k * A + C
    return np.block([[S, X], [X, S]])


def _unpack(theta: np.ndarray, nv: int, model: str) -> dict:
    i = 0
    p = {"mu": theta[i:i + nv]}; i += nv
    for comp in "ace":
        if comp in model.lower() or comp == "e":
            p[comp] = theta[i:i + nv]; i += nv
        else:
            p[comp] = np.zeros(nv)
    if nv == 2:
        for comp in "ace":
            key = "r" + comp.upper()
            if comp in model.lower() or comp == "e":
                p[key] = np.tanh(theta[i]); i += 1
            else:
                p[key] = 0.0
    return p


def fit(mz: np.ndarray, dz: np.ndarray, model: str = "ACE", start: np.ndarray | None = None) -> dict:
    """mz, dz: arrays (n_pairs, 2*nv) ordered [t1 trait1..nv, t2 trait1..nv]. model in ACE/AE/CE/E."""
    nv = mz.shape[1] // 2
    ncomp = sum(1 for c in "ace" if c in model.lower() or c == "e")
    npar = nv + ncomp * nv + (ncomp if nv == 2 else 0)
    if start is None:
        allv = np.vstack([mz[:, :nv], mz[:, nv:], dz[:, :nv], dz[:, nv:]])
        sd = np.nanstd(allv, axis=0)
        start = np.concatenate([np.nanmean(allv, 0)] + [sd / np.sqrt(ncomp)] * ncomp +
                               ([np.full(ncomp, 0.1)] if nv == 2 else []))
    assert len(start) == npar

    def nll(th):
        p = _unpack(th, nv, model)
        mu = np.concatenate([p["mu"], p["mu"]])
        return -(_loglik(mz, mu, _cov(p, 1.0, nv)) + _loglik(dz, mu, _cov(p, 0.5, nv)))

    best = None
    for s in (start, start * np.r_[np.ones(nv), np.full(npar - nv, 0.7)]):
        r = minimize(nll, s, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    p = _unpack(best.x, nv, model)
    va, vc, ve = p["a"] ** 2, p["c"] ** 2, p["e"] ** 2
    vt = va + vc + ve
    out = {"model": model, "minus2ll": 2 * best.fun, "npar": npar, "converged": bool(best.success),
           "n_mz": int(len(mz)), "n_dz": int(len(dz)), "theta": best.x}
    for j in range(nv):
        out[f"h2_{j + 1}"] = va[j] / vt[j]; out[f"c2_{j + 1}"] = vc[j] / vt[j]; out[f"e2_{j + 1}"] = ve[j] / vt[j]
    if nv == 2:
        sdt = np.sqrt(vt)
        cA = p["rA"] * p["a"][0] * p["a"][1]; cC = p["rC"] * p["c"][0] * p["c"][1]
        cE = p["rE"] * p["e"][0] * p["e"][1]
        rP = (cA + cC + cE) / (sdt[0] * sdt[1])
        out.update(rA=p["rA"], rC=p["rC"], rE=p["rE"], rP=rP,
                   rP_A=cA / (sdt[0] * sdt[1]), rP_C=cC / (sdt[0] * sdt[1]), rP_E=cE / (sdt[0] * sdt[1]))
    return out


def pair_arrays(T: pd.DataFrame, pairs: pd.DataFrame, traits: list[str], kind: str) -> np.ndarray:
    P = pairs[pairs.type == kind]
    a = T.reindex(P.id1)[traits].to_numpy(float)
    b = T.reindex(P.id2)[traits].to_numpy(float)
    X = np.hstack([a, b])
    return X[~np.all(np.isnan(X), axis=1)]


def twin_corrs(mz: np.ndarray, dz: np.ndarray, sib: np.ndarray | None = None) -> dict:
    def icc(X):
        m = ~np.isnan(X).any(1); X = X[m]
        # double-entry intraclass correlation
        v = np.r_[X[:, 0], X[:, 1]]; w = np.r_[X[:, 1], X[:, 0]]
        return np.corrcoef(v, w)[0, 1], int(m.sum())
    out = {}
    for name, X in (("MZ", mz), ("DZ", dz), ("SIB", sib)):
        if X is not None:
            out[f"r{name}"], out[f"n{name}"] = icc(X)
    out["h2_falconer"] = 2 * (out["rMZ"] - out["rDZ"])
    return out


def bootstrap(mz: np.ndarray, dz: np.ndarray, model: str, n: int, seed: int, start: np.ndarray,
              keys: list[str]) -> pd.DataFrame:
    """Pair-resampling bootstrap (within zygosity), warm-started at the full-sample optimum."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        m = mz[rng.integers(0, len(mz), len(mz))]
        d = dz[rng.integers(0, len(dz), len(dz))]
        f = fit(m, d, model, start=start)
        rows.append({k: f[k] for k in keys})
    return pd.DataFrame(rows)
