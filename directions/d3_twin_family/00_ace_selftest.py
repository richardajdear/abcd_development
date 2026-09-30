"""Self-test of ace.py: parameter recovery on simulated twin pairs.

    python directions/d3_twin_family/00_ace_selftest.py      (env abcd-spatial; ~30 s)

Simulates bivariate ACE data with known h2, c2, rA, rC, rE (a) at large n (20,000 MZ +
20,000 DZ pairs, 5 replicates) and (b) at the ABCD twin sample size (271 MZ, 432 DZ, 10%
missing, 50 replicates). The replicate mean should match the truth; the SD is the
expected sampling error at that n. Writes results/ace_selftest.tsv and fails if the
large-n mean misses by more than 0.02 (variance components, rE), 0.03 (rA) or 0.06 (rC,
poorly determined because c2_2 = 0.1: single-run SD at n = 20,000 is ~0.05).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import ace  # noqa: E402
from common import RES  # noqa: E402

TRUE = dict(h2_1=0.5, h2_2=0.3, c2_1=0.2, c2_2=0.1, rA=0.6, rC=0.0, rE=0.2)
rng = np.random.default_rng(1)


def sim(n, k, miss=0.0):
    a, c = np.sqrt([TRUE["h2_1"], TRUE["h2_2"]]), np.sqrt([TRUE["c2_1"], TRUE["c2_2"]])
    e = np.sqrt(1 - a ** 2 - c ** 2)
    cm = lambda sd, r: np.outer(sd, sd) * np.array([[1, r], [r, 1]])
    A, C, E = cm(a, TRUE["rA"]), cm(c, TRUE["rC"]), cm(e, TRUE["rE"])
    S, X = A + C + E, k * A + C
    Y = rng.multivariate_normal(np.zeros(4), np.block([[S, X], [X, S]]), n)
    Y[rng.random(Y.shape) < miss] = np.nan
    return Y


rows = []
for label, n1, n2, miss, nrep in (("n = 20,000 + 20,000 pairs, 5 reps", 20000, 20000, 0.0, 5),
                                   ("ABCD size (271 MZ, 432 DZ, 10% missing), 50 reps", 271, 432, 0.1, 50)):
    reps = [ace.fit(sim(n1, 1.0, miss), sim(n2, 0.5, miss), "ACE") for _ in range(nrep)]
    for k, v in TRUE.items():
        x = np.array([r[k] for r in reps])
        rows.append(dict(scenario=label, param=k, truth=v, estimate=x.mean(), sd=x.std()))
R = pd.DataFrame(rows)
R.to_csv(RES / "ace_selftest.tsv", sep="\t", index=False, float_format="%.4g")
print(R.round(3).to_string(index=False))
big = R[R.scenario.str.startswith("n =")]
TOL = dict(h2_1=0.02, h2_2=0.02, c2_1=0.02, c2_2=0.02, rA=0.03, rC=0.06, rE=0.02)
bad = big[(big.estimate - big.truth).abs() > big.param.map(TOL)]
assert bad.empty, f"ace.py fails large-n parameter recovery:\n{bad}"
print("self-test passed")
