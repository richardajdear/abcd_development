"""
29_c3_mechanism_inputs.py -- group-level inputs for the simplified mechanism
slide (fig5s_c3_mechanism.R), hcp_3d_ds5 fit only.

  results/c3_score_matrix.tsv   Spearman rho + spin p for every pair of the six
                                region maps CT, dCT, PLS1, PLS2, AHBA C1, AHBA C3
                                (137 parcels with AHBA coverage where C/PLS maps
                                are involved; 179 for CT-dCT)
  results/c3_tier_profile.tsv   by thirds of AHBA C3 (c3_tiers.csv, written by
                                26_gradient_scores.py): mean normative thinning
                                (um/yr) and baseline CT, and mean z-expression of
                                each cell-class / layer marker set
                                (tier_profile_3d_ds5_parcels.csv, 27_tier_profile.py)
  results/c3_q20_profile.tsv    the same for the bottom 20% / middle 60% / top 20% of C3
                                (c3_q20.csv), plus a spin test of the top-minus-bottom
                                difference in normative thinning (C3 map rotated,
                                extremes re-selected on every rotation)

The marker profile is DESCRIPTIVE: C3 is itself a weighted sum of the same
expression matrix, so marker expression differing between C3 thirds restates
what panel d's gene-level enrichment already shows. The informative quantity for
the mechanism is the thinning (and CT) difference between the thirds, whose
significance is the spin p of C3 vs dCT in c3_score_matrix.tsv.
"""
import sys
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls  # noqa: E402

RES = HERE.parent / "results"
N_SPIN, SEED = 5000, 0

M = pd.read_csv(RES / "hcp_summary_maps.csv", index_col=0)
VARS = ["CT", "dCT", "PLS1", "PLS2", "C1", "C3"]
rows = []
for a, b in combinations(VARS, 2):
    d = M[[a, b]].dropna()
    rho, p, _ = pls.spin_corr(d[a], d[b], n_perm=N_SPIN, seed=SEED, centroids=pls.HCP_CENTROIDS)
    rows.append(dict(x=a, y=b, rho=rho, p_spin=p, n=len(d)))
S = pd.DataFrame(rows)
S.to_csv(RES / "c3_score_matrix.tsv", sep="\t", index=False, float_format="%.4g")

def profile(tiers_csv: str, out_tsv: str) -> pd.DataFrame:
    T = pd.read_csv(RES / tiers_csv, index_col=0)
    T.index = T.index.str.lower()
    F = pd.read_csv(RES / "tier_profile_3d_ds5_parcels.csv", index_col=0)
    F.index = F.index.str.lower()
    M.index = M.index.str.lower()
    feat = ["Neuro-Ex", "Neuro-In", "Astro", "Micro", "Oligo", "OPC", "L1", "L2", "L3", "L4", "L5", "L6", "WM"]
    J = T.join(F[feat]).join(M[["CT", "dCT"]])
    assert len(J) == len(T) == 137 and J[feat].notna().all().all(), (len(J), J[feat].isna().sum().sum())
    prof = J.groupby("tier")[["normative_thinning_um_yr", "CT"] + feat].mean().T
    prof = prof[["C3-low", "C3-mid", "C3-high"]]
    prof["high_minus_low"] = prof["C3-high"] - prof["C3-low"]
    prof.index.name = "feature"
    prof.insert(0, "kind", ["thinning", "thickness"] + ["cell class"] * 6 + ["layer"] * 7)
    return prof, J




prof, _ = profile("c3_tiers.csv", "c3_tier_profile.tsv")
prof.to_csv(RES / "c3_tier_profile.tsv", sep="\t", float_format="%.4g")
q20, J = profile("c3_q20.csv", "c3_q20_profile.tsv")

# spin test of the top-minus-bottom 20% thinning difference: rotate the C3 map,
# re-select its bottom / top 20% on each rotation, recompute the difference
c3 = J.C3; g = J.normative_thinning_um_yr
def diff(v):
    lo, hi = v.quantile(0.2), v.quantile(0.8)
    return g[v >= hi].mean() - g[v <= lo].mean()
labs = M.loc[c3.index].index
perms = pls.cached_spin_permutations([l.replace("lh_", "lh_") for l in pd.read_csv(RES / "c3_q20.csv", index_col=0).index],
                                     n_perm=N_SPIN, seed=SEED, path=pls.HCP_CENTROIDS)
cv = c3.to_numpy()
obs = diff(c3)
null = np.array([diff(pd.Series(cv[pp], index=c3.index)) for pp in perms])
p_spin = (1 + (np.abs(null) >= abs(obs)).sum()) / (1 + N_SPIN)
q20.loc["normative_thinning_um_yr", "p_spin_high_minus_low"] = p_spin
q20.to_csv(RES / "c3_q20_profile.tsv", sep="\t", float_format="%.4g")
print(f"top - bottom 20% thinning: {obs:.2f} um/yr, spin p = {p_spin:.4f} (null sd {null.std():.2f})")

pd.set_option("display.width", 200)
print(S.round(3).to_string(index=False))
print(prof.round(3).to_string())
print(q20.round(3).to_string())
