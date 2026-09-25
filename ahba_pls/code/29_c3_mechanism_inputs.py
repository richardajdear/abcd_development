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

The marker profile is DESCRIPTIVE: C3 is itself a weighted sum of the same
expression matrix, so marker expression differing between C3 thirds restates
what panel d's gene-level enrichment already shows. The informative quantity for
the mechanism is the thinning (and CT) difference between the thirds, whose
significance is the spin p of C3 vs dCT in c3_score_matrix.tsv.
"""
import sys
from itertools import combinations
from pathlib import Path
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

T = pd.read_csv(RES / "c3_tiers.csv", index_col=0)
T.index = T.index.str.lower()
F = pd.read_csv(RES / "tier_profile_3d_ds5_parcels.csv", index_col=0)
F.index = F.index.str.lower()
M.index = M.index.str.lower()
feat = ["Neuro-Ex", "Neuro-In", "Astro", "Micro", "Oligo", "OPC", "L1", "L2", "L3", "L4", "L5", "L6", "WM"]
J = T.join(F[feat]).join(M[["CT", "dCT"]])
assert len(J) == 137 and J[feat].notna().all().all(), (len(J), J[feat].isna().sum().sum())
prof = J.groupby("tier")[["normative_thinning_um_yr", "CT"] + feat].mean().T
prof = prof[["C3-low", "C3-mid", "C3-high"]]
prof["high_minus_low"] = prof["C3-high"] - prof["C3-low"]
prof.index.name = "feature"
prof.insert(0, "kind", ["thinning", "thickness"] + ["cell class"] * 6 + ["layer"] * 7)
prof.to_csv(RES / "c3_tier_profile.tsv", sep="\t", float_format="%.4g")

pd.set_option("display.width", 200)
print(S.round(3).to_string(index=False))
print(prof.round(3).to_string())
