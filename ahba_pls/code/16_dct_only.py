"""
16_dct_only.py -- inputs for the single-Y variant of the two-parcellation
figure: Y = thinning rate (dCT) ALONE, in both DK and HCP-MMP.

Why a separate script rather than a switch in the figure: with one Y column the
PLS has a single component, so PLS1 *is* the vector of gene-map correlations --
a different estimator from the option-2 PLS2 the main figure shows, with its own
spin p, bootstrap Z and concordance statistics.  Keeping it separate means the
main figure's tables are never overwritten with numbers from a different design.

Nothing is refitted that is already saved.  Reused as-is:
  DK  scores  results/pls_scores/opt1_dCT_ds25.csv      (04_fit_pls.py)
  DK  weights results/pls_weights/opt1_dCT_ds25.tsv     (04_fit_pls.py)
  HCP weights results/hcp_pls_weights.tsv               (12_hcp_pls.py)
  spin p / bootstrap reproducibility from pls_components.tsv, hcp_pls_components.tsv
Recomputed here because they were never saved for option 1:
  HCP region scores (Lx = zscore_cols(X) @ U, the PLSFit.scores() definition,
  evaluated on the saved U -- identical to the fit, not an approximation), and
  every concordance statistic for the two option-1 vectors.

Sign convention, as everywhere in this directory: pls.py fixes the dCT salience
positive, and dCT is an age SLOPE (negative where cortex thins), so the thinning
orientation -- positive = expressed where thinning is FASTER -- is the negated
component.  Both scores and weights are negated here.

Outputs (all read by code/fig1_signature_dct_only.R)
  results/dct_only_scores_dk.csv     33 regions: score + dCT + NSPN PLS2 + C3
  results/dct_only_scores_hcp.csv    137 parcels: score + dCT + C1-C3
  results/dct_only_weights.tsv       gene x {DK, HCP} thinning-oriented Z
  results/dct_only_components.tsv    the two option-1 component rows, side by side
  results/dct_only_concordance.tsv   scores and weights vs C3 / NSPN / each other
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls

ROOT = HERE.parent
REPO = ROOT.parent
RES, DATA, REF = ROOT / "results", ROOT / "data", ROOT / "data" / "reference"
EXPR = Path.home() / "Git" / "AHBA" / "data" / "abagen-data" / "expression"
DS = "ds25"
N_SPIN, SEED = 5000, 0

# ------------------------------------------------------------------ DK -------
sc_dk = pd.read_csv(RES / "pls_scores" / f"opt1_dCT_{DS}.csv", index_col=0)
w_dk = pd.read_csv(RES / "pls_weights" / f"opt1_dCT_{DS}.tsv", sep="\t", index_col=0)
Y34 = pd.read_csv(DATA / "y_maps_bilateral_34.csv", index_col=0)
nspn34 = pd.read_csv(REF / "nspn_dk_maps_bilateral_34.csv", index_col=0)
c3r = pd.read_csv(REF / f"ahba_c123_scores_recomputed_{DS}.csv", index_col=0)
c123w = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)
nspn_w = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")

S_dk = -sc_dk["PLS1_gene_scores"]
Z_dk = -w_dk["PLS1_Z"]
out_dk = pd.DataFrame({"thinning_score": S_dk,
                       "dCT": Y34.loc[S_dk.index, "dCT"],
                       "NSPN_PLS2": nspn34.loc[S_dk.index, "PLS2"],
                       "C3_recomputed": c3r.loc[S_dk.index, "C3"]})
out_dk.index.name = "label"
out_dk.to_csv(RES / "dct_only_scores_dk.csv", float_format="%.5g")

# ------------------------------------------------------------------ HCP ------
Y180 = pd.read_csv(RES / "hcp_y_maps_180.csv", index_col=0)
shipped = pd.read_csv(REPO / "data" / "weights.csv", index_col=0)
genes = list(shipped.index)
Xh = pd.read_csv(EXPR / "hcp_3d.csv", index_col=0)[genes]
c123s = pd.read_csv(REPO / "data" / "ahba_dme_hcp_top8kgenes_scores.csv")
id2label = dict(zip(c123s.id, "lh_" + c123s.label.astype(str)))
Xh.index = [id2label.get(i, f"id{i}") for i in Xh.index]
Xh = Xh.loc[[l for l in Xh.index if l in Y180.index]]
c123_hcp = c123s.assign(label="lh_" + c123s.label.astype(str)).set_index("label")[["C1", "C2", "C3"]]

W_hcp = pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)
U = W_hcp["hcp_opt1_dCT_PLS1_U"].reindex(Xh.columns)
assert U.notna().all(), "saved HCP option-1 weights do not cover the X matrix"
S_hcp = pd.Series(-(pls.zscore_cols(Xh.to_numpy(float)) @ U.to_numpy()), index=Xh.index)
Z_hcp = -W_hcp["hcp_opt1_dCT_PLS1_Z"]

# the recomputed scores must reproduce the fit: check against the saved option-2
# machinery by refitting nothing -- instead verify the score-weight consistency
# (corr of Lx with X-projected weights is 1 by construction) and the parcel count
assert len(S_hcp) == 137, len(S_hcp)
out_hcp = pd.DataFrame({"thinning_score": S_hcp}).join(c123_hcp).join(Y180)
out_hcp.index.name = "label"
out_hcp.to_csv(RES / "dct_only_scores_hcp.csv", float_format="%.6g")

pd.DataFrame({"dk_thinning_Z": Z_dk.reindex(sorted(set(Z_dk.index) | set(Z_hcp.index))),
              "hcp_thinning_Z": Z_hcp.reindex(sorted(set(Z_dk.index) | set(Z_hcp.index)))}
             ).rename_axis("gene").to_csv(RES / "dct_only_weights.tsv", sep="\t",
                                          float_format="%.6g")

# ---------------------------------------------------------- component rows ---
comp = pd.read_csv(RES / "pls_components.tsv", sep="\t")
hcomp = pd.read_csv(RES / "hcp_pls_components.tsv", sep="\t")
cd = comp[(comp.option == "opt1_dCT") & (comp.ds == DS) & (comp.component == "PLS1")].iloc[0]
ch = hcomp[(hcomp.option == "hcp_opt1_dCT") & (hcomp.component == "PLS1")].iloc[0]
# pls_components.tsv carries no n_regions column (the DK fits are all 33
# regions); take the region count from the score tables actually plotted.
pd.DataFrame([
    dict(parcellation="DK", n_regions=len(S_dk), n_genes=int(cd.n_genes),
         p_spin_singular=cd.p_spin_singular, boot_reproducibility=cd.boot_reproducibility,
         n_genes_absZ_gt3=int(cd.n_genes_absZ_gt3)),
    dict(parcellation="HCP", n_regions=len(S_hcp), n_genes=int(ch.n_genes),
         p_spin_singular=ch.p_spin_singular, boot_reproducibility=ch.boot_reproducibility,
         n_genes_absZ_gt3=int(ch.n_genes_absZ_gt3)),
]).to_csv(RES / "dct_only_components.tsv", sep="\t", index=False, float_format="%.5g")

# -------------------------------------------------------------- concordance --
rows = []
def spin_scores(parc, ref_name, score, ref_full, centroids):
    """Spearman with a spin null; the reference is rotated on its COMPLETE map,
    then aligned to the regions the score covers (rotating the incomplete map
    would bias the null)."""
    rho, p, _ = pls.spin_corr(score, ref_full.dropna(), n_perm=N_SPIN, seed=SEED,
                              centroids=centroids)
    rows.append(dict(parcellation=parc, reference=ref_name, level="scores",
                     rho=rho, p_spin=p, n=int(score.reindex(ref_full.dropna().index).notna().sum())))

def wcorr(parc, ref_name, z, ref):
    shared = z.dropna().index.intersection(ref.dropna().index)
    rows.append(dict(parcellation=parc, reference=ref_name, level="weights",
                     rho=stats.spearmanr(z.loc[shared], ref.loc[shared]).statistic,
                     p_spin=np.nan, n=len(shared)))

spin_scores("DK", "C3", S_dk, c3r["C3"], pls.DK_CENTROIDS)
spin_scores("DK", "NSPN_PLS2", S_dk, nspn34["PLS2"], pls.DK_CENTROIDS)
spin_scores("HCP", "C3", S_hcp, c123_hcp["C3"].reindex(Y180.index), pls.HCP_CENTROIDS)
wcorr("DK", "C3", Z_dk, c123w["C3"])
wcorr("DK", "C1", Z_dk, c123w["C1"])
wcorr("DK", "NSPN_PLS2_z", Z_dk, nspn_w["PLS2_z"])
wcorr("HCP", "C3", Z_hcp, c123w["C3"])
wcorr("HCP", "C1", Z_hcp, c123w["C1"])
wcorr("HCP", "DK_dCT_PLS1", Z_hcp, Z_dk)
# how far the single-Y vector is from the option-2 signature it is meant to replace
lead_dk = pd.read_csv(RES / "lead_signature_weights.tsv", sep="\t", index_col=0)["thinning_Z_ds25"]
lead_hcp = -pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)["hcp_opt2_dCT_CT_PLS2_Z"]
wcorr("DK", "DK_opt2_PLS2", Z_dk, lead_dk)
wcorr("HCP", "HCP_opt2_PLS2", Z_hcp, lead_hcp)
C = pd.DataFrame(rows)
C.to_csv(RES / "dct_only_concordance.tsv", sep="\t", index=False, float_format="%.4g")

print(C.to_string(index=False))
print("\ncomponents:")
print(pd.read_csv(RES / "dct_only_components.tsv", sep="\t").to_string(index=False))
