"""
20_basis_vs_parcellation.py -- the ABCD dCT+CT component's score map correlates
only rho = 0.74 between its DK and HCP-MMP fits (19_resample_check.py), while
its gene weights agree at 0.90.  Which of the two differences between those fits
does that?

  * EXPRESSION BASIS -- the DK fit uses the AHBA_updated native-DK matrices
    (ds25, 12,007 genes); the HCP fit uses abagen-data hcp_3d.csv on the 7,973
    genes of the shipped C1-C3 weights.  Different build, different gene list.
  * PARCELLATION -- 33 DK regions against 137 HCP-MMP parcels.

The design separates them with a third fit: the same option-2 PLS (Y = dCT + CT)
on the same 33 DK regions but with the HCP fit's expression basis
(dk_3d.csv restricted to those 7,973 genes -- the "matched" X of 12_hcp_pls.py,
refitted here so its SCORES are available and not only its gene weights).

  native DK  vs  matched DK   -> basis effect, parcellation held fixed
  matched DK vs  HCP -> DK    -> parcellation effect, basis held fixed
  native DK  vs  HCP -> DK    -> the 0.74 being decomposed

The HCP score map is brought down to DK regions through fsaverage vertices
(code/surface.py), the same route used for the NSPN maps.

The second question this answers: NSPN PLS2 agrees with the native DK component
at rho = 0.75 and with the HCP one at 0.29.  Whichever of the two fits the
matched-basis component resembles tells us whether that gap is about the genes
or about the atlas.

Outputs
  results/basis_vs_parcellation.tsv   the decomposition
  results/basis_matched_dk_scores.csv the third fit's region scores and weights
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls, surface                                      # noqa: E402

ROOT = HERE.parent
REPO = ROOT.parent
RES, DATA, REF = ROOT / "results", ROOT / "data", ROOT / "data" / "reference"
EXPR = Path.home() / "Git" / "AHBA" / "data" / "abagen-data" / "expression"
N_BOOT, N_SPIN, SEED = 1000, 5000, 0

# ---------------------------------------------------- the matched-basis fit ---
# id -> label mapping and validation are as in 12_hcp_pls.py.
FS_DK = ["bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus", "entorhinal",
         "fusiform", "inferiorparietal", "inferiortemporal", "isthmuscingulate", "lateraloccipital",
         "lateralorbitofrontal", "lingual", "medialorbitofrontal", "middletemporal",
         "parahippocampal", "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis",
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", "precuneus",
         "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal", "superiorparietal",
         "superiortemporal", "supramarginal", "frontalpole", "temporalpole", "transversetemporal",
         "insula"]
genes = pd.read_csv(REPO / "data" / "weights.csv", index_col=0).index
Xd = pd.read_csv(EXPR / "dk_3d.csv", index_col=0)[genes]
Xd.index = ["lh_" + FS_DK[i - 1] for i in Xd.index]
_ship = pd.read_csv(REPO / "data" / "ahba_dme_dsk_scores.csv").set_index("label")
_Xz = pls.zscore_cols(Xd.to_numpy(float)); _Xz -= _Xz.mean(0, keepdims=True)
_, _, _Vt = np.linalg.svd(_Xz, full_matrices=False)
_pc1 = pd.Series(_Xz @ _Vt[0], index=Xd.index)
_sh = Xd.index.intersection(_ship.index)
assert abs(stats.spearmanr(_pc1.loc[_sh], _ship.loc[_sh, "C1"]).statistic) > 0.95, \
    "DK id -> label mapping looks wrong"

Y = pd.read_csv(DATA / "y_maps_bilateral_34.csv", index_col=0)[["dCT", "CT"]].loc[Xd.index]
fit = pls.pls_svd(Xd, Y)
bt = pls.bootstrap_weights(Xd, Y, n_boot=N_BOOT, seed=SEED)
# thinning orientation, the convention everywhere in this directory
sign = -np.sign(fit.saliences().loc["dCT", "PLS2"])
S_matched = sign * fit.scores()["PLS2_gene_scores"]
W_matched = sign * bt["Z"]["PLS2"]
pd.DataFrame({"score": S_matched}).rename_axis("label").to_csv(
    RES / "basis_matched_dk_scores.csv", float_format="%.6g")

# ------------------------------------------------------------- the three maps -
SP = pd.read_csv(RES / "design_grid_points_scores.csv")
get = lambda parc, var: (SP[(SP.parcellation == parc) & (SP.variable == var)]
                         .set_index("label").value)
S_native = get("DK", "dCT + CT")
S_hcp = get("HCP", "dCT + CT").rename(lambda s: s.replace("lh_", ""))
S_hcp2dk = surface.resample(S_hcp, "lh.HCPMMP1.annot", "lh.aparc.annot").rename(lambda s: "lh_" + s)
nspn = pd.read_csv(REF / "nspn_dk_maps_bilateral_34.csv", index_col=0)["PLS2"]
c3 = pd.read_csv(REF / "ahba_c123_scores_recomputed_ds25.csv", index_col=0)["C3"]

# Two DK versions of the NSPN map exist and they are not interchangeable at
# n = 33: the published DK-level table, and the 308-region map resampled to DK
# (they agree at r = 0.985 but that still leaves room for rho to move by ~0.2 on
# 33 points). Both are reported rather than picking one.
_t = pd.read_csv(Path.home() / "Git" / "AHBA" / "data" / "whitakervertes2016_complete_308.csv",
                 index_col=0)
_t["PLS2"] = _t["PLS2"].replace(-99, np.nan)
_lh = _t[_t.hemi == "l"].copy()
_lh["key"] = _lh.region + "_part" + _lh.n_sub_regions.astype(int).astype(str)
nspn_rs = surface.resample(_lh.set_index("key")["PLS2"], "lh.500.aparc.annot", "lh.aparc.annot",
                           src_key=lambda n: n).rename(lambda s: "lh_" + s)

# every comparison is restricted to the 33 regions the ABCD DK fits use, so the
# rows are comparable (the HCP-derived map otherwise contributes a 34th region,
# frontalpole, which has no AHBA donor coverage in the DK fits)
REGIONS = S_native.dropna().index
maps = {"native DK (AHBA_updated ds25)": S_native,
        "matched DK (abagen-data, HCP gene basis)": S_matched,
        "HCP-MMP pushed to DK": S_hcp2dk.reindex(REGIONS)}

rows = []
def cmp_(a: str, b: str, x: pd.Series, y: pd.Series, what: str, spin: bool = True):
    sh = x.dropna().index.intersection(y.dropna().index)
    rho = stats.spearmanr(x.loc[sh], y.loc[sh]).statistic
    p = np.nan
    if spin:
        _, p, _ = pls.spin_corr(x.loc[sh], y.loc[sh], n_perm=N_SPIN, seed=SEED,
                                centroids=pls.DK_CENTROIDS)
    rows.append(dict(comparison=what, a=a, b=b, rho=rho, p_spin=p, n=len(sh)))

names = list(maps)
cmp_(names[0], names[1], maps[names[0]], maps[names[1]], "basis effect (parcellation fixed)")
cmp_(names[1], names[2], maps[names[1]], maps[names[2]], "parcellation effect (basis fixed)")
cmp_(names[0], names[2], maps[names[0]], maps[names[2]], "total (the 0.74 being decomposed)")
for nm, s_ in maps.items():
    cmp_(nm, "NSPN PLS2 (published DK table)", s_, nspn, "vs NSPN PLS2")
    cmp_(nm, "NSPN PLS2 (308 resampled to DK)", s_, nspn_rs, "vs NSPN PLS2")
    cmp_(nm, "AHBA C3 (DK, ds25)", s_, c3, "vs AHBA C3")

# gene weights, for the same three fits
W_native = pd.read_csv(RES / "lead_signature_weights.tsv", sep="\t", index_col=0)["thinning_Z_ds25"]
W_hcp = -pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)["hcp_opt2_dCT_CT_PLS2_Z"]
for a, b, x, y in [(names[0], names[1], W_native, W_matched),
                   (names[1], names[2], W_matched, W_hcp),
                   (names[0], names[2], W_native, W_hcp)]:
    cmp_(a, b, x, y, "gene weights", spin=False)

T = pd.DataFrame(rows)
T.to_csv(RES / "basis_vs_parcellation.tsv", sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 200, "display.max_colwidth", 44)
print(T.round(4).to_string(index=False))
print("\nmatched-basis fit: PLS2 covariance explained "
      f"{fit.cov_explained[1]:.3f}, dCT salience {fit.saliences().loc['dCT', 'PLS2']:+.3f}")
