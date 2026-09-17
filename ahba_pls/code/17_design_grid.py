"""
17_design_grid.py -- one tidy set of tables covering THREE Y-matrix designs in
BOTH parcellations, so a single figure can show every design's agreement with
the published components (fig1_signature_designs.R).

Designs (all thinning-oriented, see below):
  dCT + CT              option 2, PLS2   -- the lead signature
  dCT alone             option 1, PLS1   -- the single-Y control (16_dct_only.py)
  CT + dCT + T1T2 + dT1T2  option 4, component SELECTED HERE as the one whose
                        gene weights correlate most strongly (|Spearman|) with
                        the AHBA C3 weights; the selection and the runner-up are
                        written to design_grid_components.tsv rather than assumed.
                        DK only: there is no HCP-MMP T1w/T2w run (the ABCD
                        release tabulates T1w/T2w in Desikan space only, and the
                        locally-derived HCP parcellation covers thickness only),
                        so the four-feature design cannot be fitted at 137
                        parcels without a new surface run.

Nothing is refitted. Everything comes from saved fits:
  DK   results/pls_scores|pls_weights/{opt1_dCT,opt2_dCT_CT,opt4_full4}_ds25.*
  HCP  results/hcp_pls_weights.tsv (+ scores recomputed from the saved U in
       16_dct_only.py for the single-Y design; the two-Y scores are saved)
Only the spin correlations are computed here.

Orientation: pls.py fixes the FIRST Y column's salience positive, which differs
per design, so each component is multiplied by -sign(salience of dCT) -- dCT is
an age slope (negative where cortex thins), so this puts every vector in the
thinning orientation: positive = expressed where thinning is FASTER.

Outputs
  results/design_grid_scores.csv       long: design, parcellation, label, score
  results/design_grid_weights.tsv      wide: gene x design|parcellation
  results/design_grid_components.tsv   per design x parcellation component stats
  results/design_grid_concordance.tsv  scores and weights vs C3 / C1 / NSPN PLS2
  results/design_grid_pairs_scores.tsv  EVERY pair of regional maps, per parcellation
  results/design_grid_pairs_weights.tsv EVERY pair of gene vectors, per parcellation
  results/design_grid_points_scores.csv the regional values behind those pairs
  results/design_grid_points_weights.tsv the gene vectors behind those pairs

NSPN PLS2 now exists in HCP-MMP space too (18_nspn_to_hcp.py), so the pair
matrices are symmetric: every regional map is available in both parcellations.
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
RES, DATA, REF = ROOT / "results", ROOT / "data", ROOT / "data" / "reference"
DS = "ds25"
N_SPIN, SEED = 5000, 0

D_TWOY, D_ONEY, D_FOUR = "dCT + CT", "dCT alone", "CT + dCT + T1T2 + dT1T2"

# ------------------------------------------------------------- references ----
Y34 = pd.read_csv(DATA / "y_maps_bilateral_34.csv", index_col=0)
Y180 = pd.read_csv(RES / "hcp_y_maps_180.csv", index_col=0)
nspn34 = pd.read_csv(REF / "nspn_dk_maps_bilateral_34.csv", index_col=0)
nspn_hcp = pd.read_csv(REF / "nspn_hcp_maps.csv", index_col=0)   # 18_nspn_to_hcp.py
pc1 = (pd.read_csv(ROOT.parent / "data" / "velmeshev_PC1_gene_loadings.csv")
       .dropna(subset=["feature_name"]).drop_duplicates("feature_name")
       .set_index("feature_name")[["PC1_herringV3", "PC1_U01V2"]])
nspn_w = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")
c3dk = pd.read_csv(REF / f"ahba_c123_scores_recomputed_{DS}.csv", index_col=0)
c123w = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)
hcp_c123 = (pd.read_csv(ROOT.parent / "data" / "ahba_dme_hcp_top8kgenes_scores.csv")
            .assign(label=lambda d: "lh_" + d.label.astype(str)).set_index("label")[["C1", "C2", "C3"]])

comp_dk = pd.read_csv(RES / "pls_components.tsv", sep="\t")
comp_hcp = pd.read_csv(RES / "hcp_pls_components.tsv", sep="\t")

def dk_fit(option: str):
    w = pd.read_csv(RES / "pls_weights" / f"{option}_{DS}.tsv", sep="\t", index_col=0)
    s = pd.read_csv(RES / "pls_scores" / f"{option}_{DS}.csv", index_col=0)
    c = comp_dk[(comp_dk.option == option) & (comp_dk.ds == DS)].set_index("component")
    return w, s, c

# ------------------------------------- component choice for the 4-feature ----
w4, s4, c4 = dk_fit("opt4_full4")
cands = [k[:-2] for k in w4.columns if k.endswith("_Z")]
rho_c3 = {}
for k in cands:
    sh = w4[f"{k}_Z"].dropna().index.intersection(c123w["C3"].dropna().index)
    rho_c3[k] = stats.spearmanr(w4.loc[sh, f"{k}_Z"], c123w.loc[sh, "C3"]).statistic
picked = max(rho_c3, key=lambda k: abs(rho_c3[k]))
runner = sorted(rho_c3, key=lambda k: -abs(rho_c3[k]))[1]
print(f"4-feature: |rho with C3| = " + ", ".join(f"{k} {abs(v):.3f}" for k, v in rho_c3.items())
      + f" -> {picked}", file=sys.stderr)

# ------------------------------------------------------------ assemble -------
w2, s2, c2 = dk_fit("opt2_dCT_CT")
w1, s1, c1 = dk_fit("opt1_dCT")
Whcp = pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)
Shcp = pd.read_csv(RES / "hcp_pls_scores.csv", index_col=0)
Shcp1 = pd.read_csv(RES / "dct_only_scores_hcp.csv", index_col=0)

def orient(c_row) -> int:
    """-sign(dCT salience): puts the component in the thinning orientation."""
    s = c_row["sal_dCT"]
    assert not np.isnan(s) and s != 0, s
    return int(-np.sign(s))

FITS = {}   # (design, parcellation) -> dict(scores, weights, comp)
FITS[(D_TWOY, "DK")] = dict(
    comp=c2.loc["PLS2"], sign=orient(c2.loc["PLS2"]),
    scores=s2["PLS2_gene_scores"], weights=w2["PLS2_Z"])
FITS[(D_ONEY, "DK")] = dict(
    comp=c1.loc["PLS1"], sign=orient(c1.loc["PLS1"]),
    scores=s1["PLS1_gene_scores"], weights=w1["PLS1_Z"])
FITS[(D_FOUR, "DK")] = dict(
    comp=c4.loc[picked], sign=orient(c4.loc[picked]),
    scores=s4[f"{picked}_gene_scores"], weights=w4[f"{picked}_Z"])
ch2 = comp_hcp[(comp_hcp.option == "hcp_opt2_dCT_CT") & (comp_hcp.component == "PLS2")].iloc[0]
ch1 = comp_hcp[comp_hcp.option == "hcp_opt1_dCT"].iloc[0]
FITS[(D_TWOY, "HCP")] = dict(
    comp=ch2, sign=orient(ch2),
    # hcp_pls_scores.csv stores the THINNING-oriented score; undo it so the one
    # orientation rule below applies to every cell identically
    scores=-Shcp["thinning_score"], weights=Whcp["hcp_opt2_dCT_CT_PLS2_Z"])
FITS[(D_ONEY, "HCP")] = dict(
    comp=ch1, sign=orient(ch1),
    # already thinning-oriented in dct_only_scores_hcp.csv, so undo before the
    # common orientation step below keeps one rule for every cell
    scores=-Shcp1["thinning_score"], weights=Whcp["hcp_opt1_dCT_PLS1_Z"])

for key, f in FITS.items():
    f["scores"] = (f["scores"] * f["sign"]).rename("score")
    f["weights"] = (f["weights"] * f["sign"]).rename("Z")

srows, wcols = [], {}
for (design, parc), f in FITS.items():
    srows.append(f["scores"].rename_axis("label").reset_index()
                 .assign(design=design, parcellation=parc))
    wcols[f"{design}|{parc}"] = f["weights"]
S = pd.concat(srows, ignore_index=True)[["design", "parcellation", "label", "score"]]
S.to_csv(RES / "design_grid_scores.csv", index=False, float_format="%.6g")
pd.DataFrame(wcols).rename_axis("gene").to_csv(RES / "design_grid_weights.tsv", sep="\t",
                                               float_format="%.6g")

crows = []
for (design, parc), f in FITS.items():
    c = f["comp"]
    crows.append(dict(design=design, parcellation=parc,
                      component=(picked if (design, parc) == (D_FOUR, "DK")
                                 else ("PLS1" if design == D_ONEY else "PLS2")),
                      n_regions=int(f["scores"].shape[0]),
                      n_genes=int(f["weights"].notna().sum()),
                      cov_explained=c["cov_explained"], p_spin=c["p_spin_singular"],
                      boot_reproducibility=c["boot_reproducibility"],
                      sal_dCT=c["sal_dCT"],
                      selection=("|rho| with C3 weights: "
                                 f"{picked} {abs(rho_c3[picked]):.2f} > {runner} {abs(rho_c3[runner]):.2f}"
                                 if (design, parc) == (D_FOUR, "DK") else "")))
pd.DataFrame(crows).to_csv(RES / "design_grid_components.tsv", sep="\t", index=False,
                           float_format="%.4g")

# ------------------------------------------------------------ concordance ----
rows = []
for (design, parc), f in FITS.items():
    sc = f["scores"]
    refs_s = ({"C3": c3dk["C3"], "NSPN_PLS2": nspn34["PLS2"]} if parc == "DK"
              else {"C3": hcp_c123["C3"].reindex(Y180.index).dropna(),
                    "NSPN_PLS2": nspn_hcp["PLS2"].dropna()})
    cent = pls.DK_CENTROIDS if parc == "DK" else pls.HCP_CENTROIDS
    for rn, rv in refs_s.items():
        rho, p, _ = pls.spin_corr(sc, rv.dropna(), n_perm=N_SPIN, seed=SEED, centroids=cent)
        rows.append(dict(design=design, parcellation=parc, level="scores", reference=rn,
                         rho=rho, p_spin=p, n=int(sc.reindex(rv.dropna().index).notna().sum())))
    for rn, rv in (("C3", c123w["C3"]), ("C1", c123w["C1"]), ("NSPN_PLS2_z", nspn_w["PLS2_z"])):
        shared = f["weights"].dropna().index.intersection(rv.dropna().index)
        rows.append(dict(design=design, parcellation=parc, level="weights", reference=rn,
                         rho=stats.spearmanr(f["weights"].loc[shared], rv.loc[shared]).statistic,
                         p_spin=np.nan, n=len(shared)))
C = pd.DataFrame(rows)
C.to_csv(RES / "design_grid_concordance.tsv", sep="\t", index=False, float_format="%.4g")

pd.set_option("display.width", 200)
print(pd.DataFrame(crows).round(3).to_string(index=False))
print("\n" + C.pivot_table(index=["design", "parcellation"], columns=["level", "reference"],
                           values="rho").round(2).to_string())
print("\nscores spin p:")
print(C[C.level == "scores"][["design", "parcellation", "reference", "rho", "p_spin", "n"]]
      .round(4).to_string(index=False))

# ================================================================ pair matrix ==
# Every pair of regional maps, and every pair of gene vectors, within each
# parcellation -- the figure plots these as two square matrices with DK below
# the diagonal and HCP-MMP above it.
S_VARS = ["dCT rate", "dCT + CT", "dCT alone", "NSPN PLS2", "AHBA C3"]
# snRNAseq PC1 = the snRNA-seq maturation axis from the transcriptional_maturation
# project (data/velmeshev_PC1_gene_loadings.csv in the repo root). Two independent
# datasets are supplied and they agree at rho ~ 0.5, so both are carried in the
# pairs table and the figure plots the Herring one.
W_VARS = ["dCT + CT", "dCT alone", "NSPN PLS2", "AHBA C3", "AHBA C1",
          "snRNAseq PC1", "snRNAseq PC1 (U01)"]

score_vecs = {
    "DK": {"dCT rate": Y34["dCT"], "dCT + CT": FITS[(D_TWOY, "DK")]["scores"],
           "dCT alone": FITS[(D_ONEY, "DK")]["scores"],
           "NSPN PLS2": nspn34["PLS2"], "AHBA C3": c3dk["C3"]},
    "HCP": {"dCT rate": Y180["dCT"], "dCT + CT": FITS[(D_TWOY, "HCP")]["scores"],
            "dCT alone": FITS[(D_ONEY, "HCP")]["scores"],
            "NSPN PLS2": nspn_hcp["PLS2"], "AHBA C3": hcp_c123["C3"].reindex(Y180.index)},
}
REFS_W = {"NSPN PLS2": nspn_w["PLS2_z"], "AHBA C3": c123w["C3"], "AHBA C1": c123w["C1"],
          "snRNAseq PC1": pc1["PC1_herringV3"].dropna(),
          "snRNAseq PC1 (U01)": pc1["PC1_U01V2"].dropna()}
weight_vecs = {
    "DK": {"dCT + CT": FITS[(D_TWOY, "DK")]["weights"],
           "dCT alone": FITS[(D_ONEY, "DK")]["weights"], **REFS_W},
    "HCP": {"dCT + CT": FITS[(D_TWOY, "HCP")]["weights"],
            "dCT alone": FITS[(D_ONEY, "HCP")]["weights"], **REFS_W},
}

# points behind the panels, so the figure fits nothing
pts = []
for parc, vecs in score_vecs.items():
    common = None
    for v in vecs.values():
        common = v.dropna().index if common is None else common.intersection(v.dropna().index)
    for name, v in vecs.items():
        pts.append(pd.DataFrame({"parcellation": parc, "label": common,
                                 "variable": name, "value": v.reindex(common).to_numpy()}))
pd.concat(pts, ignore_index=True).to_csv(RES / "design_grid_points_scores.csv",
                                         index=False, float_format="%.6g")
wcols = {f"{n}|{parc}" if n in (D_TWOY, D_ONEY, "dCT + CT", "dCT alone") else n: v
         for parc, vecs in weight_vecs.items() for n, v in vecs.items()}
pd.DataFrame(wcols).rename_axis("gene").to_csv(RES / "design_grid_points_weights.tsv",
                                               sep="\t", float_format="%.6g")

def pair_rows(vecs, vars_, parc, level, spin):
    rows = []
    cent = pls.DK_CENTROIDS if parc == "DK" else pls.HCP_CENTROIDS
    for i, a in enumerate(vars_):
        for b in vars_[i + 1:]:
            x, y = vecs[a].dropna(), vecs[b].dropna()
            sh = x.index.intersection(y.index)
            rho = stats.spearmanr(x.loc[sh], y.loc[sh]).statistic
            p = np.nan
            if spin:
                # rotate the more complete map (b) and align to a, as elsewhere
                _, p, _ = pls.spin_corr(x.loc[sh], y, n_perm=N_SPIN, seed=SEED, centroids=cent)
            rows.append(dict(parcellation=parc, level=level, var_x=a, var_y=b,
                             rho=rho, p_spin=p, n=len(sh)))
    return rows

PS = pd.DataFrame([r for parc, vecs in score_vecs.items()
                   for r in pair_rows(vecs, S_VARS, parc, "scores", True)])
PW = pd.DataFrame([r for parc, vecs in weight_vecs.items()
                   for r in pair_rows(vecs, W_VARS, parc, "weights", False)])
PS.to_csv(RES / "design_grid_pairs_scores.tsv", sep="\t", index=False, float_format="%.4g")
PW.to_csv(RES / "design_grid_pairs_weights.tsv", sep="\t", index=False, float_format="%.4g")

print("\nregional map pairs (Spearman rho; p_spin over 5,000 rotations):")
print(PS.pivot_table(index=["var_x", "var_y"], columns="parcellation",
                     values=["rho", "p_spin"]).round(3).to_string())
print("\ngene vector pairs:")
print(PW.pivot_table(index=["var_x", "var_y"], columns="parcellation", values="rho").round(3).to_string())
