"""
12_hcp_pls.py -- the option-2 PLS (dCT + CT) in the HCP-MMP parcellation, the
native space of AHBA C1-C3, and its SCZ/MDD enrichment against the DK version.

Exploratory: the HCP thickness table is derived from the release FreeSurfer
surfaces and covers 24,921 of 33,825 sessions, so the run behind this
(thickness_hcp_70_*) has 5,947 subjects against the DK run's 8,192.  The
question is only whether the finer parcellation adds gene-level signal, so the
smaller sample is accepted and stated rather than corrected.

Design, matched to the DK analysis except for the parcellation:
  X : AHBA HCP-MMP expression (hcp_3d.csv) on the 7,973 genes the shipped C1-C3
      weights were fitted on, restricted to the 137 left parcels with donor
      coverage; genes z-scored across regions.  The DK comparator uses the same
      7,973 genes from dk_3d.csv so the two are gene-matched.
  Y : bilateral (lh/rh mean) HCP dCT (age slope) and CT (intercept) from the
      run's fits/fixed.parquet.  Parcel H is dropped -- it sits on the
      FreeSurfer medial wall and has thickness 0 in ~6,350 sessions.
  Inference: 5,000 spin rotations of the complete 180-parcel map (HCP spherical
      centroids from data/hcp_centroids.csv) then restriction to the covered
      137; gene weights = Z over 1,000 region bootstraps, Procrustes-aligned.

Outputs
  results/hcp_pls_components.tsv    singular values, cov explained, spin p, saliences
  results/hcp_pls_weights.tsv       gene weights (U) and bootstrap Z per component
  results/hcp_pls_scores.csv        137-parcel scores with the C1-C3 reference scores
  results/hcp_concordance.tsv       scores and weights vs C1-C3 (and DK PLS2)
  results/hcp_vs_dk_enrichment.tsv  MAGMA gene-property, all vectors, one universe
"""
from __future__ import annotations
import subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls

ROOT = HERE.parent
REPO = ROOT.parent
RES, DATA, REF = ROOT / "results", ROOT / "data", ROOT / "data" / "reference"
GS = REF / "gene_sets"
RUN_DIR = REPO / "out" / "thickness_hcp_70_fec93121f0dd"
MAGMA = REPO / "tools" / "bin" / ("magma_mac/magma" if (REPO / "tools/bin/magma_mac/magma").exists() else "magma")
RAW = {"SCZ": REPO / "hpc/work/results/magma/SCZ.genes.raw",
       "MDD": REPO / "hpc/work/results/magma/MDD.genes.raw"}
EXPR = Path.home() / "Git" / "AHBA" / "data" / "abagen-data" / "expression"
N_PERM, N_BOOT, SEED = 5000, 1000, 0
BAD_PARCELS = ["H"]          # medial-wall parcel, thickness 0 in many sessions

# ------------------------------------------------------------------ Y --------
fixed = pd.read_parquet(RUN_DIR / "fits" / "fixed.parquet")
terms = sorted(fixed.term.unique())
slope_term = next(t for t in terms if "age" in t.lower())
icpt_term = next(t for t in terms if "intercept" in t.lower())
wide = (fixed[fixed.term.isin([slope_term, icpt_term])]
        .pivot(index="label", columns="term", values="estimate")
        .rename(columns={slope_term: "dCT", icpt_term: "CT"}))
wide = wide[~wide.index.str.contains("global", case=False)]
wide["region"] = wide.index.str.replace(r"^(lh|rh)_", "", regex=True)
Y180 = wide.groupby("region")[["dCT", "CT"]].mean()
Y180 = Y180.drop(index=[r for r in BAD_PARCELS if r in Y180.index])
Y180.index = "lh_" + Y180.index
print(f"HCP Y: {Y180.shape[0]} bilateral parcels (dropped {BAD_PARCELS})", file=sys.stderr)
# the complete map (all parcels, not just AHBA-covered) for figures
Y180.to_csv(RES / "hcp_y_maps_180.csv", float_format="%.6g")

# ------------------------------------------------------------------ X --------
c123w = pd.read_csv(REPO / "data" / "weights.csv", index_col=0)
genes = list(c123w.index)
Xh = pd.read_csv(EXPR / "hcp_3d.csv", index_col=0)[genes]
c123s = pd.read_csv(REPO / "data" / "ahba_dme_hcp_top8kgenes_scores.csv")   # id -> label
id2label = dict(zip(c123s.id, "lh_" + c123s.label.astype(str)))
Xh.index = [id2label.get(i, f"id{i}") for i in Xh.index]
Xh = Xh.loc[[l for l in Xh.index if l in Y180.index]]
regions = list(Xh.index)
print(f"HCP X: {Xh.shape[0]} covered parcels x {Xh.shape[1]} genes", file=sys.stderr)
assert Xh.shape[0] >= 130, Xh.shape

# ---- gene- AND pipeline-matched DK comparator -------------------------------
# The DK weights in results/ come from the AHBA_updated native-DK matrices, a
# different abagen build from hcp_3d.csv, so comparing them with the HCP fit
# would confound parcellation with pipeline. This refits the same option-2 PLS
# on dk_3d.csv: same abagen build, same 7,973 genes, 33 regions.
FS_DK = ["bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus", "entorhinal",
         "fusiform", "inferiorparietal", "inferiortemporal", "isthmuscingulate", "lateraloccipital",
         "lateralorbitofrontal", "lingual", "medialorbitofrontal", "middletemporal",
         "parahippocampal", "paracentral", "parsopercularis", "parsorbitalis", "parstriangularis",
         "pericalcarine", "postcentral", "posteriorcingulate", "precentral", "precuneus",
         "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal", "superiorparietal",
         "superiortemporal", "supramarginal", "frontalpole", "temporalpole", "transversetemporal",
         "insula"]
Xd = pd.read_csv(EXPR / "dk_3d.csv", index_col=0)[genes]
Xd.index = ["lh_" + FS_DK[i - 1] for i in Xd.index]
# the id -> label mapping is not documented with the matrix; validate it against
# the shipped DK component scores before using it (C1 is unmistakable).
_ship = pd.read_csv(REPO / "data" / "ahba_dme_dsk_scores.csv").set_index("label")
_Xz = pls.zscore_cols(Xd.to_numpy(float)); _Xz -= _Xz.mean(0, keepdims=True)
_, _, _Vt = np.linalg.svd(_Xz, full_matrices=False)
_pc1 = pd.Series(_Xz @ _Vt[0], index=Xd.index)
_sh = Xd.index.intersection(_ship.index)
_r = abs(stats.spearmanr(_pc1.loc[_sh], _ship.loc[_sh, "C1"]).statistic)
assert _r > 0.95, f"DK id->label mapping looks wrong: |rho(PC1, C1)| = {_r:.3f}"
Ydk = pd.read_csv(DATA / "y_maps_bilateral_34.csv", index_col=0)[["dCT", "CT"]]
fit_dk = pls.pls_svd(Xd, Ydk.loc[Xd.index])
bt_dk = pls.bootstrap_weights(Xd, Ydk.loc[Xd.index], n_boot=N_BOOT, seed=SEED)
Z_dk_matched = -bt_dk["Z"]["PLS2"]
print(f"matched DK comparator: |rho(PC1,C1)| = {_r:.3f}, "
      f"PLS2 dCT salience {fit_dk.saliences().loc['dCT', 'PLS2']:.2f}", file=sys.stderr)

# ------------------------------------------------------------------ PLS ------
rows, wcols = [], {}
for opt, ycols in (("hcp_opt2_dCT_CT", ["dCT", "CT"]), ("hcp_opt1_dCT", ["dCT"])):
    Yf = Y180[ycols]
    sp = pls.spin_test_pls(Xh, Yf, n_perm=N_PERM, seed=SEED, centroids=pls.HCP_CENTROIDS)
    bt = pls.bootstrap_weights(Xh, Yf.loc[regions], n_boot=N_BOOT, seed=SEED)
    fit = sp["fit"]
    W, sal = fit.weights(), fit.saliences()
    for j, comp in enumerate(W.columns):
        wcols[f"{opt}_{comp}_U"] = W[comp]
        wcols[f"{opt}_{comp}_Z"] = bt["Z"][comp]
        r = dict(option=opt, component=comp, n_regions=len(regions), n_genes=Xh.shape[1],
                 singular_value=fit.S[j], cov_explained=fit.cov_explained[j],
                 p_spin_singular=sp["p_singular"][j],
                 r_gene_imaging_scores=fit.lx_ly_corr[j],
                 boot_reproducibility=bt["reproducibility"][comp],
                 n_genes_absZ_gt3=int((bt["Z"][comp].abs() > 3).sum()))
        for yc in ("dCT", "CT"):
            r[f"sal_{yc}"] = sal.loc[yc, comp] if yc in sal.index else np.nan
        rows.append(r)
    if opt == "hcp_opt2_dCT_CT":
        scores = fit.scores()
comp_tab = pd.DataFrame(rows)
comp_tab.to_csv(RES / "hcp_pls_components.tsv", sep="\t", index=False, float_format="%.5g")
Wall = pd.DataFrame(wcols); Wall.index.name = "gene"
Wall.to_csv(RES / "hcp_pls_weights.tsv", sep="\t", float_format="%.6g")

# ------------------------------------------------- concordance with C1-C3 ----
c123_hcp = c123s.assign(label="lh_" + c123s.label.astype(str)).set_index("label")[["C1", "C2", "C3"]]
lead_comp = comp_tab[comp_tab.option == "hcp_opt2_dCT_CT"].set_index("component")
# thinning orientation: flip so positive = expressed where thinning is faster
Z_lead = -Wall["hcp_opt2_dCT_CT_PLS2_Z"]
S_lead = -scores["PLS2_gene"] if "PLS2_gene" in scores.columns else -scores.filter(like="PLS2").iloc[:, 0]

conc = []
for ref in ("C1", "C2", "C3"):
    rm = c123_hcp[ref].reindex(Y180.index)          # complete-map reference for the spin
    rho, p, _ = pls.spin_corr(S_lead, rm.dropna(), n_perm=N_PERM, centroids=pls.HCP_CENTROIDS)
    shared = Z_lead.index.intersection(c123w[ref].dropna().index)
    conc.append(dict(reference=ref, level="scores", rho=rho, p_spin=p, n=len(rm.dropna())))
    conc.append(dict(reference=ref, level="weights",
                     rho=stats.spearmanr(Z_lead.loc[shared], c123w[ref].loc[shared]).statistic,
                     p_spin=np.nan, n=len(shared)))
dk_lead = pd.read_csv(RES / "lead_signature_weights.tsv", sep="\t", index_col=0)["thinning_Z_ds25"]
shared = Z_lead.index.intersection(dk_lead.dropna().index)
conc.append(dict(reference="DK_PLS2_ds25", level="weights",
                 rho=stats.spearmanr(Z_lead.loc[shared], dk_lead.loc[shared]).statistic,
                 p_spin=np.nan, n=len(shared)))
C = pd.DataFrame(conc)
C.to_csv(RES / "hcp_concordance.tsv", sep="\t", index=False, float_format="%.4g")

out_scores = pd.DataFrame({"thinning_score": S_lead}).join(c123_hcp).join(Y180)
out_scores.to_csv(RES / "hcp_pls_scores.csv", float_format="%.6g")

# ------------------------------------------------------------- enrichment ----
sym2ent = (pd.read_csv(GS / "magma_SCZ_genes.tsv", sep="\t").dropna(subset=["symbol"])
           .drop_duplicates("symbol").set_index("symbol")["GENE"].astype(int))
dk_thin = {f"DK_PLS2_{ds}": pd.read_csv(RES / "lead_signature_weights.tsv", sep="\t",
                                        index_col=0)[f"thinning_Z_{ds}"]
           for ds in ("ds0", "ds25", "ds50")}
vectors = {"HCP_PLS2_thinning": Z_lead,
           "HCP_PLS1_dCTonly": -Wall["hcp_opt1_dCT_PLS1_Z"],
           "DK_PLS2_matchedX": Z_dk_matched,
           **dk_thin,
           "C3_shipped": c123w["C3"], "C1_shipped": c123w["C1"]}
cov = pd.DataFrame(vectors).dropna(how="any")
cov = cov.loc[cov.index.intersection(sym2ent.index)]
cov.index = sym2ent.loc[cov.index].values
cov = cov[~cov.index.duplicated()]
cov.index.name = "GENE"
path = RES / "magma_runs" / "hcp_vs_dk.covar"
cov.to_csv(path, sep="\t", float_format="%.6g")
print(f"shared universe: {len(cov)} genes", file=sys.stderr)

res = []
for dis, raw in RAW.items():
    prefix = RES / "magma_runs" / f"hcp_vs_dk_{dis}"
    r = subprocess.run([str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(path),
                        "--out", str(prefix)], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"MAGMA failed ({dis}):\n{r.stdout[-1200:]}")
    t = pd.read_csv(f"{prefix}.gsa.out", sep=r"\s+", comment="#")
    if "FULL_NAME" in t.columns:
        t["VARIABLE"] = t["FULL_NAME"].fillna(t["VARIABLE"])
    t.insert(0, "disorder", dis)
    res.append(t)
M = pd.concat(res, ignore_index=True)
M["se_std"] = M.SE * M.BETA_STD / M.BETA
M.to_csv(RES / "hcp_vs_dk_enrichment.tsv", sep="\t", index=False, float_format="%.5g")

pd.set_option("display.width", 200)
print("\nPLS components (HCP-MMP, 137 parcels):")
print(comp_tab[["option", "component", "cov_explained", "p_spin_singular", "sal_dCT", "sal_CT",
                "boot_reproducibility", "n_genes_absZ_gt3"]].round(3).to_string(index=False))
print("\nConcordance of the HCP thinning component with the AHBA components:")
print(C.round(3).to_string(index=False))
print(f"\nMAGMA gene-property, shared universe of {int(M.NGENES.iloc[0])} genes:")
order = ["HCP_PLS2_thinning", "DK_PLS2_matchedX", "DK_PLS2_ds0", "DK_PLS2_ds25", "DK_PLS2_ds50",
         "HCP_PLS1_dCTonly", "C3_shipped", "C1_shipped"]
print(M[M.VARIABLE.isin(order)].pivot(index="VARIABLE", columns="disorder",
                                      values=["BETA_STD", "se_std", "P"]).loc[order].round(4).to_string())

# ------------------------------------------- does HCP add beyond DK? --------
# Comparing betas across runs ignores that the two weight vectors are
# correlated (rho ~ 0.73), so the difference is tested inside MAGMA instead:
# each vector conditioned on the other, and both conditioned on C3.
cond_cov = cov[["HCP_PLS2_thinning", "DK_PLS2_matchedX", "C3_shipped"]]
cpath = RES / "magma_runs" / "hcp_vs_dk_cond.covar"
cond_cov.to_csv(cpath, sep="\t", float_format="%.6g")
cres = []
for tag, model in (("none", []),
                   ("cond_DK", ["--model", "condition=DK_PLS2_matchedX"]),
                   ("cond_HCP", ["--model", "condition=HCP_PLS2_thinning"]),
                   ("cond_C3", ["--model", "condition=C3_shipped"])):
    for dis, raw in RAW.items():
        prefix = RES / "magma_runs" / f"hcp_cond_{tag}_{dis}"
        r = subprocess.run([str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(cpath),
                            "--out", str(prefix), *model], capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"MAGMA failed ({tag} {dis}):\n{r.stdout[-1200:]}")
        t = pd.read_csv(f"{prefix}.gsa.out", sep=r"\s+", comment="#")
        if "FULL_NAME" in t.columns:
            t["VARIABLE"] = t["FULL_NAME"].fillna(t["VARIABLE"])
        t.insert(0, "disorder", dis); t.insert(0, "run", tag)
        cres.append(t)
CD = pd.concat(cres, ignore_index=True)
CD["se_std"] = CD.SE * CD.BETA_STD / CD.BETA
CD.to_csv(RES / "hcp_vs_dk_conditional.tsv", sep="\t", index=False, float_format="%.5g")
keep = CD[(CD.run == "none") | (CD.run.str.replace("cond_", "").map(
    {"DK": "DK_PLS2_matchedX", "HCP": "HCP_PLS2_thinning", "C3": "C3_shipped"}) != CD.VARIABLE)]
print("\nConditional models (same universe):")
print(keep[["run", "disorder", "VARIABLE", "BETA_STD", "se_std", "P"]].round(4).to_string(index=False))
