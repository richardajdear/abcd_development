"""
05_concordance.py -- H1: do ABCD PLS components match NSPN-PLS2 and AHBA-C3?

Two levels, kept separate:
  (a) spatial scores: Spearman rho of each component's gene-side region scores
      (X @ U, 33 regions) with the reference DK map, spin p (5000 rotations of
      the more complete reference map, then aligned to the 33 covered regions).
      References: NSPN PLS2 (34 bilateral), NSPN CT_delta/MT_delta/CT/MT maps,
      AHBA C1-C3 recomputed on the SAME DS matrix (33), shipped C1-C3 (34).
  (b) gene weights: Spearman rho of the bootstrap-Z gene weights with NSPN
      PLS1_z/PLS2_z and AHBA C1/C2/C3 weights on the shared genes; p by
      10,000 gene permutations (reported for completeness -- with thousands of
      genes any |rho| > ~0.05 is significant, so the effect size is what matters).

Also: (c) the raw Y maps vs the NSPN maps (is ABCD's thinning map the NSPN
thinning map?), and (d) the cross-option / cross-DS stability matrix of the
dCT-carrying components' gene weights.

Outputs (results/): concordance_scores.tsv, concordance_weights.tsv,
ymaps_vs_nspn.tsv, weight_stability_dCT_components.csv
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
DATA, RES, REF = ROOT / "data", ROOT / "results", ROOT / "data" / "reference"
N_SPIN, N_GPERM, SEED = 5000, 10000, 0

comp = pd.read_csv(RES / "pls_components.tsv", sep="\t")
Y34 = pd.read_csv(DATA / "y_maps_bilateral_34.csv", index_col=0)
regions33 = [l.strip() for l in open(DATA / "region_order_33.txt") if l.strip()]

nspn34 = pd.read_csv(REF / "nspn_dk_maps_bilateral_34.csv", index_col=0)
c123_ship = pd.read_csv(REF / "ahba_c123_dk_scores.csv", index_col=0)
c123_rec = {ds: pd.read_csv(REF / f"ahba_c123_scores_recomputed_{ds}.csv", index_col=0)
            for ds in ("ds0", "ds25", "ds50")}
nspn_w = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")
c123_w = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)

ref_maps_fixed = {
    "NSPN_PLS2": nspn34["PLS2"], "NSPN_CT_delta": nspn34["CT_delta"],
    "NSPN_MT_delta": nspn34["MT_delta"], "NSPN_CT": nspn34["CT"], "NSPN_MT": nspn34["MT"],
    "C1_shipped": c123_ship["C1"], "C2_shipped": c123_ship["C2"], "C3_shipped": c123_ship["C3"],
}
ref_weights = {"NSPN_PLS1_z": nspn_w["PLS1_z"], "NSPN_PLS2_z": nspn_w["PLS2_z"],
               "C1": c123_w["C1"], "C2": c123_w["C2"], "C3": c123_w["C3"]}


def gene_perm_p(a: np.ndarray, b: np.ndarray, n: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    r_obs = stats.spearmanr(a, b).statistic
    ra = stats.rankdata(a); rb = stats.rankdata(b)
    null = np.array([np.corrcoef(rng.permutation(ra), rb)[0, 1] for _ in range(n)])
    return (1 + (np.abs(null) >= abs(r_obs)).sum()) / (1 + n)


# ---- (c) raw Y maps vs NSPN maps -----------------------------------------
rows = []
for yc in Y34.columns:
    for rn in ("NSPN_CT_delta", "NSPN_MT_delta", "NSPN_CT", "NSPN_MT", "NSPN_PLS2",
               "C1_shipped", "C2_shipped", "C3_shipped"):
        r, p, _ = pls.spin_corr(Y34[yc], ref_maps_fixed[rn], n_perm=N_SPIN, seed=SEED)
        rows.append(dict(y_map=yc, reference=rn, rho=r, p_spin=p, n_regions=34))
pd.DataFrame(rows).to_csv(RES / "ymaps_vs_nspn.tsv", sep="\t", index=False, float_format="%.4g")

# ---- (a) scores and (b) weights per component ----------------------------
srows, wrows, dct_Z = [], [], {}
for _, c in comp.iterrows():
    opt, ds, k = c.option, c.ds, c.component
    sc = pd.read_csv(RES / "pls_scores" / f"{opt}_{ds}.csv", index_col=0)
    W = pd.read_csv(RES / "pls_weights" / f"{opt}_{ds}.tsv", sep="\t", index_col=0)
    Lx = sc[f"{k}_gene_scores"]
    Z = W[f"{k}_Z"].dropna()
    refs = dict(ref_maps_fixed)
    for j in ("C1", "C2", "C3"):
        refs[f"{j}_recomputed_{ds}"] = c123_rec[ds][j]
    for rn, rm in refs.items():
        r, p, _ = pls.spin_corr(Lx, rm, n_perm=N_SPIN, seed=SEED)
        srows.append(dict(option=opt, ds=ds, component=k, reference=rn, rho=r, p_spin=p,
                          n_regions=len(Lx.index.intersection(rm.dropna().index)),
                          sal_dCT=c.get("sal_dCT", np.nan), p_spin_component=c.p_spin_singular))
    for rn, rw in ref_weights.items():
        shared = Z.index.intersection(rw.dropna().index)
        a, b = Z.loc[shared].to_numpy(), rw.loc[shared].to_numpy()
        rho = stats.spearmanr(a, b).statistic
        wrows.append(dict(option=opt, ds=ds, component=k, reference=rn, rho=rho,
                          p_gene_perm=gene_perm_p(a, b, N_GPERM, SEED), n_genes=len(shared),
                          sal_dCT=c.get("sal_dCT", np.nan), p_spin_component=c.p_spin_singular))
    # collect the dCT-carrying component per option (largest |sal_dCT| among components)
    if opt != "opt5_slopePCs":
        sub = comp[(comp.option == opt) & (comp.ds == ds)]
        if k == sub.loc[sub.sal_dCT.abs().idxmax(), "component"]:
            dct_Z[f"{opt}|{ds}|{k}"] = Z
    print(opt, ds, k, file=sys.stderr)

S = pd.DataFrame(srows); Wt = pd.DataFrame(wrows)
S.to_csv(RES / "concordance_scores.tsv", sep="\t", index=False, float_format="%.4g")
Wt.to_csv(RES / "concordance_weights.tsv", sep="\t", index=False, float_format="%.4g")

# ---- (d) stability of dCT-carrying components across options / DS -------
allZ = pd.DataFrame(dct_Z)
stab = allZ.corr(method="spearman", min_periods=1000)
stab.to_csv(RES / "weight_stability_dCT_components.csv", float_format="%.3f")

# ---- console summary ------------------------------------------------------
pd.set_option("display.width", 200)
key = S[S.reference.isin(["NSPN_PLS2", "C3_recomputed_" + "ds25", "C3_shipped", "NSPN_CT_delta"])]
print("SCORES (ds25 shown for recomputed C3):")
print(key[key.ds == "ds25"][["option", "component", "sal_dCT", "reference", "rho", "p_spin"]]
      .round(3).to_string(index=False))
print("\nWEIGHTS (ds25):")
print(Wt[(Wt.ds == "ds25") & Wt.reference.isin(["NSPN_PLS2_z", "C3", "C1"])]
      [["option", "component", "sal_dCT", "reference", "rho", "n_genes"]].round(3).to_string(index=False))
print("\nY maps vs NSPN:")
ym = pd.read_csv(RES / "ymaps_vs_nspn.tsv", sep="\t")
print(ym[ym.reference.isin(["NSPN_CT_delta", "NSPN_CT", "NSPN_PLS2", "C3_shipped"])]
      .pivot(index="y_map", columns="reference", values="rho").round(2).to_string())
print("\nstability of dCT components (Spearman of Z):")
print(stab.round(2).to_string())
