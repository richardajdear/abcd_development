"""
04_fit_pls.py -- fit PLS-SVD for the five Y-matrix options at ds0/25/50.

Outputs (ahba_pls/results/):
  pls_components.tsv          one row per option x DS x component: singular value,
                              covariance explained, spin p (5000), gene-side vs
                              imaging-side score correlation, bootstrap
                              reproducibility, saliences of each Y column
  pls_weights/<opt>_<ds>.tsv  gene, U (weight), Z (bootstrap), boot_sd per component
  pls_scores/<opt>_<ds>.csv   region scores (gene side and imaging side)
  pls_null/<opt>_<ds>.npz     null singular values (for the notebook)

Sign convention (pls.py): first Y column = dCT salience positive, so positive
gene weight <=> higher expression where the thinning RATE is more positive
(i.e. LESS thinning).  Flip for presentation as "thinning-associated".
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pls

ROOT = HERE.parent
DATA, RES = ROOT / "data", ROOT / "results"
for d in ("pls_weights", "pls_scores", "pls_null"):
    (RES / d).mkdir(parents=True, exist_ok=True)

OPTIONS = {
    "opt1_dCT":          ["dCT"],
    "opt2_dCT_CT":       ["dCT", "CT"],
    "opt3_dCT_dT1T2":    ["dCT", "dT1T2"],
    "opt4_full4":        ["dCT", "CT", "dT1T2", "T1T2"],
    "opt5_slopePCs":     ["slopePC1", "slopePC2", "slopePC3"],
}
DS = ["ds0", "ds25", "ds50"]
N_PERM, N_BOOT, SEED = 5000, 1000, 0

Y34 = pd.read_csv(DATA / "y_maps_bilateral_34.csv", index_col=0)
regions33 = [l.strip() for l in open(DATA / "region_order_33.txt") if l.strip()]
assert len(regions33) == 33 and set(regions33) <= set(Y34.index)

rows = []
t0 = time.time()
for ds in DS:
    X = pd.read_parquet(DATA / f"X_{ds}.parquet").loc[regions33]
    assert list(X.index) == regions33
    for opt, ycols in OPTIONS.items():
        Yf = Y34[ycols]
        Y = Yf.loc[regions33]
        sp = pls.spin_test_pls(X, Yf, n_perm=N_PERM, seed=SEED)
        bt = pls.bootstrap_weights(X, Y, n_boot=N_BOOT, seed=SEED)
        fit = sp["fit"]
        k = len(fit.S)
        # weights
        W = fit.weights()
        out = pd.DataFrame(index=W.index)
        for j in range(k):
            c = f"PLS{j+1}"
            out[f"{c}_U"] = W[c]
            out[f"{c}_Z"] = bt["Z"][c]
            out[f"{c}_boot_sd"] = bt["U_boot_sd"][c]
        out.index.name = "gene"
        out.to_csv(RES / "pls_weights" / f"{opt}_{ds}.tsv", sep="\t", float_format="%.6g")
        fit.scores().to_csv(RES / "pls_scores" / f"{opt}_{ds}.csv", float_format="%.6g")
        np.savez_compressed(RES / "pls_null" / f"{opt}_{ds}.npz",
                            null_S=sp["null_S"], null_cov=sp["null_cov"], S=fit.S)
        sal = fit.saliences()
        for j in range(k):
            c = f"PLS{j+1}"
            r = dict(option=opt, ds=ds, n_genes=X.shape[1], component=c,
                     singular_value=fit.S[j], cov_explained=fit.cov_explained[j],
                     p_spin_singular=sp["p_singular"][j],
                     p_spin_cov_explained=sp["p_cov_explained"][j],
                     r_gene_imaging_scores=fit.lx_ly_corr[j],
                     boot_reproducibility=bt["reproducibility"][c],
                     n_genes_absZ_gt3=int((out[f"{c}_Z"].abs() > 3).sum()))
            for yc in Y34.columns:
                r[f"sal_{yc}"] = sal.loc[yc, c] if yc in sal.index else np.nan
                r[f"salZ_{yc}"] = bt["V_boot_Z"].loc[yc, c] if yc in sal.index else np.nan
            rows.append(r)
        print(f"{ds} {opt}: {time.time()-t0:.0f}s", file=sys.stderr)

comp = pd.DataFrame(rows)
comp.to_csv(RES / "pls_components.tsv", sep="\t", index=False, float_format="%.5g")
print(comp[["option", "ds", "component", "cov_explained", "p_spin_singular",
            "r_gene_imaging_scores", "boot_reproducibility", "n_genes_absZ_gt3"]]
      .to_string(index=False))
