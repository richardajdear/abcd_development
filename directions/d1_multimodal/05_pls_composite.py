"""D1 step 5a: a diffusion predictor of the thinning slope (PLS1, cross-fitted), and per-child diffusion PC scores.

    python directions/d1_multimodal/05_pls_composite.py        (env abcd-spatial)

1. Child slopes (03 BLUPs) of the 45 DTI/RSI measures predict the faster-thinning slope (-HCP thickness BLUP).
   Training uses site-residualised, z-scored slopes so PLS cannot learn the scanner drift. Folds are 5
   GroupKFold splits by family, so no child is predicted from weights that saw them or a relative. The number of
   components is chosen by inner 5-fold CV inside each training split. Output: out-of-fold prediction pls_oof
   (child level) and the fold-specific weights.
2. Per-scan composite: for each child, comp = sum_k w_k * x_scan,k / s_k with that child's fold weights (s_k = SD of
   the site-residualised child slopes used for scaling). The LMM slope of comp is therefore the PLS prediction
   built from scan-level data, so its slope reliability can be estimated like any other measure (06_*.R).
   Scans: those present in both the dMRI sample and the thickness (Figure-1) sample.
3. PCA (as 04_followups.R): site-residualised z-scored slopes of all 47 measures -> PC1..PC5 scores per child.
   Signs follow prcomp/sklearn convention and are arbitrary.

Writes (gitignored) work/pls_child.csv.gz, work/pls_scan.csv.gz; (aggregate) tables/pls_fit.tsv, pls_weights.tsv.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
W, TAB = HERE / "work", HERE / "tables"
D3W = HERE.parents[0] / "d3_twin_family" / "work"
rng = np.random.default_rng(20261027)

B = pd.read_csv(W / "blups_dmri.csv.gz")
S = pd.read_csv(W / "scan_means_dmri.csv.gz")
fam = S.drop_duplicates("subject").set_index("subject").family_id
site = B.drop_duplicates("subject").set_index("subject").site
ws = B.pivot(index="subject", columns="metric", values="slope")
DM = [m for m in ws.columns if m.startswith(("dti__", "rsi__"))]
ws = ws.dropna(subset=DM + ["smri__thk_hcp"])
y_raw = -ws["smri__thk_hcp"]


def site_resid(df: pd.DataFrame) -> pd.DataFrame:
    return df - df.groupby(site.reindex(df.index).values).transform("mean")


Xr = site_resid(ws[DM]); yr = site_resid(y_raw.to_frame())["smri__thk_hcp"]
sx = Xr.std(); Xz = (Xr - Xr.mean()) / sx; yz = (yr - yr.mean()) / yr.std()
groups = fam.reindex(ws.index).fillna(-1).astype(str).values
print(f"children with all 45 diffusion slopes + thickness: {len(ws)}; families {len(set(groups))}")


def fit_pls(X, y, g):
    """Choose n_components by inner GroupKFold CV R2, refit on all of X."""
    best, best_r2 = 1, -np.inf
    for nc in range(1, 11):
        pred = np.zeros(len(y))
        for tr, te in GroupKFold(5).split(X, y, g):
            m = PLSRegression(n_components=nc, scale=False).fit(X[tr], y[tr])
            pred[te] = m.predict(X[te]).ravel()
        r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        if r2 > best_r2:
            best, best_r2 = nc, r2
    m = PLSRegression(n_components=best, scale=False).fit(X, y)
    # predict() is affine in X; recover its exact weights (sklearn's coef_ convention varies by version)
    A = np.column_stack([np.ones(len(X)), X]); p = m.predict(X).ravel()
    sol, *_ = np.linalg.lstsq(A, p, rcond=None)
    assert np.max(np.abs(A @ sol - p)) < 1e-8, "PLS prediction is not reproduced by the recovered weights"
    m.w_eff_ = sol[1:]
    return m, best, best_r2


X, y = Xz.to_numpy(), yz.to_numpy()
oof = np.zeros(len(y)); fold = np.zeros(len(y), int); WTS, FIT = [], []
for f, (tr, te) in enumerate(GroupKFold(5).split(X, y, groups)):
    m, nc, r2_in = fit_pls(X[tr], y[tr], groups[tr])
    oof[te] = m.predict(X[te]).ravel(); fold[te] = f
    WTS.append(pd.Series(m.w_eff_, index=DM, name=f"fold{f}"))
    FIT.append(dict(fold=f, n_components=nc, inner_cv_r2=r2_in,
                    test_r=np.corrcoef(oof[te], y[te])[0, 1], n_test=len(te)))
    print(f"fold {f}: {nc} components, inner R2 {r2_in:.3f}, test r {FIT[-1]['test_r']:.3f}", flush=True)
r_oof = np.corrcoef(oof, y)[0, 1]
FIT.append(dict(fold="all (out-of-fold)", n_components=np.nan, inner_cv_r2=np.nan, test_r=r_oof, n_test=len(y)))
pd.DataFrame(FIT).to_csv(TAB / "pls_fit.tsv", sep="\t", index=False)
Wd = pd.concat(WTS, axis=1); Wd["mean"] = Wd.mean(axis=1); Wd["sd_over_folds"] = Wd.iloc[:, :5].std(axis=1)
Wd.sort_values("mean", key=np.abs, ascending=False).to_csv(TAB / "pls_weights.tsv", sep="\t")
print(f"out-of-fold r(pred, faster thinning) = {r_oof:.3f}  (R2 {r_oof**2:.3f})")
print("largest mean weights:\n", Wd["mean"].sort_values(key=np.abs, ascending=False).head(8).round(3).to_string())

# ---- per-scan composite with each child's own fold weights ---------------------------------------
fold_of = pd.Series(fold, index=ws.index)
sw = S[S.subject.isin(ws.index)].pivot_table(index=["subject", "visit"], columns="metric", values="mean_val")[DM].dropna()
coef = np.vstack([(WTS[f] / sx).to_numpy() for f in range(5)])        # per-fold weights on raw scan units
fi = fold_of.reindex(sw.index.get_level_values("subject")).to_numpy()
comp = np.einsum("ij,ij->i", sw.to_numpy(), coef[fi])
cs = pd.DataFrame({"comp": comp}, index=sw.index).reset_index()
ct = pd.read_csv(D3W / "scan_means_ct.csv")[["subject", "visit", "mean_val", "age_c", "sex", "site", "family_id"]]
cs = cs.merge(ct.rename(columns={"mean_val": "thk"}), on=["subject", "visit"], how="inner")
cs["n_visits"] = cs.groupby("subject").visit.transform("size")
cs = cs[cs.n_visits >= 2]
cs.to_csv(W / "pls_scan.csv.gz", index=False)
print(f"per-scan composite: {len(cs)} scans, {cs.subject.nunique()} children (in both dMRI and thickness samples, >= 2 scans)")

# ---- PCA scores ------------------------------------------------------------------------------------
ALL = DM + ["smri__thk_hcp", "smri__t1t2__gm"]
wa = B.pivot(index="subject", columns="metric", values="slope")[ALL].dropna()
Za = site_resid(wa); Za = (Za - Za.mean()) / Za.std()
pca = PCA(5).fit(Za.to_numpy())
pcs = pd.DataFrame(pca.transform(Za.to_numpy()), index=Za.index, columns=[f"PC{i}" for i in range(1, 6)])
print("PCA variance explained:", np.round(pca.explained_variance_ratio_, 3))
load = pd.DataFrame(pca.components_.T, index=ALL, columns=pcs.columns)
print("thickness loadings:", load.loc["smri__thk_hcp"].round(3).to_dict())
for k in ["PC1", "PC2"]:
    print(k, "top:", load[k].sort_values(key=np.abs, ascending=False).head(6).round(2).to_dict())
child = pd.DataFrame({"pls_oof": oof, "thin_fast_resid": y}, index=ws.index).join(pcs, how="outer")
child.rename_axis("subject").reset_index().to_csv(W / "pls_child.csv.gz", index=False)
load.to_csv(TAB / "pca5_loadings.tsv", sep="\t")
