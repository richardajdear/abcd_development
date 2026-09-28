"""Build per-scan transcriptional-axis scores for a single-LMM fit (C3-axis phenotypes).

Why
---
The whole-cortex slope phenotype carries a global factor (every parcel's slope
covaries positively). Subtracting each child's cortex-wide mean ("row-centring")
and decomposing the rest gives C1, C2 and then C3 as separate axes of
between-child variation. On the left hemisphere C3 is component 5 (rho with
AHBA C3 = 0.52-0.54, spin p < 0.002); diffusion map embedding gives the same
component (|r| >= 0.996). See genetic_analysis/README_HPC.md section 8.

What
----
1. Reads the per-parcel slope BLUPs of the settled HCP-MMP run and fits PCA on
   the row-centred, standardised LEFT-hemisphere slopes (179 parcels; AHBA
   C1-C3 are left-hemisphere maps). The components that best match AHBA C1, C2
   and C3 become weight vectors, signed to correlate positively with their map.
2. Defines map projections on C3, PLS2 and the normative dCT map.
3. Turns every weight vector w (on standardised slopes) into per-scan weights
   on raw thickness, w_r / sd_r, with w mean-centred (the mean-centring IS the
   row-centring: sum_r (z_r - zbar) v_r = sum_r z_r (v_r - vbar)). The same
   weight goes on each lh/rh homologue, so both hemispheres contribute.
4. Writes out/<run>_c3axis/model_table.parquet (one label per phenotype, per
   scan) plus config.yaml, ready for `Rscript R/fit_lmm.R --run-dir ...`, exactly
   as the single-LMM `global_slope` was built (README section 5.2, step 1b).
   The `global_mean` label reproduces that primary trait as a check.
5. Writes a committed, map-level validation table (no per-child data):
   genetic_analysis/c3axis/c3axis_weights_check.tsv.

Usage (repo root):
    PYTHONPATH=src python genetic_analysis/c3axis/build_c3axis_model_table.py \
        [--run out/thickness_hcp_70_aa6e91efba82]
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
AP = REPO / "ahba_pls"


def load_maps(labels_lh: list[str]) -> pd.DataFrame:
    """AHBA C1-C3 and PLS1/PLS2 (hcp_3d_ds5); map-level only, indexed by lh_ label."""
    ylab = {l.lower(): l for l in labels_lh}
    ul = pd.read_csv(AP / "data/reference/hcp_cortices/HCP-MMP1_UniqueRegionList.csv", encoding="utf-8-sig")
    id2 = {int(r.regionID): ylab.get(f"lh_{r.region}".lower()) for r in ul[ul.LR == "L"].itertuples()}
    c = pd.read_csv(REPO / "data/ahba_dme_hcp_top8kgenes_scores.csv")
    c = c[c.id.isin(id2)].assign(label=lambda d: d.id.map(id2)).dropna(subset=["label"]).set_index("label")
    pls_ = pd.read_csv(AP / "results/hcp_summary_maps.csv", index_col=0)[["PLS1", "PLS2"]]
    m = c[["C1", "C2", "C3"]].join(pls_, how="left")
    assert len(m) >= 130, len(m)
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="out/thickness_hcp_70_aa6e91efba82")
    ap.add_argument("--k", type=int, default=10)
    a = ap.parse_args()
    run = REPO / a.run
    out = run.parent / f"{run.name}_c3axis"

    # ---- slopes, LH row-centred PCA ------------------------------------------
    bl = pd.read_parquet(run / "fits/blups.parquet")
    W = bl.pivot(index="subject", columns="label", values="re_slope")
    assert W.notna().all().all()
    sd = W.std()
    lh = [c for c in W.columns if c.startswith("lh_")]
    Z = (W[lh] - W[lh].mean()) / sd[lh]
    Zc = Z.sub(Z.mean(1), axis=0)                                   # row-centring
    A = Zc.to_numpy() - Zc.to_numpy().mean(0)
    _, S, Vt = np.linalg.svd(A, full_matrices=False)
    V = pd.DataFrame(Vt[: a.k].T, index=lh, columns=[f"rc{i+1}" for i in range(a.k)])
    ev = S[: a.k] ** 2 / (S ** 2).sum()

    maps = load_maps(lh)
    fx = pd.read_parquet(run / "fits/fixed.parquet").pivot(index="label", columns="term", values="estimate")
    dct_lh = fx.loc[lh, "age_c"]
    maps["dCT"] = dct_lh.reindex(maps.index)
    sp = lambda x, y: stats.spearmanr(x, y, nan_policy="omit").statistic
    MAPCOLS = ("C1", "C2", "C3", "PLS2", "dCT")

    # ---- pick the component matching each AHBA axis --------------------------
    rows, weights = [], {}
    for target in ("C1", "C2", "C3"):
        r = pd.Series({c: sp(V[c].reindex(maps.index), maps[target]) for c in V})
        best = r.abs().idxmax()
        v = V[best] * np.sign(r[best])
        weights[f"{target.lower()}axis_rc"] = v
        rows.append(dict(phenotype=f"{target.lower()}axis_rc", source=best, var_explained=ev[int(best[2:]) - 1],
                         **{f"rho_{m}": sp(v.reindex(maps.index), maps[m]) for m in MAPCOLS}))
    c3 = rows[-1]
    if c3["source"] != "rc5" or c3["rho_C3"] < 0.45:
        print(f"WARNING: C3 component is {c3['source']} (rho {c3['rho_C3']:.2f}); expected rc5, >= 0.45. "
              "Inspect before running the genetics.", file=sys.stderr)

    # ---- map projections (weights = z-scored map over covered parcels) --------
    for name, m in (("proj_C3", maps.C3), ("proj_PLS2", maps.PLS2), ("proj_dCT", dct_lh)):
        m = m.dropna(); weights[name] = (m - m.mean()) / m.std()
        rows.append(dict(phenotype=name, source="map", var_explained=np.nan,
                         **{f"rho_{k}": sp(weights[name].reindex(maps.index), maps[k]) for k in MAPCOLS}))
    weights["global_mean"] = None
    pd.DataFrame(rows).to_csv(REPO / "genetic_analysis/c3axis/c3axis_weights_check.tsv", sep="\t",
                              index=False, float_format="%.4f")

    # ---- per-scan scores -----------------------------------------------------
    mt = pd.read_parquet(run / "model_table.parquet")
    T = mt.pivot_table(index=["subject", "visit"], columns="label", values="value")
    T = T.loc[T.notna().all(axis=1)]
    key = lambda l: l.split("_", 1)[1].lower()
    scores = {}
    for name, w in weights.items():
        if w is None:                                             # the primary single-LMM trait
            scores[name] = T.mean(axis=1); continue
        wk = {key(l): v for l, v in (w - w.mean()).items()}       # mean-centred weights = row-centring
        cols = [c for c in T.columns if key(c) in wk]
        full = np.array([wk[key(c)] / sd[c] for c in cols])      # same weight on each lh/rh homologue
        scores[name] = pd.Series(T[cols].to_numpy() @ full / len(cols), index=T.index)
    Sc = pd.DataFrame(scores)
    Sc = Sc / Sc.xs("v0", level="visit").std()                   # unit baseline SD, for the LMM
    cov_cols = ["subject", "visit", "age", "age_c", "sex", "site", "family_id",
                "age_first", "age_last", "age_span", "n_visits"]
    covs = mt.drop_duplicates(["subject", "visit"])[cov_cols]
    long = (Sc.rename_axis(columns="label").stack().rename("value").reset_index()
              .merge(covs, on=["subject", "visit"], how="left"))
    long = long.assign(metric="c3axis", release="7.0", hemi="both", region=long.label)
    assert long[["age_c", "sex", "site"]].notna().all().all()
    out.mkdir(exist_ok=True)
    long.to_parquet(out / "model_table.parquet", index=False)
    shutil.copy(run / "config.yaml", out / "config.yaml")
    print(f"wrote {out.relative_to(REPO)}/model_table.parquet: {len(long)} rows, {long.subject.nunique()} "
          f"subjects, phenotypes {sorted(long.label.unique())}")


if __name__ == "__main__":
    main()
