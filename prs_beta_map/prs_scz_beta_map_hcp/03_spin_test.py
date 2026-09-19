"""Step 3 (laptop): spin-test the SCZ-PRS beta map against AHBA C3 and the
ahba_pls lead component, on HCP-MMP.

For every results/beta_map_<cell>.tsv (358 parcel rows + cortex_mean):

  * lh beta map (179 parcels; H is excluded by the phenotype config) is the
    map that is ROTATED -- 5,000 spins on the HCP spherical centroids with a
    bijective parcel assignment (ahba_pls/code/pls.py::spin_corr, the same
    routine and seed as every spin p in ahba_pls/).  The transcriptomic maps
    (137 AHBA-covered lh parcels) are the fixed side, so the rotation is done
    on the complete map and uncovered parcels are dropped afterwards.
  * the bilateral-mean beta (lh and rh averaged per Glasser area) is spun the
    same way as the lower-noise variant.
  * context rows: beta vs the group thinning map (dCT, all 179 parcels), and
    lh vs rh beta agreement (a reliability readout for the map itself).
  * the same tests for beta_cond (beta with the whole-cortex mean slope held
    fixed: the parcel-specific part) and for r_global (each parcel's loading on
    the cortex mean).  A score that only shifts the whole cortex still draws a
    structured beta map -- the loading map -- so beta_cond vs C3 is the test of
    whether the map carries spatial information of its own, and r_global vs C3
    says how much the loading map alone already resembles C3.

References (both lh, 137 parcels, Glasser names):
  C3        data/ahba_dme_hcp_top8kgenes_scores.csv          (Dear et al. 2024)
  PLS lead  ahba_pls/results/hcp_pls_scores.csv[thinning_score]  (ahba_pls, 12_hcp_pls.py;
            oriented so positive = expressed where thinning is faster)

Writes results/spin_tests.tsv (one row per cell x map variant x reference) and
work/spin_null_<cell>.npz (the null rho vectors, for the figure).

Usage (repo root):  python prs_beta_map/prs_scz_beta_map_hcp/03_spin_test.py [--results DIR] [--n-perm 5000]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO / "ahba_pls" / "code"))
import pls  # noqa: E402  (spin_corr, HCP_CENTROIDS)

C3_FILE = REPO / "data" / "ahba_dme_hcp_top8kgenes_scores.csv"
PLS_FILE = REPO / "ahba_pls" / "results" / "hcp_pls_scores.csv"
DCT_FILE = REPO / "ahba_pls" / "results" / "hcp_y_maps_180.csv"


def references() -> dict[str, pd.Series]:
    c = pd.read_csv(C3_FILE)
    c3 = pd.Series(c.C3.to_numpy(), index="lh_" + c.label.astype(str), name="C3")
    p = pd.read_csv(PLS_FILE, index_col=0)
    lead = p["thinning_score"].rename("PLS_lead")
    dct = pd.read_csv(DCT_FILE, index_col=0)["dCT"].rename("dCT")
    return {"C3": c3, "PLS_lead": lead, "dCT": dct}


def split_map(bm: pd.DataFrame, col: str = "beta") -> dict[str, pd.Series]:
    """lh, rh-as-lh and bilateral-mean versions of one column of the beta table."""
    b = bm.set_index("label")[col]
    lh = b[b.index.str.startswith("lh_")]
    rh = b[b.index.str.startswith("rh_")]
    rh_as_lh = rh.rename(index=lambda s: "lh_" + s[3:])
    both = pd.concat([lh, rh_as_lh], axis=1).mean(axis=1)
    return {"lh": lh, "bilateral_mean": both, "rh_as_lh": rh_as_lh}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=HERE / "results")
    ap.add_argument("--work", type=Path, default=HERE / "work")
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    files = sorted(a.results.glob("beta_map_*.tsv"))
    if not files:
        sys.exit(f"no beta_map_*.tsv under {a.results} -- run step 2 on CSD3 first")
    refs = references()
    rows, nulls = [], {}
    for f in files:
        bm = pd.read_csv(f, sep="\t")
        cell = bm.cell.iloc[0]
        cortex = bm[bm.label == "cortex_mean"].iloc[0]
        parc = bm[bm.label != "cortex_mean"]
        maps = split_map(parc)
        # lh-rh agreement: the map's own reliability
        r_lr = stats.spearmanr(maps["lh"], maps["rh_as_lh"].reindex(maps["lh"].index)).statistic
        n_sig = int((parc.p < 0.05).sum())
        n_neg = int((parc.beta < 0).sum())
        # three maps per cell: the beta map; the beta map conditional on the
        # whole-cortex mean slope (parcel-specific part); and each parcel's
        # loading on the cortex mean (what a purely global effect would draw)
        for col in ("beta", "beta_cond", "r_global"):
            if col not in parc.columns:
                continue
            cmaps = split_map(parc, col)
            for variant in ("lh", "bilateral_mean"):
                beta = cmaps[variant].dropna()
                for ref_name in ("C3", "PLS_lead", "dCT"):
                    ref = refs[ref_name].dropna()
                    rho, p, null = pls.spin_corr(ref, beta, n_perm=a.n_perm, seed=a.seed,
                                                 centroids=pls.HCP_CENTROIDS)
                    nulls[f"{col}__{variant}__{ref_name}"] = null
                    rows.append(dict(cell=cell, stratum=bm.stratum.iloc[0], map=col,
                                     beta_variant=variant, reference=ref_name,
                                     n_parcels=len(ref.index.intersection(beta.index)),
                                     rho=rho, p_spin=p, n_perm=a.n_perm,
                                     cortex_beta=cortex.beta, cortex_se=cortex.se, cortex_p=cortex.p,
                                     n_children=int(cortex.n), n_parcels_p05=n_sig,
                                     n_parcels_beta_neg=n_neg, rho_lh_rh=r_lr))
                    print(f"{cell:26s} {col:9s} {variant:15s} vs {ref_name:8s} rho = {rho:+.3f}  "
                          f"p_spin = {p:.4f}  (n = {rows[-1]['n_parcels']})", file=sys.stderr)
        np.savez(a.work / f"spin_null_{cell}.npz", **nulls)
        print(f"{cell}: cortex beta {cortex.beta:+.4f} (p {cortex.p:.3g}); {n_sig}/358 parcels p<0.05, "
              f"{n_neg}/358 negative; lh-rh rho {r_lr:.3f}", file=sys.stderr)
    out = a.results / "spin_tests.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.5g")
    print(f"wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
