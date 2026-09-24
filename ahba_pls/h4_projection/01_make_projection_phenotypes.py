"""H4, step 1 (laptop): projection phenotypes -- each child's HCP-MMP slope map
dotted with a transcriptomic region map -- for GCTA GREML on CSD3.

FOLLOWUP_GENETICS.md H4, on the 7.0 HCP-MMP run (thickness_hcp_70_aa6e91efba82,
8,716 children x 358 parcels).  The construction is legacy/hpc_v3/
make_phenotypes_v3.py::projection, unchanged:

    proj[child] = sum_parcels  z(slope[child, parcel]) * (w[parcel] - mean(w))

z is across children within parcel; w is the region map applied to both
hemispheres (the AHBA maps are left-hemisphere, 137 covered parcels, so the
sum runs over 274 parcels).  Two maps:

    slope_projC3    w = AHBA C3 (data/ahba_dme_hcp_top8kgenes_scores.csv)
    slope_projPLS2  w = ahba_pls lead component, thinning-oriented
                        (ahba_pls/results/hcp_pls_scores.csv[thinning_score])

plus each residualised on global_slope (`_resid`), and `slope_cov137`: the
plain mean slope over the same 274 parcels, so a projection's h2 can be read
against an unweighted phenotype on the identical parcel set.

Characterisation (results/h4_phenotype_characterisation.tsv): correlation with
global_slope and with each other, lh/rh Spearman-Brown consistency, and the
split-half circularity check the spec asks for on the PLS map: the PLS is
refitted on the group-mean maps of a random half of the children, its lead
region scores compared with the full-sample ones, and half-B children's
projections on half-A weights compared with their projections on the full
weights.

Writes work/pheno_h4.tsv (IID + phenotype columns; gitignored, rsync to CSD3).

Usage (repo root):  python ahba_pls/h4_projection/01_make_projection_phenotypes.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO / "ahba_pls" / "code"))
import pls  # noqa: E402

RUN = REPO / "out" / "thickness_hcp_70_aa6e91efba82"
C3_FILE = REPO / "data" / "ahba_dme_hcp_top8kgenes_scores.csv"
PLS_FILE = REPO / "ahba_pls" / "results" / "hcp_pls_scores.csv"
EXPR = Path.home() / "Git" / "AHBA" / "data" / "abagen-data" / "expression" / "hcp_3d.csv"
SEED = 0


def slope_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    p = pd.read_parquet(RUN / "phenotypes" / "phenotypes.parquet",
                        columns=["subject", "label", "phenotype", "value"])
    sl = p[p.phenotype == "slope"].pivot(index="subject", columns="label", values="value")
    ct = p[p.phenotype == "intercept"].pivot(index="subject", columns="label", values="value")
    assert sl.shape == (8716, 358), sl.shape
    assert not sl.isna().any().any()
    return sl, ct


def projection(sl: pd.DataFrame, w_base: pd.Series) -> pd.Series:
    """v3 idiom: z across children per parcel, dot the centred bilateral map."""
    cols = [f"{h}_{b}" for b in w_base.index for h in ("lh", "rh")]
    w = pd.Series({f"{h}_{b}": v for b, v in w_base.items() for h in ("lh", "rh")})[cols]
    w = w - w.mean()
    z = (sl[cols] - sl[cols].mean()) / sl[cols].std()
    return z.mul(w, axis=1).sum(axis=1)


def hemi_projection(sl: pd.DataFrame, w_base: pd.Series, h: str) -> pd.Series:
    cols = [f"{h}_{b}" for b in w_base.index]
    w = pd.Series({f"{h}_{b}": v for b, v in w_base.items()})[cols]
    w = w - w.mean()
    z = (sl[cols] - sl[cols].mean()) / sl[cols].std()
    return z.mul(w, axis=1).sum(axis=1)


def resid(y: pd.Series, x: pd.Series) -> pd.Series:
    b = np.polyfit(x, y, 1)
    return y - (b[0] * x + b[1])


def sb(x, y):
    r = x.corr(y)
    return r, 2 * r / (1 + r)


def fit_lead_scores(Xh: pd.DataFrame, sl: pd.DataFrame, ct: pd.DataFrame, idx) -> pd.Series:
    """ahba_pls option-2 PLS (dCT, CT) on the group means of children `idx`;
    lead = PLS2 gene-side region scores, thinning-oriented."""
    lab = sl.columns
    base = lab.str.replace(r"^(lh|rh)_", "", regex=True)
    dct = sl.loc[idx].mean().groupby(base).mean()
    cti = ct.loc[idx].mean().groupby(base).mean()
    Y = pd.DataFrame({"dCT": dct, "CT": cti})
    Y.index = "lh_" + Y.index
    Y = Y.loc[Xh.index]
    fit = pls.pls_svd(Xh, Y)
    return -fit.scores()["PLS2_gene_scores"]


def main() -> int:
    sl, ct = slope_matrix()
    g = sl.mean(axis=1).rename("global_slope")

    c = pd.read_csv(C3_FILE)
    c3 = pd.Series(c.C3.to_numpy(), index=c.label.astype(str))
    lead = pd.read_csv(PLS_FILE, index_col=0)["thinning_score"]
    lead.index = lead.index.str.replace("^lh_", "", regex=True)
    assert set(c3.index) == set(lead.index) and len(c3) == 137
    lead = lead.reindex(c3.index)
    covered = [f"{h}_{b}" for b in c3.index for h in ("lh", "rh")]

    ph = pd.DataFrame({
        "global_slope": g,
        "slope_cov137": sl[covered].mean(axis=1),
        "slope_projC3": projection(sl, c3),
        "slope_projPLS2": projection(sl, lead),
    })
    ph["slope_projC3_resid"] = resid(ph.slope_projC3, g)
    ph["slope_projPLS2_resid"] = resid(ph.slope_projPLS2, g)
    out = HERE / "work" / "pheno_h4.tsv"
    ph.rename(index=lambda s: "NDAR_" + s.replace("sub-NDAR", "")).rename_axis("IID") \
      .to_csv(out, sep="\t", float_format="%.7g")

    # ---- characterisation --------------------------------------------------
    rows = []
    for name, w in (("slope_projC3", c3), ("slope_projPLS2", lead)):
        r, s = sb(hemi_projection(sl, w, "lh"), hemi_projection(sl, w, "rh"))
        rows.append(dict(phenotype=name, r_with_global=ph[name].corr(g),
                         r_with_projC3=ph[name].corr(ph.slope_projC3),
                         r_with_projPLS2=ph[name].corr(ph.slope_projPLS2),
                         r_lh_rh=r, spearman_brown=s))
    for name in ("slope_projC3_resid", "slope_projPLS2_resid"):
        rows.append(dict(phenotype=name, r_with_global=ph[name].corr(g),
                         r_with_projC3=ph[name].corr(ph.slope_projC3),
                         r_with_projPLS2=ph[name].corr(ph.slope_projPLS2),
                         r_lh_rh=np.nan, spearman_brown=np.nan))
    for name, cols in (("slope_cov137", covered), ("global_slope", list(sl.columns))):
        lh = sl[[x for x in cols if x.startswith("lh_")]].mean(axis=1)
        rh = sl[[x for x in cols if x.startswith("rh_")]].mean(axis=1)
        r, s = sb(lh, rh)
        rows.append(dict(phenotype=name, r_with_global=ph[name].corr(g),
                         r_with_projC3=ph[name].corr(ph.slope_projC3),
                         r_with_projPLS2=ph[name].corr(ph.slope_projPLS2),
                         r_lh_rh=r, spearman_brown=s))
    char = pd.DataFrame(rows)
    char["r_map_C3_vs_PLS2"] = c3.corr(lead, method="spearman")

    # ---- split-half circularity check on the PLS map ------------------------
    c123w = pd.read_csv(REPO / "data" / "weights.csv", index_col=0)
    Xh = pd.read_csv(EXPR, index_col=0)[list(c123w.index)]
    id2label = dict(zip(c.id, "lh_" + c.label.astype(str)))
    Xh.index = [id2label.get(i, f"id{i}") for i in Xh.index]
    Xh = Xh.loc[["lh_" + b for b in c3.index]]
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(len(sl))
    A, B = sl.index[perm[: len(sl) // 2]], sl.index[perm[len(sl) // 2:]]
    full = fit_lead_scores(Xh, sl, ct, sl.index)
    sa, sb_ = fit_lead_scores(Xh, sl, ct, A), fit_lead_scores(Xh, sl, ct, B)
    for s_ in (sa, sb_):
        if s_.corr(full) < 0:
            s_ *= -1
    full.index = full.index.str.replace("^lh_", "", regex=True)
    sa.index = sa.index.str.replace("^lh_", "", regex=True)
    sb_.index = sb_.index.str.replace("^lh_", "", regex=True)
    projB_fullw = projection(sl.loc[B], lead)
    projB_Aw = projection(sl.loc[B], sa.reindex(c3.index))
    projA_fullw = projection(sl.loc[A], lead)
    projA_Bw = projection(sl.loc[A], sb_.reindex(c3.index))
    split = pd.DataFrame([dict(
        check="split_half_PLS_lead_map", seed=SEED, n_half=len(A),
        rho_full_vs_shipped=full.corr(lead, method="spearman"),
        rho_halfA_vs_full=sa.corr(full, method="spearman"),
        rho_halfB_vs_full=sb_.corr(full, method="spearman"),
        rho_halfA_vs_halfB=sa.corr(sb_, method="spearman"),
        r_projB_halfAweights_vs_fullweights=projB_Aw.corr(projB_fullw),
        r_projA_halfBweights_vs_fullweights=projA_Bw.corr(projA_fullw))])
    char.to_csv(HERE / "results" / "h4_phenotype_characterisation.tsv", sep="\t",
                index=False, float_format="%.4f")
    split.to_csv(HERE / "results" / "h4_split_half_check.tsv", sep="\t", index=False,
                 float_format="%.5f")
    print(f"wrote {out.relative_to(REPO)}: {ph.shape[0]} children x {ph.shape[1]} columns",
          file=sys.stderr)
    print(char.round(3).to_string(index=False), file=sys.stderr)
    print(split.round(5).T.to_string(), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
