"""Build the hpc_v3 regional-subset phenotypes and their GCTA export.

Run as ``python hpc_v3/make_phenotypes_v3.py [run_dir]`` from the repo root
(needs ``PYTHONPATH=src`` or ``pip install -e .``).  Defaults to the settled
genetic run for ``$ABCD_CONFIG`` (ct_70_noglobal_mv2_genetic).

Adds four phenotypes to the five settled ones (see README.md in this
directory for the rationale and the plan they feed):

``slope_topDelta``  (primary)
    Mean of the subject's per-region slope BLUPs over the K_BILATERAL=8
    bilateral regions with the fastest *group-mean* thinning (most negative
    ``slope_total`` in ``docs/developmental_maps_noglobal.csv``, averaged
    across hemispheres).  16 of 68 labels.

``slope_topC3``  (primary)
    Same construction over the 8 bilateral regions with the highest AHBA C3
    score (``data/ahba_dme_dsk_scores.csv``; scores are left-hemisphere,
    applied to both hemispheres).  Note the sign convention: high C3 =
    association cortex = *faster* thinning (rho(slope_total, C3) = -0.546),
    so this set overlaps the topDelta set (3 of 8 regions).

``slope_projDelta``, ``slope_projC3``  (exploratory)
    Projection scores: the subject's z-scored regional slope map dotted with
    the mean-centred group map (-slope_total, so positive = faster thinning)
    or the mean-centred C3 map.  Unlike the subset means these are contrasts
    -- the global thinning component largely cancels -- so they are much less
    collinear with ``global_slope`` (r = 0.28 / 0.41 vs 0.89 / 0.92) at the
    price of lower internal consistency (lh/rh Spearman-Brown 0.63 / 0.65 vs
    0.82 / 0.79 for the subset means; global_slope is 0.87).

Selection uses group-level maps only -- no subject-level variance, no
genotypes -- so there is no selection circularity for h2 / PRS / rg.

Writes:
  <run_dir>/gcta_inputs_v3/   phenotypes_gcta.txt (9 phenotype columns),
                              covar_quant.txt, covar_categorical.txt,
                              phenotype_manifest.tsv     [NOT committed:
                              per-subject, lives under gitignored out/]
  hpc_v3/phenotypes_v3_regions.csv          region sets + maps  [committed]
  hpc_v3/phenotypes_v3_characterization.csv 9x9 phenotype correlations
  hpc_v3/phenotypes_v3_consistency.csv      lh/rh consistency + k sensitivity
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from abcd import gcta_export  # noqa: E402

K_BILATERAL = 8          # top-k bilateral base regions per subset (=16 labels)
K_SENSITIVITY = (4, 8, 12, 17)

#: manifest rows for the new columns; priorities continue the settled 0-4.
NEW_PHENOTYPES = (
    dict(name="slope_topDelta", priority=5, role="target",
         display=f"mean slope, top-{K_BILATERAL} fastest-thinning regions"),
    dict(name="slope_topC3", priority=6, role="target",
         display=f"mean slope, top-{K_BILATERAL} AHBA-C3 regions"),
    dict(name="slope_projDelta", priority=7, role="exploratory",
         display="slope map projected on group thinning map"),
    dict(name="slope_projC3", priority=8, role="exploratory",
         display="slope map projected on AHBA C3 map"),
)


def region_rankings() -> tuple[pd.Series, pd.Series]:
    """(deltaCT ranking, C3 ranking) over the 34 bilateral base regions.

    deltaCT is sorted most-negative-first (fastest thinning first); C3 is
    sorted highest-first.  Both indexed by base region name.
    """
    m = pd.read_csv(REPO / "docs/developmental_maps_noglobal.csv")
    m["base"] = m.label.str.replace(r"^(lh|rh)_", "", regex=True)
    delta = m.groupby("base")["slope_total"].mean().sort_values()
    a = pd.read_csv(REPO / "data/ahba_dme_dsk_scores.csv")
    a["base"] = a.label.str.replace(r"^lh_", "", regex=True)
    c3 = a.set_index("base")["C3"].sort_values(ascending=False)
    assert len(delta) == 34 and len(c3) == 34
    assert set(delta.index) == set(c3.index)
    return delta, c3


def labels_for(bases: list[str]) -> list[str]:
    return [f"{h}_{b}" for b in bases for h in ("lh", "rh")]


def slope_matrix(run_dir: Path) -> pd.DataFrame:
    """8192 subjects x 68 labels of the ``slope`` phenotype (BLUP + fixed)."""
    ph = pd.read_parquet(run_dir / "phenotypes" / "phenotypes.parquet")
    sl = ph[ph.phenotype == "slope"].pivot(index="subject", columns="label",
                                           values="value")
    # Guard against the synthetic fixtures in scratch/ (SETUP_CONTEXT.md s7):
    # the real settled run has exactly 8192 subjects and 68 complete regions.
    assert sl.shape == (8192, 68), (
        f"expected 8192x68 slope matrix, got {sl.shape} -- wrong run dir? "
        "Never point this at scratch/hpc_test fixtures.")
    assert not sl.isna().any().any()
    return sl


def projection(sl: pd.DataFrame, weights_base: pd.Series) -> pd.Series:
    """Z-scored slope map dotted with the mean-centred bilateral weight map."""
    w = pd.Series({f"{h}_{b}": v for b, v in weights_base.items()
                   for h in ("lh", "rh")}).reindex(sl.columns)
    assert not w.isna().any()
    w = w - w.mean()
    z = (sl - sl.mean()) / sl.std()
    return z.mul(w, axis=1).sum(axis=1)


def build_new_columns(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(new phenotype columns indexed by subject, region-set table)."""
    sl = slope_matrix(run_dir)
    delta, c3 = region_rankings()
    set_d = list(delta.index[:K_BILATERAL])
    set_c = list(c3.index[:K_BILATERAL])

    new = pd.DataFrame({
        "slope_topDelta": sl[labels_for(set_d)].mean(axis=1),
        "slope_topC3": sl[labels_for(set_c)].mean(axis=1),
        "slope_projDelta": projection(sl, -delta),   # positive = faster thinning
        "slope_projC3": projection(sl, c3),
    })
    new.index.name = "subject"

    regions = pd.DataFrame({
        "base": delta.index,
        "group_mean_slope": delta.values,
        "C3": c3.reindex(delta.index).values,
        "rank_delta": range(1, 35),
    }).assign(
        rank_C3=lambda d: d.C3.rank(ascending=False).astype(int),
        in_topDelta=lambda d: d.base.isin(set_d),
        in_topC3=lambda d: d.base.isin(set_c),
    )
    return new, regions


def consistency_table(run_dir: Path) -> pd.DataFrame:
    """lh/rh internal consistency and k-sensitivity for the subset means."""
    sl = slope_matrix(run_dir)
    delta, c3 = region_rankings()
    g = sl.mean(axis=1)
    rows = []

    def sb(x, y):
        r = x.corr(y)
        return r, 2 * r / (1 + r)

    for k in K_SENSITIVITY:
        for name, ranking in (("topDelta", delta.index), ("topC3", c3.index)):
            bases = list(ranking[:k])
            lh = sl[[f"lh_{b}" for b in bases]].mean(axis=1)
            rh = sl[[f"rh_{b}" for b in bases]].mean(axis=1)
            r, sb_ = sb(lh, rh)
            rows.append(dict(
                phenotype=name, k_bilateral=k,
                r_with_global=sl[labels_for(bases)].mean(axis=1).corr(g),
                r_lh_rh=r, spearman_brown=sb_))
    for name, wb in (("projDelta", -delta), ("projC3", c3)):
        def hemi(h):
            cols = [c for c in sl.columns if c.startswith(h)]
            w = pd.Series({f"{h}{b}": v for b, v in wb.items()}).reindex(cols)
            w = w - w.mean()
            z = (sl[cols] - sl[cols].mean()) / sl[cols].std()
            return z.mul(w, axis=1).sum(axis=1)
        r, sb_ = sb(hemi("lh_"), hemi("rh_"))
        rows.append(dict(phenotype=name, k_bilateral=34,
                         r_with_global=projection(sl, wb).corr(g),
                         r_lh_rh=r, spearman_brown=sb_))
    # global reference row
    lh = sl[[c for c in sl.columns if c.startswith("lh_")]].mean(axis=1)
    rh = sl[[c for c in sl.columns if c.startswith("rh_")]].mean(axis=1)
    r, sb_ = sb(lh, rh)
    rows.append(dict(phenotype="global_slope", k_bilateral=34,
                     r_with_global=1.0, r_lh_rh=r, spearman_brown=sb_))
    return pd.DataFrame(rows).round(4)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if args:
        run_dir = Path(args[0])
    else:
        from abcd.config import active_run_dir
        run_dir = active_run_dir()

    built = gcta_export.build(run_dir)          # settled five + covariates
    new, regions = build_new_columns(run_dir)

    new = new.reset_index()
    new["IID"] = gcta_export._to_genetics_id(new["subject"])
    phen = built["phenotypes_gcta"].merge(
        new.drop(columns="subject"), on="IID", how="left")
    # merge must not change row count or leave the new columns missing
    assert len(phen) == len(built["phenotypes_gcta"])
    assert phen[[p["name"] for p in NEW_PHENOTYPES]].notna().all().all()
    # keep FID IID first
    cols = ["FID", "IID"] + [c for c in phen.columns if c not in ("FID", "IID")]
    phen = phen[cols]

    manifest = built["phenotype_manifest"]
    manifest = pd.concat([manifest, pd.DataFrame([
        {"name": p["name"],
         "mpheno": phen.columns.get_loc(p["name"]) - 1,
         "priority": p["priority"], "role": p["role"],
         "n_nonmissing": int(phen[p["name"]].notna().sum())}
        for p in NEW_PHENOTYPES])], ignore_index=True)

    out = run_dir / "gcta_inputs_v3"
    out.mkdir(parents=True, exist_ok=True)
    phen.to_csv(out / "phenotypes_gcta.txt", sep=" ", index=False, na_rep="NA")
    built["covar_quant"].to_csv(out / "covar_quant.txt", sep=" ",
                                index=False, na_rep="NA")
    built["covar_categorical"].to_csv(out / "covar_categorical.txt", sep=" ",
                                      index=False, na_rep="NA")
    manifest.to_csv(out / "phenotype_manifest.tsv", sep="\t", index=False)
    # The v2 sbatch scripts array over manifest rows, so a manifest holding
    # only the four new phenotypes lets the unchanged v2 scripts run just
    # those (03: --array=1-4, 04: --array=1-88, 05: --array=1-4).  mpheno
    # still indexes into the full 9-column phenotypes_gcta.txt.
    new_names = [p["name"] for p in NEW_PHENOTYPES]
    manifest[manifest.name.isin(new_names)].to_csv(
        out / "phenotype_manifest_new_only.tsv", sep="\t", index=False)
    print(f"{out}: phenotypes_gcta.txt {phen.shape[0]} x {phen.shape[1]}")
    print(manifest.to_string(index=False))

    # committed summary tables (no per-subject data)
    here = REPO / "hpc_v3"
    regions.to_csv(here / "phenotypes_v3_regions.csv", index=False)
    all9 = phen.drop(columns=["FID"]).set_index("IID")
    all9.corr().round(4).to_csv(here / "phenotypes_v3_characterization.csv")
    consistency_table(run_dir).to_csv(here / "phenotypes_v3_consistency.csv",
                                      index=False)
    print(f"summary tables written to {here}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
