"""01_build_y_matrix.py — bilateral DK imaging Y-matrix for the AHBA PLS.

Inputs (settled 7.0 runs, see ahba_pls/README.md):
  out/thickness_dsk_70_139406217085/fits/fixed.parquet   -> dCT (age_c slope, mm/yr), CT ((Intercept))
  out/t1t2_ratio_dsk_70_9a62dde44370/fits/fixed.parquet  -> dT1T2 (age_c slope, /yr), T1T2 ((Intercept))
  docs/developmental_maps_noglobal.csv                   -> slopePC1..3; also slope_total / baseline_thickness
                                                            used as a cross-check of dCT / CT.

Model in both runs: value ~ sex + age_c + (1 + age_c | subject) + (1 | site), global_covariate = none,
age_c = age - mean(age) over all observations in the run (config age_centre: null -> sample mean;
manifest.json records it). The (Intercept) is therefore the fitted value at the run's mean age for the
reference sex (female) with site/subject random effects at zero.

Outputs:
  ahba_pls/data/y_maps_68.csv            68 DK regions x 7 maps
  ahba_pls/data/y_maps_bilateral_34.csv  34 regions (LH/RH mean, indexed by 'lh_' label) x 7 maps
  ahba_pls/results/y_map_lr_agreement.csv LH-vs-RH Pearson r per map
  ahba_pls/results/y_map_correlations.csv Pearson (lower) / Spearman (upper) among the 7 bilateral maps
  ahba_pls/data/Y_MAPS.md                notes
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

BASE = Path("/Users/richard/Git/abcd_development")
HERE = BASE / "ahba_pls"
RUNS = {
    "CT": BASE / "out/thickness_dsk_70_139406217085",
    "T1T2": BASE / "out/t1t2_ratio_dsk_70_9a62dde44370",
}
MAPS_CSV = BASE / "docs/developmental_maps_noglobal.csv"


def md_table(df: pd.DataFrame) -> str:
    """Minimal markdown table (avoids the tabulate dependency)."""
    cols = [df.index.name or ""] + [str(c) for c in df.columns]
    rows = [[str(i)] + [str(v) for v in r] for i, r in zip(df.index, df.values)]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def fixed_terms(run: Path, prefix: str) -> pd.DataFrame:
    f = pd.read_parquet(run / "fits/fixed.parquet")
    assert set(f.term) == {"(Intercept)", "sexM", "age_c"}, f.term.unique()
    wide = f.pivot(index="label", columns="term", values="estimate")
    assert wide.shape[0] == 68
    return pd.DataFrame({f"d{prefix}": wide["age_c"], prefix: wide["(Intercept)"]})


def main() -> None:
    ct = fixed_terms(RUNS["CT"], "CT")
    t1 = fixed_terms(RUNS["T1T2"], "T1T2")
    maps = pd.read_csv(MAPS_CSV).set_index("label")
    assert maps.shape[0] == 68

    # cross-check against the docs table (same run: slope_total / baseline_thickness)
    d_slope = (ct["dCT"] - maps.loc[ct.index, "slope_total"]).abs().max()
    d_base = (ct["CT"] - maps.loc[ct.index, "baseline_thickness"]).abs().max()
    assert d_slope < 1e-4 and d_base < 1e-4, (d_slope, d_base)

    y68 = pd.concat([ct, t1, maps[["slopePC1", "slopePC2", "slopePC3"]]], axis=1)
    y68 = y68[["dCT", "CT", "dT1T2", "T1T2", "slopePC1", "slopePC2", "slopePC3"]]
    y68.index.name = "label"
    assert y68.notna().all().all() and y68.shape == (68, 7)

    # order: lh_ labels in the docs order, then their rh_ partners
    lh = [l for l in maps.index if l.startswith("lh_")]
    rh = ["rh_" + l[3:] for l in lh]
    assert len(lh) == 34 and set(rh) <= set(y68.index)
    y68 = y68.loc[lh + rh]

    L = y68.loc[lh]
    R = y68.loc[rh].set_axis(lh)
    y34 = (L + R) / 2
    y34.index.name = "label"

    lr = pd.DataFrame(
        {m: {"pearson_r": pearsonr(L[m], R[m])[0], "spearman_rho": spearmanr(L[m], R[m])[0]} for m in y34.columns}
    ).T
    lr.index.name = "map"

    n = len(y34.columns)
    corr = pd.DataFrame(np.eye(n), index=y34.columns, columns=y34.columns)
    for i, a in enumerate(y34.columns):
        for j, b in enumerate(y34.columns):
            if i > j:
                corr.iloc[i, j] = pearsonr(y34[a], y34[b])[0]  # lower triangle: Pearson
            elif i < j:
                corr.iloc[i, j] = spearmanr(y34[a], y34[b])[0]  # upper triangle: Spearman
    corr.index.name = "pearson_lower__spearman_upper"

    (HERE / "data").mkdir(exist_ok=True)
    (HERE / "results").mkdir(exist_ok=True)
    y68.to_csv(HERE / "data/y_maps_68.csv", float_format="%.8g")
    y34.to_csv(HERE / "data/y_maps_bilateral_34.csv", float_format="%.8g")
    lr.to_csv(HERE / "results/y_map_lr_agreement.csv", float_format="%.4f")
    corr.to_csv(HERE / "results/y_map_correlations.csv", float_format="%.4f")

    centres = {k: json.load(open(r / "manifest.json")).get("age_centre") for k, r in RUNS.items()}
    with open(HERE / "data/Y_MAPS.md", "w") as fh:
        fh.write(
            "# Imaging Y maps\n\n"
            "NOTES\n"
            "- Built by `code/01_build_y_matrix.py` from the settled 7.0 runs (see README Design table).\n"
            f"- dCT/CT from `{RUNS['CT'].name}`; dT1T2/T1T2 from `{RUNS['T1T2'].name}`; slopePC1-3 from "
            "`docs/developmental_maps_noglobal.csv`.\n"
            "- `dX` = fixed-effect `age_c` slope (per year), `X` = `(Intercept)` = fitted value at the run's "
            f"mean age for females with random effects at zero. Age centres: CT run {centres['CT']}, "
            f"T1T2 run {centres['T1T2']} (sample-mean centring; config `age_centre: null`).\n"
            f"- dCT/CT agree with docs slope_total/baseline_thickness to <{max(d_slope, d_base):.1e}.\n"
            "- `y_maps_68.csv`: 68 rows (34 lh_ then 34 rh_, ggseg DK labels). "
            "`y_maps_bilateral_34.csv`: LH/RH mean, indexed by the lh_ label.\n"
            "- T1T2 run: rh_temporalpole fit is singular/non-converged (fits/diagnostics.parquet); its estimates are "
            "kept but that region's T1T2/dT1T2 values are less reliable.\n"
            "- `results/y_map_lr_agreement.csv`: LH-vs-RH Pearson r and Spearman rho per map across the 34 regions.\n"
            "- `results/y_map_correlations.csv`: bilateral map correlations, Pearson below the diagonal, Spearman above.\n\n"
            "## LH vs RH agreement\n\n" + md_table(lr.round(3)) + "\n\n"
            "## Bilateral map correlations (Pearson lower / Spearman upper)\n\n" + md_table(corr.round(2)) + "\n"
        )
    print("y68", y68.shape, "y34", y34.shape, "max |dCT diff|", f"{d_slope:.1e}", "max |CT diff|", f"{d_base:.1e}")
    print(lr.round(3).to_string())
    print(corr.round(2).to_string())


if __name__ == "__main__":
    main()
