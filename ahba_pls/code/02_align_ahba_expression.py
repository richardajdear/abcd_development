"""02_align_ahba_expression.py — AHBA DK (LH, donor-native) expression at three DS levels, aligned to Y.

Inputs:
  ~/Git/AHBA_updated/outputs/expression_levels/ahba_dk_lh_native_ds{0,25,50}.csv  (33 lh_ regions x genes)
  ahba_pls/data/y_maps_bilateral_34.csv                                            (from 01_build_y_matrix.py)

Processing:
  * genes z-scored across the 33 regions (ddof=0), per DS level independently
  * Y restricted to the 33 AHBA-covered regions and reordered to the AHBA region order;
    an explicit assertion checks X.index == Y.index for every DS level

Outputs (ahba_pls/data/):
  X_ds0.parquet, X_ds25.parquet, X_ds50.parquet   z-scored regions x genes (float64, index 'label')
  Y_33.csv                                         aligned 7-map Y (33 x 7)
  region_order_33.txt                              one lh_ label per line, the shared row order
  missing_region.txt                               the LH DK region without AHBA coverage
  INPUTS.md                                        notes on shapes, missing region, z-scoring, reload
"""
from pathlib import Path

import numpy as np
import pandas as pd

AHBA = Path("/Users/richard/Git/AHBA_updated/outputs/expression_levels")
HERE = Path("/Users/richard/Git/abcd_development/ahba_pls")
DS = [0, 25, 50]


def main() -> None:
    y34 = pd.read_csv(HERE / "data/y_maps_bilateral_34.csv", index_col="label")
    assert y34.shape == (34, 7)

    X = {ds: pd.read_csv(AHBA / f"ahba_dk_lh_native_ds{ds}.csv", index_col="label") for ds in DS}
    order = list(X[0].index)
    assert len(order) == 33 and all(list(x.index) == order for x in X.values()), "DS files differ in region order"
    assert set(order) <= set(y34.index)
    missing = sorted(set(y34.index) - set(order))
    assert len(missing) == 1, missing
    missing = missing[0]

    # gene-set nesting
    g = {ds: set(X[ds].columns) for ds in DS}
    assert g[50] <= g[25] <= g[0], "DS gene sets are not nested"
    for x in X.values():
        assert x.notna().all().all() and x.columns.is_unique

    # z-score across regions
    Xz = {}
    for ds, x in X.items():
        z = (x - x.mean(axis=0)) / x.std(axis=0, ddof=0)
        assert np.isfinite(z.values).all()
        assert np.allclose(z.mean(axis=0), 0, atol=1e-10) and np.allclose(z.std(axis=0, ddof=0), 1, atol=1e-10)
        Xz[ds] = z.astype("float64")

    Y = y34.loc[order]
    for ds in DS:
        assert list(Xz[ds].index) == list(Y.index) == order, f"row order mismatch at ds{ds}"

    out = HERE / "data"
    for ds in DS:
        Xz[ds].to_parquet(out / f"X_ds{ds}.parquet")
    Y.to_csv(out / "Y_33.csv", float_format="%.8g")
    (out / "region_order_33.txt").write_text("\n".join(order) + "\n")
    (out / "missing_region.txt").write_text(missing + "\n")

    shapes = {ds: Xz[ds].shape for ds in DS}
    (out / "INPUTS.md").write_text(
        "# Aligned PLS inputs\n\n"
        "NOTES\n"
        "- Built by `code/02_align_ahba_expression.py`; Y from `code/01_build_y_matrix.py`.\n"
        "- X: AHBA left-hemisphere Desikan-Killiany expression, donor-native parcellation "
        "(`AHBA_updated/outputs/expression_levels/ahba_dk_lh_native_ds{0,25,50}.csv`), 33 regions.\n"
        f"- Missing LH DK region (no AHBA samples): **{missing}**. Y therefore has 33 of 34 bilateral regions.\n"
        "- Genes z-scored across the 33 regions (population sd, ddof=0), independently per DS level. "
        "Z-scoring is region-wise standardisation of each gene; no region-wise centring across genes.\n"
        f"- Shapes (regions x genes): ds0 {shapes[0]}, ds25 {shapes[25]}, ds50 {shapes[50]}. "
        f"Gene sets are nested: ds50 ({len(g[50])}) ⊂ ds25 ({len(g[25])}) ⊂ ds0 ({len(g[0])}).\n"
        f"- Y_33.csv: {Y.shape}, columns {list(Y.columns)}; rows in the AHBA region order "
        "(`region_order_33.txt`), identical to every X_ds*.parquet index (asserted at build time).\n"
        "- Note the T1T2 run's rh_temporalpole fit was singular/non-converged; the bilateral lh_temporalpole "
        "T1T2/dT1T2 values inherit that (see data/Y_MAPS.md).\n\n"
        "## Reload\n\n"
        "```python\n"
        "import pandas as pd\n"
        "X = pd.read_parquet('ahba_pls/data/X_ds25.parquet')      # 33 x 12007, z-scored\n"
        "Y = pd.read_csv('ahba_pls/data/Y_33.csv', index_col='label')\n"
        "assert (X.index == Y.index).all()\n"
        "```\n"
    )
    print("missing:", missing, "| shapes:", shapes, "| Y", Y.shape)
    print("nesting ds50<ds25<ds0:", len(g[50]), len(g[25]), len(g[0]))


if __name__ == "__main__":
    main()
