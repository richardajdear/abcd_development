# Aligned PLS inputs

NOTES
- Built by `code/02_align_ahba_expression.py`; Y from `code/01_build_y_matrix.py`.
- X: AHBA left-hemisphere Desikan-Killiany expression, donor-native parcellation (`AHBA_updated/outputs/expression_levels/ahba_dk_lh_native_ds{0,25,50}.csv`), 33 regions.
- Missing LH DK region (no AHBA samples): **lh_frontalpole**. Y therefore has 33 of 34 bilateral regions.
- Genes z-scored across the 33 regions (population sd, ddof=0), independently per DS level. Z-scoring is region-wise standardisation of each gene; no region-wise centring across genes.
- Shapes (regions x genes): ds0 (33, 16009), ds25 (33, 12007), ds50 (33, 8005). Gene sets are nested: ds50 (8005) ⊂ ds25 (12007) ⊂ ds0 (16009).
- Y_33.csv: (33, 7), columns ['dCT', 'CT', 'dT1T2', 'T1T2', 'slopePC1', 'slopePC2', 'slopePC3']; rows in the AHBA region order (`region_order_33.txt`), identical to every X_ds*.parquet index (asserted at build time).
- Note the T1T2 run's rh_temporalpole fit was singular/non-converged; the bilateral lh_temporalpole T1T2/dT1T2 values inherit that (see data/Y_MAPS.md).

## Reload

```python
import pandas as pd
X = pd.read_parquet('ahba_pls/data/X_ds25.parquet')      # 33 x 12007, z-scored
Y = pd.read_csv('ahba_pls/data/Y_33.csv', index_col='label')
assert (X.index == Y.index).all()
```
