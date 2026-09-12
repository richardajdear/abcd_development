# NSPN (Whitaker & Vértes 2016, PNAS) reference signatures

## NOTES
- Built by `code/01_nspn_reference.py` then `code/01b_nspn_symbol_harmonise.py` (re-runnable; mygene response cached in `nspn_mygene_raw.json`).
- Source DK maps: `AHBA/data/whitakervertes2016_complete.csv` (68 regions, label already ggseg-style `lh_*`/`rh_*`).
  CT (mm), CT_delta (%/decade? — as in source; sign negative = thinning), MT, MT_delta, PLS2 scores.
- Source gene weights: `AHBA/data/whitakervertes2016_genes.csv` (BOM, multi-block). Only the four gene/weight columns were parsed;
  weights are the paper's bootstrapped Z (PLS1_z, PLS2_z). 26 Excel date-corrupted symbols
  (`Mar-02`×2 = MARCH2/MARC2 ambiguous, `Sep-*`, `Dec-01`, `Mar-*` = SEPT*/DEC1/MARCH*) were **dropped**.
- The AHBA DK matrices (`ahba_dk_lh_native_ds*.csv`) have 33 LH regions: **`lh_frontalpole` is the missing DK region**.
- Symbol harmonisation: NSPN uses 2013-era symbols. 7,155 of 20,710 NSPN genes were absent from the ds0 universe;
  examples: 61E3.4, A1CF, A26C1B, A2BP1, A2LD1, A3GALT2P, A4GNT, AADAC, AADACL1, AADACL2 (A2BP1→RBFOX1, AADACL1→NCEH1, ... are clearly outdated symbols).
  mygene.info (scopes symbol,alias; human) returned hits for 5,738, unique current symbol for 5,550
  (188 ambiguous, dropped). A rename OLD→NEW was applied only when NEW is in ds0, is not already an NSPN gene
  (19 collisions kept original) and is claimed by a single OLD symbol: **1,310 genes renamed**
  (full map in `nspn_symbol_map.csv`; `gene_original` column retained in the weights file).

## Files
| file | shape | content |
|---|---|---|
| nspn_dk_maps_68.csv | 68 × 5 | label-indexed CT, CT_delta, MT, MT_delta, PLS2 |
| nspn_dk_maps_bilateral_34.csv | 34 × 5 | LH/RH means, `lh_` labels |
| nspn_pls_gene_weights.csv | 20,710 × 3 | gene (harmonised), gene_original, PLS1_z, PLS2_z |
| nspn_symbol_map.csv | 5,550 rows | mygene mapping candidates with in_ds0 / collides flags |
| nspn_genes_absent_ds0.txt | 7,155 | NSPN symbols absent from ds0 before harmonisation |

## Overlap with AHBA gene universes
| universe | n genes | NSPN in universe (raw symbols) | after harmonisation |
|---|---|---|---|
| ds0 | 16,009 | 13,555 | 14,865 |
| ds25 | 12,007 | 10,195 | 11,175 |
| ds50 | 8,005 | 6,799 | 7,454 |
