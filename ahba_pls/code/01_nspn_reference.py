"""01_nspn_reference.py — NSPN (Whitaker & Vertes 2016, PNAS) reference signatures.

Inputs
  ~/Git/AHBA/data/whitakervertes2016_complete.csv   68 DK regions: CT, CT_delta, MT, MT_delta, PLS2 scores; label = ggseg DK
  ~/Git/AHBA/data/whitakervertes2016_genes.csv      multi-block CSV (BOM): PLS1/PLS2 ranked genes + bootstrapped weights
  ~/Git/AHBA_updated/outputs/expression_levels/ahba_dk_lh_native_ds{0,25,50}.csv   AHBA gene universes (column names)
Outputs (ahba_pls/data/reference/)
  nspn_dk_maps_68.csv, nspn_dk_maps_bilateral_34.csv, nspn_pls_gene_weights.csv, nspn_gene_overlap.md
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd

AHBA = Path('/Users/richard/Git/AHBA/data')
AHBA_UPD = Path('/Users/richard/Git/AHBA_updated/outputs/expression_levels')
OUT = Path('/Users/richard/Git/abcd_development/ahba_pls/data/reference'); OUT.mkdir(parents=True, exist_ok=True)

# ---- DK maps -------------------------------------------------------------
comp = pd.read_csv(AHBA / 'whitakervertes2016_complete.csv', index_col=0)
assert comp.shape[0] == 68 and comp['label'].is_unique
maps = comp.set_index('label')[['CT', 'CT_delta', 'MT', 'MT_delta', 'PLS2']].sort_index()
maps.to_csv(OUT / 'nspn_dk_maps_68.csv')
# bilateral: mean of lh_/rh_ per region, keep lh_ labels
maps['region'] = [l.split('_', 1)[1] for l in maps.index]
bil = maps.groupby('region').mean(numeric_only=True)
assert (maps.groupby('region').size() == 2).all()
bil.index = ['lh_' + r for r in bil.index]; bil.index.name = 'label'
bil.to_csv(OUT / 'nspn_dk_maps_bilateral_34.csv')

# ---- PLS gene weights ----------------------------------------------------
g = pd.read_csv(AHBA / 'whitakervertes2016_genes.csv', encoding='utf-8-sig', usecols=range(4))
g.columns = ['PLS1_gene', 'PLS1_z', 'PLS2_gene', 'PLS2_z']
p1 = g[['PLS1_gene', 'PLS1_z']].dropna().rename(columns={'PLS1_gene': 'gene'}).set_index('gene')
p2 = g[['PLS2_gene', 'PLS2_z']].dropna().rename(columns={'PLS2_gene': 'gene'}).set_index('gene')
# Excel date-corrupted symbols (e.g. 'Mar-02' for MARCH2 and MARC2, twice) — ambiguous, drop them
DATE_RE = r'^\d{1,2}-[A-Z][a-z]{2}$|^[A-Z][a-z]{2}-\d{2}$'
dropped = sorted(set(p1.index[p1.index.str.match(DATE_RE)]) | set(p2.index[p2.index.str.match(DATE_RE)]))
p1 = p1[~p1.index.isin(dropped)]; p2 = p2[~p2.index.isin(dropped)]
assert p1.index.is_unique and p2.index.is_unique
w = p1.join(p2, how='outer')
w.index.name = 'gene'
w = w.sort_values('PLS2_z', ascending=False)
w.to_csv(OUT / 'nspn_pls_gene_weights.csv')

# ---- overlap with AHBA universes --------------------------------------
res = {'dropped_date_corrupted_symbols': dropped, 'n_genes_PLS1': int(p1.shape[0]), 'n_genes_PLS2': int(p2.shape[0]), 'n_union': int(w.shape[0]),
       'n_both': int(w.dropna().shape[0])}
univ = {}
for ds in (0, 25, 50):
    cols = pd.read_csv(AHBA_UPD / f'ahba_dk_lh_native_ds{ds}.csv', nrows=0).columns[1:]
    univ[ds] = set(cols)
    res[f'n_ds{ds}_universe'] = len(cols)
    res[f'n_nspn_in_ds{ds}'] = len(set(w.index) & univ[ds])
absent0 = sorted(set(w.index) - univ[0])
res['n_absent_ds0'] = len(absent0)
res['absent_ds0_examples'] = absent0[:10]
res['absent_ds0_LOC_like'] = int(sum(a.startswith(('LOC', 'LINC', 'FAM', 'KIAA', 'C' )) and any(ch.isdigit() for ch in a) for a in absent0))
json.dump(res, open(OUT / 'nspn_gene_overlap.json', 'w'), indent=1)
pd.Series(absent0, name='gene').to_csv(OUT / 'nspn_genes_absent_ds0.txt', index=False, header=False)
print(json.dumps(res, indent=1))
print(maps.head(3)); print(bil.shape, w.shape)
