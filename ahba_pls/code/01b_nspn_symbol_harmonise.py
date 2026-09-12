"""01b_nspn_symbol_harmonise.py — update old NSPN (2016, Affymetrix-era) gene symbols to current HGNC symbols.

Queries mygene.info (POST /v3/query, scopes symbol,alias) for NSPN genes absent from the AHBA ds0 universe.
Rule: rename OLD -> NEW only if (a) the query maps to exactly one human Entrez gene, (b) NEW is in the ds0 universe,
(c) NEW is not already a gene in the NSPN list (no collisions). Everything else keeps its original symbol.
Writes nspn_pls_gene_weights.csv (harmonised; columns gene, gene_original, PLS1_z, PLS2_z) and
nspn_symbol_map.csv (all renames); raw mygene response cached in nspn_mygene_raw.json.
"""
import json
from pathlib import Path
import pandas as pd, requests

OUT = Path('/Users/richard/Git/abcd_development/ahba_pls/data/reference')
AHBA_UPD = Path('/Users/richard/Git/AHBA_updated/outputs/expression_levels')
w = pd.read_csv(OUT / 'nspn_pls_gene_weights.csv', index_col='gene')
if 'gene_original' in w.columns:  # re-run: start from originals
    w = w.set_index('gene_original').rename_axis('gene')[['PLS1_z', 'PLS2_z']]
univ = {ds: set(pd.read_csv(AHBA_UPD / f'ahba_dk_lh_native_ds{ds}.csv', nrows=0).columns[1:]) for ds in (0, 25, 50)}
absent = sorted(set(w.index) - univ[0])

cache = OUT / 'nspn_mygene_raw.json'
if cache.exists():
    hits = json.load(open(cache))
else:
    hits = []
    import time
    for i in range(0, len(absent), 250):
        for attempt in range(4):
            r = requests.post('https://mygene.info/v3/query',
                              data={'q': ','.join(absent[i:i+250]), 'scopes': 'symbol,alias',
                                    'fields': 'symbol,entrezgene', 'species': 'human'}, timeout=120)
            if r.ok: break
            time.sleep(3 * (attempt + 1))
        r.raise_for_status(); hits += r.json()
    json.dump(hits, open(cache, 'w'))
h = pd.DataFrame(hits)
h = h[h.get('notfound').isna()] if 'notfound' in h else h
h = h.dropna(subset=['symbol'])
# exactly one distinct current symbol per query
n_sym = h.groupby('query')['symbol'].nunique()
uniq = h[h['query'].isin(n_sym[n_sym == 1].index)].drop_duplicates('query').set_index('query')['symbol']
nspn_genes = set(w.index)
m = pd.DataFrame({'gene_original': uniq.index, 'gene_new': uniq.values})
m['in_ds0'] = m.gene_new.isin(univ[0])
m['collides'] = m.gene_new.isin(nspn_genes)
m['unchanged'] = m.gene_new == m.gene_original
ok = m[m.in_ds0 & ~m.collides & ~m.unchanged]
# a NEW symbol must be claimed by exactly one OLD symbol
ok = ok[~ok.gene_new.duplicated(keep=False)]
ren = dict(zip(ok.gene_original, ok.gene_new))
m.to_csv(OUT / 'nspn_symbol_map.csv', index=False)

w2 = w.copy(); w2['gene_original'] = w2.index
w2.index = [ren.get(g, g) for g in w2.index]; w2.index.name = 'gene'
assert w2.index.is_unique
w2 = w2[['gene_original', 'PLS1_z', 'PLS2_z']].sort_values('PLS2_z', ascending=False)
w2.to_csv(OUT / 'nspn_pls_gene_weights.csv')

res = {'n_absent_ds0_before': len(absent), 'n_queries_with_hit': int(h['query'].nunique()),
       'n_unique_mapping': int(len(uniq)), 'n_ambiguous_multi_hit': int((n_sym > 1).sum()),
       'n_new_in_ds0': int(m.in_ds0.sum()), 'n_collide_with_existing_nspn': int((m.in_ds0 & m.collides).sum()),
       'n_renamed': len(ren), 'examples': list(ren.items())[:8]}
for ds in (0, 25, 50):
    res[f'n_nspn_in_ds{ds}_after'] = len(set(w2.index) & univ[ds])
json.dump(res, open(OUT / 'nspn_symbol_harmonisation.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
