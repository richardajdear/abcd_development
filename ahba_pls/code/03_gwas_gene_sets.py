"""03_gwas_gene_sets.py — SCZ / MDD GWAS gene sets, MAGMA gene-level Z, and brain-expressed background universes.

Sources
  legacy/hpc/work/genesets/scz2022_ExtendedDataTable1.xlsx   Trubetskoy 2022 ED Table 1 (120 prioritised) + ST12 (685 locus pool)
  legacy/hpc/work/genesets/mdd2024_TableS21.xlsx             Adams 2025 Cell (doi 10.1016/j.cell.2024.12.002) Table S21: 308 high-confidence + 4,600 pool
  legacy/hpc/work/genesets/{scz_prioritised,mdd_highconf}.txt  the Entrez-keyed MAGMA set-annot files actually used in the hpc MAGMA run
  AHBA/data/gwas/hammerschlag2020_howard2019_trubetskoy2022.csv, howard2019_tableS9.csv   Howard 2019 MDD lists
  genetic_analysis/inputs/magma/{SCZ,MDD}.genes.out           MAGMA gene-level results (GENE = Entrez)
  data/symbol2entrez.csv                               symbol <-> Entrez map
Set definitions replicate hpc/work/make_prioritised_genesets.py and make_mdd_genesets.py exactly, but in SYMBOL space
(the hpc files are Entrez-keyed and restricted to genes in NCBI37.3 gene.loc). Symbols absent from the AHBA ds0 universe are
harmonised via mygene.info with the same conservative rule as 01b (unique current symbol, present in ds0, no collision).
"""
import json, time
from pathlib import Path
import pandas as pd, requests

HPC = Path('/Users/richard/Git/abcd_development/legacy/hpc/work')
D = Path('/Users/richard/Git/abcd_development/data')
GW = Path('/Users/richard/Git/AHBA/data/gwas')
AHBA_UPD = Path('/Users/richard/Git/AHBA_updated/outputs/expression_levels')
OUT = Path('/Users/richard/Git/abcd_development/ahba_pls/data/reference/gene_sets'); OUT.mkdir(parents=True, exist_ok=True)
univ = {ds: set(pd.read_csv(AHBA_UPD / f'ahba_dk_lh_native_ds{ds}.csv', nrows=0).columns[1:]) for ds in (0, 25, 50)}

# ---------------- symbol-space sets --------------------------------------------------------
sets, source = {}, {}
def add(name, symbols, src):
    s = sorted({x.strip() for x in symbols if isinstance(x, str) and x.strip()})
    sets[name] = s; source[name] = src

x = pd.ExcelFile(HPC / 'genesets/scz2022_ExtendedDataTable1.xlsx')
t1 = x.parse('Extended.Data.Table.1'); st12 = x.parse('ST12 all criteria')
assert len(t1) == 120 and len(st12) == 685
add('SCZ_prioritised', t1['Symbol.ID'], 'Trubetskoy2022 ED Table 1 (all 120 prioritised)')
add('SCZ_prioritised_pc', t1.loc[t1.gene_biotype == 'protein_coding', 'Symbol.ID'], 'ED Table 1, protein_coding only')
for n, c in {'SCZ_prio_finemap': 'FINEMAP.priority.gene', 'SCZ_prio_smr': 'SMR.priority.gene', 'SCZ_prio_rare': 'Rare.priority.gene'}.items():
    add(n, t1.loc[t1[c] == 1, 'Symbol.ID'], f'ED Table 1, {c} == 1')
add('SCZ_locus_pool', st12['Symbol.ID'], 'Trubetskoy2022 ST12 all criteria (685-gene candidate pool at GWS loci)')
prio = set(sets['SCZ_prioritised'])
add('SCZ_pool_not_prio', [s for s in st12['Symbol.ID'] if isinstance(s, str) and s not in prio], 'ST12 minus ED Table 1')

y = pd.ExcelFile(HPC / 'genesets/mdd2024_TableS21.xlsx')
hc = y.parse('High-confidence Gene List'); pool = y.parse('Table S21 Gene Mapping Methods')
assert len(hc) == 308 and len(pool) == 4600
add('MDD_highconf', hc['Gene'], 'Adams2025 Table S21 High-confidence Gene List (308; finemapping, expression or protein)')
for n, c in {'MDD_hc_finemap': 'Fine_mapping', 'MDD_hc_expression': 'Expression', 'MDD_hc_protein': 'Protein'}.items():
    add(n, hc.loc[hc[c].astype(str).str.lower() == 'true', 'Gene'], f'High-confidence list, {c} == TRUE')
add('MDD_pool', pool['Gene'], 'Adams2025 Table S21 all mapped genes (4,600; any of 7 methods incl. fastBAT/H-MAGMA)')
hcs = set(sets['MDD_highconf'])
add('MDD_pool_not_hc', [s for s in pool['Gene'] if isinstance(s, str) and s not in hcs], 'Table S21 pool minus high-confidence')

h = pd.read_csv(GW / 'hammerschlag2020_howard2019_trubetskoy2022.csv', encoding='utf-8-sig')
add('MDD_howard2019_32', h['MDD'].dropna(), 'AHBA/data/gwas/hammerschlag2020_howard2019_trubetskoy2022.csv MDD column (Howard 2019; 32 genes, exact derivation not recorded — legacy Dear2024 set)')
add('SCZ_trubetskoy_smr101', h['SCZ'].dropna(), 'same file, SCZ column = Trubetskoy2022 SMR P+F genes (trubetskoy2022_data.csv; legacy Dear2024 set)')
hs9 = pd.read_csv(GW / 'howard2019_tableS9.csv', encoding='utf-8-sig', skiprows=1)
add('MDD_howard2019_magma', hs9['Gene Name'].dropna(), 'Howard2019 Table S9: MAGMA genome-wide significant genes (P<2.8e-6)')

# ---------------- harmonise symbols absent from ds0 ---------------------------------------
all_syms = sorted(set().union(*sets.values()))
absent = sorted(set(all_syms) - univ[0])
cache = OUT / 'geneset_mygene_raw.json'
if cache.exists():
    hits = json.load(open(cache))
else:
    hits = []
    for i in range(0, len(absent), 250):
        for attempt in range(4):
            r = requests.post('https://mygene.info/v3/query', data={'q': ','.join(absent[i:i+250]), 'scopes': 'symbol,alias',
                              'fields': 'symbol,entrezgene', 'species': 'human'}, timeout=120)
            if r.ok: break
            time.sleep(3 * (attempt + 1))
        r.raise_for_status(); hits += r.json()
    json.dump(hits, open(cache, 'w'))
hd = pd.DataFrame(hits); hd = hd[hd['notfound'].isna()] if 'notfound' in hd else hd
hd = hd.dropna(subset=['symbol'])
nsym = hd.groupby('query')['symbol'].nunique()
uniq = hd[hd['query'].isin(nsym[nsym == 1].index)].drop_duplicates('query').set_index('query')['symbol']
uniq = uniq[uniq.isin(univ[0]) & ~uniq.isin(all_syms) & (uniq != uniq.index)]  # NEW in ds0, NEW not already a set gene
uniq = uniq[~uniq.duplicated(keep=False)]
ren = uniq.to_dict()
pd.Series(ren, name='gene_new').rename_axis('gene_original').to_csv(OUT / 'geneset_symbol_map.csv')

# ---------------- MAGMA gene-level + entrez map -----------------------------------------
s2e = pd.read_csv(D / 'symbol2entrez.csv').dropna()
s2e['entrez'] = s2e.entrez.astype(int)
e2s = s2e.drop_duplicates('entrez').set_index('entrez')['symbol']
magma = {tr: pd.read_csv(HPC / f'results/magma/{tr}.genes.out', sep=r'\s+') for tr in ('SCZ', 'MDD')}
# Entrez ids without a symbol in symbol2entrez.csv -> current symbol from mygene.info (scopes entrezgene)
unm = sorted({int(g) for m in magma.values() for g in m.GENE if int(g) not in e2s.index})
cache2 = OUT / 'magma_entrez_mygene_raw.json'
if cache2.exists():
    hits2 = json.load(open(cache2))
else:
    hits2 = []
    for i in range(0, len(unm), 500):
        for attempt in range(4):
            r = requests.post('https://mygene.info/v3/query', data={'q': ','.join(map(str, unm[i:i+500])), 'scopes': 'entrezgene',
                              'fields': 'symbol', 'species': 'human'}, timeout=120)
            if r.ok: break
            time.sleep(3 * (attempt + 1))
        r.raise_for_status(); hits2 += r.json()
    json.dump(hits2, open(cache2, 'w'))
h2 = pd.DataFrame(hits2); h2 = h2[h2['notfound'].isna()] if 'notfound' in h2 else h2
h2 = h2.dropna(subset=['symbol']).drop_duplicates('query')
e2s_mg = pd.Series(h2.symbol.values, index=h2['query'].astype(int))
for tr, m in magma.items():
    m.insert(1, 'symbol', m.GENE.map(e2s))
    m.insert(2, 'symbol_source', pd.Series('symbol2entrez', index=m.index).where(m.symbol.notna(), None))
    fill = m.symbol.isna() & m.GENE.isin(e2s_mg.index)
    m.loc[fill, 'symbol'] = m.loc[fill, 'GENE'].map(e2s_mg); m.loc[fill, 'symbol_source'] = 'mygene'
    m.to_csv(OUT / f'magma_{tr}_genes.tsv', sep='\t', index=False)
    magma[tr] = m
magma_stats = {tr: {'n_genes': int(len(m)), 'n_with_symbol': int(m.symbol.notna().sum()),
                    'n_symbol_from_symbol2entrez': int((m.symbol_source == 'symbol2entrez').sum()),
                    'n_symbol_from_mygene': int((m.symbol_source == 'mygene').sum()),
                    'n_symbol_in_ds0': int(m.symbol.isin(univ[0]).sum()),
                    'n_duplicate_symbols': int(m.symbol.dropna().duplicated().sum())} for tr, m in magma.items()}
z_genes = set(magma['SCZ'].symbol.dropna()) & set(magma['MDD'].symbol.dropna())
bg = {}
for ds in (0, 25, 50):
    b = sorted(univ[ds] & z_genes); bg[ds] = b
    Path(OUT / f'brain_background_ds{ds}.txt').write_text('\n'.join(b) + '\n')

# ---------------- MAGMA Entrez set-annot files used in hpc, with symbols -------------------
for src, dst in [('scz_prioritised.txt', 'magma_setannot_SCZ_hpc.tsv'), ('mdd_highconf.txt', 'magma_setannot_MDD_hpc.tsv')]:
    a = pd.read_csv(HPC / 'genesets' / src, sep='\t', header=None, names=['entrez', 'set'])
    a['symbol'] = a.entrez.map(e2s).fillna(a.entrez.map(e2s_mg))
    a.to_csv(OUT / dst, sep='\t', index=False)
    for name, g in a.groupby('set'):
        source[name + ' (hpc Entrez file)'] = f'{src}: {len(g)} Entrez ids, {g.symbol.notna().sum()} with symbol, {g.symbol.isin(univ[0]).sum()} in ds0'

# ---------------- write sets + summary ---------------------------------------------------
rows, long = [], []
for name, s in sets.items():
    s2 = sorted({ren.get(g, g) for g in s})
    Path(OUT / f'{name}.txt').write_text('\n'.join(s2) + '\n')
    long += [(name, g) for g in s2]
    rows.append({'set': name, 'n_genes': len(s2), 'n_renamed': sum(g in ren for g in s),
                 **{f'n_ds{ds}': len(set(s2) & univ[ds]) for ds in (0, 25, 50)},
                 **{f'n_bg{ds}': len(set(s2) & set(bg[ds])) for ds in (0, 25, 50)},
                 'source': source[name]})
summ = pd.DataFrame(rows); summ.to_csv(OUT / 'gene_sets_summary.tsv', sep='\t', index=False)
pd.DataFrame(long, columns=['set', 'gene']).to_csv(OUT / 'gene_sets_long.tsv', sep='\t', index=False)
meta = {'n_set_symbols': len(all_syms), 'n_absent_ds0_before': len(absent), 'n_renamed': len(ren),
        'rename_examples': list(ren.items())[:10], 'magma': magma_stats,
        'background': {ds: len(b) for ds, b in bg.items()}, 'hpc_entrez_sets': {k: v for k, v in source.items() if 'hpc Entrez' in k},
        'symbol2entrez_rows': int(len(s2e))}
json.dump(meta, open(OUT / 'gene_sets_meta.json', 'w'), indent=1)
print(summ.drop(columns='source').to_string(index=False)); print(json.dumps(meta, indent=1))
