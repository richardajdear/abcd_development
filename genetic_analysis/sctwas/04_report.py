"""SCTWAS checks -> $W/report/*.tsv and the generated sections of $W/README.md (between the AUTO markers).

Tables written (all under $W/report/):
  celltype_summary.tsv     per tag: features tested, eGenes, fitted, >=1 CS, genes per model, median SNPs per model,
                           s_rss distribution, share s_rss > 0.3, CS-count distribution, non-converged, errors
  trait_summary.tsv        per trait x tag x model: genes tested, lambda_GC (median z^2 / 0.455), n p < 2.5e-6
  top10.tsv                10 top genes per trait x tag x model (with GENCODE v38 symbol)
  positive_controls.tsv    SingleBrain Supp. Table 11 SCZ colocalisations (PP4 >= 0.8, major classes = Ext. Data Fig. 2;
                           Fig. 3c genes flagged), z and p in every tag, SCZ25 + PGC3
  bulk_comparison.tsv      Spearman of SCZ z vs old bulk GTEx Brain_Cortex TWAS (09_twas.csv, SCZ_EUR = PGC3 EUR)
  ext1_randomz_check.tsv   Spearman |fixed_z| vs Random_Z at the lead SNP per Ext1 feature (top_assoc)
  n_sensitivity.tsv        Ext1 n = 757 vs n = 983: weights and TWAS z agreement
"""
import glob, os, re
import numpy as np, pandas as pd
from scipy.stats import spearmanr

W = '/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas'
R = f'{W}/report'; os.makedirs(R, exist_ok=True)
TRAITS = ['SCZ25_EUR', 'PGC3_EUR', 'MDD_EUR', 'ABCD_thin_EUR', 'ABCD_base_EUR']
ORDER = ['Ext1', 'Ext1_n983', 'Ext', 'OPC', 'IN', 'OD', 'Ast', 'MG', 'Ext2', 'Ext6', 'Ext3', 'Ext4', 'Ext5', 'Ext7', 'Ext8']
LAYER = {'Ext1': 'L2/3 IT', 'Ext2': 'L4 IT', 'Ext3': 'L5/6 NP', 'Ext4': 'L5 IT', 'Ext5': 'L6b', 'Ext6': 'L6 CT',
         'Ext7': 'L6 IT Car3', 'Ext8': 'L6 IT', 'Ext': 'excitatory neurons', 'IN': 'inhibitory neurons', 'OPC': 'OPC',
         'OD': 'oligodendrocytes', 'Ast': 'astrocytes', 'MG': 'microglia', 'Ext1_n983': 'L2/3 IT (n = 983 sensitivity)'}
MAJOR = ['Ext', 'IN', 'OD', 'Ast', 'MG', 'OPC', 'End', 'MiGA']
# Fig. 3c eGenes, read off the published figure; each one's cell type matches its Supp. Table 11 PP4 >= 0.8 rows.
FIG3C = ['ANKS1B', 'CACNB2', 'CNTN4', 'CSMD1', 'YPEL4', 'CUL3', 'SLC4A10', 'SLC25A17', 'SCRIB', 'FOXP1', 'FUT9', 'SPCS1',
         'GNL3', 'KCNB1', 'PTPRU', 'AGO4', 'BDNF', 'MYT1L', 'NEGR1', 'SH3RF1', 'DPYD-AS1', 'DPYD', 'WBP1L', 'RP11-133K1.12',
         'PLCH2', 'ACTR5', 'PTPRT', 'RERE', 'PPARGC1A', 'RGS6', 'SLCO4C1', 'SNAP91', 'SOX2-OT', 'SPIRE2', 'RAB21', 'TMX1',
         'CNOT7', 'IPO11']

sym = pd.read_csv(f'{W}/ref/gencode_v38_genes.tsv.gz', sep='\t')
sym['gene'] = sym.gene_version.str.replace(r'\..*$', '', regex=True)
SYM = sym.drop_duplicates('gene').set_index('gene').symbol
tags = [t for t in ORDER if os.path.exists(f'{W}/weights/{t}.genes.tsv.gz') and os.path.exists(f'{W}/assoc/ABCD_base_EUR.{t}.top1.tsv.gz')]


def md(df, floatfmt=3):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == 'f':
            df[c] = df[c].map(lambda v: '' if pd.isna(v) else (f'{v:.2e}' if (abs(v) < 1e-3 and v != 0) else f'{v:.{floatfmt}g}'))
    h = '| ' + ' | '.join(map(str, df.columns)) + ' |\n|' + '---|' * len(df.columns) + '\n'
    return h + ''.join('| ' + ' | '.join(map(str, r)) + ' |\n' for r in df.values)


def celltype_summary():
    rows = []
    for t in tags:
        g = pd.read_csv(f'{W}/weights/{t}.genes.tsv.gz', sep='\t')
        s = pd.read_csv(f'{W}/weights/{t}.susie.tsv.gz', sep='\t')
        t1 = pd.read_csv(f'{W}/weights/{t}.top1.tsv.gz', sep='\t')
        base = t.replace('_n983', '')
        ps = pd.read_csv(f'{W}/stage1/{base}/prep_summary.tsv', sep='\t').iloc[0]
        ok = g[g.status == 'ok']
        sr = ok.s_rss
        rows.append(dict(tag=t, label=LAYER.get(t, ''), features_tested=ps.features_tested, egenes=ps.egenes,
                         egenes_with_snps=len(g), fit_errors=int((g.status != 'ok').sum()),
                         non_converged=int((ok.converged == False).sum()),
                         genes_ge1_cs=int((ok.n_cs > 0).sum()), w_flag=int(g.w_flag.sum()) if 'w_flag' in g else np.nan, susie_model_genes=s.gene.nunique(), top1_model_genes=t1.gene.nunique(),
                         median_snps_susie=float(s.groupby('gene').size().median()) if len(s) else np.nan,
                         median_cis_snps=float(g.n_snps_locus.median()),
                         cs_1=int((ok.n_cs == 1).sum()), cs_2=int((ok.n_cs == 2).sum()), cs_3plus=int((ok.n_cs >= 3).sum()),
                         s_rss_median=sr.median(), s_rss_q90=sr.quantile(.9), share_s_gt_0_3=(sr > 0.3).mean(),
                         share_s_gt_0_3_in_susie_model=(ok[ok.n_cs > 0].s_rss > 0.3).mean()))
    d = pd.DataFrame(rows); d.to_csv(f'{R}/celltype_summary.tsv', sep='\t', index=False); return d


_CACHE = {}


def load_assoc(trait, tag, model):
    """Cached; w_flag is forced to False for top1 (single marginal weight) and added if an older file lacks it."""
    k = (trait, tag, model)
    if k not in _CACHE:
        f = f'{W}/assoc/{trait}.{tag}.{model}.tsv.gz'
        a = pd.read_csv(f, sep='\t') if os.path.exists(f) else None
        if a is not None:
            a['w_flag'] = False if (model == 'top1' or 'w_flag' not in a) else a.w_flag.astype(bool)
        _CACHE[k] = a
    return None if _CACHE[k] is None else _CACHE[k].copy()


def trait_summary():
    rows, tops = [], []
    for t in TRAITS:
        for tag in tags:
            for m in ['susie', 'top1']:
                a = load_assoc(t, tag, m)
                if a is None:
                    continue
                u = a[~a.w_flag.astype(bool)]
                rows.append(dict(trait=t, tag=tag, model=m, genes=len(a), lambda_gc=np.median(a.z ** 2) / 0.455,
                                 n_p_lt_2_5e6=int((a.p < 2.5e-6).sum()), n_fdr05=int(bh(a.p).lt(0.05).sum()),
                                 n_w_flag=int(a.w_flag.sum()), lambda_gc_unflagged=np.median(u.z ** 2) / 0.455,
                                 n_p_lt_2_5e6_unflagged=int((u.p < 2.5e-6).sum())))
                x = a.nsmallest(10, 'p').assign(trait=t, tag=tag, model=m)
                x['symbol'] = x.gene.map(SYM)
                tops.append(x[['trait', 'tag', 'model', 'gene', 'symbol', 'n_model_snps', 'z', 'p', 'n_cs', 's_rss', 'w_flag']])
    d = pd.DataFrame(rows); d.to_csv(f'{R}/trait_summary.tsv', sep='\t', index=False)
    tp = pd.concat(tops); tp.to_csv(f'{R}/top10.tsv', sep='\t', index=False)
    return d, tp


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return pd.Series(np.minimum(q, 1))


def positive_controls():
    c = pd.read_csv(f'{W}/ref/singlebrain_ST11_scz_coloc.tsv', sep='\t')
    c.columns = ['disease', 'gwas', 'locus', 'gwas_lead', 'coloc_ct', 'qtl_lead', 'symbol', 'gene', 'pp4']
    c8 = c[(c.pp4 >= 0.8)]
    maj = c8[c8.coloc_ct.isin(MAJOR)]
    allc = maj.groupby('gene').coloc_ct.apply(lambda s: ','.join(sorted(set(s))))
    sub = c8[~c8.coloc_ct.isin(MAJOR)].groupby('gene').coloc_ct.apply(lambda s: ','.join(sorted(set(s))))
    genes = pd.DataFrame({'gene': sorted(set(maj.gene))})
    genes['symbol'] = genes.gene.map(maj.drop_duplicates('gene').set_index('gene').symbol)
    genes['coloc_classes'] = genes.gene.map(allc); genes['coloc_subtypes'] = genes.gene.map(sub).fillna('')
    genes['fig3c'] = genes.symbol.isin(FIG3C)
    rows = []
    for _, r in genes.iterrows():
        for tag in tags:
            for t in ['SCZ25_EUR', 'PGC3_EUR']:
                for m in ['susie', 'top1']:
                    a = load_assoc(t, tag, m)
                    hit = a[a.gene == r.gene] if a is not None else pd.DataFrame()
                    rows.append(dict(gene=r.gene, symbol=r.symbol, coloc_classes=r.coloc_classes, coloc_subtypes=r.coloc_subtypes,
                                     fig3c=r.fig3c, tag=tag, trait=t, model=m,
                                     z=hit.z.iloc[0] if len(hit) else np.nan, p=hit.p.iloc[0] if len(hit) else np.nan,
                                     w_flag=bool(hit.w_flag.iloc[0]) if len(hit) and 'w_flag' in hit else np.nan))
    d = pd.DataFrame(rows)
    base = d.tag.str.replace('_n983', '').str.replace(r'\d+$', '', regex=True)
    d['tag_in_coloc_class'] = [b in str(cc).split(',') or (tg in str(cs).split(',')) for b, tg, cc, cs in
                               zip(base, d.tag, d.coloc_classes, d.coloc_subtypes)]
    d.to_csv(f'{R}/positive_controls.tsv', sep='\t', index=False)
    return d


def bulk_comparison():
    b = pd.read_csv(f'{W}/ref/09_twas.csv')
    b = b[(b.trait == 'SCZ_EUR') & (b.tissue == 'Brain_Cortex')].set_index('ensembl').z
    rows = []
    for tag in tags:
        for t in ['PGC3_EUR', 'SCZ25_EUR']:
            for m in ['susie', 'top1']:
                a = load_assoc(t, tag, m)
                if a is None:
                    continue
                j = a.set_index('gene').z.to_frame('sc').join(b.rename('bulk'), how='inner')
                rho, p = spearmanr(j.sc, j.bulk)
                rows.append(dict(tag=tag, trait=t, model=m, shared_genes=len(j), spearman=rho, p=p,
                                 sign_agree_where_bulk_abs_z_gt_3_3=float((np.sign(j.sc) == np.sign(j.bulk))[j.bulk.abs() > 3.29].mean())))
    d = pd.DataFrame(rows); d.to_csv(f'{R}/bulk_comparison.tsv', sep='\t', index=False); return d


def ext1_randomz():
    t = pd.read_csv(f'{W}/raw/Ext1_eqtl_top_assoc.tsv.gz', sep='\t', usecols=['feature', 'fixed_z', 'Random_Z', 'qval'])
    rho_all = spearmanr(t.fixed_z.abs(), t.Random_Z)[0]
    eg = t[t.qval < 0.05]
    d = pd.DataFrame([dict(features=len(t), spearman_absfixed_vs_random_all=rho_all,
                           egenes=len(eg), spearman_absfixed_vs_random_egenes=spearmanr(eg.fixed_z.abs(), eg.Random_Z)[0],
                           share_random_z_negative=(t.Random_Z < 0).mean(), share_fixed_z_negative=(t.fixed_z < 0).mean(),
                           spearman_signedfixed_vs_random=spearmanr(t.fixed_z, t.Random_Z)[0])])
    d.to_csv(f'{R}/ext1_randomz_check.tsv', sep='\t', index=False); return d


def n_sensitivity():
    if 'Ext1_n983' not in tags or 'Ext1' not in tags:
        return None
    rows = []
    for m in ['susie', 'top1']:
        a = pd.read_csv(f'{W}/weights/Ext1.{m}.tsv.gz', sep='\t'); b = pd.read_csv(f'{W}/weights/Ext1_n983.{m}.tsv.gz', sep='\t')
        j = a.merge(b, on=['gene', 'rsid'], suffixes=('_757', '_983'))
        r = dict(model=m, genes_757=a.gene.nunique(), genes_983=b.gene.nunique(), shared_snp_weights=len(j),
                 weight_pearson=np.corrcoef(j.weight_757, j.weight_983)[0, 1] if len(j) > 2 else np.nan)
        for t in ['SCZ25_EUR', 'MDD_EUR', 'ABCD_thin_EUR']:
            x = load_assoc(t, 'Ext1', m).merge(load_assoc(t, 'Ext1_n983', m), on='gene', suffixes=('_757', '_983'))
            r[f'{t}_z_pearson'] = np.corrcoef(x.z_757, x.z_983)[0, 1]
        rows.append(r)
    d = pd.DataFrame(rows); d.to_csv(f'{R}/n_sensitivity.tsv', sep='\t', index=False); return d


if __name__ == '__main__':
    out = {}
    out['ct'] = celltype_summary(); out['tr'], out['top'] = trait_summary()
    out['pc'] = positive_controls(); out['bulk'] = bulk_comparison(); out['rz'] = ext1_randomz(); out['ns'] = n_sensitivity()
    for k, v in out.items():
        if v is not None:
            print(f'== {k}'); print(v.head(60).to_string(index=False))
    # markdown fragments for README
    with open(f'{R}/tables.md', 'w') as f:
        f.write('### Per cell type\n\n' + md(out['ct']) + '\n')
        f.write('### Per trait x cell type x model\n\n' + md(out['tr']) + '\n')
        f.write('### Comparison with the bulk GTEx Brain_Cortex TWAS (SCZ)\n\n' + md(out['bulk']) + '\n')
        f.write('### Ext1: Random_Z vs fixed_z\n\n' + md(out['rz']) + '\n')
        if out['ns'] is not None:
            f.write('### Ext1 sample-size sensitivity (n = 757 vs 983)\n\n' + md(out['ns']) + '\n')
        pc = out['pc']
        key = pc[(pc.trait == 'SCZ25_EUR') & pc.tag_in_coloc_class & pc.fig3c]
        f.write('### Positive controls: Fig. 3c genes in their colocalising cell type(s) (SCZ25_EUR)\n\n' +
                md(key[['symbol', 'coloc_classes', 'coloc_subtypes', 'tag', 'model', 'z', 'p', 'w_flag']]) + '\n')
        sel = pc[(pc.trait == 'SCZ25_EUR') & pc.symbol.isin(['CNTN4', 'SLC4A10'])]
        f.write('### CNTN4 and SLC4A10 in every cell type (SCZ25_EUR)\n\n' +
                md(sel.pivot_table(index=['symbol', 'tag'], columns='model', values=['z', 'p']).reset_index()) + '\n')
        tp = out['top'].copy()
        tp['lab'] = [f"{sy if isinstance(sy, str) else g}({z:+.1f}{'!' if fl else ''})" for sy, g, z, fl in zip(tp.symbol, tp.gene, tp.z, tp.w_flag.astype(bool))]
        for t in TRAITS:
            x = tp[tp.trait == t].groupby(['tag', 'model'], sort=False).lab.apply(' '.join).reset_index()
            f.write(f'### Top 10 genes: {t}  (symbol(z); ! = w_flag; full table report/top10.tsv)\n\n' + md(x) + '\n')
    rd = open(f'{W}/README.md').read()
    a, b = '<!-- AUTO:RESULTS -->', '<!-- /AUTO:RESULTS -->'
    body = a + '\n\n## Results (generated by 04_report.py)\n\n' + open(f'{R}/tables.md').read() + '\n' + b
    rd = rd.split(a)[0] + body + (rd.split(b)[1] if b in rd else rd.split(a)[1])
    open(f'{W}/README.md', 'w').write(rd)
