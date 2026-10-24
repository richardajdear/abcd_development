"""SCTWAS stage 1a (per cell type): SingleBrain eQTL summary statistics -> per-chromosome z inputs for SuSiE-RSS.

usage: 01_prep_eqtl.py <CT>
  * eGenes = top_assoc rows with qval < 0.05 (the authors' definition); features tested = top_assoc rows.
  * full_assoc is streamed (pyarrow) keeping eGene x kept-SNP rows (ld/snps_keep.tsv from stage 0).
  * Signed statistic = fixed_z relative to `Allele` (Random_Z is unsigned and never used as z).
    Rows with fixed_sd == 0 (SNP untested in every cohort) or non-finite fixed_z are dropped.
  * fixed_z is aligned to the panel A1: +z when Allele == a1, -z when Allele == a2, complement-matched likewise
    (ambiguous SNPs are already gone). {ref, alt} must equal {a1, a2} up to strand, and chr must agree.
Outputs: stage1/<CT>/zin_chr<k>.tsv.gz (gene_version, rsid, z), stage1/<CT>/prep_summary.tsv,
         stage1/<CT>/egenes.tsv (feature, lead variant, fixed_z, Random_Z, qval)
"""
import os, sys
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.csv as pcsv, pyarrow.compute as pc

W = '/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas'
COMP = str.maketrans('ACGT', 'TGCA')
COLS = ['feature', 'variant_id', 'chr', 'ref', 'alt', 'Allele', 'fixed_sd', 'fixed_z']


def align(d, snps, cnt):
    untested = (d.fixed_sd == 0) | ~np.isfinite(d.fixed_z) | ~np.isfinite(d.fixed_sd)
    cnt['dropped_untested'] += int(untested.sum())
    d = d[~untested].merge(snps, left_on='variant_id', right_on='rsid', how='inner')
    chr_ok = d.chr_x.str.replace('chr', '', regex=False) == d.chr_y
    A, r, a = d.Allele.str.upper(), d.ref.str.upper(), d.alt.str.upper()
    cA, cr, ca = A.str.translate(COMP), r.str.translate(COMP), a.str.translate(COMP)
    pair_ok = ((r == d.a1) & (a == d.a2)) | ((r == d.a2) & (a == d.a1)) | ((cr == d.a1) & (ca == d.a2)) | ((cr == d.a2) & (ca == d.a1))
    s = np.select([A == d.a1, A == d.a2, cA == d.a1, cA == d.a2], [1., -1., 1., -1.], np.nan)
    ok = chr_ok & pair_ok & np.isfinite(s)
    cnt['dropped_chr_mismatch'] += int((~chr_ok).sum()); cnt['dropped_allele_mismatch'] += int((chr_ok & ~ok).sum())
    return d.assign(z=s * d.fixed_z)[ok][['feature', 'rsid', 'chr_y', 'z']]


def main(ct):
    out = f'{W}/stage1/{ct}'; os.makedirs(out, exist_ok=True)
    top = pd.read_csv(f'{W}/raw/{ct}_eqtl_top_assoc.tsv.gz', sep='\t',
                      usecols=['feature', 'variant_id', 'chr', 'fixed_z', 'Random_Z', 'qval'])
    eg = top[top.qval < 0.05]
    eg.to_csv(f'{out}/egenes.tsv', sep='\t', index=False)
    egenes = pa.array(sorted(set(eg.feature)))
    snps = pd.read_csv(f'{W}/ld/snps_keep.tsv', sep='\t', dtype={'rsid': str, 'chr': str, 'a1': str, 'a2': str})
    keep_rs = pa.array(snps.rsid.values)
    print(ct, 'features tested', len(top), 'eGenes', len(eg), flush=True)

    reader = pcsv.open_csv(f'{W}/raw/{ct}_eqtl_full_assoc.tsv.gz',
                           read_options=pcsv.ReadOptions(block_size=1 << 26),
                           parse_options=pcsv.ParseOptions(delimiter='\t'),
                           convert_options=pcsv.ConvertOptions(include_columns=COLS, column_types={
                               'feature': pa.string(), 'variant_id': pa.string(), 'chr': pa.string(), 'ref': pa.string(),
                               'alt': pa.string(), 'Allele': pa.string(), 'fixed_sd': pa.float64(), 'fixed_z': pa.float64()}))
    parts, cnt = [], dict(full_rows=0, egene_rows=0, egene_rows_kept_snps=0, dropped_untested=0,
                          dropped_chr_mismatch=0, dropped_allele_mismatch=0)
    for b in reader:
        cnt['full_rows'] += b.num_rows
        b = pa.Table.from_batches([b])
        b = b.filter(pc.is_in(b['feature'], value_set=egenes)); cnt['egene_rows'] += b.num_rows
        b = b.filter(pc.is_in(b['variant_id'], value_set=keep_rs)); cnt['egene_rows_kept_snps'] += b.num_rows
        if b.num_rows:
            parts.append(align(b.to_pandas(), snps, cnt))
    d = pd.concat(parts, ignore_index=True).drop_duplicates(['feature', 'rsid'])
    for c, g in d.groupby('chr_y'):
        g[['feature', 'rsid', 'z']].rename(columns={'feature': 'gene_version'}).to_csv(
            f'{out}/zin_chr{c}.tsv.gz', sep='\t', index=False, float_format='%.6g')
    summ = dict(ct=ct, features_tested=len(top), egenes=len(eg), **cnt, final_rows=len(d), genes_with_snps=d.feature.nunique(),
                median_snps_per_gene=float(d.groupby('feature').size().median()))
    pd.DataFrame([summ]).to_csv(f'{out}/prep_summary.tsv', sep='\t', index=False)
    print(summ, flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
