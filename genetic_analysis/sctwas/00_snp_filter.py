"""SCTWAS stage 0: one-off SNP filter and GWAS z alignment (SCTWAS_INSTRUCTIONS.md, Stage 1 step 1).

Keep a 1000G EUR panel SNP when it is
  * an rsID on chr1-22, unique in the bim, single-base ACGT alleles;
  * outside the MHC (hg19 chr6:25-34 Mb; the panel bim is hg19);
  * MAF >= 0.01 in the panel (plink --freq);
  * not strand-ambiguous (A/T, C/G);
  * present in all four required GWAS with alleles that match the panel after flipping/complementing
    and a finite BETA/SE with SE > 0.
Outputs (under $W):
  ld/panel_keep.{bed,bim,fam}   panel restricted to the kept SNPs, --keep-allele-order, so bim col5 = panel A1
  ld/snps_keep.tsv              rsid chr pos a1 a2 maf (a1 = effect allele of every weight and aligned z)
  gwas/<trait>.z.tsv.gz         rsid z, z = BETA/SE aligned to panel A1 (raw sign; no trait sign flips)
  gwas/snp_filter_counts.tsv    attrition table
PGC3_EUR (optional comparison trait) is aligned to the kept set but does not enter the intersection.
"""
import os, subprocess, sys
import numpy as np, pandas as pd

REPO = '/home/rajd2/rds/hpc-work/abcd_development'
W = f'{REPO}/genetic_analysis/work/sctwas'
PLINK = f'{REPO}/legacy/hpc/work/bin/plink'
PANEL = '/rds/user/rajd2/hpc-work/magma/reference_data/g1000_eur'
AB = f'{REPO}/genetic_analysis/work/results_70tab_hcp/scan_1lmm/assoc_eur'
GWAS = {  # name: (path, kind)  -- effect allele = A1 in every file (checked, README.md)
    'ABCD_thin_EUR': (f'{AB}/global_slope_1lmm.sumstats.tsv.gz', 'abcd'),
    'ABCD_base_EUR': (f'{AB}/baseline_thickness_1lmm.sumstats.tsv.gz', 'abcd'),
    'SCZ25_EUR': (f'{REPO}/genetic_analysis/work/scores_scz2025/gwas/SCZ25_EUR.ma', 'ma'),
    'MDD_EUR': (f'{REPO}/legacy/hpc_v2/work/results_v2/prs_final/gwas/MDD_eur.ma', 'ma'),
}
OPTIONAL = {'PGC3_EUR': (f'{REPO}/legacy/hpc_v2/work/results_v2/prs_final/gwas/SCZ_eur.ma', 'ma')}
COMP = str.maketrans('ACGT', 'TGCA')


def read_gwas(path, kind):
    if kind == 'abcd':
        d = pd.read_csv(path, sep='\t', usecols=['SNP', 'A1', 'A2', 'BETA', 'SE'], dtype={'SNP': str, 'A1': str, 'A2': str})
    else:
        d = pd.read_csv(path, sep=r'\s+', usecols=['SNP', 'A1', 'A2', 'b', 'se'], dtype={'SNP': str, 'A1': str, 'A2': str})
        d = d.rename(columns={'b': 'BETA', 'se': 'SE'})
    d['A1'] = d.A1.str.upper(); d['A2'] = d.A2.str.upper()
    return d


def align(d, panel):
    """Return Series rsid -> z aligned to panel A1 (drops unmatched / bad SE / duplicated rsIDs)."""
    n0 = len(d)
    d = d[np.isfinite(d.BETA) & np.isfinite(d.SE) & (d.SE > 0)]
    d = d[~d.SNP.duplicated(keep=False)]
    m = d.merge(panel[['rsid', 'a1', 'a2']], left_on='SNP', right_on='rsid', how='inner')
    c1 = m.A1.str.translate(COMP); c2 = m.A2.str.translate(COMP)
    s = np.select([(m.A1 == m.a1) & (m.A2 == m.a2), (m.A1 == m.a2) & (m.A2 == m.a1),
                   (c1 == m.a1) & (c2 == m.a2), (c1 == m.a2) & (c2 == m.a1)], [1., -1., 1., -1.], np.nan)
    m = m.assign(z=s * m.BETA / m.SE)[np.isfinite(s)]
    return m.set_index('rsid').z, dict(rows=n0, in_panel=int(len(d.merge(panel[['rsid']], left_on='SNP', right_on='rsid'))), allele_match=len(m))


def main():
    os.makedirs(f'{W}/ld', exist_ok=True); os.makedirs(f'{W}/gwas', exist_ok=True)
    tmp = f'{W}/tmp'
    cnt = []
    bim = pd.read_csv(f'{PANEL}.bim', sep='\t', header=None, names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'],
                      dtype={'chr': str, 'rsid': str, 'a1': str, 'a2': str})
    cnt.append(('panel_bim', len(bim)))
    bim = bim[bim.chr.isin([str(i) for i in range(1, 23)]) & bim.rsid.str.match(r'^rs\d+$')]
    bim = bim[bim.a1.str.fullmatch('[ACGT]') & bim.a2.str.fullmatch('[ACGT]')]
    bim = bim[~bim.rsid.duplicated(keep=False)]
    cnt.append(('autosomal_rsid_snv_unique', len(bim)))
    bim = bim[~((bim.chr == '6') & (bim.pos >= 25_000_000) & (bim.pos <= 34_000_000))]
    cnt.append(('minus_MHC_hg19_chr6_25_34Mb', len(bim)))
    if not os.path.exists(f'{tmp}/panel_freq.frq'):
        subprocess.run([PLINK, '--bfile', PANEL, '--freq', '--keep-allele-order', '--out', f'{tmp}/panel_freq'], check=True,
                       stdout=subprocess.DEVNULL)
    frq = pd.read_csv(f'{tmp}/panel_freq.frq', sep=r'\s+', usecols=['SNP', 'MAF'], dtype={'SNP': str})
    frq = frq[~frq.SNP.duplicated(keep=False)]
    bim = bim.merge(frq, left_on='rsid', right_on='SNP').drop(columns='SNP').rename(columns={'MAF': 'maf'})
    bim = bim[bim.maf >= 0.01]
    cnt.append(('panel_MAF_ge_0.01', len(bim)))
    amb = (bim.a1 + bim.a2).isin(['AT', 'TA', 'CG', 'GC'])
    bim = bim[~amb]
    cnt.append(('minus_strand_ambiguous', len(bim)))

    Z = {}
    for name, (path, kind) in GWAS.items():
        if not os.path.exists(path):
            sys.exit(f'FATAL: missing GWAS {path}')
        z, info = align(read_gwas(path, kind), bim)
        Z[name] = z
        print(name, info, flush=True)
        cnt.append((f'gwas_{name}_rows', info['rows'])); cnt.append((f'gwas_{name}_allele_matched_in_filtered_panel', info['allele_match']))
    keep = set(bim.rsid)
    for z in Z.values():
        keep &= set(z.index)
    bim = bim[bim.rsid.isin(keep)]
    cnt.append(('KEPT_intersection_4_gwas', len(bim)))
    print('kept', len(bim), flush=True)

    bim[['rsid']].to_csv(f'{tmp}/keep.rsid', index=False, header=False)
    subprocess.run([PLINK, '--bfile', PANEL, '--extract', f'{tmp}/keep.rsid', '--keep-allele-order', '--make-bed',
                    '--out', f'{W}/ld/panel_keep'], check=True, stdout=subprocess.DEVNULL)
    nb = pd.read_csv(f'{W}/ld/panel_keep.bim', sep='\t', header=None, names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'], dtype={'rsid': str})
    chk = nb.merge(bim[['rsid', 'a1', 'a2']], on='rsid', suffixes=('', '_o'))
    assert len(nb) == len(bim) and (chk.a1 == chk.a1_o).all() and (chk.a2 == chk.a2_o).all(), 'panel allele order changed'
    out = nb[['rsid', 'chr', 'pos', 'a1', 'a2']].merge(bim[['rsid', 'maf']], on='rsid')
    out.to_csv(f'{W}/ld/snps_keep.tsv', sep='\t', index=False)

    for name, z in Z.items():
        z[z.index.isin(keep)].rename('z').reset_index().to_csv(f'{W}/gwas/{name}.z.tsv.gz', sep='\t', index=False, float_format='%.6g')
    for name, (path, kind) in OPTIONAL.items():
        if os.path.exists(path):
            z, info = align(read_gwas(path, kind), bim)
            z.rename('z').reset_index().to_csv(f'{W}/gwas/{name}.z.tsv.gz', sep='\t', index=False, float_format='%.6g')
            cnt.append((f'optional_{name}_in_kept_set', len(z)))
    pd.DataFrame(cnt, columns=['step', 'n']).to_csv(f'{W}/gwas/snp_filter_counts.tsv', sep='\t', index=False)
    print(pd.DataFrame(cnt, columns=['step', 'n']).to_string(index=False))


if __name__ == '__main__':
    main()
