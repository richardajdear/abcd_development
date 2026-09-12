# GWAS gene sets and gene-level statistics for H2

## NOTES
- Built by `code/03_gwas_gene_sets.py` (re-runnable; mygene.info responses cached in `*_mygene_raw.json`). No enrichment tests run here.
- **Set definitions replicate the hpc MAGMA run exactly** (`hpc/work/make_prioritised_genesets.py`, `make_mdd_genesets.py`), read from the
  same source tables (`hpc/work/genesets/scz2022_ExtendedDataTable1.xlsx`, `mdd2024_TableS21.xlsx`), but in **symbol space**.
  The hpc files (`scz_prioritised.txt`, `mdd_highconf.txt`) are Entrez-keyed and restricted to genes present in NCBI37.3 `gene.loc`,
  so they are smaller (e.g. 104 vs 120 prioritised; 462 vs 685 pool). Copies with symbols added: `magma_setannot_{SCZ,MDD}_hpc.tsv`.
- **MDD = Adams et al. 2025 (Cell, doi 10.1016/j.cell.2024.12.002; PGC MDD 2025).** The hpc code calls it the "2024 Cell depression GWAS";
  it is the same paper. Table S21 supplies a 308-row high-confidence list (295 unique symbols — 13 rows are duplicate symbols with
  different Ensembl ids) and a 4,600-row pool (4,180 unique symbols). The MAGMA MDD gene-level run used
  `pgc-mdd2025_no23andMe_div_v3-49-46-01` sumstats (median N_eff ≈ 1,244,628), i.e. the same study. Howard 2019 lists
  are included as legacy comparators only.
- SCZ MAGMA gene-level run used PGC3 SCZ wave 3 primary autosomal sumstats (Trubetskoy 2022; median N ≈ 168,263).
- Symbols absent from the AHBA ds0 universe: 2,681 of 4,697 set symbols; mygene.info (symbol/alias, unique hit,
  target in ds0, no collision) renamed only **11** (e.g. GPR98→ADGRV1, BAI3→ADGRB3, PCNXL3→PCNX3; `geneset_symbol_map.csv`).
  The remainder are lincRNA / clone-based names (`RP11-…`, `AC0…`) with no AHBA probe — the sets are therefore biased to protein-coding
  genes when intersected with the expression universe, as in the hpc run.
- Files inspected in `AHBA/data/gwas/` and **not** turned into sets: `trubetskoy2022.csv` / `_extended.csv` (identical to ED Table 1 / ST12 —
  the xlsx was used), `trubetskoy2022_finemap.csv` (SNP-level FINEMAP PPs, 20,765 SNPs / 653 genes — the gene-level FINEMAP flag comes from
  ED Table 1), `trubetskoy2022_allgenes_{PE,fetal,eQTLGen}.csv` (SMR gene-level results, all tested genes; the SMR-prioritised flag comes
  from ED Table 1), `FB_SCZ.genes.csv` (an older Ensembl-keyed MAGMA run, median N ≈ 69k, superseded by `hpc` SCZ.genes.out).
- **Background universe**: `brain_background_ds{0,25,50}.txt` = AHBA DS-level gene list ∩ genes with a MAGMA Z in **both** SCZ and MDD
  gene-level files: 13,763 / 10,322 / 6,931 genes. Use these for the PNAS-style permutation
  null and as the MAGMA gene-property gene universe so both tests share one universe.

## Symbol-space sets (`<set>.txt`, one symbol per line; also `gene_sets_long.tsv`)
| set | n | in ds0 | in ds25 | in ds50 | in background ds0/25/50 | definition |
|---|---|---|---|---|---|---|
| `SCZ_prioritised` | 120 | 99 | 85 | 63 | 92 / 78 / 58 | Trubetskoy2022 ED Table 1 (all 120 prioritised) |
| `SCZ_prioritised_pc` | 106 | 96 | 82 | 60 | 92 / 78 / 58 | ED Table 1, protein_coding only |
| `SCZ_prio_finemap` | 70 | 63 | 58 | 42 | 58 / 53 / 39 | ED Table 1, FINEMAP.priority.gene == 1 |
| `SCZ_prio_smr` | 55 | 41 | 31 | 23 | 38 / 28 / 20 | ED Table 1, SMR.priority.gene == 1 |
| `SCZ_prio_rare` | 7 | 7 | 6 | 5 | 7 / 6 / 5 | ED Table 1, Rare.priority.gene == 1 |
| `SCZ_locus_pool` | 682 | 422 | 333 | 238 | 394 / 311 / 221 | Trubetskoy2022 ST12 all criteria (685-gene candidate pool at GWS loci) |
| `SCZ_pool_not_prio` | 562 | 323 | 248 | 175 | 302 / 233 / 163 | ST12 minus ED Table 1 |
| `MDD_highconf` | 295 | 208 | 189 | 148 | 202 / 185 / 145 | Adams2025 Table S21 High-confidence Gene List (308; finemapping, expression or protein) |
| `MDD_hc_finemap` | 225 | 155 | 142 | 115 | 150 / 139 / 112 | High-confidence list, Fine_mapping == TRUE |
| `MDD_hc_expression` | 75 | 57 | 49 | 32 | 55 / 48 / 32 | High-confidence list, Expression == TRUE |
| `MDD_hc_protein` | 10 | 10 | 9 | 9 | 10 / 9 / 9 | High-confidence list, Protein == TRUE |
| `MDD_pool` | 4180 | 1713 | 1322 | 941 | 1626 / 1263 / 903 | Adams2025 Table S21 all mapped genes (4,600; any of 7 methods incl. fastBAT/H-MAGMA) |
| `MDD_pool_not_hc` | 3885 | 1505 | 1133 | 793 | 1424 / 1078 / 758 | Table S21 pool minus high-confidence |
| `MDD_howard2019_32` | 32 | 24 | 19 | 15 | 23 / 18 / 15 | AHBA/data/gwas/hammerschlag2020_howard2019_trubetskoy2022.csv MDD column (Howard 2019; 32 genes, exact derivation not recorded — legacy Dear2024 set) |
| `SCZ_trubetskoy_smr101` | 101 | 63 | 51 | 37 | 54 / 43 / 30 | same file, SCZ column = Trubetskoy2022 SMR P+F genes (trubetskoy2022_data.csv; legacy Dear2024 set) |
| `MDD_howard2019_magma` | 269 | 226 | 192 | 133 | 214 / 181 / 125 | Howard2019 Table S9: MAGMA genome-wide significant genes (P<2.8e-6) |

## hpc Entrez set-annot files (as used by `prioritised_gsa.sbatch`; `magma_setannot_*_hpc.tsv`)
| set | Entrez ids / with symbol / in ds0 |
|---|---|
| SCZ_locus_pool | scz_prioritised.txt: 462 Entrez ids, 462 with symbol, 406 in ds0 |
| SCZ_pool_not_prio | scz_prioritised.txt: 358 Entrez ids, 358 with symbol, 310 in ds0 |
| SCZ_prio_finemap | scz_prioritised.txt: 64 Entrez ids, 64 with symbol, 61 in ds0 |
| SCZ_prio_rare | scz_prioritised.txt: 7 Entrez ids, 7 with symbol, 7 in ds0 |
| SCZ_prio_smr | scz_prioritised.txt: 44 Entrez ids, 44 with symbol, 39 in ds0 |
| SCZ_prioritised | scz_prioritised.txt: 104 Entrez ids, 104 with symbol, 96 in ds0 |
| SCZ_prioritised_pc | scz_prioritised.txt: 104 Entrez ids, 104 with symbol, 96 in ds0 |
| MDD_hc_expression | mdd_highconf.txt: 56 Entrez ids, 56 with symbol, 53 in ds0 |
| MDD_hc_finemap | mdd_highconf.txt: 160 Entrez ids, 160 with symbol, 147 in ds0 |
| MDD_hc_protein | mdd_highconf.txt: 10 Entrez ids, 10 with symbol, 10 in ds0 |
| MDD_highconf | mdd_highconf.txt: 213 Entrez ids, 213 with symbol, 198 in ds0 |
| MDD_pool | mdd_highconf.txt: 1885 Entrez ids, 1884 with symbol, 1593 in ds0 |
| MDD_pool_not_hc | mdd_highconf.txt: 1672 Entrez ids, 1671 with symbol, 1395 in ds0 |

## MAGMA gene-level statistics (`magma_{SCZ,MDD}_genes.tsv`)
Columns: GENE (Entrez), symbol, symbol_source (`symbol2entrez` = `data/symbol2entrez.csv`; `mygene` = mygene.info entrezgene lookup),
CHR, START, STOP, NSNPS, NPARAM, N, ZSTAT, P. `ens2entrez_mygene.csv` was not needed (GENE is already Entrez).

| trait | n genes | with symbol | via symbol2entrez | via mygene | symbol in ds0 |
|---|---|---|---|---|---|
| SCZ | 18,447 | 18,321 (99.3%) | 15,555 | 2,766 | 13,789 |
| MDD | 19,040 | 18,922 (99.4%) | 16,150 | 2,772 | 14,271 |

No duplicate symbols after mapping. For the gene-property test, join PLS gene weights to these tables on `symbol`.
