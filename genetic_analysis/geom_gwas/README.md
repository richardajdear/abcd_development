# scan_geom — D1 step 9: GWAS of the area, volume and thickness+area slopes (CSD3, 2026-10-27/28)

Brief: `directions/d1_multimodal/CSD3_GWAS_INSTRUCTIONS.md` (commit 112d15f). Scripts: `genetic_analysis/geom_gwas/`.
Same GENESIS single-LMM pipeline, covariates (sex, site, baseline_age, n_visits, 10 PC-AiR PCs), kinship, sample sets,
MAF 0.01 / MAC 20 filters, GREML design and MAGMA references/annotations as `scan_1lmm` (HCP). Nothing outside
`scan_geom/` was written. Input `directions/d1_multimodal/work/gwas_pheno_geom.tsv` (8,716 children, raw BLUPs).

**Status: COMPLETE — steps 1–5 (five slopes, both arms), step 6 (EUR intercepts) and step 7 (LDSC).**
The rds `hpc-work` volume hit its 1 TB quota at 23:19 on 2026-10-27 and one intercept assoc task (`int_thk_vol_area`
chr2) failed with "Disk quota exceeded"; its partial outputs were deleted, the task was re-run on 2026-10-28 after the
storage move below, and the intercepts were then collected (22/22 chromosomes each) and gene-analysed.

**Storage (2026-10-28).** All of `genetic_analysis/work/` except git-tracked files (161 GB, 808 items) now lives on
`/rds/project/rds-Pb9UGUlrwWc/rajd2_abcd_work/`, mirrored path for path and symlinked back (absolute links), so every
script path is unchanged. Directories holding no tracked file were moved whole; in directories that hold tracked
summary tables only untracked files ≥ 1 MB were moved (git cannot track files behind a symlinked directory). Each item
was copied, checked for identical file count and bytes, then replaced by the link (plan/executor and per-item log in
`/rds/project/rds-Pb9UGUlrwWc/rajd2_logs/{plan_move.py,do_move.sbatch,moves.tsv,move_39456626.log}`).
`hpc-work` went from 1,023 GB to 873 GB used. New outputs written into a symlinked directory land on Pb9; new files in a
still-real directory (e.g. `results_70tab*/assoc/`) land on `hpc-work`.

## Gates

| gate | test | value | result |
|---|---|---|---|
| A | r(slope_thk, global_slope_1lmm), matched set | **1.0000** (n 8,596 / 8,596 matched; 0 missing in every column; 120 export rows not in the reference set) | PASS |
| A | r(int_thk, baseline_thickness_1lmm) | 1.0000 | — |
| A | laptop slope r matrix (10 pairs) | all reproduced to 2 dp | PASS |
| B | EUR per-SNP r(z) slope_thk vs global_slope_1lmm | **1.00000** over 7,470,569 SNPs; λGC 1.0041 vs 1.0041 | PASS |
| B′ | pooled per-SNP r(z), same | 1.00000 over 9,426,204 SNPs; λGC 1.0268 vs 1.0268 | PASS |
| — | GREML slope_thk vs scan_1lmm | 0.1833 ± 0.0451 vs 0.1833 ± 0.0451 | identical |
| — | MAGMA gene-ZSTAT r slope_thk vs `magma_panel` (EUR) / `magma_pooled` (pooled) global_slope_1lmm | 1.0 / 1.0 | PASS |

The EUR null model for slope_thk reproduces scan_1lmm's variance components to 7 digits (0.33392 / 0.43314), n 4,308.

## GWAS (score test; n = mean N per SNP)

| arm | phenotype | n | λGC | min P | P<5e-8 | P<1e-6 |
|---|---|---|---|---|---|---|
| EUR | slope_thk | 4,205 | 1.0041 | 6.9e-07 | 0 | 2 |
| EUR | slope_area | 4,205 | 1.0081 | 5.6e-08 | 0 | 3 |
| EUR | slope_vol | 4,205 | 1.0152 | 5.3e-07 | 0 | 1 |
| EUR | slope_thk_area | 4,205 | 1.0158 | 2.1e-07 | 0 | 5 |
| EUR | slope_thk_vol_area | 4,205 | 1.0156 | 4.1e-07 | 0 | 7 |
| pooled | slope_thk | 8,396 | 1.0268 | 4.2e-08 | 1 | 8 |
| pooled | slope_area | 8,396 | 0.9835 | 2.1e-07 | 0 | 18 |
| pooled | slope_vol | 8,396 | 1.0137 | 4.3e-07 | 0 | 2 |
| pooled | slope_thk_area | 8,396 | 1.0240 | 2.1e-08 | 2 | 12 |
| pooled | slope_thk_vol_area | 8,396 | 1.0176 | 1.1e-07 | 0 | 12 |

Intercepts (EUR arm, positive controls, step 6):

| phenotype | n | λGC | min P | P<5e-8 | P<1e-6 |
|---|---|---|---|---|---|
| int_thk | 4,205 | 1.0215 | 3.7e-07 | 0 | 5 |
| int_area | 4,205 | 0.9975 | 1.2e-09 | 1 | 3 |
| int_vol | 4,205 | 1.0091 | 1.4e-08 | 2 | 2 |
| int_thk_area | 4,205 | 1.0183 | 5.7e-07 | 0 | 4 |
| int_thk_vol_area | 4,205 | 1.0149 | 3.2e-07 | 0 | 5 |

int_thk reproduces `scan_1lmm/assoc_eur/baseline_thickness_1lmm` (λ 1.0215, same min P, MAGMA gene-Z r 1.0).
The brief's positive-control expectation is only partly met: baseline area gives a genome-wide hit (P 1.2e-09) but
**λGC 0.998, not > 1**; at n ≈ 4,200 genome-wide inflation is not a sensitive heritability readout (gene-level λ is
1.11–1.14 for every intercept, against 1.03–1.06 for the slopes). GREML was not run for intercepts (brief: slopes only).

GENESIS null-model share of variance on the kinship term (EUR, prop_kin): thk 0.44, area **0.81**, vol 0.68,
thk_area 0.57, thk_vol_area 0.59 — area's large kinship component sits beside a GREML h² of zero (below), i.e. it is
carried by close relatives (shared family environment and/or rare variation), not by common SNPs in unrelateds.

## Heritability

| phenotype | GREML h² ± se (p), n 6,011 unrelated | LDSC h² ± se (z), EUR | LDSC intercept |
|---|---|---|---|
| slope_thk | 0.183 ± 0.045 (1.3e-05) | 0.031 ± 0.093 (0.34) | 1.005 |
| slope_area | **0.000** ± 0.044 (0.5) | 0.209 ± 0.093 (2.24) | 0.993 |
| slope_vol | 0.047 ± 0.044 (0.14) | 0.142 ± 0.104 (1.36) | 1.003 |
| slope_thk_area | 0.179 ± 0.045 (2.5e-05) | 0.069 ± 0.089 (0.77) | 1.004 |
| slope_thk_vol_area | 0.124 ± 0.045 (0.0024) | 0.094 ± 0.092 (1.02) | 1.004 |

GREML is the canonical estimator (rule 11); LDSC at n ≈ 4,200 has se ≈ 0.09 and is report-only (brief step 7). The two
disagree in opposite directions for thk and area; all LDSC z < 4 (rule 13). Bivariate rg was **not run**:
`step5_reml.sbatch` has no `--reml-bivar` mode (brief: skip otherwise), and with slope_area at the h² = 0 boundary the
rg would be undefined anyway.

## MAGMA gene analysis (snp-wise=mean; EUR: 1000G EUR + 35/10 union annotation; pooled: in-sample ABCD reference)

| arm | phenotype | genes | λ(gene P) | gene-Z Spearman with slope_thk |
|---|---|---|---|---|
| EUR | slope_thk | 18,283 | 1.031 | 1 |
| EUR | slope_area | 18,283 | 1.036 | 0.028 |
| EUR | slope_vol | 18,283 | 1.063 | 0.435 |
| EUR | slope_thk_area | 18,283 | 1.031 | 0.816 |
| EUR | slope_thk_vol_area | 18,283 | 1.061 | 0.753 |
| pooled | slope_thk | 17,896 | 1.033 | 1 |
| pooled | slope_area | 17,896 | 1.082 | 0.014 |
| pooled | slope_vol | 17,896 | 1.058 | 0.440 |
| pooled | slope_thk_area | 17,896 | 1.081 | 0.830 |
| pooled | slope_thk_vol_area | 17,896 | 1.067 | 0.765 |

Intercepts (EUR): genes 18,283 each; λ(gene P) int_thk 1.135, int_area 1.115, int_vol 1.124, int_thk_area 1.126,
int_thk_vol_area 1.131; gene-Z Spearman with slope_thk −0.006 to +0.003.

Top 10 genes (gene P):

- EUR slope_thk: BAMBI 8.6e-06, TRIM6 1.6e-05, TRIM6-TRIM34 7.0e-05, PAM 9.6e-05, MYO10 2.1e-04, PPIP5K2 2.5e-04, OR52B6 4.2e-04, UBR2 4.6e-04, CSN1S1 5.4e-04, OLFML3 5.6e-04
- EUR slope_area: ZNF521 8.5e-05, METRNL 1.1e-04, IPO8 2.4e-04, DROSHA 3.6e-04, TMPRSS11D 4.2e-04, RASD1 4.4e-04, EEPD1 4.5e-04, CDH2 5.0e-04, STAR 5.1e-04, LETM2 5.2e-04
- EUR slope_vol: ZNF521 3.6e-06, GPR39 3.7e-05, ZNF57 9.1e-05, TRIM6 1.8e-04, OR52B6 2.3e-04, FAM179B 2.8e-04, XPO7 3.2e-04, CTSH 3.2e-04, ZNF469 3.3e-04, CCDC181 3.9e-04
- EUR slope_thk_area: TRIM6 3.9e-05, EXOC2 7.0e-05, TRIM6-TRIM34 1.6e-04, KHDRBS3 2.1e-04, OR52B6 2.5e-04, TCAIM 2.7e-04, FAM179B 3.2e-04, TOPAZ1 3.2e-04, OR52H1 3.6e-04, OLFML3 3.8e-04
- EUR slope_thk_vol_area: TRIM6 3.3e-05, OR52B6 1.3e-04, ZNF521 1.3e-04, TRIM6-TRIM34 1.4e-04, FAM179B 1.4e-04, OR52H1 2.2e-04, KHDRBS3 2.5e-04, EXOC2 2.8e-04, KLHL28 2.9e-04, SETMAR 3.4e-04
- EUR int_thk: EPB41L2 1.0e-04, RPL27 1.3e-04, SLC22A24 1.5e-04, PTGES3L 1.6e-04, ZNF345 1.7e-04, CST1 1.8e-04, IFI35 2.1e-04, SWI5 2.1e-04, PTGES3L-AARSD1 2.3e-04, CST4 2.6e-04
- EUR int_area: RC3H2 5.8e-05, PPP1R9A 8.6e-05, ZBTB6 9.2e-05, PDCL 1.1e-04, ZBTB26 1.4e-04, RABGAP1 1.7e-04, HMGA2 1.8e-04, HECTD1 2.1e-04, LRRC4B 2.2e-04, GPR21 2.3e-04
- EUR int_vol: UHRF1 4.9e-05, ZFP3 2.9e-04, PDCL 3.4e-04, NOS2 3.9e-04, ARRDC5 3.9e-04, CCT8L2 4.5e-04, COL27A1 4.5e-04, RC3H2 4.7e-04, SPOPL 4.8e-04, PA2G4 4.8e-04
- EUR int_thk_area: ZFP3 3.1e-05, SLC52A1 2.1e-04, ACTR1B 2.9e-04, CCT8L2 3.0e-04, COX5B 3.4e-04, ANKRD36B 3.4e-04, LAT 3.7e-04, FAM178B 4.1e-04, TMEM247 4.4e-04, AURKAIP1 4.6e-04
- EUR int_thk_vol_area: ZFP3 5.6e-05, UHRF1 1.3e-04, CCT8L2 2.6e-04, SLC52A1 3.8e-04, NOS2 4.0e-04, LAT 4.1e-04, MED15 4.8e-04, PHYKPL 4.9e-04, TMEM247 5.1e-04, SPNS1 6.2e-04
- pooled slope_thk: FZD4 1.9e-04, TBCC 3.1e-04, IL1B 3.5e-04, ZNRF3 4.0e-04, PRPH2 4.8e-04, DPF3 5.1e-04, RASSF2 5.2e-04, NAPB 6.9e-04, GZF1 7.0e-04, SBNO1 7.6e-04
- pooled slope_area: SLC8A1 7.3e-06, ZFAT 1.4e-04, TMEM237 1.4e-04, PUM2 2.3e-04, ALS2CR11 2.3e-04, MPP4 2.5e-04, FAM208A 5.5e-04, GALNT1 5.5e-04, BCL2L11 8.0e-04, CCDC102B 9.1e-04
- pooled slope_vol: ZFAT 3.3e-05, LAMB3 7.1e-05, IGDCC4 7.2e-05, DPP8 7.3e-05, LEUTX 1.5e-04, MRPL10 2.1e-04, C19orf25 2.2e-04, SCRN2 2.2e-04, OSBPL7 2.2e-04, LRRC46 2.3e-04
- pooled slope_thk_area: LAMB3 6.7e-05, RASSF2 2.2e-04, C19orf25 3.0e-04, C22orf31 4.3e-04, ADAM7 5.4e-04, FAM129B 6.7e-04, PEBP4 6.9e-04, PLXNA2 9.0e-04, CAMK1G 9.3e-04, ZFAT 1.0e-03
- pooled slope_thk_vol_area: LAMB3 2.4e-05, C19orf25 1.2e-04, ZFAT 1.9e-04, SULF2 4.5e-04, FAM129B 5.3e-04, CAMK1G 6.8e-04, ZDHHC21 8.9e-04, COL12A1 9.0e-04, DPP8 1.0e-03, IGDCC4 1.1e-03

No gene reaches Bonferroni (2.7e-06) in any cell. Gene-set / gene-property tests were deliberately not run here (local `gwas_cells/48*`).

## Compute (Slurm, icelake, vertes-sl3-cpu)

| stage | tasks | core-h | wall |
|---|---|---|---|
| null models (EUR + pooled slopes + EUR ints) | 15 | 1.6 | ~1 min each |
| GENESIS assoc (EUR slopes, pooled slopes, EUR ints) | 330 | 269 | ≤ 24 min/task; EUR slopes 21:21–21:45, pooled 22:05–22:44 |
| collect | 2 | 4.4 | ≤ 17 min |
| GREML | 5 | 0.3 | 20 s each |
| MAGMA EUR (slopes + intercepts) | 10 | ~17.5 | 10.5–11 min each |
| re-run of failed intercept task + its collect | 2 | ~1 | 13 min each |
| MAGMA pooled | 5 | 285.5 | 2 h 54 min each (23:00–01:54) |
| LDSC h² | 1 | 0.7 | 4 min |

Total ≈ 580 core-hours; slopes end to end ≈ 4.5 h wall (pooled MAGMA the long pole); intercepts finished 2026-10-28
after the quota stall.

## Deliverables (gene-level summaries only; SNP-level sumstats stay on CSD3)

```
rsync -av --prune-empty-dirs \
  --include='*/' \
  --include='magma_eur/genes/*.genes.out' --include='magma_eur/genes/*.genes.raw' --include='magma_eur/genes/*.log' \
  --include='magma_pooled/genes/*.genes.out' --include='magma_pooled/genes/*.genes.raw' --include='magma_pooled/genes/*.log' \
  --include='assoc/gwas_summary.tsv' --include='assoc_eur/gwas_summary.tsv' \
  --include='reml_imp_pooled/reml_summary.tsv' --include='ldsc_eur/ldsc_h2_summary.tsv' \
  --include='table_*.tsv' --include='README.md' \
  --exclude='*' \
  rajd2@login.hpc.cam.ac.uk:/rds/user/rajd2/hpc-work/abcd_development/genetic_analysis/work/results_70tab_hcp/scan_geom/ \
  ~/Git/abcd_development/genetic_analysis/work/results_70tab_hcp/scan_geom/
```
