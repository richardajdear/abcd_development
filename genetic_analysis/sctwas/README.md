# Cell-type-specific TWAS from SingleBrain sc-eQTLs (SuSiE-RSS weights)

Task spec: `transcriptional_maturation/gwas_cells/csd3/SCTWAS_INSTRUCTIONS.md` (commit 2f1aceb).
Scripts: `abcd_development/genetic_analysis/sctwas/` (committed). Work dir: this directory (`$W`).
Run on CSD3 starting 2026-10-24.

## Status (2026-10-25, all part A runs complete)

- **Complete (weights + assoc + report):** all 15 tags (Ext1, Ext1_n983, Ext, OPC, IN, OD, Ast, MG, Ext2–Ext8).
  361 Slurm jobs completed; the 5 cancellations were my own resubmissions. There were 0 susieR errors in any cell type.
- **`w_flag` excluded:** ABCD thinning and baseline have **0** genes at p < 2.5e-6 in every cell type and both
  models. SCZ25 keeps 17–116 such genes per cell type in susie, and λ stays about 1.9–2.6.
- **Part B (PsychAD) not started.** It needs a Synapse account and token from the user.
- **Later bug fix (`03_assoc.R`):** top1 files were first written with the SuSiE `w_flag` value. The flag is now set
  to FALSE before writing, and the existing 75 top1 files were patched to match. Weights and z were not affected.

## Deviation to note: implausible SuSiE weights (`w_flag`)

**What happened.** In the first Ext1 pass, a small number of genes had SuSiE-RSS posterior-mean weights far above the
marginal effect. On the standardised scale, |w| reached about 30, where |z_top|/√n ≈ 1.2. The large weights came as
opposite-signed effects on SNPs in near-perfect 1000G LD, the known failure mode of out-of-sample LD. Their TWAS z
values were absurd: VPS13A −134 and ODR4 −78 for ABCD thinning; MAPT-AS1 +29 and VPS13A −19 for SCZ.

**Scale.** This is not "susieR fails on many genes":

- Ext1: 17 of 4,322 modelled genes;
- Ext: 102 of 7,115;
- s_rss is unremarkable for these genes (0.1–0.2), so the s_rss > 0.3 flag does not catch them.

**What I did.** The method was not changed and no gene was dropped. I added a flag:

- w_ratio = max|w| / (|z_top|/√n), and `w_flag` = w_ratio > 3.
- The flag is in `weights/<CT>.genes.tsv.gz`, in `weights/<CT>.susie.tsv.gz` and as an extra last column of every
  `assoc/*.susie.tsv.gz`. It is always FALSE in top1.
- `report/trait_summary.tsv` gives λ and hit counts with and without flagged genes.
- **Recommendation:** exclude `w_flag` genes from the primary susie model downstream, or use top1 for them.

**Context for λ.** The old bulk GTEx Brain_Cortex TWAS (`09_twas.csv`) had λ: SCZ_EUR (PGC3) 2.02, MDD_EUR 1.88,
ABCD_thin_EUR 1.02, ABCD_base_EUR 1.11. The single-cell values below sit in the same range.

## What was done

**Stage 0. SNP filter (`00_snp_filter.py`, run once).** The 1000G EUR panel (`magma/reference_data/g1000_eur`,
503 individuals, hg19, rsIDs) was filtered to autosomal, unique, single-base rsIDs. Then:

- the MHC (hg19 chr6:25–34 Mb) was removed;
- SNPs with panel MAF < 0.01 were removed;
- strand-ambiguous A/T and C/G SNPs were removed;
- only SNPs present in **all four** required GWAS were kept. Their alleles had to match the panel after a swap or a
  complement, with finite BETA and SE > 0.

The panel was then subset with `--keep-allele-order` to `ld/panel_keep`. Panel A1 (bim column 5) is the effect allele
of every weight and every aligned z. GWAS z = BETA/SE, aligned to panel A1, raw sign, is in `gwas/<trait>.z.tsv.gz`.

**Stage 1a. eQTL inputs (`01_prep_eqtl.py`, per cell type).**

- eGenes are `top_assoc` rows with `qval < 0.05`. Features tested are all `top_assoc` rows, which matches
  "Test genes" in Supp. Table 6 (e.g. Ext1 11,893).
- `full_assoc` was streamed, keeping eGene × kept-SNP rows.
- The signed statistic is `fixed_z` relative to `Allele`. Rows with `fixed_sd == 0` were dropped.
- z was aligned to panel A1. `{ref, alt}` had to match the panel alleles up to strand, and the chromosome had to agree.

**Stage 1b. SuSiE-RSS (`02_susie_fit.R`, per cell type × chromosome).**

- **LD:** R is the 1000G EUR correlation of the gene's cis SNPs, in the same order as z. It is computed in R from
  `ld/panel_keep.bed` (`plink_bed.R`, mean-imputed), and matches `plink --r square` to 4 decimals.
- **Fit:** `susie_rss(z = z, R = R, n = n, L = 5, estimate_residual_variance = FALSE)`, other arguments at their
  defaults, with n = 757. For the Ext1 sensitivity run, n = 983 (tag `Ext1_n983`).
- **Weights:** `coef(fit)[-1]`.
- **Diagnostics:** `estimate_s_rss(z, R, n)`, `fit$converged`, and the purity of `fit$sets`.
- **Primary model (`susie`):** genes with at least one 95% credible set, keeping only the credible-set SNPs and their
  posterior-mean weights.
- **Secondary model (`top1`):** every eGene, with w = z_top / sqrt(n) at the lead |z| SNP among the filtered cis SNPs.
- **s_rss flag:** genes with s_rss > 0.3 are flagged, not dropped. The `s_rss` column is in every weights and assoc
  table.

**Stage 1c and 2 (`02b_collect.R`, `03_assoc.R`).**

- **Weights tables:** `weights/<CT>.{susie,top1}.tsv.gz`, with columns `gene gene_version rsid a1 a2 weight cs_id pip
  n_snps_locus n_cs s_rss converged`.
- **Diagnostics:** per-eGene diagnostics are in `weights/<CT>.genes.tsv.gz`.
- **LD:** the SNP lists in the weights tables, together with `ld/panel_keep`, rebuild every LD matrix. No per-gene LD
  matrices are stored.
- **Association:** z_g = wᵀz / sqrt(wᵀRw), with R from the panel over the model SNPs, and p = 2Φ(−|z_g|).
- **Output:** `assoc/<trait>.<CT>.<model>.tsv.gz` with columns `gene n_model_snps z p n_cs s_rss`.

Sign conventions are untouched. `global_slope_1lmm`: higher = slower thinning (flip locally).

## Software

| software | version / path |
|---|---|
| R | 4.5.3, `legacy/hpc/work/envs/abcdR` (env not modified) |
| susieR | **0.14.2** (CRAN), with mixsqp 0.3-54, matrixStats 1.5.0, RcppArmadillo 15.6.2-1, irlba 2.4.1, reshape 0.8.10, plyr 1.8.9, installed into `$W/Rlib` |
| data.table / Matrix | 1.18.4 / 1.7.5 (env) |
| plink | v1.9.0-b.7.7 64-bit (22 Oct 2024), `legacy/hpc/work/bin/plink` |
| Python | 3.11.15, `legacy/hpc/work/envs/abcd` (pandas 3.0.5, numpy 2.4.6, scipy 1.17.1, pyarrow 25.0.0) |

**Install note.** The relocated abcdR env still carries its pre-move prefix in `lib/R/etc/Makeconf`, and the system
gcc 8.5 cannot compile against R 4.5 headers, which need C23. Packages were therefore built with the env's conda
gcc 16.1 and a private `R_MAKEVARS_USER` file (`tmp/Makevars`) that rewrites the stale prefix. Neither the env nor
its Makeconf was changed.

## Commands

```bash
cd /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas
./download.sh                                   # Zenodo 14908182, md5-verified, raw/ -> rds-Pb9UGUlrwWc (hpc-work had 43 GB free)
sbatch ../../sctwas/00_snp_filter.sbatch        # stage 0
bash ../../sctwas/run_sctwas.sh                 # per CT: 01_prep_eqtl -> 02_susie_fit (array chr1-22) -> 03_collect_assoc; Ext1 also at n = 983
sbatch ../../sctwas/04_report.sbatch            # report/*.tsv + tables below (too slow for the login node)
```

## Effect-allele checks (GWAS)

| trait | file | effect allele | checked against |
|---|---|---|---|
| ABCD_thin_EUR | `results_70tab_hcp/scan_1lmm/assoc_eur/global_slope_1lmm.sumstats.tsv.gz` | `A1` | `R/05_collect_assoc.R`: A1 = the GENESIS counted allele, and `AF1`/`BETA` refer to it |
| ABCD_base_EUR | `.../baseline_thickness_1lmm.sumstats.tsv.gz` | `A1` | same |
| SCZ25_EUR | `scores_scz2025/gwas/SCZ25_EUR.ma` | `A1` | `step9_scz2025_normalise.sbatch`: A1 = ALT = the allele ES refers to |
| MDD_EUR | `legacy/hpc_v2/.../MDD_eur.ma` | `A1` | raw PGC MDD 2025 file: `effect_allele` equals `.ma` A1 with identical beta (rs3094315, rs3131972) |
| PGC3_EUR | `legacy/hpc_v2/.../SCZ_eur.ma` | `A1` | raw PGC3 EUR file: A1/BETA identical (rs62513865) |

## Subtype → layer mapping

Source: SingleBrain **Supplementary Table 6** (Nature Genetics 2026, s41588-026-02541-x, MOESM4 xlsx, copied to
`ref/singlebrain_supp_tables.xlsx`). The pipeline code (Zenodo 18332620) contains no subtype labels.

| SingleBrain | label (Supp. Table 6) | samples | test genes | eGenes (paper) |
|---|---|---|---|---|
| Ext1 | L2/3 IT | 935 | 11,893 | 11,118 |
| Ext2 | L4 IT (main text: "L4 IT3") | 920 | 12,048 | 10,932 |
| Ext3 | L5/6 NP | 492 | 12,567 | 7,500 |
| Ext4 | L5 IT | 928 | 11,988 | 10,886 |
| Ext5 | L6b | 545 | 12,099 | 8,287 |
| Ext6 | L6 CT | 650 | 12,309 | 6,863 |
| Ext7 | L6 IT Car3 | 780 | 11,887 | 8,855 |
| Ext8 | L6 IT | 920 | 11,867 | 10,408 |

## Positive controls (SCZ)

- **Source of the gene list:** Supp. Table 11, Trubetskoy et al. 2022 SCZ, PP4 ≥ 0.8, copied to
  `ref/singlebrain_ST11_scz_coloc.tsv`.
- **Extended Data Fig. 2:** a rendering of the major-class rows of that table, 131 genes at 81 loci.
- **Fig. 3c:** the 38 cell-type-specific eGenes, read off the figure. Each gene's cell type was checked against the
  table rows. RAB21 in OD is PP4 = 0.80.
- **Where the results are:** z and p for every listed gene in every cell type, both SCZ GWAS and both models, are in
  `report/positive_controls.tsv`.
- **Paper's effect directions** (protective allele lowers expression, so risk raises expression):
  - CNTN4 in Ext;
  - SLC4A10 in Ext and Ext2 (DPP4 locus), and of the subtypes only in Ext2 (PP4 0.87).

<!-- AUTO:RESULTS -->

## Results (generated by 04_report.py)

### Per cell type

| tag | label | features_tested | egenes | egenes_with_snps | fit_errors | non_converged | genes_ge1_cs | w_flag | susie_model_genes | top1_model_genes | median_snps_susie | median_cis_snps | cs_1 | cs_2 | cs_3plus | s_rss_median | s_rss_q90 | share_s_gt_0_3 | share_s_gt_0_3_in_susie_model |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Ext1 | L2/3 IT | 11893 | 11118 | 11116 | 0 | 16 | 4322 | 17 | 4322 | 11116 | 13 | 3.11e+03 | 3159 | 846 | 317 | 0.0711 | 0.148 | 0.00144 | 0.00185 |
| Ext1_n983 | L2/3 IT (n = 983 sensitivity) | 11893 | 11118 | 11116 | 0 | 22 | 4346 | 23 | 4346 | 11116 | 13 | 3.11e+03 | 3136 | 870 | 340 | 0.0714 | 0.149 | 0.00144 | 0.00184 |
| Ext | excitatory neurons | 12000 | 11629 | 11626 | 0 | 88 | 7115 | 102 | 7115 | 11626 | 13 | 3.1e+03 | 4170 | 1890 | 1055 | 0.0742 | 0.16 | 0.00353 | 0.0052 |
| OPC | OPC | 12029 | 9335 | 9332 | 0 | 5 | 2811 | 6 | 2811 | 9332 | 12 | 3.14e+03 | 2113 | 513 | 185 | 0.0685 | 0.144 | 0.00182 | 0.00462 |
| IN | inhibitory neurons | 12226 | 11165 | 11162 | 0 | 46 | 4938 | 51 | 4938 | 11162 | 13 | 3.12e+03 | 3468 | 1044 | 426 | 0.0722 | 0.151 | 0.00188 | 0.00324 |
| OD | oligodendrocytes | 11505 | 9898 | 9896 | 0 | 18 | 4157 | 22 | 4157 | 9896 | 13 | 3.08e+03 | 2793 | 888 | 476 | 0.068 | 0.147 | 0.00303 | 0.0065 |
| Ast | astrocytes | 12257 | 10282 | 10280 | 0 | 14 | 3684 | 17 | 3684 | 10280 | 13 | 3.11e+03 | 2712 | 703 | 269 | 0.0694 | 0.145 | 0.00165 | 0.0038 |
| MG | microglia | 12320 | 8423 | 8418 | 0 | 2 | 2194 | 3 | 2194 | 8418 | 12 | 3.12e+03 | 1740 | 332 | 122 | 0.0678 | 0.141 | 0.00131 | 0.00319 |
| Ext2 | L4 IT | 12048 | 10932 | 10930 | 0 | 24 | 3758 | 25 | 3758 | 10930 | 14 | 3.1e+03 | 2842 | 652 | 264 | 0.0705 | 0.146 | 0.00146 | 0.00319 |
| Ext6 | L6 CT | 12309 | 6863 | 6862 | 0 | 5 | 1103 | 8 | 1103 | 6862 | 14 | 3.08e+03 | 983 | 93 | 27 | 0.0677 | 0.138 | 0.00131 | 0.00363 |
| Ext3 | L5/6 NP | 12567 | 7500 | 7498 | 0 | 12 | 936 | 15 | 936 | 7498 | 13 | 3.06e+03 | 824 | 86 | 26 | 0.072 | 0.143 | 0.00107 | 0.00427 |
| Ext4 | L5 IT | 11988 | 10886 | 10884 | 0 | 21 | 3775 | 25 | 3775 | 10884 | 14 | 3.11e+03 | 2859 | 662 | 254 | 0.0697 | 0.147 | 0.00147 | 0.00265 |
| Ext5 | L6b | 12099 | 8287 | 8285 | 0 | 8 | 1204 | 10 | 1204 | 8285 | 13 | 3.09e+03 | 1051 | 125 | 28 | 0.072 | 0.146 | 0.00133 | 0.00249 |
| Ext7 | L6 IT Car3 | 11887 | 8855 | 8853 | 0 | 13 | 2019 | 16 | 2019 | 8853 | 14 | 3.08e+03 | 1670 | 255 | 94 | 0.067 | 0.14 | 0.00136 | 0.00297 |
| Ext8 | L6 IT | 11867 | 10408 | 10406 | 0 | 18 | 3102 | 19 | 3102 | 10406 | 13 | 3.11e+03 | 2420 | 519 | 163 | 0.0691 | 0.146 | 0.00192 | 0.00387 |

### Per trait x cell type x model

| trait | tag | model | genes | lambda_gc | n_p_lt_2_5e6 | n_fdr05 | n_w_flag | lambda_gc_unflagged | n_p_lt_2_5e6_unflagged |
|---|---|---|---|---|---|---|---|---|---|
| SCZ25_EUR | Ext1 | susie | 4322 | 2.57 | 82 | 462 | 17 | 2.56 | 75 |
| SCZ25_EUR | Ext1 | top1 | 11116 | 2 | 121 | 735 | 0 | 2 | 121 |
| SCZ25_EUR | Ext1_n983 | susie | 4346 | 2.61 | 82 | 466 | 23 | 2.58 | 73 |
| SCZ25_EUR | Ext1_n983 | top1 | 11116 | 2 | 121 | 735 | 0 | 2 | 121 |
| SCZ25_EUR | Ext | susie | 7115 | 2.49 | 139 | 754 | 102 | 2.47 | 116 |
| SCZ25_EUR | Ext | top1 | 11626 | 2.18 | 164 | 887 | 0 | 2.18 | 164 |
| SCZ25_EUR | OPC | susie | 2811 | 2.25 | 61 | 261 | 6 | 2.25 | 59 |
| SCZ25_EUR | OPC | top1 | 9332 | 1.88 | 104 | 524 | 0 | 1.88 | 104 |
| SCZ25_EUR | IN | susie | 4938 | 2.54 | 100 | 488 | 51 | 2.51 | 91 |
| SCZ25_EUR | IN | top1 | 11162 | 2.07 | 136 | 726 | 0 | 2.07 | 136 |
| SCZ25_EUR | OD | susie | 4157 | 2.33 | 80 | 408 | 22 | 2.33 | 75 |
| SCZ25_EUR | OD | top1 | 9896 | 2.02 | 131 | 641 | 0 | 2.02 | 131 |
| SCZ25_EUR | Ast | susie | 3684 | 2.3 | 73 | 391 | 17 | 2.29 | 70 |
| SCZ25_EUR | Ast | top1 | 10280 | 1.91 | 115 | 615 | 0 | 1.91 | 115 |
| SCZ25_EUR | MG | susie | 2194 | 2.12 | 24 | 184 | 3 | 2.12 | 23 |
| SCZ25_EUR | MG | top1 | 8418 | 1.76 | 63 | 387 | 0 | 1.76 | 63 |
| SCZ25_EUR | Ext2 | susie | 3758 | 2.58 | 74 | 403 | 25 | 2.57 | 72 |
| SCZ25_EUR | Ext2 | top1 | 10930 | 1.98 | 123 | 665 | 0 | 1.98 | 123 |
| SCZ25_EUR | Ext6 | susie | 1103 | 2.27 | 27 | 117 | 8 | 2.24 | 27 |
| SCZ25_EUR | Ext6 | top1 | 6862 | 1.71 | 72 | 325 | 0 | 1.71 | 72 |
| SCZ25_EUR | Ext3 | susie | 936 | 1.97 | 19 | 75 | 15 | 1.94 | 17 |
| SCZ25_EUR | Ext3 | top1 | 7498 | 1.83 | 50 | 329 | 0 | 1.83 | 50 |
| SCZ25_EUR | Ext4 | susie | 3775 | 2.64 | 66 | 412 | 25 | 2.63 | 60 |
| SCZ25_EUR | Ext4 | top1 | 10884 | 1.96 | 116 | 684 | 0 | 1.96 | 116 |
| SCZ25_EUR | Ext5 | susie | 1204 | 2.51 | 33 | 130 | 10 | 2.51 | 29 |
| SCZ25_EUR | Ext5 | top1 | 8285 | 1.8 | 73 | 370 | 0 | 1.8 | 73 |
| SCZ25_EUR | Ext7 | susie | 2019 | 2.44 | 44 | 203 | 16 | 2.39 | 39 |
| SCZ25_EUR | Ext7 | top1 | 8853 | 1.92 | 90 | 421 | 0 | 1.92 | 90 |
| SCZ25_EUR | Ext8 | susie | 3102 | 2.4 | 66 | 349 | 19 | 2.39 | 61 |
| SCZ25_EUR | Ext8 | top1 | 10406 | 1.92 | 111 | 595 | 0 | 1.92 | 111 |
| PGC3_EUR | Ext1 | susie | 4322 | 2.54 | 59 | 395 | 17 | 2.54 | 54 |
| PGC3_EUR | Ext1 | top1 | 11014 | 1.94 | 87 | 600 | 0 | 1.94 | 87 |
| PGC3_EUR | Ext1_n983 | susie | 4346 | 2.54 | 62 | 404 | 23 | 2.54 | 54 |
| PGC3_EUR | Ext1_n983 | top1 | 11014 | 1.94 | 87 | 600 | 0 | 1.94 | 87 |
| PGC3_EUR | Ext | susie | 7115 | 2.42 | 112 | 631 | 102 | 2.4 | 91 |
| PGC3_EUR | Ext | top1 | 11514 | 2.1 | 118 | 698 | 0 | 2.1 | 118 |
| PGC3_EUR | OPC | susie | 2811 | 2.2 | 45 | 242 | 6 | 2.2 | 43 |
| PGC3_EUR | OPC | top1 | 9238 | 1.84 | 84 | 441 | 0 | 1.84 | 84 |
| PGC3_EUR | IN | susie | 4936 | 2.43 | 75 | 422 | 51 | 2.39 | 69 |
| PGC3_EUR | IN | top1 | 11073 | 1.97 | 107 | 626 | 0 | 1.97 | 107 |
| PGC3_EUR | OD | susie | 4155 | 2.29 | 66 | 347 | 22 | 2.29 | 62 |
| PGC3_EUR | OD | top1 | 9804 | 1.91 | 103 | 553 | 0 | 1.91 | 103 |
| PGC3_EUR | Ast | susie | 3684 | 2.31 | 59 | 331 | 17 | 2.3 | 55 |
| PGC3_EUR | Ast | top1 | 10178 | 1.85 | 85 | 521 | 0 | 1.85 | 85 |
| PGC3_EUR | MG | susie | 2194 | 2.08 | 20 | 137 | 3 | 2.08 | 20 |
| PGC3_EUR | MG | top1 | 8350 | 1.7 | 52 | 293 | 0 | 1.7 | 52 |
| PGC3_EUR | Ext2 | susie | 3755 | 2.44 | 52 | 370 | 25 | 2.43 | 50 |
| PGC3_EUR | Ext2 | top1 | 10832 | 1.89 | 88 | 546 | 0 | 1.89 | 88 |
| PGC3_EUR | Ext6 | susie | 1103 | 2.12 | 22 | 111 | 8 | 2.12 | 22 |
| PGC3_EUR | Ext6 | top1 | 6782 | 1.67 | 57 | 286 | 0 | 1.67 | 57 |
| PGC3_EUR | Ext3 | susie | 936 | 1.74 | 16 | 67 | 15 | 1.74 | 14 |
| PGC3_EUR | Ext3 | top1 | 7418 | 1.66 | 41 | 262 | 0 | 1.66 | 41 |
| PGC3_EUR | Ext4 | susie | 3773 | 2.57 | 48 | 325 | 25 | 2.54 | 41 |
| PGC3_EUR | Ext4 | top1 | 10794 | 1.92 | 84 | 562 | 0 | 1.92 | 84 |
| PGC3_EUR | Ext5 | susie | 1204 | 2.49 | 21 | 127 | 10 | 2.48 | 19 |
| PGC3_EUR | Ext5 | top1 | 8212 | 1.71 | 52 | 278 | 0 | 1.71 | 52 |
| PGC3_EUR | Ext7 | susie | 2019 | 2.41 | 40 | 177 | 16 | 2.38 | 34 |
| PGC3_EUR | Ext7 | top1 | 8775 | 1.82 | 63 | 353 | 0 | 1.82 | 63 |
| PGC3_EUR | Ext8 | susie | 3101 | 2.43 | 50 | 301 | 19 | 2.42 | 45 |
| PGC3_EUR | Ext8 | top1 | 10322 | 1.87 | 75 | 501 | 0 | 1.87 | 75 |
| MDD_EUR | Ext1 | susie | 4322 | 2.27 | 49 | 335 | 17 | 2.25 | 41 |
| MDD_EUR | Ext1 | top1 | 11116 | 1.93 | 80 | 517 | 0 | 1.93 | 80 |
| MDD_EUR | Ext1_n983 | susie | 4346 | 2.25 | 51 | 329 | 23 | 2.23 | 40 |
| MDD_EUR | Ext1_n983 | top1 | 11116 | 1.93 | 80 | 517 | 0 | 1.93 | 80 |
| MDD_EUR | Ext | susie | 7115 | 2.3 | 90 | 604 | 102 | 2.25 | 60 |
| MDD_EUR | Ext | top1 | 11626 | 2.06 | 90 | 672 | 0 | 2.06 | 90 |
| MDD_EUR | OPC | susie | 2811 | 2.38 | 22 | 199 | 6 | 2.38 | 20 |
| MDD_EUR | OPC | top1 | 9332 | 1.83 | 48 | 376 | 0 | 1.83 | 48 |
| MDD_EUR | IN | susie | 4938 | 2.23 | 60 | 367 | 51 | 2.22 | 43 |
| MDD_EUR | IN | top1 | 11162 | 1.94 | 72 | 564 | 0 | 1.94 | 72 |
| MDD_EUR | OD | susie | 4157 | 2.2 | 31 | 273 | 22 | 2.2 | 25 |
| MDD_EUR | OD | top1 | 9896 | 1.87 | 49 | 421 | 0 | 1.87 | 49 |
| MDD_EUR | Ast | susie | 3684 | 2.41 | 29 | 287 | 17 | 2.39 | 26 |
| MDD_EUR | Ast | top1 | 10280 | 1.93 | 57 | 450 | 0 | 1.93 | 57 |
| MDD_EUR | MG | susie | 2194 | 2.1 | 20 | 130 | 3 | 2.1 | 19 |
| MDD_EUR | MG | top1 | 8418 | 1.83 | 41 | 293 | 0 | 1.83 | 41 |
| MDD_EUR | Ext2 | susie | 3758 | 2.34 | 36 | 312 | 25 | 2.33 | 29 |
| MDD_EUR | Ext2 | top1 | 10930 | 1.9 | 61 | 499 | 0 | 1.9 | 61 |
| MDD_EUR | Ext6 | susie | 1103 | 2.27 | 8 | 72 | 8 | 2.3 | 7 |
| MDD_EUR | Ext6 | top1 | 6862 | 1.8 | 24 | 235 | 0 | 1.8 | 24 |
| MDD_EUR | Ext3 | susie | 936 | 2.21 | 6 | 59 | 15 | 2.2 | 2 |
| MDD_EUR | Ext3 | top1 | 7498 | 1.76 | 24 | 204 | 0 | 1.76 | 24 |
| MDD_EUR | Ext4 | susie | 3775 | 2.26 | 40 | 295 | 25 | 2.23 | 35 |
| MDD_EUR | Ext4 | top1 | 10884 | 1.95 | 63 | 497 | 0 | 1.95 | 63 |
| MDD_EUR | Ext5 | susie | 1204 | 2.21 | 11 | 76 | 10 | 2.21 | 10 |
| MDD_EUR | Ext5 | top1 | 8285 | 1.71 | 38 | 235 | 0 | 1.71 | 38 |
| MDD_EUR | Ext7 | susie | 2019 | 2.26 | 23 | 166 | 16 | 2.22 | 15 |
| MDD_EUR | Ext7 | top1 | 8853 | 1.82 | 40 | 314 | 0 | 1.82 | 40 |
| MDD_EUR | Ext8 | susie | 3102 | 2.31 | 35 | 275 | 19 | 2.28 | 32 |
| MDD_EUR | Ext8 | top1 | 10406 | 1.91 | 56 | 484 | 0 | 1.91 | 56 |
| ABCD_thin_EUR | Ext1 | susie | 4322 | 0.987 | 4 | 6 | 17 | 0.975 | 0 |
| ABCD_thin_EUR | Ext1 | top1 | 11116 | 0.962 | 0 | 0 | 0 | 0.962 | 0 |
| ABCD_thin_EUR | Ext1_n983 | susie | 4346 | 0.987 | 4 | 7 | 23 | 0.976 | 0 |
| ABCD_thin_EUR | Ext1_n983 | top1 | 11116 | 0.962 | 0 | 0 | 0 | 0.962 | 0 |
| ABCD_thin_EUR | Ext | susie | 7115 | 1 | 21 | 25 | 102 | 0.994 | 0 |
| ABCD_thin_EUR | Ext | top1 | 11626 | 0.985 | 0 | 0 | 0 | 0.985 | 0 |
| ABCD_thin_EUR | OPC | susie | 2811 | 0.993 | 4 | 4 | 6 | 0.989 | 0 |
| ABCD_thin_EUR | OPC | top1 | 9332 | 0.995 | 0 | 0 | 0 | 0.995 | 0 |
| ABCD_thin_EUR | IN | susie | 4938 | 1.08 | 13 | 13 | 51 | 1.07 | 0 |
| ABCD_thin_EUR | IN | top1 | 11162 | 1.03 | 0 | 0 | 0 | 1.03 | 0 |
| ABCD_thin_EUR | OD | susie | 4157 | 1.05 | 2 | 2 | 22 | 1.05 | 0 |
| ABCD_thin_EUR | OD | top1 | 9896 | 1.02 | 0 | 0 | 0 | 1.02 | 0 |
| ABCD_thin_EUR | Ast | susie | 3684 | 1.02 | 1 | 1 | 17 | 1.02 | 0 |
| ABCD_thin_EUR | Ast | top1 | 10280 | 1.01 | 0 | 0 | 0 | 1.01 | 0 |
| ABCD_thin_EUR | MG | susie | 2194 | 1.06 | 0 | 0 | 3 | 1.06 | 0 |
| ABCD_thin_EUR | MG | top1 | 8418 | 1.01 | 0 | 0 | 0 | 1.01 | 0 |
| ABCD_thin_EUR | Ext2 | susie | 3758 | 1.04 | 6 | 6 | 25 | 1.02 | 0 |
| ABCD_thin_EUR | Ext2 | top1 | 10930 | 1.01 | 0 | 0 | 0 | 1.01 | 0 |
| ABCD_thin_EUR | Ext6 | susie | 1103 | 1.05 | 1 | 1 | 8 | 1.05 | 0 |
| ABCD_thin_EUR | Ext6 | top1 | 6862 | 1.05 | 0 | 0 | 0 | 1.05 | 0 |
| ABCD_thin_EUR | Ext3 | susie | 936 | 1.04 | 5 | 5 | 15 | 1.03 | 0 |
| ABCD_thin_EUR | Ext3 | top1 | 7498 | 0.963 | 0 | 0 | 0 | 0.963 | 0 |
| ABCD_thin_EUR | Ext4 | susie | 3775 | 1.03 | 7 | 7 | 25 | 1.02 | 0 |
| ABCD_thin_EUR | Ext4 | top1 | 10884 | 1 | 0 | 0 | 0 | 1 | 0 |
| ABCD_thin_EUR | Ext5 | susie | 1204 | 1.09 | 1 | 1 | 10 | 1.08 | 0 |
| ABCD_thin_EUR | Ext5 | top1 | 8285 | 1 | 0 | 0 | 0 | 1 | 0 |
| ABCD_thin_EUR | Ext7 | susie | 2019 | 1.08 | 4 | 4 | 16 | 1.06 | 0 |
| ABCD_thin_EUR | Ext7 | top1 | 8853 | 0.947 | 0 | 0 | 0 | 0.947 | 0 |
| ABCD_thin_EUR | Ext8 | susie | 3102 | 1.03 | 4 | 5 | 19 | 1.02 | 0 |
| ABCD_thin_EUR | Ext8 | top1 | 10406 | 0.992 | 0 | 0 | 0 | 0.992 | 0 |
| ABCD_base_EUR | Ext1 | susie | 4322 | 1.17 | 2 | 8 | 17 | 1.16 | 0 |
| ABCD_base_EUR | Ext1 | top1 | 11116 | 1.09 | 0 | 0 | 0 | 1.09 | 0 |
| ABCD_base_EUR | Ext1_n983 | susie | 4346 | 1.17 | 2 | 7 | 23 | 1.16 | 0 |
| ABCD_base_EUR | Ext1_n983 | top1 | 11116 | 1.09 | 0 | 0 | 0 | 1.09 | 0 |
| ABCD_base_EUR | Ext | susie | 7115 | 1.15 | 24 | 28 | 102 | 1.14 | 0 |
| ABCD_base_EUR | Ext | top1 | 11626 | 1.09 | 0 | 0 | 0 | 1.09 | 0 |
| ABCD_base_EUR | OPC | susie | 2811 | 1.07 | 2 | 2 | 6 | 1.07 | 0 |
| ABCD_base_EUR | OPC | top1 | 9332 | 1.07 | 0 | 0 | 0 | 1.07 | 0 |
| ABCD_base_EUR | IN | susie | 4938 | 1.1 | 14 | 20 | 51 | 1.08 | 0 |
| ABCD_base_EUR | IN | top1 | 11162 | 1.06 | 0 | 0 | 0 | 1.06 | 0 |
| ABCD_base_EUR | OD | susie | 4157 | 1.1 | 5 | 8 | 22 | 1.09 | 0 |
| ABCD_base_EUR | OD | top1 | 9896 | 1.07 | 0 | 0 | 0 | 1.07 | 0 |
| ABCD_base_EUR | Ast | susie | 3684 | 1.07 | 5 | 6 | 17 | 1.06 | 0 |
| ABCD_base_EUR | Ast | top1 | 10280 | 1.04 | 0 | 0 | 0 | 1.04 | 0 |
| ABCD_base_EUR | MG | susie | 2194 | 1.07 | 0 | 1 | 3 | 1.07 | 0 |
| ABCD_base_EUR | MG | top1 | 8418 | 1.1 | 0 | 0 | 0 | 1.1 | 0 |
| ABCD_base_EUR | Ext2 | susie | 3758 | 1.08 | 5 | 7 | 25 | 1.07 | 0 |
| ABCD_base_EUR | Ext2 | top1 | 10930 | 1.03 | 0 | 0 | 0 | 1.03 | 0 |
| ABCD_base_EUR | Ext6 | susie | 1103 | 0.977 | 1 | 2 | 8 | 0.968 | 0 |
| ABCD_base_EUR | Ext6 | top1 | 6862 | 1.04 | 0 | 0 | 0 | 1.04 | 0 |
| ABCD_base_EUR | Ext3 | susie | 936 | 1.09 | 2 | 5 | 15 | 1.07 | 0 |
| ABCD_base_EUR | Ext3 | top1 | 7498 | 1.08 | 0 | 1 | 0 | 1.08 | 0 |
| ABCD_base_EUR | Ext4 | susie | 3775 | 1.12 | 3 | 7 | 25 | 1.11 | 0 |
| ABCD_base_EUR | Ext4 | top1 | 10884 | 1.1 | 0 | 1 | 0 | 1.1 | 0 |
| ABCD_base_EUR | Ext5 | susie | 1204 | 1.11 | 1 | 1 | 10 | 1.1 | 0 |
| ABCD_base_EUR | Ext5 | top1 | 8285 | 1.05 | 0 | 1 | 0 | 1.05 | 0 |
| ABCD_base_EUR | Ext7 | susie | 2019 | 1.06 | 5 | 7 | 16 | 1.04 | 0 |
| ABCD_base_EUR | Ext7 | top1 | 8853 | 1.08 | 0 | 0 | 0 | 1.08 | 0 |
| ABCD_base_EUR | Ext8 | susie | 3102 | 1.18 | 6 | 8 | 19 | 1.18 | 0 |
| ABCD_base_EUR | Ext8 | top1 | 10406 | 1.05 | 0 | 0 | 0 | 1.05 | 0 |

### Comparison with the bulk GTEx Brain_Cortex TWAS (SCZ)

| tag | trait | model | shared_genes | spearman | p | sign_agree_where_bulk_abs_z_gt_3_3 |
|---|---|---|---|---|---|---|
| Ext1 | PGC3_EUR | susie | 3162 | 0.483 | 1.17e-184 | 0.846 |
| Ext1 | PGC3_EUR | top1 | 7213 | 0.278 | 1.45e-128 | 0.676 |
| Ext1 | SCZ25_EUR | susie | 3162 | 0.449 | 1.36e-156 | 0.84 |
| Ext1 | SCZ25_EUR | top1 | 7283 | 0.26 | 1.07e-112 | 0.678 |
| Ext1_n983 | PGC3_EUR | susie | 3179 | 0.484 | 1.56e-186 | 0.842 |
| Ext1_n983 | PGC3_EUR | top1 | 7213 | 0.278 | 1.45e-128 | 0.676 |
| Ext1_n983 | SCZ25_EUR | susie | 3179 | 0.449 | 2.80e-157 | 0.836 |
| Ext1_n983 | SCZ25_EUR | top1 | 7283 | 0.26 | 1.07e-112 | 0.678 |
| Ext | PGC3_EUR | susie | 4985 | 0.427 | 9.88e-220 | 0.797 |
| Ext | PGC3_EUR | top1 | 7505 | 0.316 | 6.42e-174 | 0.694 |
| Ext | SCZ25_EUR | susie | 4985 | 0.387 | 1.56e-177 | 0.797 |
| Ext | SCZ25_EUR | top1 | 7580 | 0.29 | 9.00e-147 | 0.689 |
| OPC | PGC3_EUR | susie | 1911 | 0.386 | 4.99e-69 | 0.8 |
| OPC | PGC3_EUR | top1 | 6021 | 0.177 | 2.03e-43 | 0.628 |
| OPC | SCZ25_EUR | susie | 1911 | 0.357 | 1.58e-58 | 0.811 |
| OPC | SCZ25_EUR | top1 | 6085 | 0.163 | 1.45e-37 | 0.616 |
| IN | PGC3_EUR | susie | 3508 | 0.399 | 2.69e-134 | 0.801 |
| IN | PGC3_EUR | top1 | 7233 | 0.257 | 2.41e-109 | 0.657 |
| IN | SCZ25_EUR | susie | 3510 | 0.365 | 4.24e-111 | 0.801 |
| IN | SCZ25_EUR | top1 | 7287 | 0.227 | 3.80e-86 | 0.65 |
| OD | PGC3_EUR | susie | 2916 | 0.385 | 1.27e-103 | 0.766 |
| OD | PGC3_EUR | top1 | 6429 | 0.232 | 1.75e-79 | 0.641 |
| OD | SCZ25_EUR | susie | 2918 | 0.354 | 5.02e-87 | 0.76 |
| OD | SCZ25_EUR | top1 | 6484 | 0.213 | 1.84e-67 | 0.651 |
| Ast | PGC3_EUR | susie | 2573 | 0.443 | 5.99e-124 | 0.802 |
| Ast | PGC3_EUR | top1 | 6661 | 0.244 | 8.52e-91 | 0.675 |
| Ast | SCZ25_EUR | susie | 2573 | 0.411 | 1.76e-105 | 0.802 |
| Ast | SCZ25_EUR | top1 | 6719 | 0.231 | 1.97e-82 | 0.643 |
| MG | PGC3_EUR | susie | 1498 | 0.287 | 6.79e-30 | 0.719 |
| MG | PGC3_EUR | top1 | 5464 | 0.118 | 1.63e-18 | 0.558 |
| MG | SCZ25_EUR | susie | 1498 | 0.263 | 3.41e-25 | 0.703 |
| MG | SCZ25_EUR | top1 | 5507 | 0.114 | 1.89e-17 | 0.537 |
| Ext2 | PGC3_EUR | susie | 2765 | 0.482 | 7.75e-161 | 0.873 |
| Ext2 | PGC3_EUR | top1 | 7072 | 0.277 | 3.44e-125 | 0.666 |
| Ext2 | SCZ25_EUR | susie | 2767 | 0.442 | 4.96e-133 | 0.866 |
| Ext2 | SCZ25_EUR | top1 | 7143 | 0.253 | 1.81e-104 | 0.656 |
| Ext6 | PGC3_EUR | susie | 798 | 0.489 | 3.80e-49 | 0.85 |
| Ext6 | PGC3_EUR | top1 | 4533 | 0.196 | 1.10e-40 | 0.645 |
| Ext6 | SCZ25_EUR | susie | 798 | 0.45 | 4.43e-41 | 0.825 |
| Ext6 | SCZ25_EUR | top1 | 4589 | 0.174 | 1.26e-32 | 0.644 |
| Ext3 | PGC3_EUR | susie | 639 | 0.405 | 1.27e-26 | 0.75 |
| Ext3 | PGC3_EUR | top1 | 4881 | 0.113 | 2.64e-15 | 0.587 |
| Ext3 | SCZ25_EUR | susie | 639 | 0.357 | 1.24e-20 | 0.75 |
| Ext3 | SCZ25_EUR | top1 | 4935 | 0.097 | 8.69e-12 | 0.579 |
| Ext4 | PGC3_EUR | susie | 2788 | 0.468 | 6.65e-152 | 0.8 |
| Ext4 | PGC3_EUR | top1 | 7054 | 0.286 | 5.10e-133 | 0.685 |
| Ext4 | SCZ25_EUR | susie | 2790 | 0.436 | 9.35e-130 | 0.772 |
| Ext4 | SCZ25_EUR | top1 | 7113 | 0.258 | 1.44e-108 | 0.671 |
| Ext5 | PGC3_EUR | susie | 881 | 0.478 | 1.34e-51 | 0.824 |
| Ext5 | PGC3_EUR | top1 | 5387 | 0.169 | 1.30e-35 | 0.62 |
| Ext5 | SCZ25_EUR | susie | 881 | 0.442 | 1.68e-43 | 0.784 |
| Ext5 | SCZ25_EUR | top1 | 5435 | 0.155 | 1.56e-30 | 0.64 |
| Ext7 | PGC3_EUR | susie | 1480 | 0.436 | 1.51e-69 | 0.728 |
| Ext7 | PGC3_EUR | top1 | 5766 | 0.185 | 2.26e-45 | 0.602 |
| Ext7 | SCZ25_EUR | susie | 1480 | 0.401 | 3.33e-58 | 0.716 |
| Ext7 | SCZ25_EUR | top1 | 5817 | 0.171 | 2.16e-39 | 0.589 |
| Ext8 | PGC3_EUR | susie | 2294 | 0.48 | 3.13e-132 | 0.852 |
| Ext8 | PGC3_EUR | top1 | 6796 | 0.264 | 1.63e-108 | 0.687 |
| Ext8 | SCZ25_EUR | susie | 2295 | 0.444 | 1.51e-111 | 0.844 |
| Ext8 | SCZ25_EUR | top1 | 6855 | 0.249 | 1.88e-97 | 0.683 |

### Ext1: Random_Z vs fixed_z

| features | spearman_absfixed_vs_random_all | egenes | spearman_absfixed_vs_random_egenes | share_random_z_negative | share_fixed_z_negative | spearman_signedfixed_vs_random |
|---|---|---|---|---|---|---|
| 11893 | 0.902 | 11118 | 0.909 | 0.0191 | 0.529 | -0.0111 |

### Ext1 sample-size sensitivity (n = 757 vs 983)

| model | genes_757 | genes_983 | shared_snp_weights | weight_pearson | SCZ25_EUR_z_pearson | MDD_EUR_z_pearson | ABCD_thin_EUR_z_pearson |
|---|---|---|---|---|---|---|---|
| susie | 4322 | 4346 | 109125 | 0.989 | 0.986 | 0.998 | 0.997 |
| top1 | 11116 | 11116 | 11116 | 1 | 1 | 1 | 1 |

### Positive controls: Fig. 3c genes in their colocalising cell type(s) (SCZ25_EUR)

| symbol | coloc_classes | coloc_subtypes | tag | model | z | p | w_flag |
|---|---|---|---|---|---|---|---|
| CUL3 | Ast | Ast1,Ast4 | Ast | susie | -4.69 | 2.80e-06 | False |
| CUL3 | Ast | Ast1,Ast4 | Ast | top1 | -4.69 | 2.72e-06 | False |
| PTPRU | Ext |  | Ext1 | susie | -6.08 | 1.18e-09 | False |
| PTPRU | Ext |  | Ext1 | top1 | -4.1 | 4.14e-05 | False |
| PTPRU | Ext |  | Ext1_n983 | susie | -6.09 | 1.14e-09 | False |
| PTPRU | Ext |  | Ext1_n983 | top1 | -4.1 | 4.14e-05 | False |
| PTPRU | Ext |  | Ext | susie | -4.2 | 2.62e-05 | False |
| PTPRU | Ext |  | Ext | top1 | -4.6 | 4.14e-06 | False |
| PTPRU | Ext |  | Ext2 | susie |  |  | nan |
| PTPRU | Ext |  | Ext2 | top1 | -0.996 | 0.319 | False |
| PTPRU | Ext |  | Ext6 | susie |  |  | nan |
| PTPRU | Ext |  | Ext6 | top1 |  |  | nan |
| PTPRU | Ext |  | Ext3 | susie |  |  | nan |
| PTPRU | Ext |  | Ext3 | top1 |  |  | nan |
| PTPRU | Ext |  | Ext4 | susie |  |  | nan |
| PTPRU | Ext |  | Ext4 | top1 | -0.413 | 0.68 | False |
| PTPRU | Ext |  | Ext5 | susie |  |  | nan |
| PTPRU | Ext |  | Ext5 | top1 | -0.97 | 0.332 | False |
| PTPRU | Ext |  | Ext7 | susie |  |  | nan |
| PTPRU | Ext |  | Ext7 | top1 | -4.37 | 1.26e-05 | False |
| PTPRU | Ext |  | Ext8 | susie | 0.0827 | 0.934 | False |
| PTPRU | Ext |  | Ext8 | top1 | 0.09 | 0.928 | False |
| SNAP91 | OD | OD2,OD3 | OD | susie | 5.83 | 5.50e-09 | False |
| SNAP91 | OD | OD2,OD3 | OD | top1 | 7.34 | 2.18e-13 | False |
| RAB21 | OD |  | OD | susie |  |  | nan |
| RAB21 | OD |  | OD | top1 | -5.35 | 8.64e-08 | False |
| IPO11 | OD |  | OD | susie |  |  | nan |
| IPO11 | OD |  | OD | top1 | 8.62 | 6.63e-18 | False |
| SLC25A17 | IN |  | IN | susie |  |  | nan |
| SLC25A17 | IN |  | IN | top1 | 5.46 | 4.85e-08 | False |
| ACTR5 | IN | IN1,IN2 | IN | susie | 7.18 | 6.75e-13 | False |
| ACTR5 | IN | IN1,IN2 | IN | top1 | 7.18 | 6.95e-13 | False |
| PPARGC1A | Ext |  | Ext1 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext1 | top1 | 2.34 | 0.0194 | False |
| PPARGC1A | Ext |  | Ext1_n983 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext1_n983 | top1 | 2.34 | 0.0194 | False |
| PPARGC1A | Ext |  | Ext | susie | -1.12 | 0.261 | False |
| PPARGC1A | Ext |  | Ext | top1 | -1.12 | 0.265 | False |
| PPARGC1A | Ext |  | Ext2 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext2 | top1 | 1.95 | 0.0518 | False |
| PPARGC1A | Ext |  | Ext6 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext6 | top1 | -0.23 | 0.818 | False |
| PPARGC1A | Ext |  | Ext3 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext3 | top1 | 0.575 | 0.565 | False |
| PPARGC1A | Ext |  | Ext4 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext4 | top1 | 2.95 | 0.00318 | False |
| PPARGC1A | Ext |  | Ext5 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext5 | top1 | 1.39 | 0.164 | False |
| PPARGC1A | Ext |  | Ext7 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext7 | top1 | -0.923 | 0.356 | False |
| PPARGC1A | Ext |  | Ext8 | susie |  |  | nan |
| PPARGC1A | Ext |  | Ext8 | top1 | 1.85 | 0.0648 | False |
| FOXP1 | OD | OD1,OD2,OD3 | OD | susie | -6.08 | 1.21e-09 | False |
| FOXP1 | OD | OD1,OD2,OD3 | OD | top1 | -6.26 | 3.95e-10 | False |
| SPCS1 | IN |  | IN | susie |  |  | nan |
| SPCS1 | IN |  | IN | top1 | -9.53 | 1.61e-21 | False |
| AGO4 | OD |  | OD | susie | 5.65 | 1.63e-08 | False |
| AGO4 | OD |  | OD | top1 | 5.43 | 5.79e-08 | False |
| TMX1 | OPC | OPC1 | OPC | susie | 4.72 | 2.34e-06 | False |
| TMX1 | OPC | OPC1 | OPC | top1 | 4.1 | 4.05e-05 | False |
| RERE | Ast | Ast1,Ast2,Ast3 | Ast | susie | 4.76 | 1.96e-06 | False |
| RERE | Ast | Ast1,Ast2,Ast3 | Ast | top1 | 5.42 | 5.80e-08 | False |
| SLC4A10 | Ext | Ext2 | Ext1 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext1 | top1 | -1.16 | 0.247 | False |
| SLC4A10 | Ext | Ext2 | Ext1_n983 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext1_n983 | top1 | -1.16 | 0.247 | False |
| SLC4A10 | Ext | Ext2 | Ext | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext | top1 | -5.13 | 2.95e-07 | False |
| SLC4A10 | Ext | Ext2 | Ext2 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext2 | top1 | -5.49 | 4.04e-08 | False |
| SLC4A10 | Ext | Ext2 | Ext6 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext6 | top1 | -1.86 | 0.0631 | False |
| SLC4A10 | Ext | Ext2 | Ext3 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext3 | top1 | 0.748 | 0.454 | False |
| SLC4A10 | Ext | Ext2 | Ext4 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext4 | top1 | -2.56 | 0.0105 | False |
| SLC4A10 | Ext | Ext2 | Ext5 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext5 | top1 |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext7 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext7 | top1 | -0.708 | 0.479 | False |
| SLC4A10 | Ext | Ext2 | Ext8 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext8 | top1 | 0.57 | 0.569 | False |
| CNTN4 | Ext |  | Ext1 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext1 | top1 | 1.09 | 0.276 | False |
| CNTN4 | Ext |  | Ext1_n983 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext1_n983 | top1 | 1.09 | 0.276 | False |
| CNTN4 | Ext |  | Ext | susie |  |  | nan |
| CNTN4 | Ext |  | Ext | top1 | 7.2 | 6.02e-13 | False |
| CNTN4 | Ext |  | Ext2 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext2 | top1 | 6.1 | 1.05e-09 | False |
| CNTN4 | Ext |  | Ext6 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext6 | top1 | 0.159 | 0.874 | False |
| CNTN4 | Ext |  | Ext3 | susie | 0.449 | 0.653 | False |
| CNTN4 | Ext |  | Ext3 | top1 | 0.449 | 0.653 | False |
| CNTN4 | Ext |  | Ext4 | susie | 5.62 | 1.90e-08 | False |
| CNTN4 | Ext |  | Ext4 | top1 | 5.89 | 3.85e-09 | False |
| CNTN4 | Ext |  | Ext5 | susie | 0.768 | 0.443 | False |
| CNTN4 | Ext |  | Ext5 | top1 | 0.574 | 0.566 | False |
| CNTN4 | Ext |  | Ext7 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext7 | top1 |  |  | nan |
| CNTN4 | Ext |  | Ext8 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext8 | top1 | 0.0747 | 0.94 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext1 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext1 | top1 | 5.62 | 1.96e-08 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext1_n983 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext1_n983 | top1 | 5.62 | 1.96e-08 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext | susie | 6.08 | 1.21e-09 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext | top1 | 5.74 | 9.65e-09 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext2 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext2 | top1 | 4.27 | 1.93e-05 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext6 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext6 | top1 |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext3 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext3 | top1 |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext4 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext4 | top1 | 1.18 | 0.237 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext5 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext5 | top1 | -1.37 | 0.171 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext7 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext7 | top1 |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext8 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext8 | top1 | -0.0632 | 0.95 | False |
| SH3RF1 | MG,MiGA | Ext7,MG2,MG3 | MG | susie | 6.4 | 1.56e-10 | False |
| SH3RF1 | MG,MiGA | Ext7,MG2,MG3 | MG | top1 | 5.72 | 1.06e-08 | False |
| SH3RF1 | MG,MiGA | Ext7,MG2,MG3 | Ext7 | susie |  |  | nan |
| SH3RF1 | MG,MiGA | Ext7,MG2,MG3 | Ext7 | top1 | -6.42 | 1.36e-10 | False |
| KCNB1 | IN | IN7 | IN | susie | 6.13 | 8.90e-10 | False |
| KCNB1 | IN | IN7 | IN | top1 | 6.02 | 1.74e-09 | False |
| GNL3 | IN | Ext4,IN4,IN6,IN7 | IN | susie | 7.41 | 1.27e-13 | False |
| GNL3 | IN | Ext4,IN4,IN6,IN7 | IN | top1 | 8.17 | 2.97e-16 | False |
| GNL3 | IN | Ext4,IN4,IN6,IN7 | Ext4 | susie |  |  | nan |
| GNL3 | IN | Ext4,IN4,IN6,IN7 | Ext4 | top1 | 9.3 | 1.47e-20 | False |
| CACNB2 | Ext |  | Ext1 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext1 | top1 | 4.17 | 3.05e-05 | False |
| CACNB2 | Ext |  | Ext1_n983 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext1_n983 | top1 | 4.17 | 3.05e-05 | False |
| CACNB2 | Ext |  | Ext | susie |  |  | nan |
| CACNB2 | Ext |  | Ext | top1 | 6.44 | 1.21e-10 | False |
| CACNB2 | Ext |  | Ext2 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext2 | top1 | -1.91 | 0.0559 | False |
| CACNB2 | Ext |  | Ext6 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext6 | top1 | -0.913 | 0.361 | False |
| CACNB2 | Ext |  | Ext3 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext3 | top1 | -1.26 | 0.207 | False |
| CACNB2 | Ext |  | Ext4 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext4 | top1 | 4.17 | 3.05e-05 | False |
| CACNB2 | Ext |  | Ext5 | susie | -1.87 | 0.062 | False |
| CACNB2 | Ext |  | Ext5 | top1 | -1.62 | 0.106 | False |
| CACNB2 | Ext |  | Ext7 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext7 | top1 | 0.827 | 0.408 | False |
| CACNB2 | Ext |  | Ext8 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext8 | top1 | 0.693 | 0.488 | False |
| WBP1L | IN | Ext2,Ext7 | IN | susie | 9.21 | 3.24e-20 | False |
| WBP1L | IN | Ext2,Ext7 | IN | top1 | 9.03 | 1.65e-19 | False |
| WBP1L | IN | Ext2,Ext7 | Ext2 | susie | 8.85 | 8.37e-19 | False |
| WBP1L | IN | Ext2,Ext7 | Ext2 | top1 | 5.17 | 2.28e-07 | False |
| WBP1L | IN | Ext2,Ext7 | Ext7 | susie |  |  | nan |
| WBP1L | IN | Ext2,Ext7 | Ext7 | top1 | 4.98 | 6.23e-07 | False |
| FUT9 | IN | Ast4,IN6,IN7,OD1,OD2 | IN | susie | 5.57 | 2.59e-08 | False |
| FUT9 | IN | Ast4,IN6,IN7,OD1,OD2 | IN | top1 | 5.76 | 8.34e-09 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1 | susie | -4.63 | 3.65e-06 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1 | top1 | -5.4 | 6.60e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1_n983 | susie | -4.66 | 3.13e-06 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1_n983 | top1 | -5.4 | 6.60e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext | susie | -5.35 | 8.78e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext | top1 | -5.4 | 6.60e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext2 | susie | -4.86 | 1.16e-06 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext2 | top1 | -5.05 | 4.31e-07 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext6 | susie | -5.66 | 1.52e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext6 | top1 | -5.78 | 7.65e-09 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext3 | susie |  |  | nan |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext3 | top1 |  |  | nan |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext4 | susie | -5.47 | 4.62e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext4 | top1 | -5.4 | 6.60e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext5 | susie |  |  | nan |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext5 | top1 | -5.4 | 6.60e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext7 | susie |  |  | nan |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext7 | top1 | -5.05 | 4.42e-07 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext8 | susie | -5.63 | 1.77e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext8 | top1 | -5.4 | 6.60e-08 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext1 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext1 | top1 | 5.15 | 2.59e-07 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext1_n983 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext1_n983 | top1 | 5.15 | 2.59e-07 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext | top1 | 5.15 | 2.59e-07 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext2 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext2 | top1 | 0.902 | 0.367 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext6 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext6 | top1 |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext3 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext3 | top1 |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext4 | susie | 4.83 | 1.40e-06 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext4 | top1 | 4.86 | 1.15e-06 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext5 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext5 | top1 | 0.95 | 0.342 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext7 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext7 | top1 |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext8 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext8 | top1 | -1.45 | 0.146 | False |
| RGS6 | Ext |  | Ext1 | susie |  |  | nan |
| RGS6 | Ext |  | Ext1 | top1 | -4.86 | 1.15e-06 | False |
| RGS6 | Ext |  | Ext1_n983 | susie |  |  | nan |
| RGS6 | Ext |  | Ext1_n983 | top1 | -4.86 | 1.15e-06 | False |
| RGS6 | Ext |  | Ext | susie | -8.95 | 3.48e-19 | False |
| RGS6 | Ext |  | Ext | top1 | -7.01 | 2.36e-12 | False |
| RGS6 | Ext |  | Ext2 | susie | -5.97 | 2.37e-09 | False |
| RGS6 | Ext |  | Ext2 | top1 | -8.05 | 8.47e-16 | False |
| RGS6 | Ext |  | Ext6 | susie | -1.2 | 0.231 | False |
| RGS6 | Ext |  | Ext6 | top1 | -1.59 | 0.112 | False |
| RGS6 | Ext |  | Ext3 | susie | -0.451 | 0.652 | False |
| RGS6 | Ext |  | Ext3 | top1 | -0.451 | 0.652 | False |
| RGS6 | Ext |  | Ext4 | susie |  |  | nan |
| RGS6 | Ext |  | Ext4 | top1 | 0.314 | 0.753 | False |
| RGS6 | Ext |  | Ext5 | susie | -0.809 | 0.419 | False |
| RGS6 | Ext |  | Ext5 | top1 | -0.451 | 0.652 | False |
| RGS6 | Ext |  | Ext7 | susie |  |  | nan |
| RGS6 | Ext |  | Ext7 | top1 | 0.0198 | 0.984 | False |
| RGS6 | Ext |  | Ext8 | susie |  |  | nan |
| RGS6 | Ext |  | Ext8 | top1 | -4.86 | 1.15e-06 | False |
| CSMD1 | OPC |  | OPC | susie | 0.499 | 0.618 | False |
| CSMD1 | OPC |  | OPC | top1 | 0.617 | 0.537 | False |
| ANKS1B | OPC | OPC1 | OPC | susie | 5.18 | 2.22e-07 | False |
| ANKS1B | OPC | OPC1 | OPC | top1 | 5.25 | 1.55e-07 | False |
| MYT1L | OPC | Ext5,OPC1 | OPC | susie | 5.91 | 3.38e-09 | False |
| MYT1L | OPC | Ext5,OPC1 | OPC | top1 | 5.76 | 8.27e-09 | False |
| MYT1L | OPC | Ext5,OPC1 | Ext5 | susie |  |  | nan |
| MYT1L | OPC | Ext5,OPC1 | Ext5 | top1 | -6.51 | 7.37e-11 | False |
| DPYD | OD | Ast4 | OD | susie | 8.79 | 1.49e-18 | False |
| DPYD | OD | Ast4 | OD | top1 | 8.78 | 1.65e-18 | False |
| PTPRT | OPC |  | OPC | susie | 4.9 | 9.81e-07 | False |
| PTPRT | OPC |  | OPC | top1 | 4.89 | 9.99e-07 | False |
| CNOT7 | IN | Ext4 | IN | susie | -1.49 | 0.135 | False |
| CNOT7 | IN | Ext4 | IN | top1 | -1.74 | 0.0816 | False |
| CNOT7 | IN | Ext4 | Ext4 | susie | 4.7 | 2.63e-06 | False |
| CNOT7 | IN | Ext4 | Ext4 | top1 | 4.67 | 3.03e-06 | False |
| SPIRE2 | MG | OPC2 | MG | susie | -5.54 | 3.11e-08 | False |
| SPIRE2 | MG | OPC2 | MG | top1 | -3.36 | 7.79e-04 | False |
| DPYD-AS1 | OD |  | OD | susie |  |  | nan |
| DPYD-AS1 | OD |  | OD | top1 | 8.78 | 1.65e-18 | False |
| SOX2-OT | IN | IN4 | IN | susie |  |  | nan |
| SOX2-OT | IN | IN4 | IN | top1 | 6.24 | 4.24e-10 | False |

### CNTN4 and SLC4A10 in every cell type (SCZ25_EUR)

| ('symbol', '') | ('tag', '') | ('p', 'susie') | ('p', 'top1') | ('z', 'susie') | ('z', 'top1') |
|---|---|---|---|---|---|
| CNTN4 | Ext |  | 6.02e-13 |  | 7.2 |
| CNTN4 | Ext1 |  | 0.276 |  | 1.09 |
| CNTN4 | Ext1_n983 |  | 0.276 |  | 1.09 |
| CNTN4 | Ext2 |  | 1.05e-09 |  | 6.1 |
| CNTN4 | Ext3 | 0.653 | 0.653 | 0.449 | 0.449 |
| CNTN4 | Ext4 | 1.90e-08 | 3.85e-09 | 5.62 | 5.89 |
| CNTN4 | Ext5 | 0.443 | 0.566 | 0.768 | 0.574 |
| CNTN4 | Ext6 |  | 0.874 |  | 0.159 |
| CNTN4 | Ext8 |  | 0.94 |  | 0.0747 |
| CNTN4 | IN | 0.0336 | 0.27 | -2.13 | -1.1 |
| CNTN4 | MG |  | 0.874 |  | -0.158 |
| CNTN4 | OD |  | 0.458 |  | 0.742 |
| CNTN4 | OPC | 0.0571 | 0.43 | 1.9 | -0.79 |
| SLC4A10 | Ast |  | 0.751 |  | -0.317 |
| SLC4A10 | Ext |  | 2.95e-07 |  | -5.13 |
| SLC4A10 | Ext1 |  | 0.247 |  | -1.16 |
| SLC4A10 | Ext1_n983 |  | 0.247 |  | -1.16 |
| SLC4A10 | Ext2 |  | 4.04e-08 |  | -5.49 |
| SLC4A10 | Ext3 |  | 0.454 |  | 0.748 |
| SLC4A10 | Ext4 |  | 0.0105 |  | -2.56 |
| SLC4A10 | Ext6 |  | 0.0631 |  | -1.86 |
| SLC4A10 | Ext7 |  | 0.479 |  | -0.708 |
| SLC4A10 | Ext8 |  | 0.569 |  | 0.57 |
| SLC4A10 | IN | 0.192 | 0.117 | -1.31 | -1.57 |
| SLC4A10 | OPC | 1.72e-04 | 1.66e-04 | -3.76 | -3.77 |

### Top 10 genes: SCZ25_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | MAPT-AS1(+29.3!) POC1B-AS1(+26.4!) POLE(+22.3!) VPS13A(-19.4!) SPCS1(-9.1) YJEFN3(-8.5) GATAD2A(-8.3) FURIN(-8.1) ARL3(-7.8) INO80E(+7.8) |
| Ext1 | top1 | YJEFN3(-8.6) GATAD2A(-8.4) RFT1(+8.2) FURIN(-8.2) ATG13(-8.1) ARL3(-7.8) MAU2(-7.7) INO80E(+7.7) INPP5D(-7.7) RPS17(-7.6) |
| Ext1_n983 | susie | MAPT-AS1(+29.4!) POC1B-AS1(+27.4!) POLE(+26.7!) PIGP(+22.2!) VPS13A(-19.5!) LRRC37A2(-10.3!) SPCS1(-9.1) YJEFN3(-8.5) GATAD2A(-8.3) FURIN(-8.1) |
| Ext1_n983 | top1 | YJEFN3(-8.6) GATAD2A(-8.4) RFT1(+8.2) FURIN(-8.2) ATG13(-8.1) ARL3(-7.8) MAU2(-7.7) INO80E(+7.7) INPP5D(-7.7) RPS17(-7.6) |
| Ext | susie | SLC17A5(+40.2!) PPP4R3B(-36.8!) VPS13A(-29.1!) CPLANE1(+24.2!) NRDC(+21.1!) MAPT-AS1(+16.6!) NPEPL1(-16.2!) SELENOI(+13.6!) EEF1AKMT1(-12.8!) MIR2052HG(-11.3!) |
| Ext | top1 | GATAD2A(-8.6) YJEFN3(-8.5) RFT1(+8.2) FURIN(-8.0) ARL3(-7.8) INO80E(+7.7) MAU2(-7.7) ATG16L1(-7.6) NEK4(+7.5) PBRM1(-7.5) |
| OPC | susie | VPS50(-67.9!) KMT2E(+8.1) INO80E(+7.7) RBM26(+7.7) TMEM219(+7.5) AS3MT(+7.1) HSPD1(+7.0) FTCDNL1(-7.0) YJEFN3(-7.0) SFMBT1(-6.8) |
| OPC | top1 | NT5DC2(+8.5) KMT2E(+8.3) INO80E(+7.8) TMEM219(+7.8) RBM26(+7.6) ARL3(-7.5) YJEFN3(-7.3) ZSCAN2(-7.3) AS3MT(+7.1) WDR73(+7.0) |
| IN | susie | TUT4(+87.1!) SELENOI(+24.0!) ANKRD36(-23.0!) AFDN(-20.4!) TTC3(-16.3!) MTHFS(-12.4!) DENND5B(-12.1!) WBP1L(+9.2) MAPT-AS1(+9.0) FURIN(-8.2) |
| IN | top1 | CACNA1C(-10.2) SPCS1(-9.5) WBP1L(+9.0) GATAD2A(-8.6) GNL3(+8.2) ARL3(-7.8) TAOK2(+7.7) FURIN(-7.7) INO80E(+7.7) DOC2A(+7.7) |
| OD | susie | TUT4(-14.7!) FOLH1(+10.0!) DPYD(+8.8) TOM1L1(+8.7!) GATAD2A(-8.4) KCTD13(-7.9) INO80E(+7.8) ASPHD1(-7.8) KMT2E(+7.6) ABCB9(+7.2) |
| OD | top1 | DPYD(+8.8) DPYD-AS1(+8.8) IPO11(+8.6) SRPK2(-8.3) ASPHD1(-7.9) KCTD13(-7.8) SEZ6L2(-7.8) FURIN(-7.7) INO80E(+7.7) GIGYF2(+7.6) |
| Ast | susie | TMC1(-12.0!) KCTD13(+8.0) TMEM219(+7.9) ARL3(-7.9) INO80E(+7.8) NDUFA6-DT(+6.7) SF3B1(+6.6) PCCB(+6.6) ARL17B(-6.5) TKT(+6.5) |
| Ast | top1 | GLT8D1(-8.8) MAU2(+8.6) TMEM219(+7.9) KCTD13(+7.9) PSMG3-AS1(+7.9) INO80E(+7.8) ARL3(-7.8) ASPHD1(+7.8) LINC01004(+7.7) PPP1R13B(-7.0) |
| MG | susie | INO80E(+7.8) ITIH4(-7.3) AS3MT(+6.9) SMDT1(+6.7) LINC01004(+6.6) PRKD3(+6.4) SH3RF1(+6.4) KANSL1-AS1(-6.2) ARL17B(-6.2) LRRC37A(-6.1) |
| MG | top1 | NPM3(-8.0) TAOK2(+7.8) INO80E(+7.7) GIGYF2(+7.6) WDR82(-7.6) TXNL1(+7.4) ITIH4(-7.3) PPM1M(+7.3) KMT2E(+7.1) FTCDNL1(-7.1) |
| Ext2 | susie | POC1B-AS1(+31.0!) WBP1L(+8.9) PBRM1(-8.6) FURIN(-8.0) GATAD2A(-7.9) ARL3(-7.9) GLT8D1(-7.6) NEK4(+7.5) LRRC37A2(-7.1) ACTR5(+6.9) |
| Ext2 | top1 | NEK4(+9.3) PBRM1(-8.2) RGS6(-8.0) ARL3(-7.8) TRIM8(+7.8) MAU2(-7.7) FURIN(-7.7) PPP1R13B(-7.6) SPCS1(-7.5) GLT8D1(-7.2) |
| Ext6 | susie | PPP1R13B(-7.4) LRRC37A(-7.3) LRRC37A2(-7.1) FTCDNL1(-6.8) PCCB(+6.6) MRPS30(-6.6) MAPK3(+6.5) ARL17B(-6.4) MAPT(+6.2) MAPT-AS1(+6.2) |
| Ext6 | top1 | SPCS1(-9.3) GLB1L2(+8.2) CACNA1C(+8.2) GLT8D1(-7.4) GATAD2A(-7.3) NSF(-7.0) LRRC37A(-7.0) STIMATE(-6.9) NDUFA6(+6.7) LINC01877(-6.6) |
| Ext3 | susie | KIF6(-9.1!) GATAD2A(-8.4) LINC02008(-8.0!) FTCDNL1(-7.0) LRRC37A(-7.0) LRRC37A2(-6.5) WBP2NL(+6.5) ARL17B(-6.2) MAPT-AS1(+6.2) PLEKHM1(+6.1) |
| Ext3 | top1 | GATAD2A(-8.5) ARL3(-7.8) DNAJA4(-7.7) HAPLN4(+7.6) GLT8D1(-7.3) FTCDNL1(-7.0) CCDC24(-6.8) MRPS30(-6.7) WBP2NL(+6.4) MAPT-AS1(+6.4) |
| Ext4 | susie | MAPT-AS1(+15.5!) AGAP9(+14.0!) CEMIP2(+13.4!) LRRC37A2(-10.5!) FURIN(-8.3) GATAD2A(-8.1) MFSD13A(+7.9!) ARL3(-7.8) INO80E(+7.8) TRIM8(+7.6) |
| Ext4 | top1 | GNL3(+9.3) NT5DC2(+8.5) FURIN(-8.2) TSPAN31(+7.8) INO80E(+7.8) ARL3(-7.8) CORO1A(+7.8) GLT8D1(-7.6) NEK4(+7.3) GATAD2A(-7.3) |
| Ext5 | susie | PRSS3(+15.6!) ARL17B(-10.4!) LRRC37A2(-9.2!) AC159540.2(-9.1!) STIMATE(-7.8) C12orf65(+7.6) GLT8D1(-7.4) ITIH4(-7.1) LRRC37A(-6.9) UBE2Q2L(+6.6) |
| Ext5 | top1 | SPCS1(-9.3) NEK4(+9.1) GNL3(+8.6) PTN(-7.9) ITIH4(-7.4) STIMATE(-7.3) GLT8D1(-7.3) FTCDNL1(-7.2) PPP1R13B(-7.0) LRRC37A(-6.8) |
| Ext7 | susie | MFSD13A(+17.1!) TUT4(-12.9!) UBR5-AS1(+10.9!) LRRC37A2(-9.1!) COP1(+8.6!) GLT8D1(-7.6) NT5DC2(+7.3) SFMBT1(-7.1) CACNA1I(+7.0) FTCDNL1(-7.0) |
| Ext7 | top1 | YJEFN3(-8.6) PHF7(-7.9) ASPHD1(+7.8) ADGRB1(-7.6) INO80E(+7.4) SFMBT1(-7.4) NEK4(+7.4) CACNA1I(+7.2) GLT8D1(-7.2) GATAD2A(-7.1) |
| Ext8 | susie | ZNF66(-19.7!) MAPT-AS1(+15.4!) KCTD21-AS1(-11.8!) KLHL23(+9.8!) ZBTB37(+8.9!) YJEFN3(-8.5) GATAD2A(-8.3) FURIN(-8.2) ARL3(-7.9) INO80E(+7.8) |
| Ext8 | top1 | YJEFN3(-8.6) GATAD2A(-8.5) FURIN(-8.2) ARL3(-7.9) INO80E(+7.8) GLT8D1(-7.6) RPS17(-7.6) NTM(-7.4) SPCS1(-7.4) NEK4(+7.4) |

### Top 10 genes: PGC3_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | MAPT-AS1(+18.6!) SPCS1(-8.7) FURIN(-8.5) YJEFN3(-8.0) GATAD2A(-8.0) INO80E(+7.6) TAOK2(+7.3) PHF21A(-7.2) NEK4(+7.1) LINC00634(+7.1) |
| Ext1 | top1 | FURIN(-8.6) YJEFN3(-8.2) RFT1(+8.0) GATAD2A(-7.8) INPP5D(-7.8) ATG13(-7.7) INO80E(+7.5) PPP1R13B(-7.3) AS3MT(+7.2) FTCDNL1(-7.2) |
| Ext1_n983 | susie | MAPT-AS1(+18.7!) SPCS1(-8.7) FURIN(-8.5) YJEFN3(-8.0) GATAD2A(-8.0) PIGP(+7.7!) INO80E(+7.6) TAOK2(+7.3) PHF21A(-7.2) NEK4(+7.1) |
| Ext1_n983 | top1 | FURIN(-8.6) YJEFN3(-8.2) RFT1(+8.0) GATAD2A(-7.8) INPP5D(-7.8) ATG13(-7.7) INO80E(+7.5) PPP1R13B(-7.3) AS3MT(+7.2) FTCDNL1(-7.2) |
| Ext | susie | NRDC(+30.9!) PPP4R3B(-28.0!) CPLANE1(+16.7!) EEF1AKMT1(-15.9!) AGAP9(+15.5!) MFSD13A(+14.0!) NPEPL1(-11.8!) MAPT-AS1(+11.8!) MIR2052HG(-11.2!) SNX14(+9.0!) |
| Ext | top1 | FURIN(-8.4) GATAD2A(-8.2) ATG16L1(-8.1) RFT1(+8.0) YJEFN3(-7.8) INO80E(+7.5) FTCDNL1(-7.4) AS3MT(+7.2) NEK4(+7.1) WBP2NL(+7.1) |
| OPC | susie | VPS50(-45.4!) KMT2E(+7.7) INO80E(+7.6) FTCDNL1(-7.5) TMEM219(+7.4) AS3MT(+7.2) RBM26(+6.8) HSPD1(+6.7) YJEFN3(-6.5) PCCB(+6.3) |
| OPC | top1 | KMT2E(+7.7) INO80E(+7.6) TMEM219(+7.5) NT5DC2(+7.4) FTCDNL1(-7.4) AS3MT(+7.2) YJEFN3(-6.9) RBM26(+6.8) ZSCAN2(-6.8) SMDT1(+6.6) |
| IN | susie | TUT4(-24.3!) SELENOI(+13.9!) AFDN(-12.9!) MTHFS(-10.4!) WBP1L(+8.2) FURIN(-8.0) INO80E(+7.5) KMT2E(+7.4) PPP1R13B(-7.3) GNL3(+7.0) |
| IN | top1 | SPCS1(-9.1) CACNA1C(-8.8) GATAD2A(-8.2) WBP1L(+8.1) GNL3(+8.0) TAOK2(+7.5) FTCDNL1(-7.5) DOC2A(+7.5) INO80E(+7.5) FURIN(-7.3) |
| OD | susie | TUT4(-9.9!) DPYD(+8.8) GATAD2A(-7.9) TOM1L1(+7.7!) INO80E(+7.6) KCTD13(-7.4) KMT2E(+7.4) ASPHD1(-7.3) AS3MT(+7.2) HSPD1(+6.8) |
| OD | top1 | DPYD(+8.8) DPYD-AS1(+8.8) GIGYF2(+8.1) IPO11(+8.1) ASPHD1(-7.7) SRPK2(-7.6) INO80E(+7.5) KCTD13(-7.5) SEZ6L2(-7.5) FTCDNL1(-7.4) |
| Ast | susie | TMC1(-14.2!) ATP5F1C(-8.1!) KCTD13(+7.8) TMEM219(+7.7) INO80E(+7.5) ARL3(-6.7) NAGA(+6.7) IRF3(-6.6) NDUFA6-DT(+6.6) DPYD(-6.6) |
| Ast | top1 | MAU2(+8.0) GLT8D1(-7.7) TMEM219(+7.7) KCTD13(+7.7) INO80E(+7.5) ASPHD1(+7.5) FTCDNL1(-7.2) PSMG3-AS1(+7.2) WBP2NL(+7.1) NT5C2(+7.1) |
| MG | susie | INO80E(+7.6) AS3MT(+7.0) TNFRSF13C(-6.8) LINC01004(+6.7) SMDT1(+6.6) ITIH4(-6.5) IMMP2L(-6.1) KANSL1-AS1(-5.7) PLEKHM1(+5.6) ARL17B(-5.6) |
| MG | top1 | GIGYF2(+8.2) TAOK2(+7.6) INO80E(+7.5) FTCDNL1(-7.4) WDR82(-7.1) AS3MT(+7.0) PPM1M(+6.9) KMT2E(+6.9) GTF3C3(+6.8) TXNL1(+6.8) |
| Ext2 | susie | MFSD13A(+9.1!) PBRM1(-8.2) WBP1L(+7.7) FURIN(-7.6) GATAD2A(-7.5) GLT8D1(-7.1) NEK4(+7.1) ARL3(-6.7) C11orf49(-6.7) SMDT1(+6.6) |
| Ext2 | top1 | NEK4(+8.8) PBRM1(-8.0) SNORC(+7.6) FTCDNL1(-7.5) PPP1R13B(-7.3) AS3MT(+7.2) FURIN(-7.2) NISCH(-7.1) WBP2NL(+7.1) IREB2(-7.0) |
| Ext6 | susie | FTCDNL1(-7.4) PPP1R13B(-7.3) LRRC37A(-6.8) LRRC37A2(-6.5) PCCB(+6.4) MAPK3(+6.2) MRPS30(-5.9) SLCO4C1(-5.7) MAPT(+5.7) MAPT-AS1(+5.7) |
| Ext6 | top1 | SPCS1(-8.7) LINC01877(-7.6) BZW1(-7.3) SPATS2L(-7.2) GLB1L2(+7.1) CACNA1C(+7.1) FTCDNL1(-7.0) GATAD2A(-6.9) GLT8D1(-6.8) NSF(-6.5) |
| Ext3 | susie | LINC02008(-13.0!) GATAD2A(-7.9) FTCDNL1(-7.4) WBP2NL(+7.1) DHFR2(-6.4!) LRRC37A(-6.4) LRRC37A2(-5.9) MAPT-AS1(+5.7) PLEKHM1(+5.6) ACTR1B(-5.6) |
| Ext3 | top1 | GATAD2A(-7.8) DNAJA4(-7.8) FTCDNL1(-7.4) SNORC(+7.2) WBP2NL(+7.1) GLT8D1(-6.9) HAPLN4(+6.8) ARL3(-6.7) PPP1R13B(-6.7) MRPS30(-6.1) |
| Ext4 | susie | AGAP9(+26.8!) MFSD13A(+19.4!) CEMIP2(+14.4!) MAPT-AS1(+11.4!) FURIN(-8.6) GATAD2A(-7.6) INO80E(+7.5) GLT8D1(-7.0) NEK4(+7.0) ARL3(-6.7) |
| Ext4 | top1 | GNL3(+8.8) FURIN(-8.6) TSPAN31(+8.0) INO80E(+7.5) CORO1A(+7.5) NT5DC2(+7.4) GLT8D1(-7.2) NISCH(-7.1) NEK4(+6.9) WBP2NL(+6.9) |
| Ext5 | susie | LRRC37A2(-9.4!) STIMATE(-7.3) WBP2NL(+7.1) GLT8D1(-7.0) C12orf65(+6.8) LRRC37A(-6.5) ITIH4(-6.3) MAPK3(+6.2) UBE2Q2L(+6.1) CATSPER2(-5.8) |
| Ext5 | top1 | SPCS1(-8.8) NEK4(+8.6) GNL3(+8.0) FTCDNL1(-7.5) PTN(-7.4) WBP2NL(+7.1) GLT8D1(-6.8) STIMATE(-6.6) SMDT1(+6.6) ITIH4(-6.5) |
| Ext7 | susie | UBR5-AS1(-15.0!) TUT4(-14.2!) MFSD13A(+13.1!) COP1(+8.8!) FTCDNL1(-7.4) GLT8D1(-7.2) WBP2NL(+7.1) NT5DC2(+7.0) C11orf49(-6.7) SMDT1(+6.6) |
| Ext7 | top1 | YJEFN3(-8.2) PHF7(-7.7) ASPHD1(+7.5) FTCDNL1(-7.2) INO80E(+7.0) GLT8D1(-6.9) WBP2NL(+6.9) GATAD2A(-6.9) NEK4(+6.9) ADGRB1(-6.8) |
| Ext8 | susie | ZNF66(-15.5!) KLHL23(+8.9!) ZBTB37(+8.9!) FURIN(-8.6) YJEFN3(-8.0) GATAD2A(-7.8) INO80E(+7.6) MFSD13A(+7.2!) WBP2NL(+7.1) GLT8D1(-7.0) |
| Ext8 | top1 | FURIN(-8.6) YJEFN3(-8.2) GATAD2A(-7.9) INO80E(+7.5) WBP2NL(+7.1) GLT8D1(-7.1) PPP1R13B(-6.9) NEK4(+6.9) SPCS1(-6.8) ARL3(-6.8) |

### Top 10 genes: MDD_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | VPS13A(-94.3!) POC1B-AS1(+41.8!) MAPT-AS1(-23.5!) POLE(-17.8!) SLC12A5(+8.5) RAB27B(+8.3) MED20(-7.8!) ODR4(+7.6!) ZCCHC7(-7.6) MED19(-7.3) |
| Ext1 | top1 | SLC12A5(+8.5) RAB27B(+8.1) ZCCHC7(-7.4) MED19(-7.3) SREBF2(-7.3) SLC25A17(-6.7) ZNF445(-6.6) ZKSCAN7-AS1(-6.6) PSEN2(+6.5) ZNF852(-6.5) |
| Ext1_n983 | susie | VPS13A(-94.8!) POC1B-AS1(+43.3!) MAPT-AS1(-23.6!) POLE(-20.8!) SLC12A5(+8.7) RAB27B(+8.3) MED20(-8.2!) CCDC32(-7.7!) ZCCHC7(-7.5) MED19(-7.3) |
| Ext1_n983 | top1 | SLC12A5(+8.5) RAB27B(+8.1) ZCCHC7(-7.4) MED19(-7.3) SREBF2(-7.3) SLC25A17(-6.7) ZNF445(-6.6) ZKSCAN7-AS1(-6.6) PSEN2(+6.5) ZNF852(-6.5) |
| Ext | susie | VPS13A(-144.5!) LINC01934(+44.8!) NBL1(+38.1!) MFSD13A(-36.2!) MIR2052HG(+27.7!) EEF1AKMT1(+26.6!) SLC17A5(+20.5!) POC1B-AS1(+19.9!) HDGFL3(+17.3!) ANAPC5(-16.0!) |
| Ext | top1 | DCC(+8.9) SLC12A5(+8.5) RAB27B(+8.2) SLC25A51(-7.5) BDNF(+7.4) GRIK3(+7.0) ZNF197(-6.9) ASXL3(+6.9) ZKSCAN7-AS1(-6.9) ZKSCAN7(-6.6) |
| OPC | susie | PHLDB1(+16.4!) SIDT1(+11.0!) NEGR1(+8.1) ASCC3(+7.5) CCDC73(-6.9) FAM120AOS(+6.7) BAG5(-6.1) SFMBT1(-5.8) SUDS3(-5.7) CYP7B1(-5.6) |
| OPC | top1 | SLC25A51(-7.9) NEGR1(+7.7) EIF3M(-7.3) CCDC73(-7.0) SLC25A17(+6.8) ASCC3(+6.8) CD59(-6.7) NDUFS3(-6.7) FAM120AOS(+6.6) NT5DC2(+6.4) |
| IN | susie | TUT4(-54.5!) POC1B-AS1(+41.7!) DENND5B(+37.3!) MTHFS(+32.1!) NIPSNAP2(-30.0!) TTC3(+25.4!) SELENOI(-19.3!) AFDN(+14.2!) AC009228.1(+11.8!) MFSD13A(-11.5!) |
| IN | top1 | SLC12A5(+8.5) RAB27B(+8.2) ST13(-7.3) SLC25A17(+7.2) KLC1(+7.0) ZKSCAN7(-6.9) ASCC3(+6.9) ZNF197(-6.7) ZNF852(-6.7) NEGR1(+6.6) |
| OD | susie | LINC01725(+15.0!) SDHA(+10.2!) MED20(-8.6!) VRK2(-7.8) ZCCHC7(-7.6) ASCC3(+7.4) MMP24OS(+6.8!) NEGR1(+6.7) FAM120AOS(+6.4) YLPM1(+6.3) |
| OD | top1 | ZCCHC7(-7.6) VRK2(-7.2) ASCC3(+7.0) ZKSCAN7(-6.8) ZKSCAN7-AS1(-6.4) FAM120AOS(+6.3) FADS1(-6.3) ECHS1(+6.1) FADS2(-6.1) YLPM1(+6.1) |
| Ast | susie | ZCCHC7(-7.7) DNAH7(-7.6!) DIAPH3(-6.6) FAM120AOS(+6.5) MYOSLID(-6.5) RASGRP1(+6.5) FADS1(-6.3) DLST(-6.3) KLC1(+6.3) NNT-AS1(+6.2!) |
| Ast | top1 | ZCCHC7(-7.9) KLC1(+7.5) CTNND1(+7.4) MYOSLID(-6.9) GRIK2(-6.8) GLT8D1(-6.7) RASGRP1(+6.6) DIAPH3(-6.6) KIAA1143(-6.5) DLST(-6.3) |
| MG | susie | FAM120AOS(+6.7) COQ8A(-6.7) KANSL1(+6.5!) PTCH1(+6.2) ITIH4(-6.1) DENND1B(+5.7) XRCC3(+5.6) SLC2A13(-5.5) METTL9(-5.4) C12orf40(-5.4) |
| MG | top1 | ZCCHC7(-7.2) PPP1R14B(+6.7) RANGAP1(+6.7) FAM120AOS(+6.6) COQ8A(-6.5) SLC25A17(+6.4) ZNF852(-6.2) MOCS1(-6.2) TXNL1(+6.1) PTCH1(+6.1) |
| Ext2 | susie | POC1B-AS1(+29.2!) LINC00662(+15.1!) MFSD13A(-12.9!) SLC12A5(+8.7) MED20(-7.9!) COP1(-7.4!) LINC00558(-6.9) PSEN2(+6.7) ZNF852(-6.7) DLST(-6.3) |
| Ext2 | top1 | SLC12A5(+8.5) RAB27B(+8.1) ZBTB5(+7.6) ZCCHC7(-7.5) LINC02796(+7.5) GRIK2(-7.4) ZKSCAN7-AS1(-6.9) LINC00558(-6.8) ASCC3(+6.8) DLST(-6.6) |
| Ext6 | susie | TBCE(-10.0!) COQ8A(-6.7) SDK1(-6.7) SGCD(-5.5) WBP2NL(+5.0) RBM6(-4.8) CCDC32(-4.8) SCRN3(+4.7) LINC02803(-4.7) GRAMD1C(+4.6) |
| Ext6 | top1 | SLC12A5(+8.5) LINC02797(+6.6) COQ8A(-6.5) SDK1(-6.4) FAM120AOS(+6.2) ZWILCH(+6.2) STIMATE(-6.1) SPPL3(-6.0) ASXL3(+5.9) TMEM106B(+5.2) |
| Ext3 | susie | KIF6(+12.8!) DHFR2(+11.9!) ASXL3(+6.1) SGCD(-5.3) LINC01709(+4.9!) SLC7A9(+4.8!) CCDC32(-4.6) GRAMD1C(+4.6) NDUFS3(+4.6) DYNC1I2(+4.5) |
| Ext3 | top1 | ZBTB16(+9.4) DCC(+8.9) SLC12A5(+8.5) SLC25A51(-7.6) PSMG3(+6.7) PHF2(-6.5) ZRANB2-AS2(-6.2) TOMM5(-6.2) ASCC3(+6.1) ASXL3(+6.0) |
| Ext4 | susie | MFSD13A(-23.6!) AGAP9(-17.3!) CEMIP2(+14.3!) MAPT-AS1(-10.5!) DCC(+8.9) SLC12A5(+8.6) RAB27B(+8.3) MED20(-8.1!) ZCCHC7(-7.6) ZKSCAN7(-6.8) |
| Ext4 | top1 | DCC(+9.0) SLC12A5(+8.5) RAB27B(+8.2) ZCCHC7(-7.4) ZKSCAN7(-6.8) ZNF852(-6.6) FURIN(-6.4) NT5DC2(+6.4) TCAIM(-6.3) FAM120AOS(+6.2) |
| Ext5 | susie | SLC12A5(+8.5) AC159540.2(+7.7!) ZCCHC7(-7.6) FAM120AOS(+6.4) ITIH4(-6.0) METTL15(-5.9) STIMATE(-5.7) ZNF502(-5.6) GLT8D1(-5.4) SGCD(-5.0) |
| Ext5 | top1 | SLC12A5(+8.5) ZCCHC7(-7.4) GRIK2(-7.3) ZKSCAN7-AS1(-6.8) FAM120AOS(+6.3) FADS3(+6.3) ADGRA1-AS1(+6.1) ITIH4(-6.0) STIMATE(-5.9) METTL15(-5.9) |
| Ext7 | susie | UBR5-AS1(-33.7!) TUT4(+13.4!) COP1(-10.3!) CYB561(+8.9!) SLC12A5(+8.5) MFSD13A(-8.5!) USP49(-8.0!) CCND3(-7.8!) ZCCHC7(-7.6) FAM120A(+7.0) |
| Ext7 | top1 | SLC12A5(+8.5) RAB27B(+8.2) SETD7(+7.6) FAM120A(+7.0) SLC25A17(-6.9) COQ8A(-6.9) LINC02796(+6.7) ZNF445(-6.4) FAM120AOS(+6.2) PTCH1(+6.1) |
| Ext8 | susie | COP1(+22.3!) ZNF66(+17.9!) MFSD13A(-11.9!) SLC12A5(+8.3) RAB27B(+8.3) ZNF445(-6.8) FURIN(-6.4) FAM120AOS(+6.4) ZNF852(-6.4) FADS1(-6.3) |
| Ext8 | top1 | SLC12A5(+8.5) RAB27B(+8.1) PSEN2(+6.6) STN1(+6.6) ZNF445(-6.5) LINC02797(+6.4) FURIN(-6.4) FAM120AOS(+6.2) FADS1(-6.2) PHF2(-6.1) |

### Top 10 genes: ABCD_thin_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | ODR4(-78.4!) VPS13A(-134.3!) PILRB(-14.3!) GATC(-6.0!) MAPT-AS1(-4.5!) RABGEF1(+4.4!) ADGRA3(-3.3) PIDD1(-3.3) B3GNTL1(-3.2) IQANK1(-3.2) |
| Ext1 | top1 | NDUFS6(-3.4) B3GNTL1(-3.4) RRAGA(-3.4) ADGRA3(-3.3) PIDD1(-3.3) TLE3(-3.3) PRKCE(-3.3) MAP1A(+3.3) DEGS1(+3.3) PITPNB(-3.2) |
| Ext1_n983 | susie | ODR4(-67.5!) VPS13A(-135.0!) PILRB(-14.6!) GATC(-6.3!) MAPT-AS1(-4.5!) RABGEF1(+4.2!) PIGP(+4.2!) ADGRA3(-3.3) PIDD1(-3.3) B3GNTL1(-3.2) |
| Ext1_n983 | top1 | NDUFS6(-3.4) B3GNTL1(-3.4) RRAGA(-3.4) ADGRA3(-3.3) PIDD1(-3.3) TLE3(-3.3) PRKCE(-3.3) MAP1A(+3.3) DEGS1(+3.3) PITPNB(-3.2) |
| Ext | susie | ODR4(-50.2!) VPS13A(-205.0!) LINC01934(-41.5!) LINC01811(+31.8!) DELE1(-20.8!) SNX14(+17.1!) CNOT9(-15.7!) SELENOI(-14.4!) POC1B-AS1(-11.7!) NRDC(-10.9!) |
| Ext | top1 | ASB7(-4.0) FAM71D(-3.8) CACHD1(-3.8) YWHAH(-3.8) TRNP1(+3.7) NPTX1(+3.4) TMOD3(-3.4) SNHG3(+3.4) ADGRA3(-3.3) SOGA3(+3.3) |
| OPC | susie | VPS50(+27.9!) VPS35L(-7.5!) SIDT1(-6.8!) PHLDB1(+4.8!) RAB27A(+3.4) CCDC146(+3.3) NAPB(+3.3) CACHD1(+3.2) MIS18BP1(-3.2) B3GNTL1(-3.1) |
| OPC | top1 | ERI2(-3.8) OGDHL(-3.8) PPP1R8(-3.7) CACHD1(+3.5) RAB11FIP5(+3.5) DTD1(+3.5) CCDC146(+3.4) TTYH1(-3.4) LINC01411(+3.4) MIS18BP1(-3.3) |
| IN | susie | NIPSNAP2(-51.1!) AC009228.1(-29.9!) SELENOI(-26.6!) TUT4(-26.6!) TTC3(-9.8!) DNAAF4(+8.4!) CCDC7(+7.9!) PILRB(-6.9!) LIX1-AS1(-6.4!) UTP4(-6.3!) |
| IN | top1 | PEX2(-3.8) SMG6(-3.6) UBR2(-3.5) ME3(+3.5) NAPB(+3.5) RBBP9(+3.4) SNHG3(+3.4) ACVR1C(-3.4) PRR13(+3.4) BUD13(+3.3) |
| OD | susie | TUT4(-13.5!) LINC01725(+4.8!) SDHA(+4.0!) MMP24OS(+3.9!) MRPL10(+3.8!) NAPB(+3.4) RAB27A(+3.4) SEMA3C(+3.4) PPA1(+3.3) RAB30(+3.2) |
| OD | top1 | FAM71D(-3.8) LPCAT3(-3.6) ZNF12(+3.5) UPF2(+3.5) EXT1(+3.5) DCLRE1A(-3.5) C19orf38(+3.5) CDC26(+3.4) MIIP(-3.3) NAPB(+3.3) |
| Ast | susie | ATP5F1C(+22.2!) PTPN13(+3.7!) NAPB(+3.4) RAB27A(+3.3) B3GNTL1(-3.1) HEATR4(+3.1) ANKRD37(-3.1) TCEA3(-3.1) ZFAS1(-3.0) OSER1(-3.0) |
| Ast | top1 | CTNNBIP1(-4.1) CACHD1(-3.8) PARG(+3.7) JRK(-3.7) RAB18(-3.7) CDC26(+3.7) NAPB(+3.6) SNAPC5(-3.5) TRNAU1AP(-3.5) LRRC47(-3.4) |
| MG | susie | KCNMA1-AS1(+3.5) KCNMA1(+3.4) RAB27A(+3.4) MAP3K20(-3.3) PARP15(-3.1) RAB30-DT(+2.9) PARP14(-2.9) AP1B1(+2.9) CNIH4(-2.9) PDE9A(-2.8) |
| MG | top1 | LBR(-3.9) UTP11(-3.7) CD46(+3.5) HERC6(-3.5) ATRNL1(-3.5) RASGEF1C(-3.4) PACSIN2(+3.4) DSCAML1(+3.3) RAB27A(+3.3) RPL23(+3.2) |
| Ext2 | susie | LINC00662(-19.7!) POC1B-AS1(-14.5!) PILRB(-5.9!) MFSD13A(+5.7!) COP1(-4.9!) C2orf92(-4.8!) ACOT6(+3.5) PRPH2(+3.4) NPTX1(+3.4) RAB27A(+3.3) |
| Ext2 | top1 | TSFM(-3.9) LYRM4(+3.8) SIPA1L2(-3.8) MAP7D1(-3.7) ACOT6(+3.6) ELOVL6(-3.5) ZNF699(-3.5) RBBP9(+3.4) MDM2(+3.4) NPTX1(+3.4) |
| Ext6 | susie | TBCE(-5.4!) PRPH2(+3.5) B3GNTL1(-3.2) DLG2(+3.2) MIS18BP1(-3.1) C2orf72(+3.0) GPR55(+3.0) MIR4713HG(-2.9) ITM2C(+2.8) PDXDC1(-2.8) |
| Ext6 | top1 | BANP(+3.7) EEF1AKMT1(-3.4) PRPH2(+3.4) PABPC4(-3.3) PARD3B(+3.3) ZNF214(+3.3) ATP8A2(-3.3) TMEM41B(+3.3) GPR55(+3.2) RPL32(-3.2) |
| Ext3 | susie | KIF6(-49.8!) LINC02008(+16.1!) DHFR2(-13.6!) DEPDC1B(+9.5!) LINC01709(-4.8!) HEATR4(+3.4) B3GNTL1(-3.3) RAB27A(+3.3) TMOD3(-3.2) PARG(+3.1) |
| Ext3 | top1 | HPS4(-3.8) MIAT(+3.8) CAV2(+3.5) B3GNTL1(-3.4) TAF9(-3.4) HABP4(-3.4) PFKM(+3.3) YEATS4(+3.3) CACNA1D(-3.3) RAB27A(+3.3) |
| Ext4 | susie | ODR4(-41.3!) CEMIP2(-32.3!) MFSD13A(+9.0!) AC159540.2(-5.9!) AGAP9(-5.6!) PILRB(-4.9!) PHACTR4(-4.7!) LMO1(-3.7) TCAIM(-3.4) TMOD3(-3.3) |
| Ext4 | top1 | TCAIM(-3.8) CACHD1(-3.8) RBBP9(+3.5) TMOD3(-3.4) NDUFS6(-3.4) TIMM23B(-3.4) RAB27A(+3.4) BCL2L2-PABPN1(+3.4) PIDD1(-3.3) SMPD1(-3.3) |
| Ext5 | susie | AC159540.2(-52.4!) FAM182B(-3.7!) PRPH2(+3.4) B3GNTL1(-3.3) KCNA1(-3.2) GPR55(+3.2) DAP3(+3.1) SORBS3(-3.0) CABYR(+2.8) CCAR2(-2.8) |
| Ext5 | top1 | PRPH2(+3.7) LRRN1(+3.5) RAB30(-3.5) NTRK2(-3.4) B3GNTL1(-3.4) KCNA1(-3.2) LINC00571(-3.2) GPR55(+3.1) ZNRF3(-3.1) ANKRD37(-3.1) |
| Ext7 | susie | COP1(+38.6!) MRPS30-DT(-16.7!) TUT4(-10.3!) MFSD13A(+8.6!) UBR5-AS1(+3.7!) YWHAH(-3.6) PCNX4(+3.5!) LINC02018(+3.2!) B3GNTL1(-3.2) HACD2(+3.0) |
| Ext7 | top1 | OSBPL6(-4.0) RARS2(-3.8) YWHAH(-3.8) ASB7(-3.8) SH3D21(-3.7) ABCD4(+3.4) RAB27A(+3.3) AK5(-3.2) RDH11(-3.2) C10orf143(+3.2) |
| Ext8 | susie | ZNF66(+22.5!) AC159540.2(-12.4!) ZBTB37(-11.7!) COP1(+6.2!) MFSD13A(+4.1!) KLHL23(-3.8!) CCDC146(+3.3) B3GNTL1(-3.2) LRRN1(+3.1) IQANK1(-3.1) |
| Ext8 | top1 | PAM(-4.2) CACHD1(-3.8) TYRO3(+3.8) CCT2(-3.7) LRPPRC(-3.7) ACOT6(+3.6) ELOVL6(-3.5) ATXN7(-3.5) ZNF780A(+3.4) B3GNTL1(-3.4) |

### Top 10 genes: ABCD_base_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | ODR4(+43.0!) GATC(-43.2!) PRKRA(-4.7!) RMC1(+4.6!) POLE(+4.4!) PILRB(+4.2!) MED20(-4.1!) ACTR1B(-4.1) KLHL23(+3.7!) SRP72(-3.6) |
| Ext1 | top1 | AC159540.2(+4.4) ACTR1B(-4.1) SLC22A25(+4.1) ANKRD36(-4.0) CRAT(-3.6) RPS29(-3.5) SRP72(-3.5) EFCAB13(-3.5) POU6F1(+3.4) COQ4(-3.4) |
| Ext1_n983 | susie | GATC(-45.8!) ODR4(+37.2!) PRKRA(-4.7!) PILRB(+4.6!) ACTR1B(-4.1) POLE(+4.1!) MED20(-4.0!) KLHL23(+3.7!) ANKRD36(-3.6) SRP72(-3.6) |
| Ext1_n983 | top1 | AC159540.2(+4.4) ACTR1B(-4.1) SLC22A25(+4.1) ANKRD36(-4.0) CRAT(-3.6) RPS29(-3.5) SRP72(-3.5) EFCAB13(-3.5) POU6F1(+3.4) COQ4(-3.4) |
| Ext | susie | ODR4(+32.7!) SNORC(-30.5!) GATC(-27.2!) LINC01811(+26.7!) CNOT9(+24.5!) NRDC(-23.4!) MIR2052HG(-22.5!) SELENOI(+19.5!) EEF1AKMT1(+13.1!) NBL1(-12.2!) |
| Ext | top1 | ACTR1B(-4.1) C2orf92(+4.1) AC159540.2(+4.1) ANKRD36(-4.0) EPB41L2(+4.0) BDNF(+3.9) BMPR1A(-3.9) ZNF345(+3.7) COQ4(-3.7) CHML(+3.6) |
| OPC | susie | VPS50(-74.5!) PHLDB1(+6.4!) ACTR1B(-3.9) AZGP1(-3.6) MAST2(+3.5) ANKRD36(-3.4) TRIOBP(+3.3) FAM184B(+3.3) EFCAB13(-3.3) DPYSL5(-3.2) |
| OPC | top1 | LMAN2L(+4.1) ACTR1B(-4.1) CPOX(+4.1) ANKRD36(-4.1) ANKRD36B(+4.0) TRUB2(+3.7) DNAH1(-3.6) DLEU7(-3.6) EFCAB13(-3.5) AGPAT4(+3.5) |
| IN | susie | NIPSNAP2(+83.7!) TUT4(+38.1!) MTHFS(-35.9!) SELENOI(+34.1!) AC009228.1(-17.0!) DHFR2(+14.6!) LINC01933(-12.2!) DENND5B(-12.1!) TTC3(-9.9!) COQ8A(-5.9!) |
| IN | top1 | BMPR1A(-4.4) ACTR1B(-4.1) ANKRD36(-4.1) COPS7B(+3.7) GPM6A(-3.7) LINC01118(+3.6) PXMP2(-3.6) SLC27A4(+3.5) SWI5(-3.5) WASHC2A(+3.5) |
| OD | susie | TUT4(+23.8!) FOLH1(-11.8!) SDHA(-7.2!) LINC01725(+7.2!) MRPL10(-6.6!) WHRN(+4.4!) PRKRA(-4.3!) ACTR1B(-4.1) TRIOBP(+3.6) EPB41L2(+3.6) |
| OD | top1 | VAT1(+4.2) ACTR1B(-4.1) C2orf92(+4.1) EPB41L2(+4.0) TBC1D24(+4.0) TRIOBP(+3.8) TMEM184B(-3.7) MED6(-3.6) HMGN3(+3.4) SGIP1(+3.3) |
| Ast | susie | ATP5F1C(+9.2!) TMC1(+8.9!) NNT-AS1(-7.0!) LINC02021(+6.4!) PRKRA(-4.8!) ACTR1B(-4.1) TRIOBP(+3.6) EPS8L2(+3.5) SRP72(-3.5) SLC27A4(+3.5) |
| Ast | top1 | ACTR1B(-4.1) ANKRD36(-4.0) ANKRD36B(+4.0) TRIOBP(+3.9) COQ4(-3.6) ANKLE2(+3.6) SLC27A4(+3.6) ETFBKMT(+3.5) SRP72(-3.5) URM1(+3.5) |
| MG | susie | ANKRD36(-4.7) UNC5C(-3.6) ANKRD36B(+3.5) TRIOBP(+3.5) DOCK5(+3.5!) DHX8(-3.5) SRP72(-3.4) ESRRB(+3.3) FAM184B(+3.1) ME3(-3.1) |
| MG | top1 | ANKRD36(-4.1) GHDC(-4.1) ANKRD36B(+4.0) CHFR(-3.8) RGS6(+3.7) USE1(-3.6) TOE1(+3.6) SRP72(-3.5) TRIOBP(+3.5) PLCB4(-3.5) |
| Ext2 | susie | COP1(-21.4!) SEM1(+10.0!) POC1B-AS1(+5.8!) PILRB(-5.7!) LINC00662(-5.2!) MED20(-4.2!) ACTR1B(-4.1) ACTR3C(+3.8) RIDA(-3.6!) AFDN(+3.6!) |
| Ext2 | top1 | ACTR1B(-4.1) ANKRD36(-4.1) ZNF540(+3.8) SWI5(-3.8) STX1A(-3.6) RPS29(-3.5) SRP72(-3.5) SDHAF1(-3.5) EFCAB13(-3.5) BUD23(-3.4) |
| Ext6 | susie | TBCE(+16.2!) ACTR1B(-3.9) ANKRD36B(+3.4) EFCAB13(-3.4) MLXIPL(-3.4) WASHC2A(+3.4) ANKRD36(-3.4) FAM184B(+3.3) ZNF721(+3.1) SUCLG2-AS1(-3.1) |
| Ext6 | top1 | AC159540.2(+4.4) ACTR1B(-4.1) ANKRD36B(+4.1) ANKRD36(-4.0) IMPDH1(+3.8) LSM6(-3.6) WASHC2A(+3.5) EFCAB13(-3.5) CCNF(-3.5) BUD13(-3.5) |
| Ext3 | susie | DHFR2(-32.0!) LINC02008(+8.1!) DEPDC1B(+4.3!) SLC7A9(+4.1!) ACTR1B(-4.1) EFCAB13(-3.4) LINC01954(-3.4!) KIF6(+3.3!) FAM184B(+3.2) POLE(+3.1) |
| Ext3 | top1 | LINC01118(+4.6) ACTR1B(-4.1) ANKRD36(-4.0) ANKRD36B(+4.0) CACNA1D(-3.9) GOLGA3(-3.9) CRIPT(-3.9) NAA35(-3.7) EXOSC1(-3.4) RPRML(+3.3) |
| Ext4 | susie | ODR4(+21.0!) AGAP9(+9.6!) MFSD13A(+6.0!) PRKRA(-4.7!) PHACTR4(+4.6!) CEMIP2(+4.4!) ACTR1B(-4.1) ANKRD36(-3.5) SLC27A4(+3.5) SLU7(+3.4) |
| Ext4 | top1 | PRKCE(-4.6) ACTR1B(-4.1) ANKRD36(-4.0) SPESP1(-3.8) PLPP5(+3.8) LAP3(-3.6) LINC01118(+3.6) PCYOX1L(-3.6) CHML(+3.6) SLC27A4(+3.6) |
| Ext5 | susie | AC159540.2(-22.6!) ANKRD36(-3.6) SMPDL3B(+3.5) ANKRD36B(+3.4) XKR8(+3.4) WASHC2A(+3.3) SETD1A(-3.2) USP13(-3.1) POLE(+3.0) FAM184B(+2.9) |
| Ext5 | top1 | DAPK1(-4.6) ACTR1B(-4.1) ANKRD36(-4.1) HECTD3(-4.0) ANKRD36B(+4.0) ZNF569(+3.8) DCBLD1(+3.7) C1orf216(-3.6) ZNF2(-3.5) MLXIPL(-3.4) |
| Ext7 | susie | COP1(+45.6!) UBR5-AS1(+29.3!) TUT4(+18.7!) VPS13A(-11.6!) PILRB(+6.2!) PCNX4(+4.2!) ACTR1B(-3.9) MRPS30-DT(-3.7!) LINC01954(-3.6!) ANKRD36(-3.6) |
| Ext7 | top1 | AC159540.2(+4.4) ACTR1B(-4.1) ANKRD36(-4.0) ANKRD36B(+4.0) TTLL1(+3.9) EIF2B3(+3.8) TESK2(+3.6) SRP72(-3.5) PIK3R3(-3.5) TMEM131(-3.5) |
| Ext8 | susie | ZBTB37(+19.2!) AC159540.2(-10.3!) GATC(-7.9!) COP1(-6.5!) ZNF66(-5.3!) KLHL23(-5.2!) PRKRA(-4.5!) ACTR1B(-4.1) ANKRD36(-3.4) ZNF721(+3.3) |
| Ext8 | top1 | AC159540.2(+4.4) ACTR1B(-4.1) ANKRD36(-4.0) RPL27(+3.9) IFT140(-3.9) MIER2(-3.7) RPSA(-3.5) SLC9A1(+3.5) NUDT21(+3.5) MRPS7(+3.4) |


<!-- /AUTO:RESULTS -->

## Positive-control notes

- **CNTN4, Ext:** top1 z = +7.20 (p 6e-13), in the expected direction (risk allele raises expression). The susie
  model fits no credible set for CNTN4 in Ext, so it is absent from `assoc/*.Ext.susie`.
- **SLC4A10, Ext and Ext2:** top1 z = −5.13 (Ext, p 3e-7) and −5.49 (Ext2, p 4e-8). These are strong but opposite
  to the paper's direction (protective rs4572609-G lowers expression). The paper's lead and fine-mapped SNP rs4572609
  is C/G, so the required strand-ambiguity filter removed it. top1 is therefore tagging a different SNP whose LD phase
  with rs4572609 decides the sign. SuSiE fits no credible set for SLC4A10 in either cell type.
- **CNTN4, Ext2:** top1 z = +6.10 (p 1e-9), the same direction as in Ext.
- **Fig. 3c gene × its colocalising cell type(s), SCZ25:** subtype rows count when the subtype or its parent class colocalises.
  - susie: 44 rows are modelled, 34 have p < 0.05, and 23 have z > 0.
  - top1: 110 rows, 68 have p < 0.05, and 61 have z > 0.
  - The TWAS sign need not match the coloc direction when the lead SNP was filtered out (see SLC4A10 above).

## Where the ExN outputs are (CSD3)

Base: `/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/`

| tag | cell type | status | files |
|---|---|---|---|
| Ext1 | L2/3 IT (priority) | done | `assoc/<trait>.Ext1.{susie,top1}.tsv.gz`, `weights/Ext1.{susie,top1,genes}.tsv.gz` |
| Ext1_n983 | L2/3 IT, n = 983 sensitivity | done | `assoc/<trait>.Ext1_n983.*`, `weights/Ext1_n983.*` |
| Ext | all excitatory | done | `assoc/<trait>.Ext.*`, `weights/Ext.*` |
| Ext2 | L4 IT | done | `assoc/<trait>.Ext2.*`, `weights/Ext2.*` |
| Ext3 | L5/6 NP | done | `assoc/<trait>.Ext3.*`, `weights/Ext3.*` |
| Ext4 | L5 IT | done | `assoc/<trait>.Ext4.*`, `weights/Ext4.*` |
| Ext5 | L6b | done | `assoc/<trait>.Ext5.*`, `weights/Ext5.*` |
| Ext6 | L6 CT | done | `assoc/<trait>.Ext6.*`, `weights/Ext6.*` |
| Ext7 | L6 IT Car3 | done | `assoc/<trait>.Ext7.*`, `weights/Ext7.*` |
| Ext8 | L6 IT | done | `assoc/<trait>.Ext8.*`, `weights/Ext8.*` |

`<trait>` is one of `SCZ25_EUR`, `PGC3_EUR`, `MDD_EUR`, `ABCD_thin_EUR` and `ABCD_base_EUR`.

## Deliverables / rsync

Run this on the laptop. Weights contain SNP lists only; no LD matrices are stored.

ExN only (Ext, Ext1, Ext1_n983, Ext2–Ext8):

```bash
rsync -av --include='*.Ext*' --exclude='*.tsv.gz' rajd2@login-icelake.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/assoc/ ~/Git/transcriptional_maturation/gwas_cells/work/sctwas/assoc/
rsync -av --include='Ext*' --exclude='*' rajd2@login-icelake.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/weights/ ~/Git/transcriptional_maturation/gwas_cells/work/sctwas/weights/
```

Everything:

```bash
rsync -av --relative rajd2@login-icelake.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/./{assoc,weights,report,README.md} ~/Git/transcriptional_maturation/gwas_cells/work/sctwas/
```
