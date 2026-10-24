# Cell-type-specific TWAS from SingleBrain sc-eQTLs (SuSiE-RSS weights)

Task spec: `transcriptional_maturation/gwas_cells/csd3/SCTWAS_INSTRUCTIONS.md` (commit 2f1aceb).
Scripts: `abcd_development/genetic_analysis/sctwas/` (committed). Work dir: this directory (`$W`).
Run on CSD3 starting 2026-10-24.

## Status (2026-10-24, interim push)

- **Complete (weights + assoc):** Ext1 (L2/3 IT), Ext1_n983, Ext, OPC, IN, OD. Tables below cover only the tags whose
  assoc files existed when `04_report.py` last ran.
- **Still running** on SL3 (fits or chained assoc): Ast, MG, Ext2, Ext3, Ext4, Ext5, Ext6, Ext7, Ext8. Each finished
  tag gets its files in `assoc/` and `weights/` automatically. Rerun `04_report.py` to refresh the tables.
- **Part B (PsychAD) not started.** It needs a Synapse account and token from the user.

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
../../../legacy/hpc/work/envs/abcd/bin/python ../../sctwas/04_report.py   # report/*.tsv + tables below
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
| OPC | OPC | 12029 | 9335 | 9332 | 0 | 5 | 2811 |  | 2811 | 9332 | 12 | 3.14e+03 | 2113 | 513 | 185 | 0.0685 | 0.144 | 0.00182 | 0.00462 |
| IN | inhibitory neurons | 12226 | 11165 | 11162 | 0 | 46 | 4938 | 51 | 4938 | 11162 | 13 | 3.12e+03 | 3468 | 1044 | 426 | 0.0722 | 0.151 | 0.00188 | 0.00324 |
| OD | oligodendrocytes | 11505 | 9898 | 9896 | 0 | 18 | 4157 | 22 | 4157 | 9896 | 13 | 3.08e+03 | 2793 | 888 | 476 | 0.068 | 0.147 | 0.00303 | 0.0065 |

### Per trait x cell type x model

| trait | tag | model | genes | lambda_gc | n_p_lt_2_5e6 | n_fdr05 | n_w_flag | lambda_gc_unflagged | n_p_lt_2_5e6_unflagged |
|---|---|---|---|---|---|---|---|---|---|
| SCZ25_EUR | Ext1 | susie | 4322 | 2.57 | 82 | 462 | 17 | 2.56 | 75 |
| SCZ25_EUR | Ext1 | top1 | 11116 | 2 | 121 | 735 | 0 | 2 | 121 |
| SCZ25_EUR | Ext1_n983 | susie | 4346 | 2.61 | 82 | 466 | 23 | 2.58 | 73 |
| SCZ25_EUR | Ext1_n983 | top1 | 11116 | 2 | 121 | 735 | 0 | 2 | 121 |
| SCZ25_EUR | Ext | susie | 7115 | 2.49 | 139 | 754 | 102 | 2.47 | 116 |
| SCZ25_EUR | Ext | top1 | 11626 | 2.18 | 164 | 887 | 0 | 2.18 | 164 |
| SCZ25_EUR | OPC | susie | 2811 | 2.25 | 61 | 261 | 0 | 2.25 | 61 |
| SCZ25_EUR | OPC | top1 | 9332 | 1.88 | 104 | 524 | 0 | 1.88 | 104 |
| SCZ25_EUR | IN | susie | 4938 | 2.54 | 100 | 488 | 51 | 2.51 | 91 |
| SCZ25_EUR | IN | top1 | 11162 | 2.07 | 136 | 726 | 0 | 2.07 | 136 |
| SCZ25_EUR | OD | susie | 4157 | 2.33 | 80 | 408 | 22 | 2.33 | 75 |
| SCZ25_EUR | OD | top1 | 9896 | 2.02 | 131 | 641 | 0 | 2.02 | 131 |
| PGC3_EUR | Ext1 | susie | 4322 | 2.54 | 59 | 395 | 17 | 2.54 | 54 |
| PGC3_EUR | Ext1 | top1 | 11014 | 1.94 | 87 | 600 | 0 | 1.94 | 87 |
| PGC3_EUR | Ext1_n983 | susie | 4346 | 2.54 | 62 | 404 | 23 | 2.54 | 54 |
| PGC3_EUR | Ext1_n983 | top1 | 11014 | 1.94 | 87 | 600 | 0 | 1.94 | 87 |
| PGC3_EUR | Ext | susie | 7115 | 2.42 | 112 | 631 | 102 | 2.4 | 91 |
| PGC3_EUR | Ext | top1 | 11514 | 2.1 | 118 | 698 | 0 | 2.1 | 118 |
| PGC3_EUR | OPC | susie | 2811 | 2.2 | 45 | 242 | 0 | 2.2 | 45 |
| PGC3_EUR | OPC | top1 | 9238 | 1.84 | 84 | 441 | 0 | 1.84 | 84 |
| PGC3_EUR | IN | susie | 4936 | 2.43 | 75 | 422 | 51 | 2.39 | 69 |
| PGC3_EUR | IN | top1 | 11073 | 1.97 | 107 | 626 | 0 | 1.97 | 107 |
| PGC3_EUR | OD | susie | 4155 | 2.29 | 66 | 347 | 22 | 2.29 | 62 |
| PGC3_EUR | OD | top1 | 9804 | 1.91 | 103 | 553 | 0 | 1.91 | 103 |
| MDD_EUR | Ext1 | susie | 4322 | 2.27 | 49 | 335 | 17 | 2.25 | 41 |
| MDD_EUR | Ext1 | top1 | 11116 | 1.93 | 80 | 517 | 0 | 1.93 | 80 |
| MDD_EUR | Ext1_n983 | susie | 4346 | 2.25 | 51 | 329 | 23 | 2.23 | 40 |
| MDD_EUR | Ext1_n983 | top1 | 11116 | 1.93 | 80 | 517 | 0 | 1.93 | 80 |
| MDD_EUR | Ext | susie | 7115 | 2.3 | 90 | 604 | 102 | 2.25 | 60 |
| MDD_EUR | Ext | top1 | 11626 | 2.06 | 90 | 672 | 0 | 2.06 | 90 |
| MDD_EUR | OPC | susie | 2811 | 2.38 | 22 | 199 | 0 | 2.38 | 22 |
| MDD_EUR | OPC | top1 | 9332 | 1.83 | 48 | 376 | 0 | 1.83 | 48 |
| MDD_EUR | IN | susie | 4938 | 2.23 | 60 | 367 | 51 | 2.22 | 43 |
| MDD_EUR | IN | top1 | 11162 | 1.94 | 72 | 564 | 0 | 1.94 | 72 |
| MDD_EUR | OD | susie | 4157 | 2.2 | 31 | 273 | 22 | 2.2 | 25 |
| MDD_EUR | OD | top1 | 9896 | 1.87 | 49 | 421 | 0 | 1.87 | 49 |
| ABCD_thin_EUR | Ext1 | susie | 4322 | 0.987 | 4 | 6 | 17 | 0.975 | 0 |
| ABCD_thin_EUR | Ext1 | top1 | 11116 | 0.962 | 0 | 0 | 0 | 0.962 | 0 |
| ABCD_thin_EUR | Ext1_n983 | susie | 4346 | 0.987 | 4 | 7 | 23 | 0.976 | 0 |
| ABCD_thin_EUR | Ext1_n983 | top1 | 11116 | 0.962 | 0 | 0 | 0 | 0.962 | 0 |
| ABCD_thin_EUR | Ext | susie | 7115 | 1 | 21 | 25 | 102 | 0.994 | 0 |
| ABCD_thin_EUR | Ext | top1 | 11626 | 0.985 | 0 | 0 | 0 | 0.985 | 0 |
| ABCD_thin_EUR | OPC | susie | 2811 | 0.993 | 4 | 4 | 0 | 0.993 | 4 |
| ABCD_thin_EUR | OPC | top1 | 9332 | 0.995 | 0 | 0 | 0 | 0.995 | 0 |
| ABCD_thin_EUR | IN | susie | 4938 | 1.08 | 13 | 13 | 51 | 1.07 | 0 |
| ABCD_thin_EUR | IN | top1 | 11162 | 1.03 | 0 | 0 | 0 | 1.03 | 0 |
| ABCD_thin_EUR | OD | susie | 4157 | 1.05 | 2 | 2 | 22 | 1.05 | 0 |
| ABCD_thin_EUR | OD | top1 | 9896 | 1.02 | 0 | 0 | 0 | 1.02 | 0 |
| ABCD_base_EUR | Ext1 | susie | 4322 | 1.17 | 2 | 8 | 17 | 1.16 | 0 |
| ABCD_base_EUR | Ext1 | top1 | 11116 | 1.09 | 0 | 0 | 0 | 1.09 | 0 |
| ABCD_base_EUR | Ext1_n983 | susie | 4346 | 1.17 | 2 | 7 | 23 | 1.16 | 0 |
| ABCD_base_EUR | Ext1_n983 | top1 | 11116 | 1.09 | 0 | 0 | 0 | 1.09 | 0 |
| ABCD_base_EUR | Ext | susie | 7115 | 1.15 | 24 | 28 | 102 | 1.14 | 0 |
| ABCD_base_EUR | Ext | top1 | 11626 | 1.09 | 0 | 0 | 0 | 1.09 | 0 |
| ABCD_base_EUR | OPC | susie | 2811 | 1.07 | 2 | 2 | 0 | 1.07 | 2 |
| ABCD_base_EUR | OPC | top1 | 9332 | 1.07 | 0 | 0 | 0 | 1.07 | 0 |
| ABCD_base_EUR | IN | susie | 4938 | 1.1 | 14 | 20 | 51 | 1.08 | 0 |
| ABCD_base_EUR | IN | top1 | 11162 | 1.06 | 0 | 0 | 0 | 1.06 | 0 |
| ABCD_base_EUR | OD | susie | 4157 | 1.1 | 5 | 8 | 22 | 1.09 | 0 |
| ABCD_base_EUR | OD | top1 | 9896 | 1.07 | 0 | 0 | 0 | 1.07 | 0 |

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
| PTPRU | Ext |  | Ext1 | susie | -6.08 | 1.18e-09 | False |
| PTPRU | Ext |  | Ext1 | top1 | -4.1 | 4.14e-05 | False |
| PTPRU | Ext |  | Ext1_n983 | susie | -6.09 | 1.14e-09 | False |
| PTPRU | Ext |  | Ext1_n983 | top1 | -4.1 | 4.14e-05 | False |
| PTPRU | Ext |  | Ext | susie | -4.2 | 2.62e-05 | False |
| PTPRU | Ext |  | Ext | top1 | -4.6 | 4.14e-06 | False |
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
| FOXP1 | OD | OD1,OD2,OD3 | OD | susie | -6.08 | 1.21e-09 | False |
| FOXP1 | OD | OD1,OD2,OD3 | OD | top1 | -6.26 | 3.95e-10 | False |
| SPCS1 | IN |  | IN | susie |  |  | nan |
| SPCS1 | IN |  | IN | top1 | -9.53 | 1.61e-21 | False |
| AGO4 | OD |  | OD | susie | 5.65 | 1.63e-08 | False |
| AGO4 | OD |  | OD | top1 | 5.43 | 5.79e-08 | False |
| TMX1 | OPC | OPC1 | OPC | susie | 4.72 | 2.34e-06 | False |
| TMX1 | OPC | OPC1 | OPC | top1 | 4.1 | 4.05e-05 | False |
| SLC4A10 | Ext | Ext2 | Ext1 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext1 | top1 | -1.16 | 0.247 | False |
| SLC4A10 | Ext | Ext2 | Ext1_n983 | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext1_n983 | top1 | -1.16 | 0.247 | False |
| SLC4A10 | Ext | Ext2 | Ext | susie |  |  | nan |
| SLC4A10 | Ext | Ext2 | Ext | top1 | -5.13 | 2.95e-07 | False |
| CNTN4 | Ext |  | Ext1 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext1 | top1 | 1.09 | 0.276 | False |
| CNTN4 | Ext |  | Ext1_n983 | susie |  |  | nan |
| CNTN4 | Ext |  | Ext1_n983 | top1 | 1.09 | 0.276 | False |
| CNTN4 | Ext |  | Ext | susie |  |  | nan |
| CNTN4 | Ext |  | Ext | top1 | 7.2 | 6.02e-13 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext1 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext1 | top1 | 5.62 | 1.96e-08 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext1_n983 | susie |  |  | nan |
| PLCH2 | Ext | Ext2,IN7 | Ext1_n983 | top1 | 5.62 | 1.96e-08 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext | susie | 6.08 | 1.21e-09 | False |
| PLCH2 | Ext | Ext2,IN7 | Ext | top1 | 5.74 | 9.65e-09 | False |
| KCNB1 | IN | IN7 | IN | susie | 6.13 | 8.90e-10 | False |
| KCNB1 | IN | IN7 | IN | top1 | 6.02 | 1.74e-09 | False |
| GNL3 | IN | Ext4,IN4,IN6,IN7 | IN | susie | 7.41 | 1.27e-13 | False |
| GNL3 | IN | Ext4,IN4,IN6,IN7 | IN | top1 | 8.17 | 2.97e-16 | False |
| CACNB2 | Ext |  | Ext1 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext1 | top1 | 4.17 | 3.05e-05 | False |
| CACNB2 | Ext |  | Ext1_n983 | susie |  |  | nan |
| CACNB2 | Ext |  | Ext1_n983 | top1 | 4.17 | 3.05e-05 | False |
| CACNB2 | Ext |  | Ext | susie |  |  | nan |
| CACNB2 | Ext |  | Ext | top1 | 6.44 | 1.21e-10 | False |
| WBP1L | IN | Ext2,Ext7 | IN | susie | 9.21 | 3.24e-20 | False |
| WBP1L | IN | Ext2,Ext7 | IN | top1 | 9.03 | 1.65e-19 | False |
| FUT9 | IN | Ast4,IN6,IN7,OD1,OD2 | IN | susie | 5.57 | 2.59e-08 | False |
| FUT9 | IN | Ast4,IN6,IN7,OD1,OD2 | IN | top1 | 5.76 | 8.34e-09 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1 | susie | -4.63 | 3.65e-06 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1 | top1 | -5.4 | 6.60e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1_n983 | susie | -4.66 | 3.13e-06 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext1_n983 | top1 | -5.4 | 6.60e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext | susie | -5.35 | 8.78e-08 | False |
| SLCO4C1 | Ext | Ext1,Ext4,Ext5,Ext6,Ext7,Ext8 | Ext | top1 | -5.4 | 6.60e-08 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext1 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext1 | top1 | 5.15 | 2.59e-07 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext1_n983 | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext1_n983 | top1 | 5.15 | 2.59e-07 | False |
| SCRIB | Ext | Ext1,Ext2 | Ext | susie |  |  | nan |
| SCRIB | Ext | Ext1,Ext2 | Ext | top1 | 5.15 | 2.59e-07 | False |
| RGS6 | Ext |  | Ext1 | susie |  |  | nan |
| RGS6 | Ext |  | Ext1 | top1 | -4.86 | 1.15e-06 | False |
| RGS6 | Ext |  | Ext1_n983 | susie |  |  | nan |
| RGS6 | Ext |  | Ext1_n983 | top1 | -4.86 | 1.15e-06 | False |
| RGS6 | Ext |  | Ext | susie | -8.95 | 3.48e-19 | False |
| RGS6 | Ext |  | Ext | top1 | -7.01 | 2.36e-12 | False |
| CSMD1 | OPC |  | OPC | susie | 0.499 | 0.618 | False |
| CSMD1 | OPC |  | OPC | top1 | 0.617 | 0.537 | False |
| ANKS1B | OPC | OPC1 | OPC | susie | 5.18 | 2.22e-07 | False |
| ANKS1B | OPC | OPC1 | OPC | top1 | 5.25 | 1.55e-07 | False |
| MYT1L | OPC | Ext5,OPC1 | OPC | susie | 5.91 | 3.38e-09 | False |
| MYT1L | OPC | Ext5,OPC1 | OPC | top1 | 5.76 | 8.27e-09 | False |
| DPYD | OD | Ast4 | OD | susie | 8.79 | 1.49e-18 | False |
| DPYD | OD | Ast4 | OD | top1 | 8.78 | 1.65e-18 | False |
| PTPRT | OPC |  | OPC | susie | 4.9 | 9.81e-07 | False |
| PTPRT | OPC |  | OPC | top1 | 4.89 | 9.99e-07 | False |
| CNOT7 | IN | Ext4 | IN | susie | -1.49 | 0.135 | False |
| CNOT7 | IN | Ext4 | IN | top1 | -1.74 | 0.0816 | False |
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
| CNTN4 | IN | 0.0336 | 0.27 | -2.13 | -1.1 |
| CNTN4 | OD |  | 0.458 |  | 0.742 |
| CNTN4 | OPC | 0.0571 | 0.43 | 1.9 | -0.79 |
| SLC4A10 | Ext |  | 2.95e-07 |  | -5.13 |
| SLC4A10 | Ext1 |  | 0.247 |  | -1.16 |
| SLC4A10 | Ext1_n983 |  | 0.247 |  | -1.16 |
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
| OPC | susie | VPS50(-67.9) KMT2E(+8.1) INO80E(+7.7) RBM26(+7.7) TMEM219(+7.5) AS3MT(+7.1) HSPD1(+7.0) FTCDNL1(-7.0) YJEFN3(-7.0) SFMBT1(-6.8) |
| OPC | top1 | NT5DC2(+8.5) KMT2E(+8.3) INO80E(+7.8) TMEM219(+7.8) RBM26(+7.6) ARL3(-7.5) YJEFN3(-7.3) ZSCAN2(-7.3) AS3MT(+7.1) WDR73(+7.0) |
| IN | susie | TUT4(+87.1!) SELENOI(+24.0!) ANKRD36(-23.0!) AFDN(-20.4!) TTC3(-16.3!) MTHFS(-12.4!) DENND5B(-12.1!) WBP1L(+9.2) MAPT-AS1(+9.0) FURIN(-8.2) |
| IN | top1 | CACNA1C(-10.2) SPCS1(-9.5) WBP1L(+9.0) GATAD2A(-8.6) GNL3(+8.2) ARL3(-7.8) TAOK2(+7.7) FURIN(-7.7) INO80E(+7.7) DOC2A(+7.7) |
| OD | susie | TUT4(-14.7!) FOLH1(+10.0!) DPYD(+8.8) TOM1L1(+8.7!) GATAD2A(-8.4) KCTD13(-7.9) INO80E(+7.8) ASPHD1(-7.8) KMT2E(+7.6) ABCB9(+7.2) |
| OD | top1 | DPYD(+8.8) DPYD-AS1(+8.8) IPO11(+8.6) SRPK2(-8.3) ASPHD1(-7.9) KCTD13(-7.8) SEZ6L2(-7.8) FURIN(-7.7) INO80E(+7.7) GIGYF2(+7.6) |

### Top 10 genes: PGC3_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | MAPT-AS1(+18.6!) SPCS1(-8.7) FURIN(-8.5) YJEFN3(-8.0) GATAD2A(-8.0) INO80E(+7.6) TAOK2(+7.3) PHF21A(-7.2) NEK4(+7.1) LINC00634(+7.1) |
| Ext1 | top1 | FURIN(-8.6) YJEFN3(-8.2) RFT1(+8.0) GATAD2A(-7.8) INPP5D(-7.8) ATG13(-7.7) INO80E(+7.5) PPP1R13B(-7.3) AS3MT(+7.2) FTCDNL1(-7.2) |
| Ext1_n983 | susie | MAPT-AS1(+18.7!) SPCS1(-8.7) FURIN(-8.5) YJEFN3(-8.0) GATAD2A(-8.0) PIGP(+7.7!) INO80E(+7.6) TAOK2(+7.3) PHF21A(-7.2) NEK4(+7.1) |
| Ext1_n983 | top1 | FURIN(-8.6) YJEFN3(-8.2) RFT1(+8.0) GATAD2A(-7.8) INPP5D(-7.8) ATG13(-7.7) INO80E(+7.5) PPP1R13B(-7.3) AS3MT(+7.2) FTCDNL1(-7.2) |
| Ext | susie | NRDC(+30.9!) PPP4R3B(-28.0!) CPLANE1(+16.7!) EEF1AKMT1(-15.9!) AGAP9(+15.5!) MFSD13A(+14.0!) NPEPL1(-11.8!) MAPT-AS1(+11.8!) MIR2052HG(-11.2!) SNX14(+9.0!) |
| Ext | top1 | FURIN(-8.4) GATAD2A(-8.2) ATG16L1(-8.1) RFT1(+8.0) YJEFN3(-7.8) INO80E(+7.5) FTCDNL1(-7.4) AS3MT(+7.2) NEK4(+7.1) WBP2NL(+7.1) |
| OPC | susie | VPS50(-45.4) KMT2E(+7.7) INO80E(+7.6) FTCDNL1(-7.5) TMEM219(+7.4) AS3MT(+7.2) RBM26(+6.8) HSPD1(+6.7) YJEFN3(-6.5) PCCB(+6.3) |
| OPC | top1 | KMT2E(+7.7) INO80E(+7.6) TMEM219(+7.5) NT5DC2(+7.4) FTCDNL1(-7.4) AS3MT(+7.2) YJEFN3(-6.9) RBM26(+6.8) ZSCAN2(-6.8) SMDT1(+6.6) |
| IN | susie | TUT4(-24.3!) SELENOI(+13.9!) AFDN(-12.9!) MTHFS(-10.4!) WBP1L(+8.2) FURIN(-8.0) INO80E(+7.5) KMT2E(+7.4) PPP1R13B(-7.3) GNL3(+7.0) |
| IN | top1 | SPCS1(-9.1) CACNA1C(-8.8) GATAD2A(-8.2) WBP1L(+8.1) GNL3(+8.0) TAOK2(+7.5) FTCDNL1(-7.5) DOC2A(+7.5) INO80E(+7.5) FURIN(-7.3) |
| OD | susie | TUT4(-9.9!) DPYD(+8.8) GATAD2A(-7.9) TOM1L1(+7.7!) INO80E(+7.6) KCTD13(-7.4) KMT2E(+7.4) ASPHD1(-7.3) AS3MT(+7.2) HSPD1(+6.8) |
| OD | top1 | DPYD(+8.8) DPYD-AS1(+8.8) GIGYF2(+8.1) IPO11(+8.1) ASPHD1(-7.7) SRPK2(-7.6) INO80E(+7.5) KCTD13(-7.5) SEZ6L2(-7.5) FTCDNL1(-7.4) |

### Top 10 genes: MDD_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | VPS13A(-94.3!) POC1B-AS1(+41.8!) MAPT-AS1(-23.5!) POLE(-17.8!) SLC12A5(+8.5) RAB27B(+8.3) MED20(-7.8!) ODR4(+7.6!) ZCCHC7(-7.6) MED19(-7.3) |
| Ext1 | top1 | SLC12A5(+8.5) RAB27B(+8.1) ZCCHC7(-7.4) MED19(-7.3) SREBF2(-7.3) SLC25A17(-6.7) ZNF445(-6.6) ZKSCAN7-AS1(-6.6) PSEN2(+6.5) ZNF852(-6.5) |
| Ext1_n983 | susie | VPS13A(-94.8!) POC1B-AS1(+43.3!) MAPT-AS1(-23.6!) POLE(-20.8!) SLC12A5(+8.7) RAB27B(+8.3) MED20(-8.2!) CCDC32(-7.7!) ZCCHC7(-7.5) MED19(-7.3) |
| Ext1_n983 | top1 | SLC12A5(+8.5) RAB27B(+8.1) ZCCHC7(-7.4) MED19(-7.3) SREBF2(-7.3) SLC25A17(-6.7) ZNF445(-6.6) ZKSCAN7-AS1(-6.6) PSEN2(+6.5) ZNF852(-6.5) |
| Ext | susie | VPS13A(-144.5!) LINC01934(+44.8!) NBL1(+38.1!) MFSD13A(-36.2!) MIR2052HG(+27.7!) EEF1AKMT1(+26.6!) SLC17A5(+20.5!) POC1B-AS1(+19.9!) HDGFL3(+17.3!) ANAPC5(-16.0!) |
| Ext | top1 | DCC(+8.9) SLC12A5(+8.5) RAB27B(+8.2) SLC25A51(-7.5) BDNF(+7.4) GRIK3(+7.0) ZNF197(-6.9) ASXL3(+6.9) ZKSCAN7-AS1(-6.9) ZKSCAN7(-6.6) |
| OPC | susie | PHLDB1(+16.4) SIDT1(+11.0) NEGR1(+8.1) ASCC3(+7.5) CCDC73(-6.9) FAM120AOS(+6.7) BAG5(-6.1) SFMBT1(-5.8) SUDS3(-5.7) CYP7B1(-5.6) |
| OPC | top1 | SLC25A51(-7.9) NEGR1(+7.7) EIF3M(-7.3) CCDC73(-7.0) SLC25A17(+6.8) ASCC3(+6.8) CD59(-6.7) NDUFS3(-6.7) FAM120AOS(+6.6) NT5DC2(+6.4) |
| IN | susie | TUT4(-54.5!) POC1B-AS1(+41.7!) DENND5B(+37.3!) MTHFS(+32.1!) NIPSNAP2(-30.0!) TTC3(+25.4!) SELENOI(-19.3!) AFDN(+14.2!) AC009228.1(+11.8!) MFSD13A(-11.5!) |
| IN | top1 | SLC12A5(+8.5) RAB27B(+8.2) ST13(-7.3) SLC25A17(+7.2) KLC1(+7.0) ZKSCAN7(-6.9) ASCC3(+6.9) ZNF197(-6.7) ZNF852(-6.7) NEGR1(+6.6) |
| OD | susie | LINC01725(+15.0!) SDHA(+10.2!) MED20(-8.6!) VRK2(-7.8) ZCCHC7(-7.6) ASCC3(+7.4) MMP24OS(+6.8!) NEGR1(+6.7) FAM120AOS(+6.4) YLPM1(+6.3) |
| OD | top1 | ZCCHC7(-7.6) VRK2(-7.2) ASCC3(+7.0) ZKSCAN7(-6.8) ZKSCAN7-AS1(-6.4) FAM120AOS(+6.3) FADS1(-6.3) ECHS1(+6.1) FADS2(-6.1) YLPM1(+6.1) |

### Top 10 genes: ABCD_thin_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | ODR4(-78.4!) VPS13A(-134.3!) PILRB(-14.3!) GATC(-6.0!) MAPT-AS1(-4.5!) RABGEF1(+4.4!) ADGRA3(-3.3) PIDD1(-3.3) B3GNTL1(-3.2) IQANK1(-3.2) |
| Ext1 | top1 | NDUFS6(-3.4) B3GNTL1(-3.4) RRAGA(-3.4) ADGRA3(-3.3) PIDD1(-3.3) TLE3(-3.3) PRKCE(-3.3) MAP1A(+3.3) DEGS1(+3.3) PITPNB(-3.2) |
| Ext1_n983 | susie | ODR4(-67.5!) VPS13A(-135.0!) PILRB(-14.6!) GATC(-6.3!) MAPT-AS1(-4.5!) RABGEF1(+4.2!) PIGP(+4.2!) ADGRA3(-3.3) PIDD1(-3.3) B3GNTL1(-3.2) |
| Ext1_n983 | top1 | NDUFS6(-3.4) B3GNTL1(-3.4) RRAGA(-3.4) ADGRA3(-3.3) PIDD1(-3.3) TLE3(-3.3) PRKCE(-3.3) MAP1A(+3.3) DEGS1(+3.3) PITPNB(-3.2) |
| Ext | susie | ODR4(-50.2!) VPS13A(-205.0!) LINC01934(-41.5!) LINC01811(+31.8!) DELE1(-20.8!) SNX14(+17.1!) CNOT9(-15.7!) SELENOI(-14.4!) POC1B-AS1(-11.7!) NRDC(-10.9!) |
| Ext | top1 | ASB7(-4.0) FAM71D(-3.8) CACHD1(-3.8) YWHAH(-3.8) TRNP1(+3.7) NPTX1(+3.4) TMOD3(-3.4) SNHG3(+3.4) ADGRA3(-3.3) SOGA3(+3.3) |
| OPC | susie | VPS50(+27.9) VPS35L(-7.5) SIDT1(-6.8) PHLDB1(+4.8) RAB27A(+3.4) CCDC146(+3.3) NAPB(+3.3) CACHD1(+3.2) MIS18BP1(-3.2) B3GNTL1(-3.1) |
| OPC | top1 | ERI2(-3.8) OGDHL(-3.8) PPP1R8(-3.7) CACHD1(+3.5) RAB11FIP5(+3.5) DTD1(+3.5) CCDC146(+3.4) TTYH1(-3.4) LINC01411(+3.4) MIS18BP1(-3.3) |
| IN | susie | NIPSNAP2(-51.1!) AC009228.1(-29.9!) SELENOI(-26.6!) TUT4(-26.6!) TTC3(-9.8!) DNAAF4(+8.4!) CCDC7(+7.9!) PILRB(-6.9!) LIX1-AS1(-6.4!) UTP4(-6.3!) |
| IN | top1 | PEX2(-3.8) SMG6(-3.6) UBR2(-3.5) ME3(+3.5) NAPB(+3.5) RBBP9(+3.4) SNHG3(+3.4) ACVR1C(-3.4) PRR13(+3.4) BUD13(+3.3) |
| OD | susie | TUT4(-13.5!) LINC01725(+4.8!) SDHA(+4.0!) MMP24OS(+3.9!) MRPL10(+3.8!) NAPB(+3.4) RAB27A(+3.4) SEMA3C(+3.4) PPA1(+3.3) RAB30(+3.2) |
| OD | top1 | FAM71D(-3.8) LPCAT3(-3.6) ZNF12(+3.5) UPF2(+3.5) EXT1(+3.5) DCLRE1A(-3.5) C19orf38(+3.5) CDC26(+3.4) MIIP(-3.3) NAPB(+3.3) |

### Top 10 genes: ABCD_base_EUR  (symbol(z); ! = w_flag; full table report/top10.tsv)

| tag | model | lab |
|---|---|---|
| Ext1 | susie | ODR4(+43.0!) GATC(-43.2!) PRKRA(-4.7!) RMC1(+4.6!) POLE(+4.4!) PILRB(+4.2!) MED20(-4.1!) ACTR1B(-4.1) KLHL23(+3.7!) SRP72(-3.6) |
| Ext1 | top1 | AC159540.2(+4.4) ACTR1B(-4.1) SLC22A25(+4.1) ANKRD36(-4.0) CRAT(-3.6) RPS29(-3.5) SRP72(-3.5) EFCAB13(-3.5) POU6F1(+3.4) COQ4(-3.4) |
| Ext1_n983 | susie | GATC(-45.8!) ODR4(+37.2!) PRKRA(-4.7!) PILRB(+4.6!) ACTR1B(-4.1) POLE(+4.1!) MED20(-4.0!) KLHL23(+3.7!) ANKRD36(-3.6) SRP72(-3.6) |
| Ext1_n983 | top1 | AC159540.2(+4.4) ACTR1B(-4.1) SLC22A25(+4.1) ANKRD36(-4.0) CRAT(-3.6) RPS29(-3.5) SRP72(-3.5) EFCAB13(-3.5) POU6F1(+3.4) COQ4(-3.4) |
| Ext | susie | ODR4(+32.7!) SNORC(-30.5!) GATC(-27.2!) LINC01811(+26.7!) CNOT9(+24.5!) NRDC(-23.4!) MIR2052HG(-22.5!) SELENOI(+19.5!) EEF1AKMT1(+13.1!) NBL1(-12.2!) |
| Ext | top1 | ACTR1B(-4.1) C2orf92(+4.1) AC159540.2(+4.1) ANKRD36(-4.0) EPB41L2(+4.0) BDNF(+3.9) BMPR1A(-3.9) ZNF345(+3.7) COQ4(-3.7) CHML(+3.6) |
| OPC | susie | VPS50(-74.5) PHLDB1(+6.4) ACTR1B(-3.9) AZGP1(-3.6) MAST2(+3.5) ANKRD36(-3.4) TRIOBP(+3.3) FAM184B(+3.3) EFCAB13(-3.3) DPYSL5(-3.2) |
| OPC | top1 | LMAN2L(+4.1) ACTR1B(-4.1) CPOX(+4.1) ANKRD36(-4.1) ANKRD36B(+4.0) TRUB2(+3.7) DNAH1(-3.6) DLEU7(-3.6) EFCAB13(-3.5) AGPAT4(+3.5) |
| IN | susie | NIPSNAP2(+83.7!) TUT4(+38.1!) MTHFS(-35.9!) SELENOI(+34.1!) AC009228.1(-17.0!) DHFR2(+14.6!) LINC01933(-12.2!) DENND5B(-12.1!) TTC3(-9.9!) COQ8A(-5.9!) |
| IN | top1 | BMPR1A(-4.4) ACTR1B(-4.1) ANKRD36(-4.1) COPS7B(+3.7) GPM6A(-3.7) LINC01118(+3.6) PXMP2(-3.6) SLC27A4(+3.5) SWI5(-3.5) WASHC2A(+3.5) |
| OD | susie | TUT4(+23.8!) FOLH1(-11.8!) SDHA(-7.2!) LINC01725(+7.2!) MRPL10(-6.6!) WHRN(+4.4!) PRKRA(-4.3!) ACTR1B(-4.1) TRIOBP(+3.6) EPB41L2(+3.6) |
| OD | top1 | VAT1(+4.2) ACTR1B(-4.1) C2orf92(+4.1) EPB41L2(+4.0) TBC1D24(+4.0) TRIOBP(+3.8) TMEM184B(-3.7) MED6(-3.6) HMGN3(+3.4) SGIP1(+3.3) |


<!-- /AUTO:RESULTS -->

## Positive-control notes (interim)

- **CNTN4, Ext:** top1 z = +7.20 (p 6e-13), in the expected direction (risk allele raises expression). The susie
  model fits no credible set for CNTN4 in Ext, so it is absent from `assoc/*.Ext.susie`.
- **SLC4A10, Ext:** top1 z = −5.13 (p 3e-7), opposite to the paper's direction. The paper's lead and fine-mapped SNP
  rs4572609 is C/G, so the strand-ambiguity filter removed it. top1 is therefore tagging a different SNP. The Ext2
  result, where the paper finds subtype-only colocalisation, is pending.
- **Fig. 3c gene × own-class rows (susie, SCZ25):** 24 are modelled, 21 have p < 0.05, and 14 have z > 0.

## Where the ExN outputs are (CSD3)

Base: `/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/`

| tag | cell type | status | files |
|---|---|---|---|
| Ext1 | L2/3 IT (priority) | done | `assoc/<trait>.Ext1.{susie,top1}.tsv.gz`, `weights/Ext1.{susie,top1,genes}.tsv.gz` |
| Ext1_n983 | L2/3 IT, n = 983 sensitivity | done | `assoc/<trait>.Ext1_n983.*`, `weights/Ext1_n983.*` |
| Ext | all excitatory | done | `assoc/<trait>.Ext.*`, `weights/Ext.*` |
| Ext2–Ext8 | L4 IT, L5/6 NP, L5 IT, L6b, L6 CT, L6 IT Car3, L6 IT | running | same pattern, appear when each tag's assoc job finishes |

`<trait>` is one of `SCZ25_EUR`, `PGC3_EUR`, `MDD_EUR`, `ABCD_thin_EUR` and `ABCD_base_EUR`.

## Deliverables / rsync

Run this on the laptop. Weights contain SNP lists only; no LD matrices are stored.

ExN only (Ext, Ext1, Ext1_n983, and Ext2–8 as they finish):

```bash
rsync -av --include='*.Ext*' --exclude='*.tsv.gz' rajd2@login-icelake.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/assoc/ ~/Git/transcriptional_maturation/gwas_cells/work/sctwas/assoc/
rsync -av --include='Ext*' --exclude='*' rajd2@login-icelake.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/weights/ ~/Git/transcriptional_maturation/gwas_cells/work/sctwas/weights/
```

Everything:

```bash
rsync -av --relative rajd2@login-icelake.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas/./{assoc,weights,report,README.md} ~/Git/transcriptional_maturation/gwas_cells/work/sctwas/
```
