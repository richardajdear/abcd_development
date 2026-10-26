# MAGMA windows: CSD3 report (2026-10-26)

**Status.** All 20 main window × trait gene analyses completed. `w35_10` reproduces the existing annotation byte for byte,
and the ABCD_thin, ABCD_base and MDD gene results exactly (r = 1, max |ΔZ| = 0).

**Open items**
1. **SCZ25_EUR does not reproduce exactly** (Pearson 0.9999, Spearman 1.0000, max |ΔZ| 1.19). Inputs, annotation, SNP
   counts, NPARAM and N are identical, and the two MAGMA logs differ only in paths. Only 25 of 18,546 genes differ, all
   with old Z > 3.5 (22 of them Z > 4.4). Many have rounded P values (2e-9, 5e-10), which suggests MAGMA's sampling
   fallback for very small P. That fallback is stochastic, so it would explain the gap. SCZ has many more such genes than
   MDD (61 vs 6 with P < 1e-6). Rerun job 39250552 (`check_rerun/SCZ25_EUR`) tests this: if the rerun's Z values move
   for the same genes, the gap is sampling noise. The rank-based uses (Spearman, top-1,000 sets) are unaffected
   (Spearman = 1).
2. **`w35_10_multi` for SCZ25_EUR and MDD_EUR timed out** at 4 h. They were resubmitted with 12 h (job 39250551). The
   two ABCD multi runs are complete.

Gene files (`genes/<window>/<trait>.genes.{out,raw}` and `.log`) are on CSD3 at
`abcd_development/genetic_analysis/work/magma_windows/`, which git ignores. Fetch them with the rsync command at the
end of this report. Scripts: `abcd_development/genetic_analysis/magma_windows/`.

---


Task: `transcriptional_maturation/gwas_cells/csd3/MAGMA_WINDOWS_INSTRUCTIONS.md`. Scripts: `genetic_analysis/magma_windows/` (common.sh, annotate.sbatch, genes.sbatch, collect.py).

MAGMA: `MAGMA version: v1.10 (linux/s)`. LD ref g1000_eur, NCBI37.3 gene loc, union-3 snploc, default gene model (snp-wise=mean).

## Reproduction (w35_10 vs existing)

- Annotation identical to `magma_scz2025/ref/ncbi37.window35-10.union3.genes.annot`: **yes** (cmp).

| trait | n_old | n_new | n_common | pearson | spearman | max_abs_dZ |
| --- | --- | --- | --- | --- | --- | --- |
| ABCD_thin_EUR | 18283 | 18283 | 18283 | 1 | 1 | 0 |
| ABCD_base_EUR | 18283 | 18283 | 18283 | 1 | 1 | 0 |
| SCZ25_EUR | 18546 | 18546 | 18546 | 0.9999 | 1 | 1.194 |
| MDD_EUR | 19132 | 19132 | 19132 | 1 | 1 | 0 |

## Per window (annotation)

| window | genes | snps_annot | med_snps_per_gene | p90_snps_per_gene | med_genes_per_snp | mean_genes_per_snp | frac_snps_multi_gene |
| --- | --- | --- | --- | --- | --- | --- | --- |
| w0_0 | 19273 | 26048234 | 515 | 3168 | 1 | 1.016 | 0.013 |
| w10_10 | 19359 | 31078640 | 968 | 3597 | 1 | 1.125 | 0.103 |
| w35_10 | 19367 | 35053290 | 1554 | 4175 | 1 | 1.301 | 0.213 |
| w50_50 | 19378 | 39859920 | 2852 | 5488 | 1 | 1.737 | 0.39 |
| w100_100 | 19379 | 44818564 | 5171 | 7950 | 2 | 2.497 | 0.557 |

## Per window x trait

λGC = median χ²₁(P) / 0.4549. Overlap and Spearman over genes tested in both windows. For w35_10_multi, P = P_MULTI.

| window | trait | genes_tested | lambdaGC | n_P_lt_2_5e6 | spearman_vs_w35_10 | top1000_overlap_vs_w35_10 |
| --- | --- | --- | --- | --- | --- | --- |
| w0_0 | ABCD_thin_EUR | 17882 | 1.002 | 0 | 0.8155 | 680 |
| w0_0 | ABCD_base_EUR | 17882 | 1.134 | 0 | 0.8191 | 673 |
| w0_0 | SCZ25_EUR | 18318 | 3.916 | 580 | 0.8703 | 802 |
| w0_0 | MDD_EUR | 18734 | 3.654 | 353 | 0.8657 | 784 |
| w10_10 | ABCD_thin_EUR | 18213 | 1.041 | 0 | 0.9197 | 778 |
| w10_10 | ABCD_base_EUR | 18213 | 1.143 | 0 | 0.9219 | 774 |
| w10_10 | SCZ25_EUR | 18499 | 4.383 | 706 | 0.9514 | 898 |
| w10_10 | MDD_EUR | 19074 | 4.09 | 425 | 0.9479 | 870 |
| w35_10 | ABCD_thin_EUR | 18283 | 1.031 | 0 | 1 | 1000 |
| w35_10 | ABCD_base_EUR | 18283 | 1.135 | 0 | 1 | 1000 |
| w35_10 | SCZ25_EUR | 18546 | 4.645 | 776 | 1 | 1000 |
| w35_10 | MDD_EUR | 19132 | 4.359 | 438 | 1 | 1000 |
| w50_50 | ABCD_thin_EUR | 18353 | 1.02 | 0 | 0.9037 | 743 |
| w50_50 | ABCD_base_EUR | 18353 | 1.136 | 0 | 0.9093 | 769 |
| w50_50 | SCZ25_EUR | 18585 | 5.256 | 961 | 0.9433 | 890 |
| w50_50 | MDD_EUR | 19206 | 4.976 | 541 | 0.9401 | 859 |
| w100_100 | ABCD_thin_EUR | 18431 | 1.023 | 0 | 0.7759 | 581 |
| w100_100 | ABCD_base_EUR | 18431 | 1.165 | 0 | 0.7878 | 629 |
| w100_100 | SCZ25_EUR | 18621 | 6.253 | 1199 | 0.8701 | 816 |
| w100_100 | MDD_EUR | 19266 | 5.889 | 704 | 0.8591 | 742 |
| w35_10_multi | ABCD_thin_EUR | 18283 | 1.007 | 0 | 0.8599 | 675 |
| w35_10_multi | ABCD_base_EUR | 18283 | 1.206 | 0 | 0.8601 | 709 |

### Top 10 genes by ZSTAT

| window | trait | top10 |
| --- | --- | --- |
| w0_0 | ABCD_thin_EUR | LILRB4 (3.92), CDO1 (3.87), PPIP5K2 (3.81), CSN1S1 (3.55), PRELID2 (3.46), MYO10 (3.36), UBR2 (3.35), PRKACG (3.26), PAM (3.24), CLDND2 (3.24) |
| w0_0 | ABCD_base_EUR | ZNF345 (3.82), AKR1A1 (3.81), C10orf111 (3.78), EPB41L2 (3.78), RPL27 (3.76), ZNF829 (3.59), CST1 (3.56), MYL4 (3.50), MSRB1 (3.50), TMEM33 (3.47) |
| w0_0 | SCZ25_EUR | DPYD (9.80), TCF4 (9.43), PPP1R16B (8.74), CACNB2 (8.74), HIST1H3J (8.68), CACNA1C (8.56), SGCD (8.51), TSSK6 (8.44), HIST1H2BN (8.21), VWA7 (8.21) |
| w0_0 | MDD_EUR | DCC (10.82), SGCZ (9.80), ASTN2 (9.22), SOX5 (9.00), SORCS3 (8.75), SDK1 (8.69), TCF4 (8.17), EYS (8.06), PCLO (7.92), NICN1 (7.91) |
| w10_10 | ABCD_thin_EUR | PPIP5K2 (3.77), OLFML3 (3.55), MYO10 (3.49), BAMBI (3.47), PAM (3.42), PRELID2 (3.39), TRIM6 (3.39), PRPH2 (3.36), CSN1S1 (3.35), OR52B6 (3.29) |
| w10_10 | ABCD_base_EUR | RPL27 (3.87), ZNF345 (3.81), ZNF829 (3.78), EPB41L2 (3.74), MSRB1 (3.71), IFI35 (3.66), MYL4 (3.65), SWI5 (3.53), AKR1A1 (3.50), RUNDC1 (3.49) |
| w10_10 | SCZ25_EUR | DPYD (10.12), TCF4 (9.66), SGCD (8.65), HLA-B (8.65), CACNB2 (8.63), CACNA1C (8.58), NOTCH4 (8.22), TRIM31 (8.20), GATAD2A (8.18), LSM2 (8.13) |
| w10_10 | MDD_EUR | DCC (10.50), SGCZ (9.87), ASTN2 (9.24), SORCS3 (8.88), SDK1 (8.86), SOX5 (8.86), TCF4 (8.33), EYS (8.09), PCLO (7.98), GPX1 (7.93) |
| w35_10 | ABCD_thin_EUR | BAMBI (4.30), TRIM6 (4.15), TRIM6-TRIM34 (3.81), PAM (3.73), MYO10 (3.52), PPIP5K2 (3.48), OR52B6 (3.34), UBR2 (3.31), CSN1S1 (3.27), OLFML3 (3.26) |
| w35_10 | ABCD_base_EUR | EPB41L2 (3.71), RPL27 (3.66), SLC22A24 (3.61), PTGES3L (3.60), ZNF345 (3.59), CST1 (3.56), IFI35 (3.53), SWI5 (3.53), PTGES3L-AARSD1 (3.50), CST4 (3.47) |
| w35_10 | SCZ25_EUR | DPYD (10.37), TCF4 (8.81), MAD1L1 (8.68), SGCD (8.66), SLC17A3 (8.65), CACNB2 (8.56), NOTCH4 (8.51), IER3 (8.46), PPP1R16B (8.43), BTN2A2 (8.41) |
| w35_10 | MDD_EUR | DCC (10.71), SGCZ (9.91), ASTN2 (9.09), SOX5 (8.88), SDK1 (8.87), DRD2 (8.82), SORCS3 (8.71), TCF4 (8.35), PCDH9 (8.25), EYS (8.14) |
| w50_50 | ABCD_thin_EUR | OR52B6 (4.10), TRIM6-TRIM34 (4.07), TRIM34 (3.88), TRIM6 (3.88), BAMBI (3.74), MYO10 (3.67), PAM (3.62), UBR2 (3.50), PRPH2 (3.44), PPIP5K2 (3.39) |
| w50_50 | ABCD_base_EUR | ZAP70 (3.80), EPB41L2 (3.72), CST1 (3.69), ZNF345 (3.68), TMEM131 (3.56), CPOX (3.55), FAM73B (3.52), ZNF829 (3.47), AKR1A1 (3.43), PAICS (3.42) |
| w50_50 | SCZ25_EUR | DPYD (10.26), NOTCH4 (9.27), TCF4 (9.17), IGSF9B (8.75), CACNB2 (8.68), HIST1H4C (8.62), SGCD (8.58), NT5DC2 (8.52), HIST1H2BK (8.51), PITPNM2 (8.47) |
| w50_50 | MDD_EUR | DCC (10.23), SGCZ (9.94), ASTN2 (9.18), DRD2 (8.83), SORCS3 (8.82), TCF4 (8.81), SDK1 (8.80), SOX5 (8.57), NCAM1 (8.54), NRXN1 (8.45) |
| w100_100 | ABCD_thin_EUR | OR52B6 (3.96), TRIM6 (3.83), LOC101927763 (3.67), TRIM34 (3.50), TRIM6-TRIM34 (3.49), C5orf30 (3.49), UBR2 (3.47), PAM (3.45), OR52H1 (3.37), ACOT1 (3.34) |
| w100_100 | ABCD_base_EUR | ZAP70 (4.00), COX5B (3.90), ACTR1B (3.80), SLC39A10 (3.69), SMLR1 (3.67), SLC22A25 (3.65), EPB41L2 (3.63), TMEM131 (3.62), FAM73B (3.56), CST4 (3.53) |
| w100_100 | SCZ25_EUR | DPYD (10.66), PPT2 (10.08), C6orf10 (10.06), AGPAT1 (10.04), EGFL8 (9.91), PBX2 (9.84), NOTCH4 (9.53), GPSM3 (9.50), RNF5 (9.39), DRD2 (9.35) |
| w100_100 | MDD_EUR | DCC (10.29), SGCZ (10.05), ASTN2 (9.45), TCF4 (9.36), NCAM1 (9.26), SORCS3 (9.12), TTC12 (8.90), SDK1 (8.70), SOX5 (8.69), DRD2 (8.53) |
| w35_10_multi | ABCD_thin_EUR | BAMBI (4.55), MYO10 (3.92), TRIM6 (3.87), PPIP5K2 (3.73), TRIM6-TRIM34 (3.69), PAM (3.67), UBR2 (3.66), WDR26 (3.47), C5orf30 (3.46), CSN1S1 (3.45) |
| w35_10_multi | ABCD_base_EUR | RPL27 (3.90), PTGES3L (3.90), EPB41L2 (3.82), PTGES3L-AARSD1 (3.80), ZNF345 (3.79), ZNF790 (3.78), IFI35 (3.73), CPOX (3.72), NASP (3.71), MIR205HG (3.69) |

## Files

`genes/<window>/<trait>.genes.{out,raw}` and `<trait>.log` (MAGMA names the log `<out>.log`). Annotations in `ref/` (not copied back; large).

Rsync to local (run on the local machine):

```bash
rsync -av --include='*/' --include='*.genes.out' --include='*.genes.raw' --include='*.log' --include='README.md' --exclude='*' rajd2@login.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/magma_windows/ ~/Git/transcriptional_maturation/gwas_cells/work/magma_windows/ --prune-empty-dirs
```
