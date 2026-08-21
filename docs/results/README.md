# Results bundle — for figure design

Everything here is **aggregate**. It is the input a figure/report agent needs,
without any individual-level data. See "What is deliberately absent" below.

Narrative and interpretation: [`../../hpc/README_HPC.md`](../../hpc/README_HPC.md)
§8.15 (heritability) and §8.16 (GWAS / rg / MAGMA / PRS).
Method rationale: [`../multianc/README.md`](../multianc/README.md).

## Naming convention

Every filename says which analysis arm produced it. **The arm matters more than
the numbers** — several results differ between arms, and the ancestry-matched
arm is the one to believe (README_HPC §8.16).

| suffix | meaning |
|---|---|
| `00_published_EUR_array` | the original 4.0 EUR result, for comparison |
| `allanc_*` | 7.0 **array** genotypes (456k SNPs), multi-ancestry |
| `imp_pooled*` | 7.0 **imputed** genotypes, pooled multi-ancestry |
| `imp_eur_*` | 7.0 imputed, **EUR-stratified** (ancestry-matched) |
| `*_eur` (gwas/ldsc/magma) | EUR-stratified arm |

## Tables (`tables/`)

**Heritability** — `reml_*.tsv`, columns `phenotype h2 se pval n`.
The headline is `reml_imp_pooled.tsv`: `global_slope` h² = 0.137 ± 0.046,
p = 0.0012. `reml_imp_eur_maf001.tsv` is the validation arm that reproduces the
published 0.575 (it gives 0.474 ± 0.142). `reml_stratified.tsv` has the
per-ancestry estimates — **the non-EUR strata are uninformative, SE 0.23–0.45,
not null**; plot them with their error bars or not at all.

**GWAS** — `gwas_*_summary.tsv` (λ_GC, hit counts) plus thinned per-SNP
statistics in `gwas_thinned/` for Manhattan and QQ plots.

**Genetic correlation** — `ldsc_*_rg.tsv`. Use `ldsc_imp_eur_rg.tsv` as primary:
the multi-ancestry arm applies **European LD scores to a 32 %-non-European
sample** and its intercept (1.0995) shows the misspecification. Every slope row
is flagged `underpowered` in both arms — that flag should reach the figure.

**MAGMA** — `magma_*.tsv` (gene-set covariates incl. AHBA C1–C3) and
`magma_prio_*_{scz,mdd}.tsv` (prioritised sets). Same LD caveat: prefer `_eur`.
**Do not plot the AHBA C1− result as a finding** — it is p = 5.1e-04 in the
multi-ancestry arm and p = 0.20 in the ancestry-matched one (§8.16.4).

**PRS** — `prs_imputed_multianc.tsv`, one row per
disorder × threshold × phenotype × stratum, with `beta se p p_adj r2_partial`
and `n`, `n_families`, `n_snps`. **This is the per-threshold data** for a
threshold-sweep plot. `stratum` is `EUR` or `full`; **`EUR` is the primary arm by
design** (European-discovery scores transfer poorly), and the headline result is
significant in `full` but not `EUR`. A figure showing only `full` would
misrepresent it.

**Diagnostics** — `grm_diagnostics_by_k.tsv`, `selfreport_by_cluster_by_k.tsv`,
`unrelated_composition_by_k.tsv`, `impqual_ksweep_chr22.tsv` back the four
figures in `../multianc/fig/`, which can be restyled from these.

## What is deliberately absent, and why

**Per-subject polygenic scores are NOT here.** `results/prs_imp/*.profile`
contains one row per ABCD participant, keyed by participant ID. That is
individual-level genetic data on minors, covered by the ABCD Data Use
Certification, and **this repository is public**. It is excluded by `.gitignore`
and must not be force-added. The aggregate association results in
`prs_imputed_multianc.tsv` carry every number a figure needs.

Also absent, for size rather than governance: full GWAS summary statistics
(~700 MB each), genotype filesets, and GRMs. All are regenerable from
`hpc/` — see README_HPC §9 and the job IDs in §8.16.

## If you need the individual-level data

Do not route it through this repository. It stays on CSD3 under
`/rds/user/rajd2/hpc-work/abcd_development/hpc/work/results/prs_imp/`.
