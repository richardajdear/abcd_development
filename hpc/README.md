# HPC pipeline: heritability, GWAS, and gene-set enrichment

**Status: written, not executed.** CSD3 refused key authentication during
development (`Permission denied (publickey,keyboard-interactive,hostbased)`)
while an RCS storage recovery was in progress. Every script here is complete
and parameterised; none has been run against 5.1 in this form. Treat the first
run as a dry run and check the assertions in `00_check_inputs.sh` before
launching anything expensive.

These scripts assume the CSD3 conventions from the previous analysis:
account `VERTES-SL2-CPU`, partition `cclake`, working tree
`/home/rajd2/rds/hpc-work/ABCD/`. Override via `hpc/config.sh`.

## What runs where, and why

The local pipeline (`src/abcd/`, `R/fit_lmm.R`) produces subject-level
phenotypes. Everything downstream of that is genotype work: it needs the
imputed genotypes, several GB of LD reference, and GCTA/MAGMA binaries, all of
which live on the cluster. Nothing here needs a GPU; all of it is
CPU-and-IO-bound, and the GRM step is the only one that needs real memory.

## Order of operations

1. `00_check_inputs.sh` — verify genotypes, phenotype file, and binaries exist and that phenotype IDs intersect the `.fam` file. Cheap, run interactively.
2. `01_grm.sbatch` — build the genetic relatedness matrix (GCTA), then prune to unrelated at 0.05. ~2 h, 32 GB.
3. `02_reml.sbatch` — SNP heritability of each phenotype (GCTA REML). Array job over phenotypes.
4. `03_gwas.sbatch` — per-phenotype GWAS (GCTA fastGWA mixed model). Array job.
5. `04_magma.sbatch` — SNP→gene aggregation, then competitive gene-set tests against AHBA components and the snRNA-seq PC1 gene sets.

Submit with dependencies so a failure stops the chain:

```bash
cd hpc
jid1=$(sbatch --parsable 01_grm.sbatch)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 02_reml.sbatch)
jid3=$(sbatch --parsable --dependency=afterok:$jid1 03_gwas.sbatch)
sbatch --dependency=afterok:$jid3 04_magma.sbatch
```

## The phenotype decision that matters most

Use the `ct_genetics.yaml` config, **not** `ct_baseline.yaml`, for anything
here. The two differ in whether the mixed model carries a `family_id` random
intercept. With it, family variance is partialled out of the subject BLUPs and
the DZ/sibling correlation of the intercept phenotype goes to −0.62, which is
biologically impossible and would destroy any heritability estimate. Twin
correlations on the whole-cortex mean (207 MZ / 617 DZ pairs):

| phenotype | family RE | MZ *r* | DZ *r* | Falconer *h²* |
|---|---|---|---|---|
| intercept | yes | 0.14 | −0.62 | — |
| intercept | no | 0.88 | 0.47 | 0.83 |
| slope | yes | — | — | 0.54 |
| slope | no | — | — | 0.58 |

The random-effects structure *is* part of the phenotype definition, not a
fitting detail.

## Expected power, stated up front

On release 5.1 the subject-level slope phenotype has reliability 0.15 (2 visits)
to 0.19 (3 visits) per region, and 0.12 for the whole-cortex mean. Averaging
regions does not help because the dominant noise is per-occasion, not
per-region. A slope GWAS on 5.1 is therefore attenuated roughly 5–17× relative
to a perfectly measured phenotype and is **not** expected to yield genome-wide
hits; run it to establish the pipeline, not to discover loci. The intercept
phenotype (reliability 0.86–0.90 without the family term) is the sanity check —
it should recover the known thickness heritability of ~0.8.

This is the main argument for release 7.0: the added third and fourth timepoints
raise slope reliability, and it is reliability, not sample size, that currently
binds.

## Files this pipeline needs that are not in the repo

| what | where it was | notes |
|---|---|---|
| imputed genotypes | `hpc-work/ABCD/genotypes/` | PLINK bed/bim/fam |
| SCZ sumstats | `gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.fixed.tsv` | SNP col `ID`, p col `PVAL` |
| MDD sumstats (2025) | `gwas/pgc-mdd2025_no23andMe_div_v3-49-46-01_formatted.tsv` | SNP col `rsid`, p col `p_value` |
| MDD sumstats (2018) | `gwas/PGC_UKB_23andMe_depression_10000.txt` | SNP col `MarkerName`, p col `P` |
| MAGMA SNP→gene annot | `magma/snp_gene_annotation/` | `ncbi37.window35-10` etc. |
| LD reference | `magma/g1000_eur/` | 1000G EUR |

## Interpreting `04_magma.sbatch` output

The gene-set test is **competitive**, not self-contained: it asks whether the
AHBA/snRNA-seq gene sets are more associated with the disorder than other genes
of similar size, density, and MAC — those covariates are conditioned internally
(see the `CONDITIONED_INTERNAL` header line in any `.gsa.out`). A
self-contained test would be confounded by gene length alone.

Note that this is the route to the SCZ/MDD hypothesis that does *not* depend on
the map-level result. The spatial analysis found no correlation between the
group-mean developmental change map and AHBA C1–C3 (`figures/abcd_spatial_nulls.png`),
but that is a statement about *where* cortex changes on average, not about
which genes carry disorder risk in a phenotype's association signal. Keep the
two conclusions separate.
