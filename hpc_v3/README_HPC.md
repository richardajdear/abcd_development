# hpc_v3 — Experiment A on the cluster: progress and results

**Run log for [`README.md`](README.md) §2, started 2026-09-08 on CSD3
(account VERTES-SL3-CPU, partition icelake).** Documentation contract is
`hpc_v2/README_HPC.md` §10: every result carries its job ID, its script, the
**N** it ran on, the file it landed in, and whether it replaces or supplements
an existing number. Failures and dead ends are recorded with their diagnosis.

**Scope decision carried from the instruction that opened this run:**
**C+T only for PRS at this stage.** PRS-CS and SBayesR exist for the settled
five (`hpc_v2/README_HPC.md` §13.5) and are deliberately not run for the new
phenotypes yet. This is a power/scope choice, not a verdict on the methods —
and it means the §13.5 caveat (the SCZ→`global_slope` C+T association does not
survive PRS-CS) still hangs over every C+T number below.

---

## Status

| step | what | job | state |
|:---|:---|:---|:---|
| 0 | build + align the v3 GCTA export | local | **done** |
| 2 | GENESIS null models, pooled | 35089435 | **done** (n = 8,082) |
| 3 | GENESIS association, pooled (4 × 22) | 35089436 | **done** (88/88) |
| 2e | GENESIS null models, EUR arm | 35089468 | **done** (n = 4,116) |
| 3e | GENESIS association, EUR arm (4 × 22) | 35089470 | **done** (88/88) |
| 3c | collection, both arms | 35089437 + 35089472 **FAILED** → 35093677 | **done** (§6.1) |
| 4 | Zaitlen REML, κ ≥ 0.0250 | 35089439 | **done** |
| 4 | Zaitlen REML, κ ≥ 0.0884 | 35089440 | **done** |
| 6 | C+T PRS, population + within-family | 35089628 | **done** |
| 6b | paired Δβ, min-p permutation, random-region null | 35089862 | **done** |
| 5+7 | LDSC (9 phenotypes) and MAGMA (4), EUR arm | 35089930 cancelled → 35093678 | **done** |
| — | rebuild `h2_vs_scz_relevance.csv` | local | **done**, 9 rows |

**Experiment A is complete, and the answer is negative on every readout.**
Restricting the thinning phenotype to the fastest-thinning or the
highest-AHBA-C3 cortex does not improve heritability (§8, §7), the disorder
association (§5), or the transcriptional coupling (§10). `slope_topC3` is
weaker than ~80 % of *random* 8-region sets (§5.2). Full verdict: §11.

New scripts written for this run, all committed:

| file | what |
|:---|:---|
| `align_export_v3.py` | fixes the three export mismatches in §1; asserts the settled five reproduce v2's `$PHENO` |
| `01_prs_ct.sbatch` | C+T PRS, population + within-family, over the combined score directory |
| `03_prs_paired.sbatch` + `prs_paired_delta.py` | paired Δβ against `global_slope`, min-p permutation, random-region null |
| `02_ldsc_magma.sbatch` | LDSC (9 phenotypes) and MAGMA (4) on the EUR sumstats, into fresh output directories |
| `regen_h2_relevance.py` | rebuilds `h2_vs_scz_relevance.csv` from source files; reproduces all five committed rows to 0 |

**No `hpc_v2` or `hpc` pipeline code was modified.** Every step above either
runs an existing script with environment overrides, or is a new script that
consumes existing outputs — which is what `README.md` §2 asks for.

---

## 1. Step 0 — the export, and three silent mismatches it had to be fixed for

`hpc_v3/make_phenotypes_v3.py` ran locally against
`out/thickness_dsk_70_139406217085` under `ABCD_CONFIG=ct_70_noglobal_mv2_genetic`
and reproduced all three committed summary tables **byte-identically**
(`git status hpc_v3/` clean afterwards), so the phenotype definitions in
`README.md` §1 are the ones analysed here.

Its `gcta_inputs_v3/` output is **not** what the cluster consumes, though.
Comparing it against `hpc/work/pheno_allanc_prs/` — the `$PHENO` every v2 step
actually ran on — found three differences, none of which would have raised an
error:

| # | difference | what it would have done |
|:---|:---|:---|
| 1 | IID spelled `NDAR_INV005V6D2C`; every `.fam`, `.grm.id` and `.profile` on CSD3 spells it `sub-005V6D2C` | GENESIS normalises to the 8-char token and is fine; **GCTA matches FID+IID literally**, so step 05 would have returned a zero-subject REML |
| 2 | phenotypes in native units | v1 z-scored every column (`hpc/work/standardise_pheno.py`) because `global_slope` has Var = 1.06e-06. Every β in this project is **per phenotype SD**; a native-unit column would have made the paired Δβ against `global_slope` a comparison of two different scales |
| 3 | 8,192 rows (all phenotyped) vs the 8,082 phenotyped-**and-genotyped** analysis sample | standardising over 8,192 puts `global_slope` on a slightly different scale than v2's published number |

Also, the fresh export carries **no ancestry PCs**: `gcta_export.build`'s
release-adapter lookup for `genetic_pcs` raises and the exception is swallowed
into `n_pcs = 0`. GENESIS substitutes PC-AiR PCs regardless, but `prs_assoc.R`
and GCTA both read `covar_quant.txt` and would have run unadjusted for
structure — silently, with only a warning in a log.

**Fix:** `hpc_v3/align_export_v3.py` (new, committed) adopts the reference
export's FID/IID spelling, subject set and row order, z-scores over the 8,082,
and copies v2's covariate files rather than re-deriving them (their
`baseline_age`, `n_visits`, `sex` and `site` values were verified to agree
**exactly** with the fresh export for all 8,082 shared subjects, so the copy
adds the PCs and changes nothing else).

**Verification that the chain is the same one v2 ran.** After alignment the
five settled columns reproduce v2's `$PHENO` to
**max |v3 − v2| = 8.9e-16** — the script asserts this and refuses to write
otherwise. So the four new columns sit in a file that is a column-wise superset
of the one behind every published v2 number.

Native-unit scale factors (multiply a standardised β by these to get native
units), computed over the 8,082:

| phenotype | mean | SD (native) |
|:---|---:|---:|
| `baseline_thickness` | 2.71951 | 0.0627813 |
| `global_slope` | −0.0189205 | 0.00120025 |
| `slope_PC3` | 0.015452 | 15.3242 |
| `slope_PC2` | −0.013762 | 16.4818 |
| `slope_PC1` | −0.0192588 | 16.3192 |
| **`slope_topDelta`** | **−0.0252729** | **0.0017015** |
| **`slope_topC3`** | **−0.0227597** | **0.00162507** |
| **`slope_projDelta`** | 8.29874e-06 | 0.0624419 |
| **`slope_projC3`** | −0.00161796 | 12.8579 |

The two subset means thin faster than the cortical average in native units
(−0.0253 and −0.0228 mm/yr against −0.0189), which is what selecting on the
group thinning map and on C3 is supposed to do, and a first check that the
region sets are not mis-joined.

**Landed in:** `hpc_v2/work/pheno_v3/` (gitignored; per-subject).
Preflight `hpc_v2/00_check_inputs.sh` passed with
`8082 of 8082 phenotyped subjects present in the array .fam`,
`1250 multi-member families`, `manifest covers 4 phenotypes within 9 exported
columns`.

## 2. Steps 1–2 — kinship reused, null models fitted

Step 01 (GDS) and 02 (kinship) were **not** re-run, per `README.md` §2 and
`SETUP_CONTEXT.md` §6: `work/results_v2/kinship/` (pooled) and
`kinship_eur/` (EUR arm) were reused unchanged.

**Job 35089435**, `hpc_v2/03_null_model.sbatch --array=1-4`, 56 s per task.
`nullmodel/<pheno>_null_summary.tsv`:

| phenotype | n | varcomp_kin | varcomp_resid | prop_kin |
|:---|---:|---:|---:|---:|
| `slope_topDelta` | 8082 | 0.3518 | 0.5666 | 0.383 |
| `slope_topC3` | 8082 | 0.3708 | 0.5259 | 0.413 |
| `slope_projDelta` | 8082 | 0.3320 | 0.6455 | 0.340 |
| `slope_projC3` | 8082 | 0.3693 | 0.5535 | 0.400 |
| *(`global_slope`, v2, for reference)* | *8082* | *0.3400* | *0.5415* | *0.386* |

**n = 8,082 exceeds v1's 5,649**, which is the check `hpc_v2/README_HPC.md` §6
asks for: relatives were kept. `prop_kin` is the *total* familial share
(pedigree, not SNP h²) and sits in the same 0.34–0.41 band as the settled
phenotypes — the projections are not obviously less familial than the subset
means despite their lower internal consistency.

## 3. Step 4 — Zaitlen two-GRM heritability

The κ ≥ 0.0884 arm (job 35089440) landed within minutes; the κ ≥ 0.0250 arm
(35089439) took up to 1 h 50 m per task. **Both arms and the full comparison
`README.md` §2 asks for are in [§8](#8-step-4-completed--zaitlen-reml-at-both-thresholds)**,
written once the slower arm finished rather than twice.

Headline, so it is not buried: **no new phenotype is more SNP-heritable than
`global_slope` at either kinship threshold**, and the 31–51 % threshold
instability `SETUP_CONTEXT.md` §7 records for `global_slope` reproduces on all
four.

## 4. Method notes fixed before any result was read

Recorded here because both are the kind of silent defect
`hpc_v2/README_HPC.md` §10 asks for, and both were found while wiring the run
rather than by a result looking wrong.

### 4.1 Every published PRS β in this project was fitted without an age covariate

`tools/prs_assoc.R`'s model line reads

```r
terms <- c("scale(PRS)", "sex", "age_c", pc_cols)
terms <- terms[terms %in% names(sub) | terms == "scale(PRS)"]
```

`age_c` is not a column of `covar_quant.txt` — `abcd.gcta_export` renames it
`baseline_age` on export — so the filter silently drops it and the fitted model
is `scale(y) ~ scale(PRS) + sex + PC1..PC10 + (1 | family_id)`. This affects
v1's `prs_association.tsv`, v2's control-disorder tables and §12.8, not just
this run.

**It is almost certainly immaterial** — the phenotype is a *slope* already
adjusted for age in the mixed model that produced it, and `baseline_age` is a
covariate of the export rather than of the trait. But it is not what the code
claims to fit. **Handled by running both:** the primary arm here reproduces the
model exactly as fitted (a comparison against v2's numbers is worthless if the
model differs), and `--with-age` refits with `baseline_age` and `n_visits`
added. Both tables are reported below; if they disagree, the disagreement is
the finding.

### 4.2 Fresh output directories for LDSC and MAGMA

`hpc/05_ldsc_rg.sbatch` and `hpc/04_magma.sbatch` write
`ldsc_rg_summary.tsv` and `magma_gsa_summary.tsv` covering **exactly the
phenotypes in `$MANIFEST`**. Pointing them at `ldsc_eur/` and `magma_eur/`
with a four-phenotype manifest would have replaced the published v2 tables
with four-row versions — no error, no backup. `02_ldsc_magma.sbatch` therefore
writes to `ldsc_eur_v3/` and `magma_eur_v3/` and symlinks the disorder-side
intermediates (munged sumstats, `.genes.raw`) in so they are reused rather than
recomputed.

The GENESIS collector does not have this problem: `R/05_collect_assoc.R`
updates `gwas_summary.tsv` in place per phenotype, so the pooled and EUR
summaries gain rows rather than losing them (verified — the five settled rows
are still present alongside the new ones).

## 5. Step 6 — C+T PRS. **HEADLINE: the SNR-concentration hypothesis fails.**

**Jobs 35089628** (`hpc_v3/01_prs_ct.sbatch`, population + within-family via
the unchanged `tools/prs_assoc.R` and `hpc_v2/R/06_prs_family.R`) and
**35089862** (`hpc_v3/03_prs_paired.sbatch` → `prs_paired_delta.py`, the paired
Δβ, the min-p permutation and the random-region null). n = 8,082 pooled,
4,116 EUR (the anchor-set definition, §1). Nothing was re-scored — a
polygenic score applies consortium weights to ABCD genotypes and never touches
our phenotypes, so the existing C+T `score_*.profile` files are exactly the
ones a new phenotype needs. **C+T only, per the instruction opening this run.**

**Reproduction checks passed before anything new was read:**

| quantity | published | this run |
|:---|---:|---:|
| SCZ p<0.5 → `global_slope`, pooled, population β (p) | −0.041 (4.3e-03) `hpc/README_HPC.md` | **−0.0413 (4.29e-03)** |
| SCZ p<0.5 → `global_slope`, EUR, β_between (p) | −0.0460 (3.9e-03) `SETUP_CONTEXT.md` §7 | **−0.0460 (3.88e-03)** |

### 5.1 The primary readout — paired Δβ against `global_slope`

`README.md` §1 is explicit that a bigger point estimate is not a result. Δβ is
computed **exactly**: both phenotypes are regressed on the same design in the
same subjects, so the difference of the PRS coefficients *is* the PRS
coefficient of a regression on the difference score, and its family-clustered
sandwich SE is the right SE for it. The 2,000-replicate family bootstrap
`README.md` §1 asks for is reported alongside and agrees throughout.

**Δβ > 0 means the new phenotype is *weaker* than the global mean.**
SCZ, p<0.5:

| stratum | phenotype | β | p | **Δβ** | SE | p(Δβ) | bootstrap 95 % CI |
|:---|:---|---:|---:|---:|---:|---:|:---|
| EUR | `global_slope` | −0.0472 | 0.0020 | — | — | — | — |
| EUR | `slope_topDelta` | −0.0432 | 0.0053 | +0.0040 | 0.0076 | 0.60 | [−0.0104, +0.0193] |
| EUR | `slope_topC3` | −0.0374 | 0.018 | +0.0098 | 0.0068 | 0.15 | [−0.0032, +0.0229] |
| EUR | `slope_projDelta` | −0.0185 | 0.24 | +0.0287 | 0.0180 | 0.11 | [−0.0063, +0.0626] |
| EUR | `slope_projC3` | −0.0077 | 0.64 | **+0.0395** | 0.0171 | **0.021** | [+0.0066, +0.0734] |
| pooled | `global_slope` | −0.0436 | 0.0025 | — | — | — | — |
| pooled | `slope_topDelta` | −0.0439 | 0.0029 | −0.0003 | 0.0071 | 0.96 | [−0.0148, +0.0137] |
| pooled | `slope_topC3` | −0.0392 | 0.0079 | +0.0044 | 0.0063 | 0.48 | [−0.0078, +0.0165] |
| pooled | `slope_projDelta` | −0.0082 | 0.59 | **+0.0354** | 0.0175 | **0.043** | [−0.0003, +0.0675] |
| pooled | `slope_projC3` | −0.0093 | 0.54 | **+0.0342** | 0.0164 | **0.037** | [+0.0026, +0.0664] |

**No Δβ is negative beyond noise.** (Checked against
`prs_paired_delta.tsv` across all 320 rows: 81 are negative, all trivially —
max |Δβ| = 0.0065, every p ≥ 0.41, including this table's pooled
`slope_topDelta` at −0.0003. An earlier draft of this section said "every Δβ
is positive", which its own table contradicts.) The two subset means are statistically
indistinguishable from `global_slope` (Δβ p = 0.15–0.96, bootstrap intervals
straddling zero); the two projections — the traits that are genuinely different
from the global mean (r = 0.28/0.41) — **significantly lose** the association.
The same pattern holds for MDD (Δβ p = 0.14–0.99 for the subset means,
0.046–0.048 for the projections in the pooled arm).

**This answers Experiment A's question in the negative.** Restricting to the
fastest-thinning or the highest-AHBA-C3 cortex does not sharpen the disorder
association. `README.md` §1 anticipated exactly this reading: the subset means
share ~80–84 % of their variance with `global_slope`, and what little is
distinctive about them does not carry disorder signal.

### 5.2 The random-region null — and this is the stronger result

The paired test asks whether a subset beats the global mean. It does not ask
whether *this* subset beats an arbitrary one. 1,000 random 8-bilateral-region
subset means, built from the same BLUPs and pushed through the same
regression (percentile = fraction of random sets with a *more negative* β, so
**lower = stronger**):

| stratum | phenotype | SCZ p<0.5 | SCZ p<0.01 | MDD p<0.5 | MDD p<0.01 |
|:---|:---|---:|---:|---:|---:|
| EUR | `global_slope` | 19.3 | 37.6 | 45.2 | 40.7 |
| EUR | `slope_topDelta` | 43.3 | 65.5 | 44.7 | 21.8 |
| EUR | **`slope_topC3`** | **78.0** | **80.4** | **75.6** | **80.3** |
| pooled | `global_slope` | 32.8 | 45.4 | 35.6 | 35.9 |
| pooled | `slope_topDelta` | 31.5 | 53.8 | 34.8 | 25.7 |
| pooled | **`slope_topC3`** | **53.8** | 49.3 | **78.1** | **81.0** |

`slope_topC3` is **weaker than roughly 80 % of randomly chosen 8-region sets**
in the EUR arm, consistently across both disorders and both thresholds.
`slope_topDelta` sits near the middle — indistinguishable from an arbitrary
subset. `global_slope` is at or below the null median, which is what averaging
all 68 regions should buy.

This is the same shape as `hpc/README_HPC.md` §2's finding that a
**PRS-selected** top-10-region phenotype "sits at the median of random
10-region sets", and it extends it: the selection here was made from *group
maps only*, with no PRS and no genotypes involved, so it is not the same
circularity — and it still does not beat chance. The C3-selected set does
worse than chance.

### 5.3 Multiplicity done properly, and the control panel

`hpc_v2/README_HPC.md` §7: Bonferroni over 8 nested C+T thresholds is the wrong
correction. Family-block permutation of min-p (2,000 permutations; whole
families exchanged with families of the same size, so within-family correlation
survives and only the score link is broken):

| stratum | disorder | phenotype | min p | **p_perm** | p_Bonferroni |
|:---|:---|:---|---:|---:|---:|
| EUR | SCZ | `global_slope` | 0.0020 | **0.0070** | 0.016 |
| EUR | SCZ | `slope_topDelta` | 0.0053 | **0.027** | 0.043 |
| EUR | SCZ | `slope_topC3` | 0.018 | 0.065 | 0.14 |
| pooled | SCZ | `global_slope` | 0.0025 | **0.013** | 0.020 |
| pooled | SCZ | `slope_topDelta` | 0.0029 | **0.010** | 0.023 |
| pooled | SCZ | `slope_topC3` | 0.0065 | **0.033** | 0.052 |

Permutation is consistently less conservative than Bonferroni, as §7 predicts.
`global_slope` remains the strongest phenotype in the EUR arm; no new phenotype
improves on it.

**The control panel does not discriminate — again.** At best-threshold, EUR:
`ALZ → slope_topC3` β = −0.0346, p = 0.023 against `SCZ → slope_topC3`
β = −0.0374, p = 0.016. This is §13.6's finding reproducing on the new
phenotypes, and for the same reason: a control panel can only demonstrate
specificity when the target effect is itself detectable.

**One apparent control failure that is not one.** The min-p permutation flags
`ASD → global_slope` in the pooled arm (min p = 0.0018, p_perm = 0.013). It
comes entirely from the p<1e-5 score — **69 SNPs** — and its sign is
**+0.036**, the *opposite* direction to SCZ. Min-p is direction-agnostic, so it
fires on a signed-opposite association. This is not evidence of a non-specific
confound in the SCZ direction; recorded so the next reader does not have to
re-derive it.

### 5.4 Within-family (Fulker) decomposition

`hpc_v2/R/06_prs_family.R`, 686 informative EUR pairs / 1,339 pooled. SCZ
p<0.5:

| stratum | phenotype | β_W (p) | β_B (p) | p_diff |
|:---|:---|---:|---:|---:|
| EUR | `global_slope` | −0.0569 (0.30) | −0.0460 (0.0039) | 0.85 |
| EUR | `slope_topDelta` | −0.0138 (0.80) | −0.0451 (0.0054) | 0.59 |
| EUR | `slope_topC3` | −0.0565 (0.31) | −0.0358 (0.027) | 0.72 |
| pooled | `global_slope` | +0.0070 (0.89) | −0.0456 (0.0025) | 0.31 |
| pooled | `slope_topDelta` | +0.0287 (0.57) | −0.0481 (0.0017) | 0.14 |
| pooled | `slope_topC3` | −0.0060 (0.90) | −0.0409 (0.0077) | 0.50 |

**No `p_diff` is significant**, so this run gives no evidence that the
population estimates for the new phenotypes are inflated by between-family
confounding — and, at ~18 % power, no evidence against it either. The script
prints that caveat itself and §2.3 of `hpc_v2/README_HPC.md` governs how to
read it. A null β_W here is uninformative, not negative.

### 5.5 The missing age covariate (§4.1) is immaterial

Both arms were run. Over all 320 (disorder × threshold × stratum × phenotype)
rows, adding `baseline_age` and `n_visits`:

| quantity | max |primary − with-age| |
|:---|---:|
| β, new phenotype | 0.0055 |
| β, `global_slope` | 0.0010 |
| Δβ | 0.0052 |
| p(Δβ) | 0.157 → the largest single p-value shift is 0.021 → 0.025 |

No conclusion in §5.1–5.4 changes. **The §4.1 defect is real but harmless** —
worth fixing in `tools/prs_assoc.R` so the code fits what it claims, not worth
revisiting any published number over. Files:
`prs_paired_delta_withage.tsv`, `prs_minp_permutation_withage.tsv`,
`prs_random_region_null_withage.tsv`.

## 6. Step 3 — the association scans, and a collector failure worth recording

**Jobs 35089436** (pooled, 88 tasks) and **35089470** (EUR, 88 tasks), both
`hpc_v2/04_assoc.sbatch` unchanged. **176/176 tasks COMPLETED.**

### 6.1 Defect: both collection jobs died on temp-directory space

`run_all.sh`'s final collection passes (**35089437** pooled, **35089472** EUR)
each failed after 5 seconds:

```
gzip: stdout: No space left on device
Error in fread(cmd = sprintf("gzip -dc %s", shQuote(f))) :
  External command failed with exit code 1. This can happen when the disk is
  full in the temporary directory ('/tmp/RtmpaNR80v').
```

`data.table::fread(cmd=)` decompresses each chromosome file through R's
`tempdir()`, which on these nodes is a small node-local `/tmp` shared with
every other job on the node. **Nothing was wrong with the scan** — all 176
tasks completed and all 22 chromosome files were present for all four
phenotypes in both arms. The diagnosis matters because the failure looks like
a data problem and is not one.

Two knock-on effects, both handled:

- 35089472's failure left **35089930** (LDSC/MAGMA) unrunnable on its
  `afterok` dependency. Cancelled and resubmitted behind the fixed collector.
- The inline `chr == 22` collects inside `04_assoc.sbatch` had already written
  **partial** `<pheno>.sumstats.tsv.gz` files (1, 3 and 6 chromosomes) and
  partial `gwas_summary.tsv` rows. A non-empty-file check would have passed on
  those. **`n_chr` in `gwas_summary.tsv` is the check that catches it** — it is
  22 for every row below.

**Fix: environmental, not a code change.** `hpc_v3/04_collect.sbatch`
(**35093677**, 20 m 32 s) reruns the collector with `TMPDIR` on `/rds`
(258 G free) and 48 G of memory. The collector is idempotent and rebuilds from
whatever chromosome files exist, which is exactly what this needed.

### 6.2 GWAS results

Pooled arm, n = 8,082 (mean per-variant N 7,893), 9,413,281 variants:

| phenotype | λ_GC | p<5e-8 | p<1e-5 | min p | n_chr |
|:---|---:|---:|---:|---:|---:|
| `baseline_thickness` | 1.0363 | 0 | 302 | 8.2e-08 | 22 |
| `global_slope` | 1.0139 | 0 | 117 | 5.1e-07 | 22 |
| `slope_PC1` / `PC2` / `PC3` | 1.0142 / 1.0358 / 1.0048 | 0 / 1 / 2 | 50 / 121 / 107 | — | 22 |
| **`slope_topDelta`** | **1.0202** | **0** | 91 | 2.1e-07 | 22 |
| **`slope_topC3`** | **1.0135** | **0** | 193 | 5.7e-08 | 22 |
| **`slope_projDelta`** | **1.0342** | **0** | 144 | 8.0e-08 | 22 |
| **`slope_projC3`** | **1.0366** | **0** | 184 | 7.2e-08 | 22 |

EUR arm, n = 4,116 (mean N 4,017), 7,469,278 variants:

| phenotype | λ_GC | p<5e-8 | p<1e-5 | min p | n_chr |
|:---|---:|---:|---:|---:|---:|
| `global_slope` | 1.0001 | 0 | 151 | 1.1e-07 | 22 |
| **`slope_topDelta`** | **0.9952** | **0** | 77 | 1.1e-07 | 22 |
| **`slope_topC3`** | **0.9997** | **0** | 41 | 1.1e-06 | 22 |
| **`slope_projDelta`** | **0.9978** | **0** | 38 | 5.3e-08 | 22 |
| **`slope_projC3`** | **1.0018** | **0** | 25 | 4.1e-07 | 22 |

**Zero genome-wide hits on any new phenotype in either arm, at λ_GC between
0.995 and 1.037.** The λ values are the point: `SETUP_CONTEXT.md` §7 requires
them before a hit count means anything, and v1's 20 "hits" at λ = 1.107 were
stratification. These scans are calibrated, so the zeros are real zeros.
`slope_topC3` has the largest p<1e-5 count in the pooled arm (193 vs
`global_slope`'s 117) at a *lower* λ (1.0135 vs 1.0139) — suggestive, but
0 of them clear 5e-8 and the LDSC h² below says the excess is not
heritability.

## 7. Step 5 — LDSC. **The new phenotypes are LESS heritable, and rg is uninformative.**

**Job 35093678** (`hpc_v3/02_ldsc_magma.sbatch`, 3 h 32 m), LDSC on the EUR
arm, all nine phenotypes in one call. Output: `ldsc_eur_v3/`.

**Reproduction check:** all five settled rows reproduce their published
`h2_obs`, `h2_obs_se`, `h2_z` and `rg` **exactly** (verified programmatically —
max |diff| = 0 on every column; see §9).

| phenotype | LDSC h²_obs | SE | **h² z** | rg with SCZ | SE |
|:---|---:|---:|---:|---:|---:|
| `baseline_thickness` | 0.4671 | 0.104 | **4.49** | +0.026 | 0.052 |
| `slope_PC2` | 0.3773 | 0.109 | **3.46** | −0.052 | 0.060 |
| `global_slope` | 0.1238 | 0.107 | 1.16 | −0.149 | 0.109 |
| `slope_PC1` | 0.0676 | 0.103 | 0.66 | +0.094 | 0.150 |
| **`slope_projDelta`** | 0.0700 | 0.106 | **0.66** | +0.061 | 0.137 |
| **`slope_projC3`** | 0.0595 | 0.096 | **0.62** | −0.194 | 0.233 |
| **`slope_topC3`** | 0.0388 | 0.098 | **0.40** | −0.309 | **0.453** |
| **`slope_topDelta`** | 0.0341 | 0.114 | **0.30** | −0.072 | 0.242 |
| `slope_PC3` | −0.2398 | 0.094 | −2.55 | n/a | — |

**Every new phenotype has a LOWER LDSC h² z than `global_slope`'s already
inadequate 1.16.** The two subset means are the worst of the nine positive
ones (z = 0.30 and 0.40). This agrees with the GENESIS REML (§3, §8) and is
the second independent way Experiment A comes out negative.

**The rg verdict, as `README.md` §2 item 3 requires:** *rg is uninformative,
not null, for all four.* `slope_topC3` has the largest |rg| in the whole
project (−0.309) and it means nothing — its SE is **0.453**, an interval
running from −1 to +0.58. Reporting that number without its SE would be the
single most misleading thing this run could produce.

**A concrete demonstration of the go/no-go concern.** MDD × `slope_topDelta`
did not merely return a wide interval — **LDSC crashed**:

```
FloatingPointError: invalid value encountered in sqrt
  regressions.py:705  np.multiply(hsq1.tot_delete_values, hsq2.tot_delete_values)
```

rg is gencov / √(h²₁·h²₂); at h² z = 0.30 the jackknife h² goes negative in
some blocks and the denominator has no square root. This is not a bug to route
around — it is LDSC declining to normalise by a heritability that is not
distinguishable from zero. **It is direct, mechanical evidence for
`SETUP_CONTEXT.md` §5.1's gate**, and it fired on a phenotype from this very
experiment.

## 8. Step 4 (completed) — Zaitlen REML at both thresholds

**Job 35089439** (κ ≥ 0.0250, 42 m – 1 h 50 m per task) and **35089440**
(κ ≥ 0.0884). n = 8,082 throughout. `README.md` §2 requires both, because
`global_slope` swung 0.183 → 0.115 between them.

| phenotype | h²_SNP at κ ≥ 0.0250 | h²_SNP at κ ≥ 0.0884 | swing |
|:---|---:|---:|---:|
| `baseline_thickness` | 0.5393 ± 0.055 | 0.5480 ± 0.052 | +0.9 % |
| `global_slope` | 0.1830 ± 0.049 | 0.1149 ± 0.056 | **−37 %** |
| `slope_PC1` | 0.1779 ± 0.051 | 0.1394 ± 0.056 | −22 % |
| `slope_PC2` | 0.2296 ± 0.050 | 0.1925 ± 0.058 | −16 % |
| `slope_PC3` | 0.0961 ± 0.052 | 0.0992 ± 0.056 | +3 % |
| **`slope_topDelta`** | **0.1535 ± 0.042** | **0.0850 ± 0.056** | **−45 %** |
| **`slope_topC3`** | **0.1792 ± 0.052** | **0.1234 ± 0.056** | **−31 %** |
| **`slope_projDelta`** | **0.1645 ± 0.039** | **0.0806 ± 0.056** | **−51 %** |
| **`slope_projC3`** | **0.1945 ± 0.041** | **0.1206 ± 0.057** | **−38 %** |

**No new phenotype beats `global_slope` at either threshold.** `slope_topC3`
tracks it almost exactly (0.179 vs 0.183; 0.123 vs 0.115); `slope_topDelta` is
below it at both. The threshold instability §7 of `SETUP_CONTEXT.md` records
for `global_slope` **reproduces on all four new phenotypes**, at 31–51 %, so
it is a property of the kinship threshold and not of any one trait. Only
`baseline_thickness` — the one well-powered phenotype — is stable.

**Recorded so it is not mistaken for a result:** `slope_topDelta` and
`slope_topC3` both printed `variance component 1 is dropped from the model`
mid-fit at κ ≥ 0.0250 and restarted from EM-REML priors. Both converged
normally afterwards; the values above are the converged ones. It is an
intermediate state of GCTA's AI-REML, not a failure.

## 9. Step 5 item 5 — the h²-vs-relevance table, out of sample

`hpc_v3/regen_h2_relevance.py` rebuilds `h2_vs_scz_relevance.csv` from the
source files. The column provenance was undocumented and had to be
reverse-engineered; the script's five rebuilt settled rows reproduce the
committed file to **max |diff| = 0** on every column, which confirms it:

| column | source |
|:---|:---|
| `h2_obs`, `h2_obs_se`, `h2_z`, `rg` | SCZ rows of LDSC EUR `ldsc_rg_summary.tsv` |
| `p_between` | C+T SCZ p<0.5, EUR, `beta_between`'s p in `prs_withinfamily*.tsv` |
| `p_v2` | COVAR p in MAGMA `SCZ_vs_<pheno>.gsa.out`, EUR arm |
| `prs_evidence` | −log10(`p_between`) |

Now nine rows, sorted by h² z:

| phenotype | h² z | \|rg\| | PRS evidence (−log10 p) |
|:---|---:|---:|---:|
| **`slope_topDelta`** | **0.30** | 0.072 | **2.26** |
| **`slope_topC3`** | **0.40** | 0.309 * | 1.56 |
| **`slope_projC3`** | **0.62** | 0.194 | 0.16 |
| `slope_PC1` | 0.66 | 0.094 | 1.07 |
| **`slope_projDelta`** | **0.66** | 0.061 | 0.56 |
| `global_slope` | 1.16 | 0.149 | 2.41 |
| `slope_PC2` | 3.46 | 0.052 | 0.63 |
| `baseline_thickness` | 4.49 | 0.026 | 0.24 |

\* SE 0.453 — uninformative, see §7.

**The anti-ranking weakens substantially out of sample.** `SETUP_CONTEXT.md`
§4(c) reported Spearman(h² z, PRS −log10 p) = **−0.80** over the four
positively-heritable phenotypes, while stating correctly that at n = 4 this is
not a hypothesis test. Adding Experiment A's four:

| statistic | original 4 | all 8 |
|:---|---:|---:|
| Spearman(h² z, \|rg\|) | −0.800 | −0.647 |
| Spearman(h² z, PRS −log10 p) | **−0.800** | **−0.287** |

The rg version roughly holds; **the PRS version — the one that matters, since
§7 shows the rg column is noise — largely collapses.** The two projections
break it: they are *both* low-h² *and* low-association, which the anti-ranking
does not allow. The subset means do fit the pattern (lowest h² z, high PRS
evidence), but they fit it by having low h² while merely *matching*
`global_slope`'s association, not exceeding it.

**Reading.** The anti-ranking was always a description of five points, and its
real content was a warning: do not optimise heritability, because in these data
heritability and disorder-relevance are not aligned. That warning survives —
`baseline_thickness` still has 15× the h² z and a tenth the PRS evidence. What
does not survive is any suggestion that *low* heritability predicts relevance.
It does not: `slope_projC3` has h² z = 0.62 and the weakest association in the
table.

## 10. Step 7 — MAGMA, both directions

`magma_eur_v3/`, MAGMA run through v1's own `04_magma.sbatch` with
`MAGMA_DIR`/`GWAS_DIR` redirected. Written to a fresh directory (§4.2).

**Reverse direction — disorder gene Z on phenotype gene Z** (the project's
actual hypothesis), 18,261 genes for SCZ / 18,228 for MDD:

| disorder | phenotype | BETA | SE | P |
|:---|:---|---:|---:|---:|
| SCZ | `slope_topDelta` | +0.0144 | 0.0115 | 0.21 |
| SCZ | **`slope_topC3`** | +0.0180 | 0.0115 | **0.117** |
| SCZ | `slope_projDelta` | +0.0167 | 0.0114 | 0.14 |
| SCZ | `slope_projC3` | +0.0134 | 0.0114 | 0.24 |
| SCZ | *(`global_slope`, v2 reference)* | *+0.0161* | *0.0115* | *0.159* |
| MDD | `slope_topDelta` | −0.0180 | 0.0114 | 0.115 |
| MDD | `slope_topC3` | −0.0013 | 0.0114 | 0.91 |
| MDD | `slope_projDelta` | +0.0039 | 0.0113 | 0.73 |
| MDD | `slope_projC3` | −0.0150 | 0.0113 | 0.18 |

Nothing significant, and nothing meaningfully different from `global_slope`'s
own p = 0.159. `slope_topC3`'s 0.117 is the best of the four and is not a
finding.

**Forward direction — AHBA C1–C3 gene-property**, 6,941 genes, all six
signed half-components as covariates:

Every one of the 24 tests is null (p = 0.157 to 0.990). The row that matters:

| phenotype | C3+ BETA | SE | P |
|:---|---:|---:|---:|
| `slope_topC3` | **+0.0084** | 0.0701 | **0.905** |

**The phenotype constructed from the highest-AHBA-C3 regions shows no C3
enrichment in its own GWAS signal.** Selecting cortex by a transcriptional map
does not make that map's genes show up in the genetics of the selected trait.
That is the cleanest single statement Experiment A produces about the
transcriptional hypothesis at this N.

---

## 11. Verdict — Experiment A's definition of done

`README.md` §2 lists six things. All six are delivered; here they are, with
the answer.

**1. Per new phenotype: λ_GC, GENESIS h² (both GRM thresholds), LDSC h² (+z).**
§6.2, §8, §7. All four scans calibrated (λ 0.995–1.037), **zero genome-wide
hits**, h²_SNP 0.085–0.195 depending on threshold, LDSC h² z **0.30–0.66 —
every one below `global_slope`'s 1.16**.

**2. SCZ/MDD PRS β with ASD/AD controls, population + within-family, paired Δβ
vs `global_slope`, permutation-corrected p.** §5. **No Δβ is negative beyond
noise** (81/320 rows < 0, max |Δβ| 0.0065, all p ≥ 0.41 — see §5.1; none
indicates improvement on the global mean). Subset means indistinguishable from
`global_slope`; projections significantly worse. Min-p permutation:
`global_slope` p_perm = 0.007 (EUR), best new phenotype 0.027. Controls do not
discriminate — ALZ → `slope_topC3` p = 0.023 against SCZ's 0.016.

**3. LDSC rg vs SCZ/MDD iff h² z supports it; otherwise the z and the
statement.** It does not support it. §7 gives the z values and says so, and
records that MDD × `slope_topDelta` crashed LDSC's rg denominator outright.

**4. MAGMA gene-property both directions.** §10. All null. The C3-selected
phenotype shows **no C3 enrichment** (p = 0.905).

**5. Four new rows on the h²-vs-relevance table.** §9. Written.
**The anti-ranking's first out-of-sample test weakens it**: Spearman(h² z,
PRS −log10 p) falls from −0.80 (n = 4) to −0.287 (n = 8).

**6. An explicit statement of what the result licenses.** Below.

### What this licenses, and what it does not

**Licensed.** *Restricting the thinning-rate phenotype to the fastest-thinning
or the highest-AHBA-C3 cortex does not improve any genetic readout.* Not
heritability (GENESIS h² at two thresholds, LDSC h² z), not calibration, not
the disorder association (paired Δβ against `global_slope` in the same
subjects, with a family-cluster bootstrap and a permutation-corrected min-p),
not the transcriptional coupling (MAGMA gene-property), not the reverse
gene-based test. Four independent readouts, one direction.

Stronger, and the result worth carrying forward: **`slope_topC3` is weaker
than ~80 % of random 8-region sets** (§5.2) — consistently across two
disorders, two thresholds and both ancestry strata. The C3 selection is not
merely uninformative, it is worse than arbitrary. Because the selection used
group maps only — no genotypes, no PRS — this is not the circularity that
`hpc/README_HPC.md` §2 diagnosed for a PRS-selected region set. It is a
cleaner version of the same negative.

**NOT licensed.**

- **This is not evidence against the transcriptional hypothesis.** Every test
  above is bounded by the same ceiling §1 of `SETUP_CONTEXT.md` describes:
  `global_slope`'s h² z of 1.16. A phenotype cannot show more disorder signal
  than its heritability can carry, and the new phenotypes have *less*
  heritability. What has been shown is that **region subsetting is not the
  lever** — not that the C3 axis is unrelated to thinning.
- **The subset means are not independent findings.** r = 0.89–0.92 with
  `global_slope`; only the paired difference is a result, and it is null.
- **`slope_topC3`'s rg of −0.309 must never be quoted alone.** SE 0.453.
- **The Δβ nulls do not establish equivalence.** They establish that the
  difference is not detectable at this N, which is a weaker claim.
- **The within-family arm confirms nothing.** ~18 % power; a null β_W is
  uninformative by construction (§5.4).
- **C+T only.** §13.5 of `hpc_v2/README_HPC.md` records that the SCZ →
  `global_slope` C+T association **does not survive PRS-CS** (p = 0.195). Every
  β in §5 lives under that caveat. Since the new phenotypes do not beat
  `global_slope` under C+T, running PRS-CS on them was not worth the compute at
  this stage — but nothing here rehabilitates the C+T result.

### What this means for Experiments B–D

`README.md` §4 reason 3 was: *"If `slope_topC3` carries a stronger SCZ
association than `global_slope`, it becomes the MTAG target."* **It does not,
so it does not.** `global_slope` remains the primary trait, and the A→B→C→D
order paid for itself exactly as argued — this cost roughly one v2-sized pass
and removed a candidate that would otherwise have been a plausible target for
C's 68-region compute.

**The gate in §5.1 of `SETUP_CONTEXT.md` looks worse, not better, after A.**
It asks whether per-trait LDSC h² z supports stable rg. Experiment A adds four
data points at z = 0.30–0.66 — *below* `global_slope`'s 1.16, and the regional
slopes B and C would use will be noisier still, being single regions rather
than 8- or 68-region averages. One of the four already crashed LDSC's rg
denominator. **B's go/no-go should be expected to fail**, and B should be run
to establish that on the record rather than skipped — its rg matrix is cheap
(the 9 sumstats now all exist; only the ~10 pilot regional GWAS remain).

Nothing here changes the conclusion `hpc_v2/README_HPC.md` §13.6 reached:
*every failure traces to `global_slope`'s heritability, and raising that number
is worth more than any further re-analysis of the current phenotype.*
Experiment A is a well-powered demonstration that **region selection is not a
way to raise it.**

### Files

**Committed.** Scripts: `align_export_v3.py`, `prs_paired_delta.py`,
`regen_h2_relevance.py`, `01_prs_ct.sbatch`, `02_ldsc_magma.sbatch`,
`03_prs_paired.sbatch`, `04_collect.sbatch`. Evidence: this file, the updated
`h2_vs_scz_relevance.csv` (nine rows; the five settled ones byte-identical),
and `hpc_v3/prs_tables/` — the eight PRS tables behind §5.

`prs_tables/` sits outside the gitignored `work/` deliberately. Its rows are
aggregates (disorder × threshold × phenotype × stratum: β, SE, p, n), carry no
per-subject values, and are the only record of §5 that survives the scratch
tree. This is the same exception `hpc_v2` already makes for
`results_v2/prs_family/prs_withinfamily.tsv`; `SETUP_CONTEXT.md` §8 forbids
**per-subject** genotype-derived files, which these are not.

Also committed, under `hpc_v2/work/results_v2/` and matching that tree's
existing `*_summary.tsv` rule: the four new rows in each of `assoc/` and
`assoc_eur/`'s `gwas_summary.tsv`, the eight new `nullmodel*/…_null_summary.tsv`,
`ldsc_eur_v3/ldsc_rg_summary.tsv`, `magma_eur_v3/magma_summary.tsv`, and both
`reml_zaitlen_pcrel*/reml_zaitlen_summary.tsv` — regenerated over all nine
phenotypes with `hpc_v2/work/setup/collect_reml_pcrel.sh`, the settled five
reproducing exactly.

Those REML summaries carry a column §8 does not: **pedigree h²**. It is
`0.217–0.301` for the new phenotypes against `global_slope`'s `0.252–0.273`,
i.e. no more familial than the trait they were carved out of — the same
conclusion the SNP component gives. (Two of them return GCTA's boundary
`Pval = 5.0e-01` in the κ ≥ 0.0250 arm; read those LRTs as uninformative.)

**Not committed** (gitignored, per `SETUP_CONTEXT.md` §8): `hpc_v3/work/`,
`hpc_v2/work/pheno_v3/`, the per-chromosome and combined sumstats under
`assoc/` and `assoc_eur/`, and the bulk of `ldsc_eur_v3/` and `magma_eur_v3/`.

### Job ledger

| job | script | shape | elapsed |
|:---|:---|:---|:---|
| 35089435 | `hpc_v2/03_null_model.sbatch` | 4 | 56 s |
| 35089436 | `hpc_v2/04_assoc.sbatch` pooled | 88 | ~2 h wall |
| 35089468 | `03_null_model` EUR | 4 | ~40 s |
| 35089470 | `04_assoc` EUR | 88 | ~2 h wall |
| 35089437, 35089472 | `run_all.sh` collect | 2 | **FAILED** (§6.1) |
| 35089439 | `reml_zaitlen_pcrelate` κ≥0.0250 | 4 | 42 m – 1 h 50 m |
| 35089440 | `reml_zaitlen_pcrelate_2nd` κ≥0.0884 | 4 | ~2 m |
| 35089628 | `hpc_v3/01_prs_ct.sbatch` | 1 | ~9 m |
| 35089862 | `hpc_v3/03_prs_paired.sbatch` | 1 | ~11 m |
| 35089930 | `02_ldsc_magma` (first attempt) | 1 | **cancelled** (§4.2) |
| 35093677 | `hpc_v3/04_collect.sbatch` | 1 | 20 m 32 s |
| 35093678 | `hpc_v3/02_ldsc_magma.sbatch` | 1 | 3 h 32 m |

---

## 12. MAGMA locus pools, and the collection slide — 2026-09-09

Two additions after §11, both driven by the slide
(`hpc_v3/slide_results_v3.py`, `make_results_v3.py`).

### 12.1 The slide's PRS source was the wrong arm

`make_results_v3.py` read `hpc/work/results/prs/prs_association.tsv` — the
**array-genotype** arm (n = 3,725 EUR / 4,126 pooled). Every published PRS
number in this project comes from the **imputed** arm. The difference is not
cosmetic: `global_slope` × SCZ × EUR is −0.0385 (p = 0.015) under the array arm
and −0.0468 (p = 0.0022) under the imputed one. Repointed at
`hpc_v3/prs_tables/prs_association_v3.tsv`.

That table was the only committed *imputed* PRS association in the repo —
`hpc/work/results/prs_imp/` was gitignored, because `.gitignore` names
`prs/prs_association.tsv` as a one-off exception and `prs_imp/` was never
added. It is now (one line), so v1's own imputed PRS is reproducible from the
repo rather than only from `$V1_ROOT`.

### 12.2 The prioritised gene-set panel did not exist for the new phenotypes

`02_ldsc_magma.sbatch` produced the `.genes.raw` files but not the
**prioritised / high-confidence set tests** (`SCZ_locus_pool`, `MDD_pool`),
which live in a separate harness (`hpc/work/prioritised_gsa.sbatch`, invoked
for the settled five by `hpc_v2/work/setup/eur_prio_gsa.sbatch`). This is the
panel carrying the project's one surviving strand (`hpc_v2/README_HPC.md`
§13.6), so it mattered.

**Job 35110779** (`hpc_v3/05_prio_gsa.sbatch`, 3 m 17 s), into
`magma_prio_eur_v3/` — a fresh directory for the §4.2 reason: the collector
globs `$OUT/*_$TAG.gsa.out` and rewrites `<TAG>_gsa_summary.tsv` from what it
finds, so pointing it at `magma_prio_eur/` would have rewritten a published v2
table.

**One documented pipeline modification** (permitted by `README.md` §2, which
says to document it here): `prioritised_gsa.sbatch` had its phenotype list
hard-coded as `for name in baseline_thickness global_slope slope_PC1 ...`. It
is now `NAMES=${NAMES:-<that exact list>}`, in the same spirit as the
`REPO`/`MAGMA_DIR`/`SETS`/`TAG`/`COND_SET` overrides the script already had.
Every existing caller is unaffected.

**Results** (marginal model, EUR arm, gene-set β ± SE):

| phenotype | SCZ locus pool | p | MDD pool | p |
|:---|---:|---:|---:|---:|
| `baseline_thickness` | 0.1874 ± 0.051 | **1.1e-04** | 0.0209 ± 0.026 | 0.21 |
| `global_slope` | 0.1244 ± 0.051 | **0.0074** | 0.0293 ± 0.026 | 0.13 |
| **`slope_topDelta`** | **0.1253 ± 0.051** | **0.0069** | **0.0502 ± 0.026** | **0.027** |
| **`slope_topC3`** | 0.0918 ± 0.051 | 0.035 | 0.0242 ± 0.026 | 0.18 |
| **`slope_projDelta`** | 0.0077 ± 0.051 | 0.44 | 0.0402 ± 0.026 | 0.064 |
| **`slope_projC3`** | 0.0208 ± 0.051 | 0.34 | −0.0299 ± 0.026 | 0.87 |

**This is the one readout on which a new phenotype does not simply lose.**

- **SCZ locus pool: `slope_topDelta` matches `global_slope` almost exactly**
  (β 0.1253 vs 0.1244, p 0.0069 vs 0.0074). It ties; it does not beat it.
  `slope_topC3` is *lower* (0.0918, p = 0.035), consistent with everything in
  §5 and §7.
- **MDD pool: `slope_topDelta` is β = 0.0502, p = 0.027, against
  `global_slope`'s 0.0293, p = 0.13** — nominally significant where the anchor
  is not, and the only place in this experiment where a subset phenotype
  exceeds the trait it was carved from.

**Do not promote that MDD result.** Four reasons, and they are the same
reasons the rest of this document gives for not promoting anything:

1. It is **uncorrected across 12 gene-set tests** (6 phenotypes × 2 disorders)
   in this panel alone; 0.027 × 12 = 0.33.
2. `hpc_v2/README_HPC.md` §12.9 already records that **MDD gene-set results
   weakened by 2–4× in p** between v1 and v2 — this is the arm with the known
   instability.
3. The **MAGMA test is unsigned**, so it cannot say thinning is *faster* with
   higher MDD risk. §5's signed PRS Δβ for `slope_topDelta` × MDD is
   **−0.0001 (p = 0.99)** — no difference from `global_slope` whatsoever.
4. `slope_topDelta`'s **LDSC h² z is 0.30**, the lowest of the nine.

Taken with §5 and §7, the honest summary is unchanged: **`slope_topDelta` ties
`global_slope` on gene-set enrichment and loses or ties everywhere else.**
The §11 verdict stands.

### 12.3 The slide

`slide_results_v3.png`, regenerated. Eight rows: the two anchors and the four
Experiment A phenotypes through v2, plus the two anchors through v1's own
pipeline, so the pipelines sit side by side.

Encoding, and the reason there is no fourth channel: **colour = disorder**,
**fill + shape = ancestry stratum**, **row = phenotype × pipeline**. Pipeline
is deliberately not a colour or a fill — both are taken, and a third
overlapping channel is how a forest plot stops being readable. The two-hue
palette was checked with the `dataviz` skill's validator (light mode, white
surface): lightness band, chroma floor, CVD separation (adjacent-pair ΔE 12.9
deutan / 10.7 tritan), normal-vision floor (ΔE 19.9) and contrast all **PASS**.
Every row is directly labelled, so identity is never colour-alone.

Three things the first render got wrong, fixed:

- **h² drew with the marker the legend defines as "EUR stratum."** Both REML
  arms are the **pooled** multi-ancestry sample (n = 8,082 v2 / 5,649 v1), so
  the rows now carry the pooled marker.
- **v1's pooled PRS row duplicated v2's exactly.** A polygenic score never
  touches our GWAS, so for the pooled stratum the two pipelines are the *same
  estimate on the same 8,082 subjects* — drawn twice it reads as a replication
  and is only a duplicate. `make_results_v3.py` now **asserts** the two agree
  (to 1e-9; they do) and draws v1 only for EUR, which is what genuinely
  differs: a PC-distance cut (n = 5,361) against the anchor set (n = 4,116).
- **`p=%.3f` rounded 0.0014 to "0.001"**, an order of magnitude better than it
  is. Now two significant figures.

`results_v3_summary.csv` also gains a `pipeline` column, a `caveat` column (so
the figure prints the two conditions under which the v1/v2 rows may be read as
like-for-like, rather than the caveats living only in prose), and a
`not_estimable` status distinct from `pending` — MDD × `slope_topDelta` rg is
not a missing input, it is LDSC declining to divide by a heritability
indistinguishable from zero (§7), and the slide draws it as an × rather than a
grey placeholder.
