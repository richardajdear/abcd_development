# Genes driving adolescent cortical development: pipeline rebuild and first results on ABCD 5.1

**Status.** Refactor complete and tested on release 5.1. Modelling framework
validated. One hypothesis test run to completion with a negative result and a
positive control that passes. Genetic analysis written but not executed
(cluster authentication unavailable). 79 tests pass.

---

## 1. What was asked and what was built

The scientific goal is to identify genes driving adolescent cortical
development, under two hypotheses: that such genes are enriched for SCZ and MDD
GWAS signal, and that they connect to two transcriptional axes we already have
— AHBA component C3 (Dear et al., *Nat Neurosci* 2024) and PC1 of the
single-nucleus maturation data.

The starting point was `thesis_abcd.qmd` from the 5.1 analysis. It has been
replaced by a modular pipeline in
[`transcriptional_maturation`](https://github.com/richardajdear/transcriptional_maturation),
branch `abcd-longitudinal`. The design brief was: config-driven so a 7.0 adapter
drops in without touching analysis code; metric-agnostic so functional and
myelin-proxy metrics run through the same path; and R for the statistics where R
is better.

### Language recommendation, as requested

**Both, split on a clean seam.** Python for data assembly, QC, spatial
statistics and gene work (`src/abcd/`); R for model fitting (`R/`). The seam is
Parquet: Python writes a tidy long table, R fits and writes fixed effects and
variance components back as Parquet, Python reads them for phenotypes and maps.

The reason is not preference. `lme4` handles crossed random effects with
correlated slopes, boundary-singular fits, and REML properly, and reports
singularity honestly; `statsmodels`' mixed-model support does not cover crossed
effects at all. Conversely Python's spatial-null and gene-annotation ecosystem
is better, and the assembly code is ordinary dataframe work where pandas is
fine. Nothing crosses the seam except Parquet files, so neither side needs to
know how the other works.

---

## 2. Critique of the 5.1 analysis

Reading `thesis_abcd.qmd` first, since several findings below are consequences
of fixing what it did.

**What it got right.** The scientific framing — longitudinal change as the
imaging-genetics phenotype rather than baseline — is the right and novel move,
and the reason this project is worth doing. The AHBA/MAGMA infrastructure is
sound and reused here.

**Problems found, in descending order of consequence.**

1. **No family random effect, and no analysis of what that implies.** ABCD has
   ~5,989 families for 6,937 subjects including ~800 MZ pairs. The old model
   ignored family entirely. Including it turns out to be *wrong for phenotype
   extraction* but for a subtle reason that had not been worked out either way —
   see §4. Getting this wrong in either direction invalidates the heritability
   analysis, so it needed resolving rather than defaulting.
2. **Difference scores rather than modelled slopes** in parts of the analysis.
   A two-timepoint difference has no residual degrees of freedom, so it cannot
   report its own reliability, and its variance is dominated by measurement
   error. This turned out to matter enormously: slope reliability on 5.1 is
   ~0.15, meaning ~85% of apparent between-subject variation in rate is noise.
   The old analysis had no way to know that.
3. **Naive permutation for spatial inference.** Permuting region labels destroys
   spatial autocorrelation and so builds a null of rough maps that no real map
   resembles. Measured false-positive rate at nominal 0.05: **0.24 to 0.58**
   depending on map smoothness (§5). Several borderline results in the old
   analysis would not survive a proper null.
4. **Gene projection through three components only.** Mapping a cortical map to
   genes via C1/C2/C3 loadings confines the answer to a 3-dimensional subspace
   of gene space no matter what the map is. Any negative result obtained that way
   is uninformative (§6).
5. **Global-mean covariate not decomposed.** A single whole-cortex covariate
   conflates "this subject has thin cortex" with "this scan measured thin",
   which are respectively biology and motion artefact.
6. **Release-specific paths and column names inline.** Would have required
   rewriting analysis code for 7.0.

---

## 3. Pipeline architecture

```
src/abcd/
  config.py       RunConfig: every choice is a field; hash -> run_id
  paths.py        release-agnostic path resolution (ABCD_ROOT, ABCD_AHBA_DIR)
  io.py           release adapters (5.1 now, 7.0 later) + parcellations
  qc.py           QC as code, with a ledger of every exclusion
  assemble.py     -> out/<run_id>/model_table.parquet  (tidy long)
  phenotype.py    BLUPs + per-subject reliability -> phenotypes.parquet
  spatial.py      spin test (bijective), variogram surrogates, symmetry
  genemaps.py     AHBA components, full expression matrix, gene omnibus
  gcta_export.py  GCTA/PLINK-format phenotypes and covariates
  magma_export.py MAGMA gene-covariate files
R/
  model_spec.R    formulas BUILT from config -- no hand-written duplicates
  fit_lmm.R       per-region fits, parallel -> fixed/varcomp/diagnostics
  compare_models.R nested model ladder
hpc/              GRM -> REML -> fastGWA -> MAGMA (written, not run)
notebooks/        01 longitudinal model (R), 02 maps and genes (Python)
```

Two properties worth stating because they were design goals rather than
accidents:

**Config-driven and release-agnostic.** Everything that could change is a
`RunConfig` field: release, metric, parcellation, age basis, random-effects
structure, whether the global covariate is included, whether family is modelled.
The config hashes to a `run_id`, so two runs differing in any choice cannot
overwrite each other, and `out/<run_id>/config.yaml` records exactly what
produced the numbers. Adding 7.0 means writing one adapter class in `io.py`.

**Metric-agnostic.** `metric: thickness` is a config field; the T1w/T2w ratio
already runs through the same path (`out/t1t2_ratio_dsk_51_*`). Adding a
functional metric means adding a table mapping, not new analysis code.

### Sample flow

| stage | scans | subjects |
|:---|---:|---:|
| imaging data present | 22,854 | 11,802 |
| passed structural QC (Euler + non-Philips) | 20,162 | 11,162 |
| complete covariates | 20,162 | 11,162 |
| **≥2 visits (analysis sample)** | **15,937** | **6,937** |
| 3 visits (subset) | 6,189 | 2,063 |

Every exclusion is written to `out/<run_id>/qc_ledger.csv` with its reason.

---

## 4. The longitudinal model

Full exposition in `notebooks/01_longitudinal_model.qmd`. The specification:

```
value ~ age_c + sex + global_between_c + global_within
        + (1 + age_c | subject) + (1 | site)   [+ (1 | family_id)]
```

Model comparison over a nested ladder (ΔAIC vs the best model, DK parcellation,
`lh_precentral`):

| specification | ΔAIC | slope reliability |
|:---|---:|---:|
| 1 intercept only | 4462.2 | — |
| 2 slope, uncorrelated | 4449.3 | 0.56 |
| 3 slope, correlated | 4427.7 | 0.83 |
| 4 + site | 4209.4 | 0.84 |
| 5 + site + family | 4105.1 | 0.85 |
| **6 + global (split within/between)** | **26.0** | — |
| 7 quadratic age | 7.3 | 0.999 |
| 8 spline age (df = 3) | 7.3 | 0.999 |
| 9 sex × age | 16.2 | 0.999 |

Reading it: the correlated random slope is strongly preferred over uncorrelated
(ΔAIC 21) and over intercept-only (ΔAIC 34). The **split global covariate is by
far the largest single improvement** (4105 → 26), confirming that shared
scan-level variance was dominating the residual. Quadratic and spline age are
marginally preferred over linear (ΔAIC ~19 and ~19 relative to model 6) but
their fitted trajectories are nearly identical to linear over 5.1's short span
(slope correlation 0.999), so **linear is retained** — the added flexibility
changes the phenotype not at all while making it harder to interpret. This
decision should be revisited for 7.0, where the span roughly doubles.

### Variance decomposition

Across 68 DK regions, median share of variance:

| component | median share |
|:---|---:|
| subject intercept | 53.5% |
| family | 18.2% |
| site | 4.8% |
| **subject slope** | **0.2%** |
| residual | 20.2% |

The subject slope carries a fraction of a percent of total variance. This is the
central quantitative fact about the project. It is not a modelling failure —
thickness at a visit is mostly *who you are*, and two years of change is a small
perturbation on that — but it sets the ceiling on everything downstream.

### The family random effect is a phenotype decision {#family}

This is the most consequential finding of the refactor, and it is the opposite
of what I expected.

A family random intercept absorbs variance shared between relatives. Variance
shared between relatives *is* the genetic signal. So subject BLUPs from a model
including family are deviations from the family mean — a phenotype with
heritability deliberately removed.

Twin correlations on the whole-cortex mean (207 MZ, 617 DZ/sibling pairs):

| phenotype | family RE | MZ *r* | DZ/sib *r* | Falconer *h²* |
|:---|:---|---:|---:|---:|
| intercept | included | 0.14 | **−0.62** | — |
| intercept | omitted | 0.88 | 0.47 | **0.83** |
| slope | included | — | — | 0.54 |
| slope | omitted | — | — | 0.58 |

A DZ correlation of −0.62 is biologically impossible; it is the arithmetic
signature of subtracting a family mean from two siblings. Without the family
term, intercept *h²* = 0.83, matching published twin heritability of cortical
thickness. Intercept reliability also rises from ~0.65 to 0.88–0.91.

The pipeline therefore ships two configs that produce *different phenotypes on
purpose*: `ct_genetics.yaml` (no family term — use for GWAS) and
`ct_baseline.yaml` (family term — use for descriptive inference, where omitting
it understates standard errors).

### Reliability, and why 7.0 is the answer

Slope reliability is $\tau_s^2/(\tau_s^2 + \sigma^2/S_{xx})$ where $S_{xx}$ is
the within-subject spread of age. Measured medians across DK regions: **0.12
(2 visits), 0.16 (3 visits)**; means 0.15 and 0.19. Intercept: 0.88 / 0.91.

Two consequences follow from the formula and both are counterintuitive:

- **Span beats visit count.** $S_{xx}$ grows quadratically with follow-up span,
  only linearly with number of visits. Three visits in two years are worth less
  than two visits eight years apart.
- **Averaging regions does not help.** The dominant noise is per-occasion
  (motion, positioning) and therefore shared across regions. The whole-cortex
  mean slope has *lower* reliability (0.12) than the better individual regions.

At *r* ≈ 0.15, observed effect sizes are attenuated by $\sqrt{r}$ ≈ 0.39 and the
sample size for equal power scales as $1/r$ ≈ 6.7×. **A slope GWAS on 5.1 is a
pipeline test, not a discovery attempt.** The first thing to check on 7.0 is
whether reliability rises as the $S_{xx}$ expression predicts; that single
number determines whether the genetic analysis is viable.

---

## 5. Spatial inference

`src/abcd/spatial.py`. Nulls calibrated on synthetic smooth maps with
known-zero true correlation; entries are false-positive rates at nominal 0.05:

| map smoothness | naive label permutation | variogram surrogates | spin test |
|---:|---:|---:|---:|
| 0.10 | 0.235 | 0.180 | 0.165 |
| 0.20 | 0.270 | 0.170 | 0.105 |
| 0.35 | 0.430 | 0.200 | 0.115 |
| 0.60 | 0.575 | 0.325 | **0.090** |

The naive null degrades badly with smoothness — at 0.60 it is wrong more often
than right. The spin test stays at or below nominal. Two implementation points:

- The rotation-to-parcel assignment is solved as a **bijection**
  (`linear_sum_assignment`), not by independent nearest-neighbour lookup. The
  common shortcut reuses some parcels and drops others, which is not a
  permutation and distorts the null's variance. Verified: 180/180 unique targets
  per rotation.
- Parcel centroids for DK and HCP-MMP are generated once from the fsaverage
  annotation and committed (`data/*_centroids.csv`), so downstream analysis
  needs no surface software. Label ordering was validated empirically, not
  trusted: HCP parcels were confirmed spatially contiguous on the
  vertex subset (dispersion ratio 0.121 vs random) and all 137 AHBA regions map
  onto the imaging parcellation.

---

## 6. Linking maps to genes: a negative result, and one provisional positive

Full detail in `notebooks/02_maps_and_genes.qmd`.

### The map is reliable, so a null result means something

Split-half over subjects for the HCP-space developmental change map: **ρ =
0.975**. Left–right symmetry 0.935. The map is essentially noise-free.

### Map-level test against AHBA components

| map | component | ρ | p_spin | p_naive |
|:---|:---|---:|---:|---:|
| **baseline thickness** | **C1** | **−0.58** | **0.001** | 0.001 |
| baseline thickness | C2 | 0.17 | 0.270 | 0.041 |
| absolute change | C1 | −0.14 | 0.566 | 0.093 |
| absolute change | C3 | −0.18 | 0.209 | 0.048 |
| conditioned change | C2 | 0.22 | 0.136 | 0.017 |

Note the `p_naive` column: three rows would have been "significant" under the
naive null and are not under the spin null. That is §5 in action.

The developmental change maps show no reliable relationship with C1–C3. The
positive control does (baseline vs C1, ρ = −0.58, p = 0.001), which establishes
that the machinery detects real map-to-map structure.

### Removing the projection bottleneck

The gene-level test used to route map → C1/C2/C3 → loadings → gene ranking. The
loading vectors are near-orthogonal (pairwise |cos| ≤ 0.104), so they span a
genuine 3-dimensional subspace — and every map's gene ranking lies inside it,
determined by three numbers regardless of the map. A negative result obtained
that way says nothing.

The fix used the full expression matrix
(`AHBA/data/abagen-data/expression/hcp_3d_ds5.csv`, 137 regions × 7,973 genes,
your own optimal processing). Its rows are keyed by integer code with no names,
so **row order was validated rather than assumed**: components reconstructed as
z(expression) @ loadings must track the published scores (ρ = 0.998 / 0.974 /
0.912), and `genemaps._assert_row_order` raises otherwise. Row misalignment is
the failure mode that invalidates every number downstream while raising no
error.

### The omnibus test, and why per-gene counts mislead

The raw count of genes at p < 0.05 for the change map was 800 against 399
expected — apparently 2× enrichment. It is not. Genes share expression
gradients, so one alignment with a dominant gradient produces thousands of
"significant" genes. The correct null is over *maps*: rotate the map, recompute
the whole gene correlation vector, compare mean |ρ| to its rotation
distribution. Gene–gene dependence is then identical in observed and null.

| map | mean \|ρ\| | null | p |
|:---|---:|---:|---:|
| **baseline thickness** | 0.318 | 0.098 | **0.001** |
| absolute change | 0.123 | 0.125 | 0.443 |
| conditioned change | 0.121 | 0.095 | 0.169 |
| slope SD (individual differences) | 0.100 | 0.137 | 0.817 |
| **slope SD / residual SD** | 0.188 | 0.083 | **0.001** |
| residual SD (noise) alone | 0.237 | 0.171 | 0.118 |
| parcel size (confound) alone | — | — | 0.177 |

**The group-mean developmental change map is not gene-linked.** With the full
matrix, a reliable map, and a valid null, this is a result to trust rather than a
failure to detect.

**Regional slope reliability (slope SD / residual SD) is gene-linked**, and it
survived every confound test constructed: p = 0.009 partialling residual noise,
0.018 partialling parcel size, 0.003 partialling baseline thickness, and neither
confound is itself gene-linked. **But** that map's left–right symmetry is 0.35
where real cortical maps here exceed 0.9, which suggests it is substantially
estimation noise in the variance components — and variance components are
precisely what is estimated worst when most subjects have two visits. **Treat as
provisional pending 7.0.**

### What this does and does not say about the SCZ/MDD hypothesis

It does not test it. Where cortex thins *on average* is a group-level anatomical
fact; whether disorder risk genes are enriched in the genetic architecture of
*individual differences* in development is a question about association signal.
That route is GWAS of the subject-level phenotype followed by competitive
gene-set testing — the `hpc/` pipeline — and it is independent of everything in
this section. The negative spatial result does not bear on it.

---

## 7. Genetic analysis: written, not run

CSD3 refused key authentication throughout (`Permission denied
(publickey,keyboard-interactive,hostbased)`, with an RCS storage recovery banner).
Per instruction, `hpc/` contains complete runnable scripts rather than results:
`00_check_inputs.sh` (preflight), `01_grm.sbatch`, `02_reml.sbatch`,
`03_gwas.sbatch`, `04_magma.sbatch`, chained with sbatch dependencies. Account
`VERTES-SL2-CPU`, partition `cclake`, paths from the previous analysis.

Decisions documented there because they change the answer, not the runtime:

- Use `ct_genetics.yaml`, not `ct_baseline.yaml` — see [§4](#family).
- Both unrelated (`--grm-cutoff 0.05`) and full-GRM REML are run; if
  *h²*(full) ≫ *h²*(unrelated), shared environment is inflating the estimate.
  fastGWA retains related subjects via a sparse GRM rather than discarding ~⅓ of
  the cohort.
- Gene sets go in as **continuous covariates** (`--gene-covar`), not hard sets,
  so there is no arbitrary cutoff; each signed component is split into ± halves
  because a component's two poles are different biology and one signed covariate
  averages opposing enrichments toward zero.
- Expected power is stated up front so the first run is not misread.

Export modules are verified locally: 6,937 subjects export in GCTA format with
IDs normalised to the genetics convention (ABCD writes `NDAR_INV…` in genetics
tables and `sub-NDARINV…` in imaging tables; the preflight checks the `.fam`
intersection because that join fails silently). Gene covariates: 7,643 AHBA
genes with zero identifier loss, 17,618 snRNA-seq PC1 genes (99.9% mapped via
NCBI `gene_info`, synonyms used only where the primary symbol failed).

---

## 8. Honest summary of where this leaves the hypothesis

**Established.** A clean, tested, config-driven pipeline that runs end to end on
5.1 and needs one adapter class for 7.0. A phenotype whose reliability is known
rather than assumed. Spatial and gene-level inference with calibrated nulls and a
passing positive control.

**Resolved by measurement, against expectation.** The family random effect must
be *omitted* for genetic phenotypes. Splitting the global covariate matters more
than any other single specification choice. Follow-up span, not visit count,
drives slope reliability.

**Negative, and trustworthy.** The group-mean cortical development map is not
transcriptionally patterned along AHBA C1–C3, at either map or gene level, with
a reliable map and a valid null.

**Provisional.** Regional slope reliability is gene-linked and survives confound
partialling, but fails a symmetry check that suggests noise dominance.

**Blocked.** SCZ/MDD enrichment. Needs the cluster.

**The binding constraint is reliability, not sample size.** At slope reliability
0.15, release 5.1 cannot answer the question, and no amount of additional
subjects at two visits each would fix it. Release 7.0's added timepoints attack
exactly the term ($S_{xx}$) that the algebra says matters. The first analysis on
7.0 should be the two-line reliability calculation, because it determines
whether the rest is worth running.

---

## 9. Next steps, in order

1. **Wire the 7.0 adapter** — one class in `io.py`; the config system and
   everything downstream is unchanged.
2. **Recompute slope reliability on 7.0** before anything else. Prediction:
   roughly triples with four visits over eight years. This is the go/no-go.
3. **Revisit linear vs spline age** over the longer span (`R/compare_models.R`).
   If the trajectory curves, the slope phenotype must be redefined.
4. **Re-test the provisional slope-reliability gene result**, checking whether
   its left–right symmetry rises with better variance-component estimates.
5. **Run `hpc/`** once cluster access returns: intercept first as a positive
   control (expect *h²* ≈ 0.8), then slope.
6. **Consider metrics with better scan-rescan reliability** — the T1w/T2w ratio
   path already runs. Given §4, reliability matters more than effect size.

## Reproducing

```bash
export ABCD_ROOT=/path/to/ABCD ABCD_AHBA_DIR=/path/to/AHBA/data/abagen-data/expression
python -m abcd.assemble configs/ct_genetics.yaml
Rscript R/fit_lmm.R --run-dir out/<run_id> --cores 8
python -m abcd.phenotype out/<run_id>
python -m pytest tests/ -q                     # 79 tests
python tools/check_notebook_chunks.py notebooks/01_longitudinal_model.qmd
```

Figures in `figures/`: sample flow, model comparison, model diagnostics,
phenotype reliability, spatial nulls, gene-level maps.
