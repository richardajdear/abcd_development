# Adolescent cortical development in ABCD 7.0

**Status.** Modelling framework settled and validated on release 7.0. Phenotype
selection for imaging genetics complete. Genetic analyses (GCTA GRM-based
heritability, GWAS) specified and ready to dispatch to the HPC cluster.

Supersedes `REPORT_5.1_legacy.md` (the 5.1 draft) and all earlier versions of
this document. **Two substantive conclusions in the previous version of this
report are overturned, and two reported numbers were wrong** — see §11 for all
four corrections and why they happened.

---

## Contents

1. [Scientific question and the logic of the design](#1-scientific-question-and-the-logic-of-the-design)
2. [Model specification](#2-model-specification)
3. [The framework reproduces across releases](#3-the-framework-reproduces-across-releases)
4. [Design decisions, and what settled them](#4-design-decisions-and-what-settled-them)
   - 4.1 [The global covariate must not be included](#41-the-global-covariate-must-not-be-included)
   - 4.2 [The family random effect destroys the genetic signal](#42-the-family-random-effect-destroys-the-genetic-signal)
   - 4.3 [`min_visits: 3` costs power; use ≥2](#43-min_visits-3-costs-power-use-2)
5. [The developmental phenotype: regional maps](#5-the-developmental-phenotype-regional-maps)
6. [Heritability](#6-heritability)
7. [Transcriptional comparison: AHBA components](#7-transcriptional-comparison-ahba-components)
8. [Structural covariance of developmental slopes](#8-structural-covariance-of-developmental-slopes)
9. [Site and scanner effects](#9-site-and-scanner-effects)
10. [Phenotype priority for GWAS](#10-phenotype-priority-for-gwas)
11. [Corrections to the previous report](#11-corrections-to-the-previous-report)
12. [Limitations](#12-limitations)
13. [What runs next, on the cluster](#13-what-runs-next-on-the-cluster)
14. [Reproducing this report](#14-reproducing-this-report)

---

## 1. Scientific question and the logic of the design

The goal is to identify genes driving adolescent cortical development. The
hypothesis is that such genes are (i) enriched for GWAS signal from disorders
with adolescent onset — schizophrenia and major depression — and (ii) connected
to prior transcriptional results, specifically AHBA component C3 from earlier
work in this group and the leading component of an independent snRNA-seq
analysis.

The design choice that distinguishes this work from prior imaging genetics is
the phenotype. Previous studies have almost universally run GWAS on
*cross-sectional* brain measures — thickness or volume at one timepoint. Here
the phenotype is the *rate of change*: a per-subject slope of cortical thickness
on age, estimated from repeated scans. Release 7.0 is what makes this feasible;
it supplies enough subjects with three or more imaging visits that
subject-specific slopes carry usable signal (§4.3).

This matters because the genetic architecture of a developmental *rate* need not
resemble that of a static measure. If adolescent cortical thinning is under
genetic control, the relevant variants act on the trajectory, and a
cross-sectional phenotype averages over exactly the variation of interest.

---

## 2. Model specification

For each of 68 Desikan-Killiany cortical regions, thickness is modelled across
visits with a linear mixed model fitted in R (`lme4`), with random intercepts
*and* random slopes per subject:

```
value ~ age_c + sex + (age_c | subject)
```

`age_c` is age in years centred on the sample mean. The random-effect term
`(age_c | subject)` yields per-subject BLUPs for both intercept (thickness at
mean age) and slope (annual rate of change, mm/yr). The slope BLUPs are the
developmental phenotype; everything downstream is a function of them.

Three specification choices are deliberate and each is defended empirically in
§4: no global-thickness covariate, no family random effect, and a minimum of two
visits per subject. The settled configuration is
`configs/ct_70_noglobal_mv2_genetic.yaml`.

Modelling runs in R because `lme4`'s handling of correlated random effects and
its BLUP extraction are more dependable than the Python equivalents; all
downstream spatial, transcriptional and genetic analysis is Python. The two
communicate through parquet files in `out/<run_id>/`.

---

## 3. The framework reproduces across releases

Group-level developmental maps — the regional mean of subject slopes, 68 DK
regions — are compared between releases on the *same* specification
(`global_covariate: none`, `min_visits: 2`), so that a release difference is
not confounded with the specification change of §4.1.

| comparison | Spearman ρ | Pearson r |
|---|---|---|
| release 5.1 vs 7.0, matched specification | 0.937 | 0.975 |
| visit filter within 7.0: ≥2 vs ≥3 | 0.998 | 0.999 |
| family effect within 7.0: on vs off | **1.000** | **1.000** |

The release upgrade and the adapter rewrite did not disturb the cortical
pattern, so downstream differences are attributable to sample composition rather
than to pipeline changes. Mean cortical rate is −0.0185 mm/yr on 5.1 and −0.0189
on 7.0.

Two notes on reading this table. First, these are *group-level* maps — the
`age_c` fixed effect per region — not subject BLUPs, whose group mean is zero by
construction. Second, the family effect leaves the group map exactly unchanged
(ρ = 1.000) while destroying the genetic signal (§4.2): a specification can be
invisible at the group level and fatal at the subject level, which is why the
family-effect problem went undetected in the earlier draft.

An earlier version of this report gave the release correlation as 0.977 without
naming the coefficient. The Pearson value is 0.975 and the Spearman value is
0.937; the table above reports both, computed on matched specifications.

Mean whole-cortex thinning is −0.019 mm/yr on 7.0 under the settled
specification.

---

## 4. Design decisions, and what settled them

Three specification questions turned out to matter more than any analysis
choice downstream. Each was settled by a measurement, and in two cases the
measurement contradicted the earlier draft.

### 4.1 The global covariate must not be included

The 5.1 draft conditioned each regional model on the subject's mean thickness
across the cortex, intending to isolate regionally specific development. This
is the wrong choice for a genetic phenotype, and the cost is large.

Conditioning on the global mean removes the component of development that is
shared across the cortex — but that shared component is where most of the
heritable signal lives. Removing it does not sharpen the regional signal; it
deletes the signal.

| specification | h² (global mean slope) | 95% CI |
|---|---|---|
| **no global covariate** *(settled)* | **0.584** | 0.354 – 0.821 |
| global mean as covariate | 0.385 | 0.129 – 0.634 |

The same comparison at `min_visits: 3` gives 0.547 without versus 0.221 with.
Conditioning on the global mean costs roughly a third of the heritability at
≥2 visits and more than half at ≥3.

This decision also changes what the phenotype *is*, and hence its relationship
to everything downstream. The transcriptional comparison in §7 only works on the
unconditioned phenotype — which is the substance of the first correction in §11.

### 4.2 The family random effect destroys the genetic signal

**Any run intended for heritability or GWAS must set `family_effect: false`.**

Adding `(1 | family_id)` partitions between-family variance into its own random
effect. The subject-level BLUPs that remain are *within-family deviations* — and
the between-family component, which is exactly what genetic relatedness
explains, has been removed. With two phenotyped members per family the family
effect centres them at zero, so their deviations must sum to zero and
anti-correlate by construction.

Measured on 7.0 with genotype-confirmed zygosity (260 MZ, 871 DZ-or-sibling
pairs), whole-cortex thickness:

| phenotype | `(1\|family_id)` | r(MZ) | r(DZ) | Falconer h² |
|---|---|---|---|---|
| baseline thickness *(positive control)* | with | 0.581 | **−0.173** | **1.508** |
| baseline thickness *(positive control)* | without | 0.876 | 0.478 | 0.796 |
| developmental slope *(target)* | with | 0.474 | 0.148 | 0.651 |
| developmental slope *(target)* | without | 0.502 | 0.210 | 0.584 |

The positive control is what exposes the problem. Cortical thickness is strongly
familial, so a *negative* DZ correlation is not a weak result but a structurally
impossible one, and the h² it implies — 1.508 — is outside the parameter space.
Note the trap: on the *target* phenotype the same bug returns h² = 0.651, which
is merely high rather than obviously impossible. Had we looked only at the
slope, the artefact would have read as a positive finding.

`gcta_export.py` therefore guards on the model specification
(`FamilyEffectConflict`) rather than on the h² value, with regression tests in
`tests/test_gcta_export.py`.

Relatedness must be handled where it belongs: in the GRM for GCTA, or via a
mixed-model GWAS.

Dropping the family effect costs almost nothing in the slope estimates
themselves — BLUPs correlate r = 0.997 with and without it, at every visit
filter — so there is no tradeoff to weigh. The family effect is simply wrong
for this purpose.

![Family effect on heritability]({{artifact:art_92617adf-e522-4f53-9efd-62951ecbc602}})

### 4.3 `min_visits: 3` costs power; use ≥2

Per-subject slope reliability rises with follow-up length, but retaining fewer
subjects costs more than the reliability gain returns. Effective N — the
quantity GWAS power scales with — peaks at the permissive filter:

| run | subjects | mean visits | mean slope reliability | effective N | vs 5.1 |
|---|---|---|---|---|---|
| 5.1, ≥2 visits | 6,937 | 2.30 | 0.165 | 1,145 | 1.00 |
| **7.0, ≥2 visits** | **8,192** | **2.86** | **0.179** | **1,463** | **1.28** |
| 7.0, ≥3 visits | 5,195 | 3.35 | 0.206 | 1,068 | 0.93 |
| 7.0, 4 visits | 1,830 | 4.00 | 0.249 | 456 | 0.40 |

Effective N = subjects × mean reliability. The ≥3-visit filter yields *less*
power than 5.1 did (0.93×) despite the larger release; ≥2 visits on 7.0 buys a
28% gain. Note that 7.0's advantage is longer follow-up (2.86 versus 2.30 mean
visits) as much as more subjects.

Crossing the visit filter with the family effect over all eight runs confirms
that the two decisions are independent — the reliability ordering is identical
with and without the family effect, and effective N is maximised at ≥2 visits
in both cases.

![Reliability design grid]({{artifact:art_949ca80d-6a2c-4113-855c-6ae97c1d489a}})

![Release power tradeoff]({{artifact:art_0dd79d35-d9cc-47a4-94a9-b8e244c49274}})

Reliability is modest in absolute terms — mean 0.166 across regions, ranging
0.033 to 0.273. This is the fundamental constraint on the whole enterprise and
§12 returns to it.

---

## 5. The developmental phenotype: regional maps

Under the settled specification, all 68 regions thin over adolescence — there is
no region of net thickening. Mean rate is −0.019 mm/yr; the range across regions
is −0.031 (left frontal pole) to −0.002 mm/yr (left entorhinal).

By lobe, thinning is fastest in parietal cortex (−0.022 mm/yr) and frontal
(−0.021), slowest in temporal (−0.015) and occipital (−0.017).

Between-subject variability in rate (τ, the SD of subject slopes) tracks the mean
rate closely and peaks in the same frontal-pole region (0.011 mm/yr) where it is
lowest in lingual cortex (0.004).

![Developmental maps]({{artifact:art_73b07822-69cb-441d-a80b-e8629a646bd8}})

### Components of the slope map

Principal components of the subject × region slope matrix give three
interpretable axes. Variance explained is 24.1%, 5.5% and 4.5% (see §8 for an
important note on the variance convention).

- **PC1** contrasts parietal and superior-frontal cortex (loading most
  negatively: precuneus −1.46, inferior parietal −1.36) against cingulate and
  insula (caudal anterior cingulate +1.82, entorhinal +1.86). Lobe means run
  from parietal −0.97 to insula +1.40.
- **PC2** and **PC3** carry less variance but, as §7 shows, far stronger
  transcriptional coupling.

Regional loadings are in `docs/regional_slope_pc_loadings.csv`; the full map
table with every column used below is `docs/developmental_maps_noglobal.csv`.

---

## 6. Heritability

All heritability here is **Falconer's estimator** — twice the difference between
MZ and DZ correlations — computed on genotype-confirmed zygosity. Release 7.0
does not ship pi-hat, but the 5.1 genotype file covers 3,670 of the 7.0 subjects,
giving 260 MZ and 871 DZ-or-sibling pairs with phenotypes. §12 states the bias
this carries; GCTA on the cluster is the proper replacement.

**The developmental slope phenotype is heritable: h² = 0.584** (95% CI
0.354–0.821) for whole-cortex mean slope.

The positive control passes: baseline thickness gives h² = 0.796 (CI
0.662–0.929), comfortably inside the parameter space and consistent with the
published twin literature on cortical thickness. That control is what
distinguishes the settled specification from the family-effect run, where the
same control returned an impossible 1.508 (§4.2).

Candidate phenotypes:

| phenotype | h² (full sample) |
|---|---|
| baseline thickness *(control)* | 0.796 |
| slope PC1 | 0.593 |
| global mean slope | 0.584 |
| slope PC3 | 0.497 |
| slope PC2 | 0.441 |

![Heritability]({{artifact:art_d03d2023-cf1b-4e61-be6e-f683761ff1ff}})

### Per-region heritability, and why not to trust it region by region

Regional h² has median 0.321 (IQR 0.191–0.488), peaks at 0.864 in left
precentral gyrus, and is negative — an out-of-range estimate, i.e. r(MZ) < r(DZ)
— in 2 of 68 regions. By lobe it is highest parietal (0.450) and frontal (0.365),
lowest cingulate (0.170) and insula (0.198).

**The regional map is not stable enough to rank individual regions.** Split-half
Spearman correlation of the per-region h² map is 0.32 across 20 family-level
splits. Any claim of the form "region X is the most heritable" would not survive
resampling, and §10 shows what happens when a phenotype is built on such a
claim.

Two negative results constrain everything downstream:

1. Per-region h² is **uncorrelated with regional thinning rate** (ρ = 0.018).
   Where the cortex thins fastest is not where thinning is most heritable.
2. Per-region h² shows **no association with any AHBA component** that survives
   the spin test (§7).

Heritability and transcriptional relevance are therefore separate axes. No single
regional criterion optimises both, which is precisely why the phenotype decision
in §10 required an explicit tradeoff rather than a maximisation.

---

## 7. Transcriptional comparison: AHBA components

Maps are compared to the three AHBA developmental components (C1–C3) from prior
work in this group, using spatial permutation ("spin") tests with 5,000
rotations. The lowest attainable p-value is therefore 1/5001 ≈ 0.0002, and
p-values at that floor are reported as bounds, not as zero.

The comparison runs in the DK parcellation on the 34 bilateral regions where the
AHBA components are defined. This coarsens the transcriptional side
substantially — effect sizes should be read as a lower bound, and the spin test
has few points to work with.

### The replication holds

**Absolute thinning rate versus AHBA C3: ρ = −0.546, p_spin = 0.0016.** This
closely matches the thesis 5.1 result (−0.550, p = 0.004). The phenotype refit
without the global covariate recovers the original association; the sign, the
magnitude and the significance all reproduce across a release change and a full
pipeline rewrite.

### The strongest associations in the project are on the components

| map | component | ρ | p_spin |
|---|---|---|---|
| **slope PC3** | **C2** | **+0.854** | **< 0.0002** |
| **slope PC2** | **C1** | **−0.812** | **< 0.0002** |
| between-subject SD of rate | C1 | −0.505 | 0.027 |
| slope PC1 | C1 | −0.511 | 0.035 |
| absolute thinning rate | C3 | −0.546 | 0.0016 |
| per-region h² | *(all)* | — | n.s. |

Slope PC2 and PC3 track AHBA C1 and C2 at effect sizes substantially exceeding
the mean-rate/C3 replication. These are the strongest spatial associations
anywhere in this analysis, and they reordered the phenotype priority in §10:
these components carry the transcriptional signal the project is built around,
even though they are not the most heritable candidates.

All 18 map × component pairs, with both spin and naive p-values, are in
`docs/ahba_vs_maps_noglobal.csv`.

![AHBA comparison]({{artifact:art_435c9e2f-1dfb-42b2-a774-ea403c03d437}})

---

## 8. Structural covariance of developmental slopes

The 68 × 68 correlation matrix of subject slopes across regions — structural
covariance of *development*, rather than of static thickness.

The matrix is almost entirely positive: only 13 of 2,278 region pairs are
negative, and the most negative is −0.076. Mean off-diagonal r is 0.205 (median
0.186, max 0.757). Subjects who thin faster in one region thin faster nearly
everywhere, which is the same fact that makes the global covariate so damaging
(§4.1).

Coupling decomposes by anatomy:

| pair type | mean r |
|---|---|
| homotopic (same region, opposite hemisphere) | 0.455 |
| within-lobe | 0.344 |
| within-hemisphere | 0.226 |
| between-lobe | 0.182 |

Region identity dominates hemisphere: the homotopic mean (0.455) far exceeds the
within-hemisphere mean (0.226). Within-lobe coupling is strongest in parietal
(0.472) and occipital (0.458), weakest in cingulate (0.171).

![Structural covariance]({{artifact:art_0d2955df-b5cb-4c75-a383-dd322983bfe8}})

### The SC-PCs and the slope PCs are the same decomposition

Eigendecomposing the structural covariance matrix and running PCA on the
subject × region slope matrix return **identical** components: loading
correlation r = 1.0 with maximum absolute difference 0.0 for PC1–PC3, and
subject scores that match exactly. This is expected — PCA on standardised
columns *is* the eigendecomposition of their correlation matrix — but the two
are implemented by independent code paths, so the identity is a useful
cross-check. It is a coding verification, not independent evidence: the two
routes are not separate findings and must not be reported as such.

`test_sc_pcs_equal_slope_pcs` pins this to exact equality rather than a
correlation threshold, so that a refactor dropping the per-region
standardisation would fail rather than silently produce two plausible but
different maps.

### Variance convention — read this before quoting a percentage

`covariance.py` exposes two conventions and they differ by more than threefold:

| component | λ / Σλ *(standard)* | λ² / Σλ² *(as reported in the thesis)* |
|---|---|---|
| PC1 | **24.1%** | 79.8% |
| PC2 | 5.5% | 4.1% |
| PC3 | 4.5% | 2.8% |

The thesis figure of roughly 80% of variance in PC1 is the **squared**
convention. The standard convention gives 24.1%. Any statement about variance
explained must name which is meant; this report uses the standard convention
throughout except where explicitly stated.

---

## 9. Site and scanner effects

The analysis sample spans 18 acquisition sites and 31 distinct scanner serials,
and subjects sometimes change scanner mid-study, so site and scanner are
candidate confounds for a longitudinal phenotype.

They are present but small. Between-site ICC for the global slope phenotype is
0.021 across 18 sites; baseline scanner serial gives 0.023 across 29 serials;
manufacturer gives 0.002. All are significant given n = 8,192 but account for
about 2% of variance.

Scanner changes are common — 26.5% of subjects switch serial, 19.0% switch
model, 0.7% switch manufacturer — yet their effect on the phenotype is modest:
variance ratio 1.35 for switchers versus non-switchers, Cohen's d = 0.14.

Critically, heritability is insensitive to all of this:

| phenotype variant | h² | 95% CI |
|---|---|---|
| global slope, raw *(settled)* | 0.584 | 0.354 – 0.821 |
| site-mean removed | 0.599 | 0.367 – 0.835 |
| site + scanner-mean removed | 0.604 | 0.371 – 0.843 |
| non-switchers only (n = 6,025) | 0.556 | 0.279 – 0.804 |

Removing site and scanner means moves h² by +0.02, well inside the confidence
interval, and restricting to non-switchers costs sample size without changing
the estimate. **No site or scanner correction is applied to the primary
phenotype**, and the GWAS should not need scanner covariates beyond the standard
set — though including them is harmless.

![Site and scanner supplement]({{artifact:art_371877ab-0f2b-4ddd-9a39-fc2706f9a1ff}})

![Site ICC map]({{artifact:art_d490090f-758e-4b26-8bcf-096d1d3da742}})

---

## 10. Phenotype priority for GWAS

Which phenotype should go to GWAS first? This was re-derived from scratch on the
settled specification; the earlier ordering was computed on the global-adjusted
design and does not transfer.

**Method.** Any phenotype whose *definition* consumes the twin data must be
scored out of sample, or its heritability is inflated by selection. Across 40
family-level splits, PCA loadings are fitted on the training half and
*projected* onto the held-out half, so no phenotype definition ever sees the
twins it is evaluated on. A "top-10 most heritable regions" phenotype is
included precisely because it violates this in the obvious way.

| phenotype | held-out h² (SD) | vs global mean | best AHBA | definition |
|---|---|---|---|---|
| top-10 h² regions | 0.649 (0.15) | +0.055, p = 0.017 | none | **unstable** |
| slope PC1 | 0.611 (0.15) | +0.017, p = 0.51 | C1, −0.511 | deterministic |
| global mean slope | 0.593 (0.13) | — | C3, −0.546 | deterministic |
| slope PC3 | 0.511 (0.15) | −0.082, p = 0.016 | C2, **+0.854** | deterministic |
| slope PC2 | 0.443 (0.13) | −0.151, p < 0.001 | C1, **−0.812** | deterministic |

![Phenotype priority]({{artifact:art_b67e7df5-0c5d-4f2d-b105-f51d3830c6bf}})

### Why the most heritable phenotype is rejected

Selecting the top-10 most heritable regions on a training half gives the
*highest* held-out h² of any candidate, 0.649. On the earlier global-adjusted
design this phenotype collapsed out of sample (0.363), and that collapse was the
stated reason for rejecting it. On the settled design it does not collapse, so
that objection no longer holds.

The correct objection is different, and still disqualifying. The optimism gap
remains 0.10 (in-sample 0.808 versus held-out 0.706), and — decisively — only
**6.1 of 10** selected regions are re-selected across splits. The phenotype is
heritable, but its *definition* is not reproducible. For a GWAS that must be
pre-registered and replicated in an independent cohort, a phenotype whose
constituent regions change with the sample is not usable, however well it
scores.

### Why transcriptional signal decides the ranking

Held-out h² spans only 0.443 to 0.649 across all five candidates, and the
largest paired difference is 0.055. Slope PC1 versus the global mean gives
p = 0.51 — statistically indistinguishable across 40 paired splits. AHBA effect
sizes, by contrast, span 0.51 to 0.85.

Heritability is nearly flat across the choice set; transcriptional coupling is
not. Since the scientific question is about *genes*, the criterion that
discriminates should drive the ranking.

### Recommended order

1. **Global mean slope.** Highest h² among deterministic phenotypes, replicates
   the thesis C3 association, zero selection risk. The safe primary.
2. **Slope PC3** (C2, ρ = +0.854) and 3. **slope PC2** (C1, ρ = −0.812). By far
   the strongest transcriptional coupling in the analysis. They cost 0.08–0.15
   in h², which GWAS power can absorb more easily than it can absorb a weak
   brain–transcriptome link.
4. **Slope PC1.** h² indistinguishable from the global mean but weaker AHBA
   association, so it adds little the mean does not already provide.

**Top-h² regions: not recommended at any position.**

Tables: `docs/gwas_phenotype_priority.csv`, `docs/heldout_h2_by_split.csv`,
`docs/heldout_h2_paired_tests.csv`, `docs/topH2_selection_optimism.csv`.

---

## 11. Corrections to the previous report

Two conclusions in the previous version of this document were wrong. Both
followed from the same root cause — the global-thickness covariate — and both are
recorded here rather than quietly edited, because the pattern is instructive.

**Correction 1: the maps do resemble the AHBA components.** The previous report
concluded that developmental maps showed no transcriptional association,
contradicting the thesis. That null was computed on the global-adjusted
phenotype. Refitting without the global covariate recovers the thesis result
almost exactly (ρ = −0.546 versus −0.550) and additionally reveals two much
stronger associations on the slope components (§7). The reported null was an
artefact of the specification, not a failure to replicate.

**Correction 2: the headline heritability was h² ≈ 0.39, now 0.584.** The
earlier figure came from the global-adjusted run *with* the family effect
enabled. Both choices depress it. The settled specification gives 0.584, and its
positive control passes where the earlier one returned an impossible 1.508
(§4.2).

**Correction 3: a reported mean-visit count was computed by the wrong
weighting.** `docs/fit_summary.csv` and `docs/reliability_grid.csv` reported
mean visits averaged over observation rows rather than over subjects, giving
3.06 where the subject-level mean is 2.86. Both tables have been corrected.
Effective N never used that column and is unaffected.

**Correction 4: the cross-release map correlation was unsourced.** The previous
report gave ρ = 0.977 without naming the coefficient, and no run on disk could
reproduce it — no 5.1 run had been fitted under the no-global specification. That
run has now been fitted, and the matched comparison gives Spearman ρ = 0.937 /
Pearson r = 0.975 (§3, `docs/release_map_reproducibility.csv`). The old figure
appears to have been a Pearson correlation reported as if it were the rank
correlation used elsewhere in the document.

The general lesson is that all four defects were caught by *controls and
cross-checks*, not by inspecting the headline numbers, which looked plausible
throughout. The positive control caught the family effect; the release
comparison caught the visit weighting; the thesis replication caught the global
covariate.

---

## 12. Limitations

**The twin estimator is biased upward.** Falconer's h² requires MZ and DZ
classes. Release 7.0 does not ship pi-hat, so zygosity comes from the 5.1
genotype file, which covers 3,670 of the 7.0 subjects. More importantly, the
"DZ" class is diluted with full siblings, who share less than the 0.5 the
estimator assumes; this inflates the MZ–DZ difference and therefore h². Every
heritability figure in this report should be read as an upper bound. **GCTA on
the GRM is the proper replacement and is the first cluster job** (§13).

**Slope reliability is low.** Mean per-subject slope reliability is 0.166. This
is inherent to estimating an individual trajectory from two to four points, and
it attenuates every correlation reported here — the AHBA associations included.
It also caps GWAS power: effective N is roughly 1,463, not 8,192.

**The transcriptional comparison is coarse.** AHBA components are compared on 34
bilateral DK regions. This avoids resampling the imaging data but coarsens the
transcriptional side considerably, so effect sizes are a lower bound and the spin
test has few points. Reported ρ values are conservative.

**Per-region heritability is not interpretable region by region.** Split-half
ρ = 0.32 (§6). Use it as a map, never as a ranking.

**Structural covariance is dominated by a single global factor.** With mean
off-diagonal r = 0.205 and only 13 of 2,278 pairs negative, the matrix is close
to a one-factor structure. PC2 and PC3 explain 5.5% and 4.5% of variance — small
in absolute terms, even though they carry the strongest transcriptional
associations.

**No fMRI or multi-metric integration yet.** The framework accepts other
metrics (see `configs/t1t2_baseline.yaml`), but everything here is cortical
thickness.

---

## 13. What runs next, on the cluster

1. **GCTA GRM heritability** on the four recommended phenotypes, to replace the
   Falconer estimates with unbiased SNP heritability. `gcta_export.py` writes
   the phenotype and covariate files; it refuses to export from a run with the
   family effect enabled.
2. **GWAS** on the global mean slope first, then slope PC3 and PC2, using a
   mixed model to handle relatedness rather than pruning relatives.
3. **MAGMA gene-level enrichment**, then test for enrichment against SCZ and MDD
   GWAS, which is the primary hypothesis. `magma_export.py` prepares inputs.
4. **Compare the resulting gene sets** to AHBA C1–C3 and to the snRNA-seq
   leading component.

The h²-versus-transcription tradeoff in §10 means step 3 should be run on both
the global mean *and* the components — they may implicate different genes, and
which one connects to the disorder GWAS is the empirical question the project
exists to answer.

---

## 14. Reproducing this report

Settled configuration: `configs/ct_70_noglobal_mv2_genetic.yaml`
(`global_covariate: none`, `min_visits: 2`, `family_effect: false`).

| section | table | figure |
|---|---|---|
| §3 | `release_map_reproducibility.csv` | — |
| §4.3 | `reliability_grid.csv`, `handoff_release_comparison.csv` | `reliability_design_grid.png`, `release_power_tradeoff.png` |
| §4.1–4.2 | `h2_global_slope_by_design.csv`, `family_effect_heritability_noglobal.csv`, `family_effect_blup_agreement.csv` | `family_effect_heritability.png` |
| §5 | `developmental_maps_noglobal.csv`, `regional_slope_pc_loadings.csv` | `developmental_maps_noglobal.png` |
| §6 | `h2_candidate_phenotypes.csv`, `h2_positive_control.csv`, `h2_map_split_half.csv`, `regional_h2.csv` | `heritability_noglobal.png` |
| §7 | `ahba_vs_maps_noglobal.csv` | `ahba_vs_maps_noglobal.png` |
| §8 | `sc_matrix_summary.csv`, `sc_matrix_lobe_ordered.csv`, `sc_vs_slope_pc_identity.csv` | `structural_covariance.png` |
| §9 | `site_scanner_icc.csv`, `site_scanner_icc_by_region.csv`, `scanner_switching_summary.csv`, `h2_site_scanner_sensitivity.csv` | `site_scanner_supplement.png`, `site_icc_map.png` |
| §10 | `gwas_phenotype_priority.csv`, `heldout_h2_by_split.csv`, `heldout_h2_paired_tests.csv`, `topH2_selection_optimism.csv` | `gwas_phenotype_priority.png` |

Tables are in `docs/`, figures in `docs/figures/`. Test suite: 118 tests,
`pytest` from the repository root.
