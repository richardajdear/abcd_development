# D3 — twin and family designs

Implements direction **D3** of [`docs/DIRECTIONS.md`](../../docs/DIRECTIONS.md): the three
designs that SNP-based genetics cannot supply at this sample size.

1. **Twin ACE models** (univariate and bivariate) give the heritability of the thinning
   rate independent of SNP tagging, and its genetic correlation (rA) with co-developing
   traits: puberty, the T1w/T2w slope, CBCL symptom change and cognitive gain.
2. **Within-family polygenic scores** separate the direct genetic effect of a score from
   stratification, assortative mating and genetic nurture.
3. **Co-twin control**: does the faster-thinning twin of an MZ pair have more symptoms at
   15–17 (DIRECTIONS.md prediction T6)?

Everything here is from ABCD 7.0, run locally on 2026-09-30. Figure:
[`figures/fig_d3_twin_family.png`](figures/fig_d3_twin_family.png).

## Results

Numbers are from `results/`; the figure reads the same tables.

**Sample** (`pair_counts.tsv`). The release has 419 MZ, 659 DZ and 846 non-twin sibling
pairs with genetically inferred zygosity. Of these, **271 MZ, 432 DZ and 579 sibling
pairs** have both members imaged (≥ 2 scans), after keeping one pair per birth event and
one sibling pair per family. Only 61 of the 432 DZ pairs (14 %) are opposite-sex, a
property of the release rather than of the pair selection (the rate is 13.5 % across all
659 DZ pairs).

**1. The thinning rate is heritable, and the whole of its familial resemblance is genetic**
(`ace_univariate.tsv`, `twin_correlations.tsv`)

| trait | rMZ | rDZ | A (ACE) [95 % CI] | C |
|:--|--:|--:|:--|--:|
| thinning rate ΔCT | 0.43 | 0.19 | **0.46 [0.37, 0.54]** | 0.00 |
| ΔCT, scan-quality adjusted | 0.43 | 0.20 | 0.46 [0.36, 0.54] | 0.00 |
| cortical thickness (positive control) | 0.89 | 0.50 | 0.83 [0.67, 0.91] | 0.08 |
| T1w/T2w slope (DK) | 0.53 | 0.11 | 0.51 [0.40, 0.62] | 0.00 |
| puberty timing | 0.93 | 0.54 | 0.76 [0.61, 0.91] | 0.16 |
| CBCL depressive change | 0.57 | 0.24 | 0.57 [0.46, 0.67] | 0.00 |

The twin estimate for ΔCT is 2.5 times its GREML SNP h² (0.18), and the gap for thickness
is larger still (0.83 vs 0.23). The sibling correlation (0.11) is lower than the DZ
correlation (0.19). Siblings share half their genes like DZ twins, but they are scanned at
different ages and on different days, so their slopes span different developmental
windows. Cognitive gain scores (NIH Toolbox 00A→06A) are too noisy to decompose: fluid
and total gain have rMZ < rDZ.

**2. Genes for earlier puberty are genes for faster thinning; genes for thinning are not
genes for symptom change** (`ace_bivariate.tsv`, AE model; C for ΔCT is at its bound of 0)

| ΔCT with … | rA [95 % CI] | rE [95 % CI] | rP |
|:--|:--|:--|--:|
| puberty timing | **−0.24 [−0.33, −0.15]** | 0.00 [−0.13, 0.12] | −0.15 |
| puberty tempo | −0.06 [−0.15, 0.04] | −0.02 | −0.04 |
| cortical thickness | 0.18 [0.08, 0.28] | −0.01 | 0.12 |
| T1w/T2w slope | 0.08 [−0.12, 0.25] | −0.03 | 0.02 |
| CBCL depressive change | 0.08 [−0.09, 0.24] | −0.11 [−0.22, 0.02] | −0.01 |
| CBCL internalising change | 0.11 [−0.05, 0.28] | **−0.16 [−0.27, −0.03]** | −0.02 |
| CBCL externalising / p-factor change | 0.04 / 0.07 (all CIs span 0) | −0.05 / −0.10 | 0.00 |

Higher ΔCT means slower thinning, so a negative rA with puberty timing means that
children genetically inclined to earlier puberty are genetically inclined to faster
thinning. The genetic term accounts for essentially all of the phenotypic correlation
(rP = −0.15, of which A contributes −0.15). The trait is parent-reported, so rater bias
cannot produce a genetic correlation with an MRI measure. This is the first evidence in
the project for DIRECTIONS.md prediction **T2** (the same liability shifts another
maturational clock), and it holds under scan-quality adjustment (−0.25 [−0.35, −0.16]).
It agrees with the earlier phenotypic result (genetic_analysis/c3axis, README_HPC §8.6:
earlier timing → faster global thinning). Puberty *tempo* (the PDS slope) shows no
genetic overlap with the thinning rate.

The symptom pattern is the reverse. No rA with symptom change excludes zero, and every
point estimate has the opposite sign to the population association. The negative
(faster-thinning ↔ more symptoms) covariance sits in the **unique environment** term,
which is nominally significant for internalising change. Unique environment includes
measurement error, which cannot correlate with a parent's report, and any non-shared
cause.

The thinning rate and the T1w/T2w slope are both heritable (0.46, 0.51) but share little
genetics (rA 0.08, CI −0.12 to 0.25). This is a first, weak strike against the reading
that the thinning rate *is* the myelination process (DIRECTIONS.md §2.2). The T1w/T2w
slope is DK-only and uncorrected for bias field, so the test is not decisive.

**3. Within-family PRS: not assessable at this size** (`within_family_prs.tsv`)

1,087 families (2,214 children) have ≥ 2 genotyped, imaged children once MZ co-twins are
removed. The population betas reproduce Figure 1f (SCZ 2025 SBayesRC −0.033, p = .001;
MDD SBayesRC −0.027, p = .012). The within-family SE is 0.034–0.040, so the effect
detectable at 80 % power is **2.9–4.2 times the population effect**. For MDD, the
within-family estimates (−0.026 PRS-CS, −0.023 SBayesRC) equal the population ones.
For SCZ they fall on the other side of zero (+0.025, +0.030; tested against the
population effect, p = .11 and .07). Neither is a finding; a 95 % CI of the
within/population ratio spans roughly −3 to +1.4 for SCZ and −2 to +4 for MDD. EA scores
are cluster-only, so the EA question named in DIRECTIONS.md is not assessed here.

**4. Co-twin control: the slope–symptom association does not shrink within MZ pairs**
(`cotwin_control.tsv`)

| β of ΔCT on CBCL at 15–17 given baseline (per SD) | population | within MZ (250 pairs) | within DZ (418) | within SIB (546) |
|:--|--:|--:|--:|--:|
| depressive | −0.032 (p .003) | −0.077 (p .17) | −0.055 | +0.003 |
| internalising | −0.017 | −0.092 (p .047) | −0.063 | +0.027 |
| externalising | −0.022 (p .026) | −0.103 (p .050) | +0.020 | −0.027 |
| p-factor | −0.020 (p .047) | −0.069 (p .13) | −0.052 | +0.005 |

The population column reproduces Figure 1g. Within MZ pairs, where genes and shared
family environment are fixed, every coefficient has the population sign and is 2.4–5.5
times larger, and 2 of 4 CIs exclude zero. So the association is **not explained by
familial confounding**; if genes or family explained it, the within-MZ estimates would
fall towards zero. The CIs are wide (SE ≈ 0.05), so the design cannot show the
within-MZ effect equals the population one, nor that it is larger. Within sibling pairs the estimates are null. Siblings are scanned at different ages, so
their slopes and symptom windows cover different parts of adolescence, which the pair
design does not align. Adjusting for scan quality (mean and slope of log topological
defects), baseline thickness, or both shifts the population estimates by at most 0.004
and the within-MZ estimates by at most 0.009 (`adjust` column). DZ and sibling estimates
move by up to 0.011 and 0.021, within their SEs.

Read together with result 2, the thinning-rate–symptom link runs through factors that
are not shared by co-twins: non-shared environment, or a direct effect. It does not
run through the genetics that make the thinning rate heritable. That fits the project's
existing null for the PRS → thinning → symptom chain at the individual level.

## Measurement caveats specific to twins (`cotwin_measurement.tsv`)

- **Twins are scanned in the same session** (median age gap at a matched visit: 0 days;
  siblings 476 days). Scan-level residuals from each child's own line correlate 0.15 in
  MZ and 0.18 in DZ pairs, and 0.00 in siblings. This is shared session error, equal in
  MZ and DZ pairs, so in an ACE model it loads on C, not A. C for ΔCT is estimated at 0
  even so.
- **Image quality is twin-similar**: log topological defects correlate 0.31 (MZ) vs 0.16
  (DZ). Because MZ > DZ, quality could inflate A. Adjusting for it (ΔCT_rq) leaves A at
  0.46, rA with puberty at −0.25, and the co-twin estimates unchanged. The defect count is
  the only quality measure available locally; head motion (`mr_y_qc__mot`) is on the
  fetch list and should replace or join it.
- **The LMM slope reliability understates co-twin agreement.** The mean per-child BLUP
  reliability is 0.30, yet rMZ is 0.43 after covariates. Under independent errors rMZ
  cannot exceed the reliability, so twin errors are not independent: the shared session
  error above is one reason. This is also why a classical attenuation argument for the
  MZ-difference design does not apply.
- CBCL is parent-reported, with one parent usually rating both twins. Rater effects
  appear as C, or as contrast effects (rMZ > 2 rDZ, as for depressive change), in the
  univariate symptom models. They cannot create rA or rE with the MRI slope.
- The DZ sample is 86 % same-sex, and sex is regressed out of the means only
  (no sex-limitation model).

## Pipeline

From the repo root; Python in `abcd-spatial` (`PYTHONPATH=src` for 03), R in `r`.

| step | script | writes (results/ = committed aggregates; work/ = individual-level, gitignored) |
|:--|:--|:--|
| 0 | `00_ace_selftest.py` | `ace_selftest.tsv` — parameter recovery of `ace.py` on simulated twins (passes) |
| 1 | `01_scan_means.py` | `work/scan_means_{ct,t1t2}.csv` |
| 2 | `Rscript 02_fit_slopes.R` | `work/blups_{ct,t1t2}.csv`; `lmm_variance_components.tsv`, `slope_reliability.tsv` |
| 3 | `03_build_traits_pairs.py` | `work/traits.parquet`, `work/pairs.csv`; `pair_counts.tsv`, `trait_availability.tsv` |
| 4 | `04_twin_ace.py [n_boot=500]` | `twin_correlations.tsv`, `ace_univariate.tsv`, `ace_bivariate.tsv` (~22 min) |
| 5 | `05_within_family_prs.py` | `within_family_prs.tsv`, `within_family_sample.tsv` |
| 6 | `06_cotwin_control.py` | `cotwin_control.tsv`, `cotwin_measurement.tsv` |
| 7 | `Rscript 07_openmx_check.R` | optional, **not run**: needs an R with OpenMx (see below) |
| fig | `Rscript fig_d3_twin_family.R` | `figures/fig_d3_twin_family.png` |

Design choices worth knowing:

- **ΔCT is the Figure-1 trait**: the single LMM on the per-scan HCP-MMP cortical mean,
  with no family effect, since a family effect would shrink co-twin BLUPs towards each
  other. Step 3 asserts r > 0.999 with `global_slope_1lmm` from the genetics export.
- **Covariates leave the mean model before the twin fit** (sex, site, age at first scan,
  span, visit count), because twins share them and they would otherwise count as C.
- **Twin models are written out in `ace.py`** (normal-theory FIML, the likelihood OpenMx
  maximises), because conda has no OpenMx build for osx-arm64. A CRAN source build
  failed in the conda R toolchain (no clang wrapper for `mvtnorm`, no cmake for
  `RcppParallel`). `00_ace_selftest.py` shows unbiased recovery (large-n means within
  0.01 for variance components, 0.01 for rA). `07_openmx_check.R` compares the fits
  against OpenMx where it is installed, e.g. on CSD3.
- **AE is the primary bivariate model**: C for ΔCT is at its bound. Rows where either
  trait's A < 0.02 are flagged `rA_identified = False`; their rA sits at ±1 in the ACE
  fit (fluid/total cognitive gain).
- **Bootstrap CIs** resample pairs within zygosity: 500 for univariate, 250 for bivariate.

## Next steps this suggests

1. **Puberty × thinning in the one-stage LMM (D2).** The genetic overlap is with timing,
   not tempo. That is the pattern a phase-advance (T1) reading predicts, and D2 is the
   direct test.
2. **Add head motion** (`mr_y_qc__mot`) to the quality sensitivity once fetched, and
   re-run steps 3–6.
3. **Within-family EA**: rsync the EA profiles from CSD3 and add them to
   `common.SCORES`. Expect the same power limit, since MDE ≈ 3× the population effect.
4. The co-twin result makes non-shared environment the place to look for what links
   thinning to symptoms (D5: adversity and life events, now vendored as
   `mh_{y,p}_ple`, `le_l_adi`).
