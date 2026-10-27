# D1 — multimodal slopes (myelination vs neuropil; composite slope)

See `docs/DIRECTIONS.md` §2.2, §2.4 and D1. Each modality gets the Figure 1 treatment: a whole-cortex mean per scan and
one single LMM `mean ~ age_c + sex + (1 + age_c | subject) + (1 | site)`. The per-child slope then goes into the
Figure 1 PRS model, the D3 twin pipeline and, eventually, a cross-modal composite.

| step | file | status |
|---|---|---|
| 1 | `01_t1t2_prs.R`: PRS → T1w/T2w slope, plus the quality of the measure (DK; HCP-MMP T1w/T2w does not exist yet) | done 2026-10-27 |
| — | `CSD3_INSTRUCTIONS.md`: HCP-MMP and DK GWC, T1/T2 intensity at GM mid-depth and in superficial WM, T1w/T2w ratio, and DTI/RSI if processed maps exist | GWC, T1 gm50, T1 wm1 done on CSD3 2026-10-27 (`CSD3_REPORT.md`); T2, T1w/T2w, DTI held for a decision |
| 2 | `02_dmri_scan_means.py`: per-scan whole-cortex means for 45 DK DTI/RSI release tables | done 2026-10-27 |
| 3 | `03_fit_compare.R`: single-LMM slopes, quality, coupling with thinning / T1w/T2w, PRS | done 2026-10-27 |
| 4 | `04_followups.R`: PRS adjusted for the thinning slope; PCA of the 47 slopes | done 2026-10-27 |
| 5–6 | `05`–`08`: PLS composite, diffusion PCs × PRS, single-factor test (twin rA, tetrad) | done 2026-10-27 |
| — | area, volume, GWC and T1 gm50/wm1 (CSD3 tables, rsync pending) → same pipeline; twin h² / rA for every new slope | to do |

## Step 1 result (`tables/t1t2_*.tsv`)
**Gate.** The base model reproduces `genetic_analysis/fig1_inputs/hcp70_prs_key_arms.tsv` for ΔCT in all 8 cells
(max |Δβ| 1e-11; identical n).

**The ABCD T1w/T2w measure is dominated by site** (`t1t2_measure_qc.tsv`):

| | thickness (HCP) | T1w/T2w (DK) |
|---|--:|--:|
| scan-level site ICC | 0.044 | 0.193 |
| R² of site on the child slope | 0.030 | 0.261 |
| children with slope \|z\| > 8 | 0 | 10 |
| mean per-child slope reliability, 2 / 3 / 4 visits | 0.17 / 0.32 / 0.38 | 0.26 / 0.43 / 0.50 |
| r(faster thinning, faster T1w/T2w rise), all / within site | | −0.002 / 0.013 |

The T1w/T2w "reliability" is higher than thickness's only because the LMM has site on the intercept and not on the
slope. Scanner drift within a site is therefore consistent within a child, and it counts as reliable slope variance.

**PRS → T1w/T2w slope: nothing consistent** (β per SD, z-scored outcome; positive = faster rise / faster thinning):

| score | arm | thinning (base) | T1w/T2w slope, base | T1w/T2w slope, + site |
|---|---|--:|--:|--:|
| SCZ PRS-CS | EUR | +0.036 (p 0.013) | +0.034 (0.031) | +0.020 (0.14) |
| SCZ PRS-CS | pooled | +0.029 (0.004) | +0.001 (0.93) | −0.005 (0.60) |
| MDD PRS-CS | EUR | +0.029 (0.043) | −0.009 (0.58) | +0.002 (0.85) |
| MDD PRS-CS | pooled | +0.024 (0.019) | −0.025 (0.027) | −0.004 (0.65) |

SBayesRC gives the same pattern. The two nominal T1w/T2w cells take opposite signs, do not replicate across arms and
disappear once site is added. The thinning effects do not move with site. The 80 % MDE is 0.035 SD, so an effect as
large as the PRS → thinning effect (0.024–0.036) cannot be excluded either. The test is underpowered as well as null.

**A second, unplanned observation:** a higher SCZ score goes with a lower T1w/T2w **level** (intercept). Pooled
PRS-CS gives −0.028 (p 0.010), and −0.029 (p 0.008) with site. SBayesRC pooled gives −0.030 (p 0.006). The EUR arm has
the same sign (−0.020 to −0.027, p 0.09–0.20). MDD is null. This was not pre-specified, so treat it as a lead for the
HCP-MMP/GWC data, not as a finding.

**Reading.** The DK release T1w/T2w cannot decide myelination vs neuropil: about a quarter of its slope variance is
site, and it does not covary with thinning even phenotypically. The D3 rA of 0.08 inherits the same problem.
Re-test with the CSD3 measures: GWC from the FreeSurfer surfaces, and T1w/T2w built from a T2 registered to each
subject.

## Steps 2–4 result: DTI and RSI (DK, release tables; `tables/slope_summary.tsv`, `dmri_prs*.tsv`, `slope_pca.tsv`)
**Sample.** Every one of the 45 tables gives 26,257 scans from 8,596 children. The filters are the release dMRI include
flag, no Philips, and at least 2 scans. No LMM was singular.

**Age trends go the expected way** (SD of the between-child intercept per year):
- White-matter FA rises (+0.24 to +0.27).
- MD and TD fall in grey and white matter.
- RSI restricted signal rises in white matter (RND +0.33, RNI +0.55) and grey matter (RNI +0.50).
- Grey-matter free water (FNI) rises slightly (+0.08).
- The contrast terms (`gwc`) of FA and RND rise (+0.45 to +0.60).

The %/yr column is not meaningful for `gwc` measures whose mean is near 0.

**Every diffusion slope carries much more site variance than thickness.** The R² of site on child slopes is 0.07–0.40,
against 0.03 for thickness. The `gwc` contrast measures are the least site-sensitive as a group: median 0.09 vs 0.18
for the `gm`/`wm` measures, range 0.07–0.28. Thirteen of the fifteen are at most 0.13; the exceptions are the LD
contrasts, `is__ld__gwc` 0.14 and `fs__ld__gwc` 0.28. `is__md__gwc` has a scan-level site ICC of 0.05 and slope R²
0.07, the only one close to thickness. Slope reliabilities (3 visits:
0.12–0.41) are of the same order as thickness (0.32) but, as for T1w/T2w, partly reflect site drift.

**Coupling with faster thinning (within-site r of child slopes):**

| direction | measures (within-site r) |
|---|---|
| strongest positive | grey-matter MD +0.35 and TD +0.34 (inner shell); grey-matter free water (FNI) +0.33; grey-matter LD +0.28 |
| negative | MD, LD and TD contrast (`gwc`) −0.27 to −0.30; grey-matter hindered water (HNI/HNT) −0.23 |
| white-matter myelination markers | FA −0.13, RND −0.13, RNT −0.14 (opposite to the myelination prediction) |
| T1w/T2w slope | \|r\| ≤ 0.11 with every diffusion slope |

The grey-matter MD/free-water coupling has the signature of partial volume. Thicker cortex goes with lower grey-matter
MD (intercept r −0.40) and lower FNI (−0.48), so a thinner ribbon puts more CSF into the sampled voxels.

**PCA.** The 47 site-residualised slopes span about 5 axes (variance 0.33, 0.17, 0.12, 0.10, 0.08):
- PC1: white-matter maturation (FA, RND, RNT vs TD).
- PC2: grey-matter MD / free water.
- PC3: full-shell diffusivities.
- PC4: hindered-water fractions.

Thickness loads at most 0.12 on any of them, so the thinning slope is largely its own axis.

**PRS (Figure-1 model, + site; `dmri_prs.tsv`).**
- Nominal p < 0.05 in 58 of 376 cells against 19 expected. They are concentrated in the EUR arm (36 of 180 vs 14 of 180
  pooled) and in MDD.
- **Nothing survives BH-FDR over the 360 diffusion tests (min q 0.23).** The tests are highly correlated, so the
  excess is less than it looks.
- Thickness is still the only measure significant in all 8 score × arm cells.
- Leading cells (PRS-CS, EUR arm): MDD → grey-matter FNI +0.045 (p 0.002); MDD → FA contrast +0.045 (p 0.002);
  MDD → HND contrast +0.044 (p 0.003); MDD → grey-matter TD/MD +0.043/+0.042 (p 0.004).

**Adjusting for the thinning slope (`dmri_prs_adj_thin.tsv`) splits these into two groups:**
- **Grey-matter MD/TD/FNI** shrink by a quarter (MDD EUR, e.g. FNI +0.045 → +0.036) to two-thirds (SCZ and pooled
  arms, e.g. SCZ EUR +0.023 → +0.010). They are partly the thinning effect measured again, plausibly through partial
  volume.
- **Contrast measures (FA, RND and HND `gwc`)** do not change: MDD EUR FA contrast +0.045 → +0.047, SCZ EUR RND contrast
  +0.033 → +0.036 (p 0.016). Higher scores go with a faster rise of white-vs-grey diffusion contrast at the boundary,
  independently of thinning. They are nominal, EUR-arm-led, and not FDR-significant, so treat them as candidates for
  the HCP-MMP / CSD3 data and the twin pipeline, not findings.

**Reading so far.**
1. The release diffusion measures do not support the myelination account of the PRS-linked thinning. The white-matter
   myelination markers do not rise faster in faster-thinning children (r −0.13), and they carry no PRS signal.
2. What tracks thinning is grey-matter water mobility, which is most simply partial volume of a thinner ribbon. It is
   also consistent with neuropil loss; the two cannot be separated with these data.
3. A cross-modal composite will not raise the thinning slope's reliability: thinning shares at most about 12 % of its
   between-child variance with any diffusion axis (r ≤ 0.35), and the largest share is the partial-volume one.

## Step 5: a constructed diffusion predictor of thinning, and PRS on the diffusion axes
Scripts `05_pls_composite.py`, `06_pls_reliability_prs.R` and `07_pc_prs_adj.R`; tables `pls_fit.tsv`,
`pls_weights.tsv`, `pls_reliability.tsv`, `pls_prs.tsv`, `pca5_loadings.tsv` and `pc_prs_adj_thin.tsv`.

**PLS1 predictor.** The 45 site-residualised diffusion slopes predict the faster-thinning slope. It is cross-fitted
with 5 family-grouped folds and an inner CV choice of components; every fold chose 10, the cap.
- Out-of-fold r = 0.43 (R² 0.18), against 0.35 for the best single measure.
- The largest weights are grey-matter FA/TD/LD, superficial-WM free and hindered water, and hindered-water contrast.

**Reliability does not improve.** Each child's out-of-fold weights were applied scan by scan, and each series fitted
with the Figure-1 LMM. Same 8,366 children:

| series | mean slope reliability (2 / 3 / 4 visits) | site R² on slope | r with thickness slope |
|---|---|---|---|
| thickness | 0.296 (0.175 / 0.314 / 0.384) | 0.030 | 1 |
| diffusion composite | 0.178 (0.103 / 0.188 / 0.235) | 0.099 | 0.35 |
| 75 % thickness + 25 % composite | 0.296 | 0.025 | 0.93 |
| 50 / 50 | 0.259 | 0.045 | 0.75 |
| 25 % thickness + 75 % composite | 0.212 | 0.075 | 0.54 |

The composite is noisier from scan to scan than thickness, so any weight on it adds error at least as fast as signal.
It is also mostly the partial-volume axis (below), not an independent measurement of the same process.

**PRS (Figure-1 model + site; PRS-CS):**
- The composite slope is null (+0.005 to +0.018).
- The child-level out-of-fold prediction is nominal only for MDD EUR (+0.033).
- The 75/25 mix tracks thickness (+0.021 to +0.033), adding nothing.

**Diffusion PCs** (site-residualised, z-scored slopes of all 47 measures; signs arbitrary). Correlation with faster
thinning:

| PC | r with faster thinning |
|---|--:|
| dPC1 | −0.20 |
| dPC2 | +0.34 |
| dPC3 | −0.09 |
| dPC4 | −0.07 |
| dPC5 | −0.04 |

- **dPC1 (white-matter maturation: FA, RND, RNT vs TD): no PRS association** (\|β\| ≤ 0.015 in every cell).
- **dPC2 (grey-matter MD/TD/free water, the partial-volume axis):** MDD EUR +0.048 (p 0.002), MDD pooled +0.023 (p 0.04),
  SCZ pooled +0.022 (p 0.049); SBayesRC is weaker for SCZ. With the thinning slope as a covariate these fall to
  +0.037 (p 0.011), +0.015 and +0.012 (both n.s.). So the dPC2 association is largely the thinning association
  measured again, apart from part of the MDD EUR effect.
- **dPC4 (superficial-WM hindered water and grey/white contrast of MD/FA/RND):** SCZ only.
  - PRS-CS EUR +0.041 (p 0.010), pooled +0.023 (p 0.043); SBayesRC +0.041 (p 0.010) and +0.023 (p 0.043). The two
    scores correlate 0.96, so the near-identical values are expected.
  - Adjusting for thinning makes it slightly *larger* (EUR +0.047, p 0.003; pooled +0.026, p 0.021).
  - MDD is null.
  - Loading signs (this PCA): higher dPC4 = slower fall of superficial-WM hindered (extracellular) water, slower
    WM FA rise, and a faster rise of grey/white MD/FA/RND contrast.
  - Twenty score × arm tests on the PCs, so p 0.003 roughly meets Bonferroni (0.0025) only in the thinning-adjusted
    EUR fit. Treat it as the lead candidate for the twin pipeline and for HCP-MMP regional maps, not a finding.

## Step 6: is there one process behind thinning and diffusion PC4? (`08_single_factor.py`; `tables/single_factor_*.tsv`)
**The four predictions.** If one latent process drives both measures:
1. the child slopes correlate;
2. a composite is more reliable than either part;
3. PRS effects are proportional across scores (the tetrad = 0);
4. the twin rA is high.

Models: OLS with family-clustered SE and the Figure-1 covariates + site; PRS-CS.

| test | thinning × dPC4 | thinning × dPC2 (grey-matter MD / free water) |
|---|---|---|
| phenotypic r of child slopes | −0.07 (−0.10 in the twin sample) | +0.34 |
| twin h² (AE) of the diffusion axis | 0.67 [0.56, 0.80] | 0.50 [0.41, 0.59] |
| twin rA with faster thinning | **−0.09 [−0.25, +0.05]** | **+0.47 [+0.32, +0.61]** |
| dPC/thinning PRS ratio, SCZ vs MDD | EUR 0.86 vs 0.41; pooled 0.63 vs −0.34 | — |
| tetrad = 0? | EUR p 0.63; pooled p 0.12. Not rejected, but the bootstrap CIs are wide. | — |

**Composites, z(thinning) + w · z(dPC4)** (PRS β, z):

| score / arm | w = 0 (thinning) | 0.25 | 0.5 | 1 | 2 | dPC4 alone |
|---|---|---|---|---|---|---|
| SCZ EUR | +0.040 (2.9) | +0.049 (3.5) | +0.054 (3.8) | +0.056 (3.8) | +0.051 (3.3) | +0.035 (2.1) |
| SCZ pooled | +0.030 (3.0) | +0.034 (3.4) | +0.036 (3.5) | +0.036 (3.4) | +0.031 (2.8) | +0.019 (1.7) |
| SCZ non-EUR | +0.023 (1.6) | +0.025 (1.7) | +0.025 (1.7) | +0.023 (1.5) | +0.019 (1.2) | +0.009 (0.6) |
| MDD EUR | +0.027 (1.8) | +0.030 (2.0) | +0.030 (2.0) | +0.028 (1.9) | +0.023 (1.5) | +0.011 (0.7) |
| MDD pooled | +0.024 (2.2) | +0.021 (2.0) | +0.018 (1.7) | +0.011 (1.1) | +0.003 (0.3) | −0.008 (−0.7) |

**Reading.**
- **Thinning and dPC4 are not one process.** They are uncorrelated phenotypically and genetically. Measurement error
  goes to E in the twin model, so it cannot explain the near-zero rA.
- **Why the SCZ composite still gains z (2.9 → 3.8 in EUR).** It adds two independent signals; that is aggregation, not
  reliability.
- **The gain is also optimistic.** dPC4 was singled out partly for its SCZ association in these same EUR / pooled
  children. In the only held-out subsample (non-EUR, n ≈ 4,000), dPC4 alone is z 0.6 and the composite gain is 1.6 → 1.7.
- **The composite does not help MDD.** dPC4 carries no MDD signal, so adding it dilutes the MDD effect.
- **dPC2 is the axis that does share genetics with thinning** (rA 0.47). Its PRS association disappears once thinning
  is adjusted (step 5), and adding it to thickness lowers reliability (PLS composite). Together with the intercept
  pattern (thicker cortex ↔ lower grey-matter MD and free water), this fits dPC2 being the same cortical geometry
  seen through diffusion (partial volume), not a second proxy of a common biological process.
- **So at whole-cortex resolution, no measure tried so far shares genetic liability with thinning other than through
  thickness itself.**
  - T1w/T2w: rA 0.08 (D3).
  - dPC4: rA −0.09.
  - dPC2 (partial volume): rA 0.47.

  Either the disease-relevant process is specific to thinning, or a shared process is regional and diluted in
  whole-cortex means. The regional HCP-MMP measures from CSD3 test the second possibility.

## Step 7: scanner terms, and the HCP-MMP structural measures (`09`–`11`; `tables/scanner_terms_*.tsv`, `twin_smri.tsv`)
**Scanner changes within child.** In the Figure-1 sample, 27 % of children change scanner unit and 63.5 % change
software version between visits; there are 66 unit × software epochs. Epochs are partly confounded with age
(within-epoch age SD 0.88 y vs 2.45 y), so an epoch term can absorb real change.

Models:
- M0: Figure 1.
- M1: M0 + (1 | scanner unit).
- M2: M0 + (1 | unit × software epoch).

The structural measures use exactly the Figure-1 scans (26,946 scans, 8,716 children; same DAIRC surfaces). y is
z-scored by measure, and PRS models are Figure-1 + site.

| measure | reliability M0 → M2 | site R² of slopes M0 → M2 | r with faster thinning (within site) | twin h² (M0) | rA with faster thinning (M0) |
|---|---|---|---|---|---|
| thickness | 0.303 → 0.298 | 0.030 → 0.013 | — | 0.45 | — |
| area | 0.357 → 0.416 | 0.081 → 0.040 | −0.20 | 0.73 | **−0.36 [−0.51, −0.20]** |
| GM volume | 0.449 → 0.485 | 0.030 → 0.026 | −0.74 | 0.69 | −0.72 [−0.77, −0.66] |
| GWC | 0.180 → 0.194 | 0.104 → 0.015 | −0.32 | 0.56 | **−0.27 [−0.40, −0.13]** |
| T1 GM mid-depth | 0.198 → 0.119 | 0.198 → 0.056 | +0.10 | 0.48 | +0.09 [−0.07, +0.25] |
| T1 WM 1 mm | 0.179 → 0.086 | 0.153 → 0.030 | −0.19 | 0.32 | −0.08 [−0.30, +0.25] |
| T1w/T2w (DK) | 0.408 → 0.254 | 0.260 → 0.060 | 0.01 | 0.51 | −0.07 [−0.25, +0.10] |
| MD contrast (`is__md__gwc`) | 0.321 → 0.319 | 0.067 → 0.006 | −0.30 | 0.61 | −0.44 [−0.56, −0.32] |
| GM MD (`is__md__gm`) | 0.319 → 0.317 | 0.161 → 0.009 | +0.35 | 0.59 | +0.42 [+0.28, +0.55] |

Negative rA means children genetically inclined to faster thinning also tend to lose area / volume / GWC faster.

- **Thickness is robust to scanner.** The epoch term leaves its reliability and its PRS effects unchanged: SCZ EUR
  z 2.6 → 2.8, pooled 2.8 → 2.8; MDD z 2.1–2.9 throughout.
- **For the intensity and diffusion measures, much of the apparently reliable slope variance was scanner drift.**
  The epoch term cuts site R² of slopes to ≤ 0.06 and reliability falls: T1 wm1 0.18 → 0.09, FA wm 0.17 → 0.13,
  RND wm 0.12 → 0.06, T1w/T2w 0.41 → 0.25. Some fixed age slopes shift too (GM MD −0.04 → +0.03 SD/yr), which is the
  age–epoch confounding.
- **Grey-matter volume is the one measure that competes with thickness.** It is more reliable (0.45–0.49) and its SCZ
  signal is as strong or slightly stronger (higher score = faster volume loss):

  | score / arm | thickness (faster thinning) | GM volume (M0 → M2) |
  |---|---|---|
  | SCZ PRS-CS EUR | +0.036 (z 2.6) | −0.038 (z −2.8) → −0.043 (z −3.1) |
  | SCZ PRS-CS pooled | +0.028 (2.8) | −0.025 (−2.6) → −0.028 (−2.8) |
  | SCZ SBayesRC pooled | +0.031 (3.1) | −0.029 (−3.0) → −0.032 (−3.3) |
  | MDD PRS-CS EUR / pooled | +0.030 (2.1) / +0.028 (2.7) | −0.022 (−1.5) / −0.019 (−2.0) |

  Volume's extra reliability comes from area, which carries only a weak score effect in the same direction (SCZ
  −0.015 to −0.021, z −1.1 to −1.8; MDD EUR −0.025 to −0.032). So the PRS z barely moves even though reliability
  rises by half.
- **Thinning shares genetics with area loss** (rA −0.36), and modestly with GWC decline (rA −0.27, the direction the
  myelination account predicts). It does not share genetics with any intensity measure: T1 GM, T1 WM, T1w/T2w, and
  diffusion dPC4 (step 6) all have rA ≈ 0.
- **Caveat on GWC.** GWC samples grey matter at a fixed fraction of thickness and white matter at 1 mm, so thinner
  cortex may lower it mechanically. Its MDD association also runs the other way: MDD EUR +0.035, z 2.4 = slower GWC
  fall.
- **The MDD score lowers the T1 GM-intensity slope** (EUR −0.036, z −2.7; pooled −0.022, z −2.2; same under M2 and
  SBayesRC). This is independent of thinning genetically (rA +0.09). It is a lead, not a finding.

**Reading.** In these data, the cross-measure shared genetic liability with thinning is in **surface geometry**
(thickness, area, volume) and in boundary measures that depend on it (GWC, MD contrast, GM MD), not in tissue
intensity or myelin proxies. Volume is a reasonable co-primary phenotype for the SCZ score (equal or slightly
stronger). Thickness remains better for MDD and is the only measure significant in every cell, under every scanner
model.
