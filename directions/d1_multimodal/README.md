# D1 — multimodal slopes (myelination vs neuropil; composite slope)

See `docs/DIRECTIONS.md` §2.2, §2.4 and D1. Each modality gets the Figure 1 treatment: a whole-cortex mean per scan and
one single LMM `mean ~ age_c + sex + (1 + age_c | subject) + (1 | site)`. The per-child slope then goes into the
Figure 1 PRS model, the D3 twin pipeline and, eventually, a cross-modal composite.

| step | file | status |
|---|---|---|
| 1 | `01_t1t2_prs.R`: PRS → T1w/T2w slope, plus the quality of the measure (DK; HCP-MMP T1w/T2w does not exist yet) | done 2026-10-27 |
| — | `CSD3_INSTRUCTIONS.md`: HCP-MMP and DK GWC, T1/T2 intensity at GM mid-depth and in superficial WM, T1w/T2w ratio, and DTI/RSI if processed maps exist | brief written |
| 2 | `02_dmri_scan_means.py`: per-scan whole-cortex means for 45 DK DTI/RSI release tables | done 2026-10-27 |
| 3 | `03_fit_compare.R`: single-LMM slopes, quality, coupling with thinning / T1w/T2w, PRS | done 2026-10-27 |
| 4 | `04_followups.R`: PRS adjusted for the thinning slope; PCA of the 47 slopes | done 2026-10-27 |
| — | area, volume, GWC (CSD3) → same pipeline; twin h² / rA for the lead diffusion measures; composite via multivariate LMM | to do |

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
