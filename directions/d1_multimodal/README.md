# D1 — multimodal slopes (myelination vs neuropil; composite slope)

See `docs/DIRECTIONS.md` §2.2, §2.4 and D1. Each modality gets the Figure 1 treatment: a whole-cortex mean per scan and
one single LMM `mean ~ age_c + sex + (1 + age_c | subject) + (1 | site)`. The per-child slope then goes into the
Figure 1 PRS model, the D3 twin pipeline and, eventually, a cross-modal composite.

| step | file | status |
|---|---|---|
| 1 | `01_t1t2_prs.R`: PRS → T1w/T2w slope, plus the quality of the measure (DK; HCP-MMP T1w/T2w does not exist yet) | done 2026-10-27 |
| — | `CSD3_INSTRUCTIONS.md`: HCP-MMP and DK GWC, T1/T2 intensity at GM mid-depth and in superficial WM, T1w/T2w ratio, and DTI/RSI if processed maps exist | GWC, T1 gm50, T1 wm1 done on CSD3 2026-10-27 (`CSD3_REPORT.md`); T2, T1w/T2w, DTI held for a decision |
| 2 | slopes for area, volume, GWC, DTI/RSI → PRS, twin h² and rA with ΔCT | waiting for data |
| 3 | composite slope: multivariate LMM with correlated random slopes **and** correlated residuals; reliability of the common factor vs ΔCT alone | waiting for data |

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
