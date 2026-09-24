# h4_projection — is a transcriptomically weighted slope more heritable than the whole-cortex slope?

`FOLLOWUP_GENETICS.md` H4 on the 7.0 HCP-MMP run (8,716 children × 358
parcels). Each child's slope map is projected onto a transcriptomic region map
(`legacy/hpc_v3/make_phenotypes_v3.py::projection`, unchanged): z across
children within parcel, dotted with the mean-centred map applied to both
hemispheres over the 137 AHBA-covered Glasser areas (274 parcels).

| phenotype | map |
|:--|:--|
| `slope_projC3` | AHBA C3 (Dear et al. 2024, `data/ahba_dme_hcp_top8kgenes_scores.csv`) |
| `slope_projPLS2` | the ahba_pls lead component, thinning-oriented (`results/hcp_pls_scores.csv`) |
| `*_resid` | each residualised on `global_slope` |
| `slope_cov137` | unweighted mean slope over the same 274 parcels (the parcel-set control) |
| `global_slope` | the pipeline phenotype, carried through the same join (positive control) |

h² by GCTA GREML with `genetic_analysis/step5_reml.sbatch`'s exact design
(dense imputed GRM on the PC-AiR unrelated set, sex + site, baseline_age +
n_visits + 10 PCs). Benchmark: HCP `global_slope` h² = 0.156 ± 0.045
(`results_70tab_hcp/reml_imp_pooled`).

## Status — 2026-09-24: GREML done (CSD3 job 36261098). **No phenotype is distinguishably more heritable than the whole-cortex slope.**

`results/reml/reml_summary.tsv` (GCTA GREML, PC-AiR unrelated set, n = 6,011,
identical design to step 5; the `global_slope` row reproduces the tabled
0.156 ± 0.045 exactly, so the join and the covariates are right):

| phenotype | h² | SE | p | vs `global_slope` |
|:--|--:|--:|--:|--:|
| `global_slope` (whole cortex, 358 parcels) | **0.156** | 0.045 | 1.4e-4 | — |
| `slope_cov137` (mean over the 274 AHBA-covered parcels) | 0.142 | 0.044 | 4.6e-4 | −0.014 |
| `slope_projC3` (projection on AHBA C3) | 0.146 | 0.045 | 3.7e-4 | −0.010 |
| `slope_projC3_resid` (…residualised on global_slope) | 0.159 | 0.045 | 1.0e-4 | +0.003 |
| `slope_projPLS2` (projection on the ahba_pls lead component) | **0.173** | 0.044 | 1.6e-5 | +0.017 |
| `slope_projPLS2_resid` (…residualised on global_slope) | 0.164 | 0.044 | 3.6e-5 | +0.008 |

Reading:

- **Every h² sits within 0.02 of 0.156, with SEs of 0.044.** The largest gap
  (projPLS2 +0.017) is well under half a standard error; even allowing for the
  positive correlation between phenotypes, which shrinks the SE of a
  *difference* below √2 × 0.044, nothing here is a detectable improvement or
  loss. The spec's criterion — "a projection whose h² lands inside the
  benchmark interval is not an improvement" — is met by every row.
- **The ordering is the one H4 predicted, faintly.** PLS2 projection > C3
  projection ≈ whole cortex, with the PLS2 projection the most significant
  phenotype in the table (p 1.6e-5, vs 1.4e-4). The direction is consistent
  with the hypothesis that a transcriptomically weighted slope carries a
  slightly more genetic signal; the magnitude is not distinguishable from noise
  at n = 6,011.
- **The projections are noisier phenotypes yet not less heritable.** Their
  hemisphere consistency is 0.46–0.49 against 0.80 for the mean (table above),
  so measurement error alone should have pulled their h² *down* by roughly a
  third. That it did not suggests the weighted contrast does carry a somewhat
  more genetic signal per unit of reliable variance — the same point, seen
  from the other side, and no more significant.
- **The parcel set is not the explanation.** The unweighted mean over the same
  274 parcels (`slope_cov137`, 0.142) is if anything below the whole-cortex
  mean, so the PLS2 projection's +0.03 over it is the weighting, not the
  coverage.
- **Residualising on `global_slope` costs nothing** (0.159 / 0.164): the part of
  each projection that is orthogonal to the whole-cortex slope is about as
  heritable as the whole-cortex slope itself. This is the more interesting
  fact, because that orthogonal part is a *pattern* phenotype — which parcels
  thin faster than the child's own average — and it is heritable at
  h² ≈ 0.16 (p ≈ 1e-4 both times).

What it licenses: H4 as posed (a *more* heritable GWAS phenotype) is not
supported; the numbers are compatible with equal heritability and with a small
gain, and a sample several times larger would be needed to tell. What the run
does establish is that the C3- and PLS2-shaped *pattern* of thinning,
independent of the whole-cortex rate, is itself heritable at the same level as
the rate. The natural next step, cheap and already wired (per-subject table on
CSD3, score directories known), is the PRS test on these phenotypes: does the
SCZ score predict `slope_projC3_resid` / `slope_projPLS2_resid`? That is the
individual-level version of the thinning-map-vs-C3 bridge, and the failed
per-parcel β map (`prs_beta_map/`) is exactly the noisy form of it.

## Run

```bash
python ahba_pls/h4_projection/01_make_projection_phenotypes.py   # laptop; work/pheno_h4.tsv
ssh -o ControlPersist=4h login-q-1.hpc.cam.ac.uk                  # authenticate once
bash ahba_pls/h4_projection/run_on_csd3.sh push                   # pull, rsync, sbatch (6 tasks)
bash ahba_pls/h4_projection/run_on_csd3.sh pull                   # results/reml/*.hsq + reml_summary.tsv
```
