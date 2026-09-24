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

## Status — 2026-09-24: phenotypes built and characterised; GREML pending on CSD3

`results/h4_phenotype_characterisation.tsv`:

| phenotype | r with global_slope | r with projC3 | lh/rh r (Spearman–Brown) |
|:--|--:|--:|--:|
| `slope_projC3` | 0.44 | 1 | 0.49 (0.66) |
| `slope_projPLS2` | 0.19 | 0.69 | 0.46 (0.63) |
| `slope_cov137` | 0.99 | 0.45 | 0.80 (0.89) |
| `global_slope` | 1 | 0.44 | 0.80 (0.89) |

The projections are contrasts, so they are noisier phenotypes than the mean
(hemisphere consistency 0.46–0.49 against 0.80). The C3 and PLS2 region maps
correlate ρ = 0.70. Split-half circularity check on the PLS map
(`results/h4_split_half_check.tsv`): the lead component refitted on the group
means of a random half of the children correlates ρ = 0.9999 with the
full-sample one, and half-B children's projections on half-A weights correlate
r = 0.99998 with their full-weight projections — group means over 4,358
children are stable enough that fitting the PLS on the same children's means
cannot inflate h².

## Run

```bash
python ahba_pls/h4_projection/01_make_projection_phenotypes.py   # laptop; work/pheno_h4.tsv
ssh -o ControlPersist=4h login-q-1.hpc.cam.ac.uk                  # authenticate once
bash ahba_pls/h4_projection/run_on_csd3.sh push                   # pull, rsync, sbatch (6 tasks)
bash ahba_pls/h4_projection/run_on_csd3.sh pull                   # results/reml/*.hsq + reml_summary.tsv
```
