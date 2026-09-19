# prs_scz_beta_map_hcp — the SCZ polygenic score, parcel by parcel, against the transcriptomic maps

The integrating analysis from the status review: regress each HCP-MMP parcel's
thinning slope on the primary SCZ score, spin-test the resulting β map against
AHBA C3 and the ahba_pls lead component. This is the figure that connects the
individual-level genetics (`genetic_analysis/`) to the transcriptomic arm
(`ahba_pls/`). HCP-MMP throughout: 358 parcels (H excluded), 179 left-hemisphere
parcels rotated in the spin, 137 of them AHBA-covered.

## Status — 2026-09-19

| step | where | state |
|:--|:--|:--|
| 1 per-child, per-parcel slopes exported (8,716 × 358; `work/`, gitignored) | laptop | **done** |
| 2 per-parcel mixed-model association with the score | CSD3 | **scripted and smoke-tested, not run** — the per-subject scores live only on CSD3 (README_HPC rule 16) and the login needs an interactive MFA |
| 3 spin tests (5,000 rotations, `ahba_pls/code/pls.py::spin_corr`, HCP centroids) | laptop | scripted, tested |
| 4 figure + caption | laptop | scripted, rendered on a **synthetic** score to check the code path; no real figure yet |

The whole laptop→cluster→laptop loop was exercised end to end with a synthetic
score against the DK export's covariates (same 8,716 children); nothing from that
run is kept in this folder.

## The cells

The "primary SCZ score" is read as the step-9 primary cells of
`genetic_analysis/README_HPC.md` (2025 multi-ancestry SCZ GWAS, SBayesRC):

| cell | score | arm | tabled HCP `global_slope` β to reproduce |
|:--|:--|:--|--:|
| `SCZ25_META_SBayesRC_zanc` (**primary**) | META weights, z within ancestry cluster | pooled, n 8,596 | −0.0276 (p_adj 0.008) |
| `SCZ25_EUR_SBayesRC_raw` | EUR weights, raw | EUR anchor, n 4,308 | −0.0418 (p_adj 0.004) |

Each `beta_map_<cell>.tsv` carries a `cortex_mean` row: the mean of the 358 parcel
slopes is the per-region `global_slope`, so that row must match the tabled β
before the map is read.

## Model

`tools/prs_assoc.R`'s model, unchanged, on each parcel's standardised slope:

    scale(slope_parcel) ~ scale(PRS) + sex + age_c + PC1..PC10 + (1 | family_id)

plus a second fit per parcel with `scale(cortex_mean)` as a covariate
(`beta_cond`): the parcel-specific association with the whole-cortex effect held
fixed. A score that only shifts the whole cortex still draws a structured β map
— each parcel's loading on the mean (`r_global`) — and the synthetic test showed
exactly that (a pure noise-plus-global score gave β vs C3 ρ = −0.24). So the
figure reports both, and the loading map's own correlation with C3 is the null
to beat.

## Run

```bash
# laptop
python prs_beta_map/prs_scz_beta_map_hcp/01_export_parcel_slopes.py
ssh -o ControlPersist=4h login-q-1.hpc.cam.ac.uk          # authenticate once, leave open
bash prs_beta_map/prs_scz_beta_map_hcp/run_on_csd3.sh push  # rsync slopes + scripts, sbatch (2 cells, ~10 min)
bash prs_beta_map/prs_scz_beta_map_hcp/run_on_csd3.sh status
bash prs_beta_map/prs_scz_beta_map_hcp/run_on_csd3.sh pull  # results/beta_map_*.tsv (358 rows each; no per-subject data)
python prs_beta_map/prs_scz_beta_map_hcp/03_spin_test.py
python prs_beta_map/prs_scz_beta_map_hcp/04_figure.py            # primary cell
python prs_beta_map/prs_scz_beta_map_hcp/04_figure.py --cell SCZ25_EUR_SBayesRC_raw
```

`python` is `~/mambaforge/envs/abcd/bin/python` (pyarrow, scipy, matplotlib);
the R side needs data.table, lme4, lmerTest (no optparse).

## Outputs

- `results/beta_map_<cell>.tsv` — label, beta, se, t, p, n, beta_cond…, r_global
- `results/spin_tests.tsv` — cell × map (beta / beta_cond / r_global) × variant (lh / bilateral mean) × reference (C3 / PLS_lead / dCT): ρ, p_spin, plus cortex β and lh–rh agreement
- `results/fig_prs_beta_map_<cell>.png` + `.md` caption (numbers read from the tables)

## Figure design

dataviz skill: diverging blue ↔ red with a neutral grey midpoint for the signed
maps (palette reference pair); one series in the scatters (ink points, slot-1
blue fit line); C3 and the lead component drawn on the same flat polygons
(`ahba_pls/data/hcp_polygons.csv`, ggsegGlasser) with uncovered parcels grey;
spin null of ρ as an inset strip; 180 mm wide, 7/6/5.5 pt type as
`genetic_analysis/fig1_draft.py`.
