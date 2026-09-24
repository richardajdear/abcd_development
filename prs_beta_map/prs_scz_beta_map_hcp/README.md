# prs_scz_beta_map_hcp — the SCZ polygenic score, parcel by parcel, against the transcriptomic maps

The integrating analysis from the status review: regress each HCP-MMP parcel's
thinning slope on the primary SCZ score, spin-test the resulting β map against
AHBA C3 and the ahba_pls lead component. This is the figure that connects the
individual-level genetics (`genetic_analysis/`) to the transcriptomic arm
(`ahba_pls/`). HCP-MMP throughout: 358 parcels (H excluded), 179 left-hemisphere
parcels rotated in the spin, 137 of them AHBA-covered.

## Status — 2026-09-19: run complete, both cells. **Null.**

The β map of the SCZ score does not resemble C3 or the ahba_pls lead
component. Figures `results/fig_prs_beta_map_<cell>.png`, captions alongside,
all numbers in `results/spin_tests.tsv` (5,000 spins).

| cell | whole-cortex β (reproduces table) | β vs C3 | β vs lead | β, cortex mean held fixed: vs C3 / vs lead | loading map vs C3 | lh–rh ρ |
|:--|--:|--:|--:|--:|--:|--:|
| pooled, `SCZ25_META_SBayesRC_zanc` (n 8,596) | −0.0276 (p 0.0079) ✓ | −0.08, p_spin 0.48 | +0.12, 0.33 | +0.02, 0.85 / +0.16, 0.16 | **+0.36, 0.021** | 0.29 |
| EUR, `SCZ25_EUR_SBayesRC_raw` (n 4,308) | −0.0418 (p 0.0043) ✓ | +0.03, 0.80 | +0.20, 0.067 | +0.14, 0.12 / +0.26, **0.008** | **+0.36, 0.018** | 0.33 |

Reading:

- **The parcel β map is mostly noise at this n.** Homologous left and right
  parcels agree at ρ ≈ 0.3, and 292–298 of 358 βs are negative: the score's
  effect is a near-uniform shift of the whole cortex (the tabled global β),
  with little reproducible spatial structure on top. 73 (pooled) / 50 (EUR)
  parcels reach p < 0.05, against 18 expected.
- **Raw β vs C3 and vs the lead component: null in both arms**, lh and
  bilateral-mean variants alike (|ρ| ≤ 0.20, all p_spin ≥ 0.07).
- **The one nominal hit** — EUR, parcel-specific β (cortex mean held fixed) vs
  the lead component, ρ +0.26, p_spin 0.008 — is one of 36 tests, the lh-only
  variant (bilateral mean ρ +0.15, p 0.22), and absent in the larger pooled arm
  (+0.16, p 0.16). It is recorded, not led with.
- **The loading map is C3-like** (each parcel's correlation with the
  whole-cortex mean slope vs C3: ρ +0.36, p_spin ≈ 0.02, in both arms and in
  the synthetic test). Any score that only shifts the whole cortex would
  therefore draw a weakly C3-shaped β map; here even that is not visible in
  the raw β (ρ −0.08), because the map's noise swamps it.

What it licenses: the SCZ score → faster thinning association is a
whole-cortex effect; the parcel-level genetics do not (yet) connect to the
transcriptomic arm. The imaging→genes bridge stays the group-level one
(thinning map vs C3, ρ −0.55 p_spin 0.002), and the individual-level bridge is
the whole-cortex PRS β. A more powerful version of this test would project
each child's 358 slopes onto C3 / the lead component (one phenotype per child,
`ahba_pls/FOLLOWUP_GENETICS.md` H4) and regress that on the score, instead of
358 noisy βs.

Cluster run: CSD3 job 35810927 (two tasks, ~8 min each; log in
`slurm/prs_beta_map_35810927_*.log`); first submission 35810839 failed because
the GENESIS R's data.table cannot open .gz files (fixed with a gzip pipe).

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
