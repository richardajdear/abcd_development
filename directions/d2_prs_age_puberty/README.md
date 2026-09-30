# D2 — timing versus rate: one-stage PRS × age and PRS × puberty models

*Implements `docs/DIRECTIONS.md` D2 (prediction T1). 2026-09-30. Laptop only; nothing fetched.*

## Question

Figure 1 shows that polygenic risk for schizophrenia and depression predicts a
*faster* per-child rate of cortical thinning but not thickness at 12.8. Two
developmental readings fit that: (i) **rate scaling** — risk carriers thin faster
throughout 9–17; (ii) **a timing shift** — risk carriers are further along the
same trajectory (phase advance), so their slope differs only because they sit
at a different point of a curved curve. The tempo framing in DIRECTIONS §2
leans on (ii). The two are separable here because the normative trajectory is
curved within the ABCD age range and because a phase advance has to show up in
the *level* as well as the slope.

## What was run

`01_build_tables.py` → `02_fit_models.R` → `03_figure.R`, from the repo root
(env `abcd-spatial` for Python, `r` for R). Individual-level tables go to the
gitignored `genetic_analysis/work/results_70tab_hcp/d2_prs_age_puberty/`; only
`tables/*.tsv` and the figure are committed.

- **Trait.** Per-scan mean HCP-MMP thickness (unweighted mean of the 358
  parcels, as `fig1_prep_1lmm.R`), the settled specification
  (`out/thickness_hcp_70_aa6e91efba82`). 26,597 scans of the 8,596 genotyped
  children; 13,547 scans of the 4,308 EUR-arm children. Age centred at 12.797.
- **Scores.** SCZ 2025 and MDD, SBayesRC primary and PRS-CS as sensitivity;
  pooled arm = multi-ancestry weights z-scored within ancestry cluster
  (`_zanc`), EUR arm = EUR weights standardised within the EUR arm. Same
  score files as Figure 1. EA, ALZ, ASD scores are not on this machine (see
  *Gaps*).
- **One-stage LMM.** `mean_ct ~ PRS × f(age) + sex + PC1–10 + (1 + age_c |
  child) + (1 | site) + (1 | family)`, REML, bobyqa. The PRS main effect is the
  thickness difference at 12.8; `PRS:age_c` is the effect on the rate.
  Age bases: linear; hinge at 13 (`PRS:h13` = change in the rate effect after
  13 — the T1 test); quadratic; three segments with knots at 12 and 14
  (segment rate effects and SEs as linear combinations of the fixed effects).
- **Puberty at the scan.** PDS mean (youth report; parent as sensitivity) at
  the same session, centred, added as `PRS × pds_c` alongside `PRS × age_c`
  (and instead of it), on the 24,747 scans with PDS.
- **Pubertal timing as a child trait.** `PDS ~ age_c × sex + (1 + age_c |
  child) + (1 | site)` over all 59,845 yearly youth-report waves; the child
  intercept is PDS at 12.8 ("stage"; higher = earlier puberty), the child slope
  "tempo". Then `stage ~ PRS`, `tempo ~ PRS`, and the Figure 1 slope on PRS with
  and without stage + tempo, all with sex, PC1–10, `(1 | family)`.

## Results (pooled arm, SBayesRC; `tables/onestage_terms.tsv`, `prs_rate_by_age.tsv`)

**1. The rate effect replicates in one stage and there is no level effect.**
SCZ: `PRS:age_c` = −0.28 µm/yr per SD (SE 0.08, p = 7 × 10⁻⁴); MDD −0.17 (SE
0.08, p = 0.04). The SD of the child slope (√τ₁) is 4.2 µm/yr, so these are
−0.067 and −0.040 SD of the *unshrunk* slope per SD of score — about twice the
Figure 1 βs (−0.034, −0.025), which are on the shrunk-BLUP scale; the two are
consistent. The PRS main effect (thickness at 12.8) is −0.12 µm (SE 0.76) for SCZ
and −0.28 µm (SE 0.77) for MDD. EUR arm: SCZ −0.25 µm/yr (p = 0.02), MDD −0.17
(p = 0.12), levels +1.4 / −1.5 µm (n.s.). PRS-CS gives the same numbers to within
0.03 µm/yr.

**2. The normative trajectory is curved, so a timing shift was testable.** The
fixed rate is −18.4 µm/yr before 12, −23.9 in 12–14, −17.1 after 14 (quadratic
term +0.093 µm/yr², p = 0.005; hinge at 13 +0.56, p = 0.049). Thinning of the
cortical mean peaks in 12–14 within this sample.

**3. The PRS effect on the rate does not depend on age (T1 fails).** Change in
the rate effect after 13: SCZ +0.16 µm/yr/SD (SE 0.29, p = 0.57), MDD +0.05
(SE 0.28, p = 0.86); quadratic PRS × age² term p = 0.30 / 0.86. By segment (SCZ):
<12 −0.69 (SE 0.24, p = 0.005), 12–14 +0.29 (SE 0.34, p = 0.40), ≥14 −0.49 (SE
0.28, p = 0.08); MDD −0.28 / +0.03 / −0.27, none significant. A phase advance
across a rate peak would make the effect *negative before the peak and
positive after it*; the observed effect is same-signed at both ends and, if
anything, weaker in the middle. The middle-segment dip is not significant
(SCZ `PRS:h12` p = 0.06) and is not reproduced by MDD.

**4. A phase advance is bounded to about one month per SD.** On a trajectory
thinning at 20 µm/yr, being Δ years ahead thins the cortex by 20Δ µm at every
age. The 95 % upper bound on the level effect is 1.6 µm (SCZ), 1.8 µm (MDD), so
Δ < 1.1 and < 1.0 months per SD of score. The rate effect over 9–17 (−0.28
µm/yr × 8 yr = −2.2 µm per SD) is *larger* than any level difference at the
centre — it is divergence that accumulates through the window, not a head
start.

**5. Pubertal stage does not carry the PRS effect on the rate.** On the
24,747 scans with youth PDS, adding `PRS × pds_c` leaves `PRS × age_c` intact
(SCZ −0.24 → −0.30 µm/yr, MDD −0.18 → −0.34) and the puberty interaction itself
is null (SCZ +0.16, SE 0.44; MDD +0.30, SE 0.45). Only when the age interaction
is *removed* does `PRS × pds_c` pick up the effect (SCZ −0.55, p = 0.03; MDD
−0.53, p = 0.04), because pubertal stage tracks age. Parent-report PDS gives the
same picture. Pubertal stage itself matters for thickness: −7.3 µm per PDS
unit at fixed age (p < 10⁻⁵⁰).

**6. Pubertal timing relates to both PRS and thinning, but is not the route
(`tables/puberty_tempo.tsv`).** MDD PRS predicts earlier puberty: stage at
12.8 +0.055 SD/SD (p = 8 × 10⁻⁷), tempo −0.065 (p = 6 × 10⁻¹⁰ — the linear PDS
slope is *lower* in advanced children because PDS saturates at 4, so a low
slope here means "already near ceiling", not "slow"; stage and slope BLUPs
correlate −0.49). SCZ PRS: stage +0.015 (p = 0.17), tempo −0.031 (p = 0.003).
Earlier puberty goes with faster thinning (stage → slope −0.068 SD/SD, p = 6 ×
10⁻⁹) and slightly thinner cortex at 12.8 (−0.025, p = 0.04). But adjusting the
Figure 1 slope for stage and tempo leaves the PRS effect unchanged: SCZ −0.034 →
−0.034, MDD −0.025 → −0.023. Puberty is a parallel correlate, not a mediator.

## Reading

Within 9–17 the polygenic effect is a **uniform scaling of the thinning rate**,
not a shift in developmental timing, and not a pubertal effect. This closes
the "phase advance" version of the tempo hypothesis in DIRECTIONS §2.1 (T1) and
leaves the "amplitude" version: risk carriers lose more cortex per year across
the whole window, with the deficit accumulating (≈2 µm per SD over eight years
on a 2.8 mm cortex). What that is at tissue level is D1's question; whether it
continues past 17 or is a phase advance on a longer timescale (thinning that
would otherwise happen after 17 brought into the window) ABCD cannot yet say —
the 8-year follow-up will.

## Gaps

- **EA, ALZ ± APOE, ASD** scores are on CSD3 only. The T3 test (EA opposite in
  sign on every readout) needs them. To run the same models, pull the
  profiles into the same relative locations and re-run `01`–`03`:
  ```
  rsync -av --include='*/' --include='score_*.profile' --exclude='*' \
    rajd2@login.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/ \
    genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/
  ```
  then add the new profile paths to `SCORES` in `01_build_tables.py`.
- Pubertal "tempo" from a linear PDS slope is confounded with stage by the
  PDS ceiling; a growth-curve (logistic) age-at-midpoint would be the clean
  tempo measure if it is ever needed. Stage at 12.8 is the measure used here.
- lmerTest Satterthwaite p-values; the ~50 convergence warnings are the usual
  large-n `checkConv` messages, not singular fits (all variance components are
  positive; `tables/variance_components.tsv`).

## Files

| file | content |
|:--|:--|
| `01_build_tables.py` | scan-level and child-level tables (gitignored output) |
| `02_fit_models.R` | all models; writes `tables/` |
| `03_figure.R` | `figures/fig_d2_prs_age_puberty.png`, reading only `tables/` |
| `tables/onestage_terms.tsv` | every fixed-effect term, every model, both arms and methods |
| `tables/prs_rate_by_age.tsv` | PRS effect on the rate by age segment, with SE and per-slope-SD version |
| `tables/variance_components.tsv` | τ₀, τ₁, site, family, residual per model; PDS LMM |
| `tables/puberty_tempo.tsv` | PRS → pubertal stage / tempo; stage / tempo → thinning; PRS → thinning given puberty |
| `tables/trajectories.tsv` | model-implied mean thickness 9–17 at PRS −1 / 0 / +1 SD |
