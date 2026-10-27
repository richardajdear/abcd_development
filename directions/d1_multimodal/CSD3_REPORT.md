# D1 on CSD3: report (measures 1–5 done; 6–7 held)

2026-10-27. Brief: `CSD3_INSTRUCTIONS.md`. Discovery and pilot: `CSD3_DISCOVERY.md`. Code:
`src/abcd/surf_metrics.py`, `tools/surf_metrics.sbatch`, `tools/surf_metrics_collect.py`. Helper scripts in this
directory: `csd3_discovery.py`, `csd3_pilot.sbatch`, `csd3_pilot_gates.py`, `csd3_report_stats.py`.
Slurm jobs: discovery 39326028, pilot 39327812, full array 39327815 (169 × 200 sessions, all COMPLETED), collect
39331883.
Measures 4–5 (§ below): `tools/t2reg.sbatch` (registration + T2w bias field), `csd3_t2reg_speed.sbatch`,
`csd3_t2reg_compare.py`, `csd3_bias_check.py`. Registration array 39352586 (144 × 200) plus fixes 39376453 and
39377622; T2 sampling 39380802 (169 × 200); collect 39380803.

## What was produced
| # | measure | stem (HCP / DK twin) | status |
|---|---|---|---|
| 1 | GWC (`surf/?h.w-g.pct.mgh`) | `mr_y_smri__gwc__hcp` / `__dsk` | **done** |
| 2 | T1w (`nu.mgz`) at GM mid-depth | `mr_y_smri__t1__gm50__hcp` / `__dsk` | **done** |
| 3 | T1w (`nu.mgz`) 1 mm below white | `mr_y_smri__t1__wm1__hcp` / `__dsk` | **done** |
| 4 | T2w (fMRIPrep `desc-preproc_T2w`, N4-corrected here) at GM mid-depth and 1 mm below white | `mr_y_smri__t2__gm50__hcp`, `mr_y_smri__t2__wm1__hcp` / `__dsk` | **done** (28,689 sessions) |
| 5 | T1w/T2w at GM mid-depth (fMRIPrep T1w / corrected T2w) | `mr_y_smri__t1t2__gm50__hcp` / `__dsk` | **done** (28,689 sessions) |
| 6–7 | DTI FA/MD, RSI | — | **skipped per brief.** Tensor maps exist (24,792 sessions, qsiprep T1w space), but there is no registration to the DAIRC subject. RSI does not exist on CSD3. |

Area and GM volume were not recomputed. `processed/hcp/mr_y_smri__{area,vol}__hcp.tsv` (33,795 rows each) and
`processed/dsk_local/mr_y_smri__{area,vol}__dsk.tsv` are present.

Method: as in the brief; details in `PROVENANCE_surf_metrics.txt`. The intensity source for T1 is `mri/nu.mgz`. Parcel
values are vertex-area-weighted means; hemisphere and whole-cortex values are area-weighted over all non-medial-wall
vertices. The HCP annotation is resolved as in `abcd.hcp_stats`: backfill for 8,874 sessions, the FS `label/` copy
for 24,921. None was regenerated.

## Sessions processed
| measure | ok | failed |
|---|--:|---|
| gwc | 33,794 | 1: `sub-2ZLBBMTU/ses-04A`, missing input (`surf/?h.w-g.pct.mgh` absent; that recon has no `recon-all.done` CMDARGS either) |
| t1__gm50 | 33,795 | 0 |
| t1__wm1 | 33,795 | 0 |

Per visit (ok): ses-00A 11,749 · ses-02A 8,093 · ses-04A 6,348 (gwc 6,347) · ses-06A 7,605.
No non-finite samples. The white-surface normal agrees with the white→pial direction at a median 99.97 % of
vertices per hemisphere (worst hemisphere 98.9 %).

## Gates (pilot, 20 sessions; details in `CSD3_DISCOVERY.md`)
- **A** GWC vs FreeSurfer `w-g.pct.stats`: r = 1.000000, max |Δ| = 0.0001 (DK; also HCP vs rr480's HCP segstats). PASS.
- **B** HCP vertex counts = `.log` `nverts`: 7,200/7,200 exact. PASS.
- **C** vs `mri_vol2surf` (trilinear): gm50 vs the same midpoint surface r = 1.00000; wm1 vs `--projdist -1`
  r = 0.993–0.995. PASS. (`--projfrac 0.5` gives r ≈ 0.95 because it is a different point.)
- **D** DK t1__gm50 vs release `mr_y_smri__t1__gm__dsk`: per-region r median 0.980, min 0.896; whole-cortex r = 0.994.
  PASS. No release DTI/RSI tables on CSD3.
- **E** 1.6 s/session for all three measures → ~15 core-hours in total.

## Whole-cortex mean by visit (median [IQR])
| measure | ses-00A | ses-02A | ses-04A | ses-06A | expected |
|---|---|---|---|---|---|
| gwc__hcp | 33.06 [31.56–37.01] | 32.51 [31.05–36.69] | 31.84 [30.47–35.97] | 31.44 [30.08–35.96] | falls ✔ |
| gwc__dsk | 33.21 [31.72–37.18] | 32.66 [31.20–36.86] | 32.00 [30.63–36.14] | 31.59 [30.23–36.13] | falls ✔ |
| t1__gm50__hcp | 71.01 [66.64–72.27] | 71.53 [67.05–72.75] | 72.12 [67.89–73.30] | 72.53 [67.97–73.77] | rises ✔ |
| t1__gm50__dsk | 70.91 [66.52–72.17] | 71.43 [66.93–72.65] | 72.01 [67.76–73.19] | 72.42 [67.86–73.66] | rises ✔ |
| t1__wm1__hcp | 103.46 [102.69–103.99] | 103.79 [103.10–104.30] | 104.16 [103.55–104.62] | 104.40 [103.77–104.86] | rises ✔ |
| t1__wm1__dsk | 103.44 [102.66–103.97] | 103.76 [103.07–104.28] | 104.12 [103.51–104.59] | 104.37 [103.74–104.83] | rises ✔ |
| thk__hcp (ref.) | 2.72 [2.67–2.77] | 2.68 [2.63–2.74] | 2.64 [2.59–2.69] | 2.61 [2.56–2.66] | falls ✔ |

## Scanner sensitivity: the main caveat
Share of variance in the session-level whole-cortex mean explained by site, and by manufacturer (one-way R²,
`ab_g_dyn__design_site` / `__design_mr__manufact`, all sessions):

| measure | R² site | R² manufacturer |
|---|--:|--:|
| gwc (hcp / dsk) | 0.740 / 0.741 | 0.731 / 0.732 |
| t1__gm50 (hcp / dsk) | 0.827 / 0.828 | 0.814 / 0.815 |
| t1__wm1 (hcp / dsk) | 0.476 / 0.480 | 0.452 / 0.456 |
| thk__hcp (reference) | 0.050 | 0.037 |

The **level** of every intensity measure is dominated by scanner. Almost all of the site variance is manufacturer,
which is why the per-visit IQRs above are so wide. This is much stronger than the site share already found for the
release T1w/T2w (`README.md`: scan-level site ICC 0.19). `nu.mgz` is bias-corrected and globally rescaled but not
contrast-normalised, so manufacturer differences in T1 contrast pass through. GWC is a ratio and should partly cancel
this, but it still has R² 0.74. For the laptop LMMs this argues for:
- manufacturer (or site) in the model, and a check of slopes **within manufacturer**, especially for children who
  change scanner between visits;
- treating the T1 level as unusable across scanners, so that only within-child change carries signal.

This is a descriptive sanity figure; no model was run on CSD3.

## Output files (CSD3; nothing committed)
Release dir: `/home/rajd2/rds/hpc-work/abcd_development/abcd-data-release-7.0/`, which resolves to
`/rds/user/rajd2/hpc-work/...` (same path as in the brief).

| file | rows × cols | size |
|---|---|--:|
| `processed/hcp/mr_y_smri__gwc__hcp.tsv` | 33,794 × 365 | 98.8 MB |
| `processed/hcp/mr_y_smri__t1__gm50__hcp.tsv` | 33,795 × 365 | 98.9 MB |
| `processed/hcp/mr_y_smri__t1__wm1__hcp.tsv` | 33,795 × 365 | 109.8 MB |
| `processed/dsk_local/mr_y_smri__gwc__dsk.tsv` | 33,794 × 73 | 19.9 MB |
| `processed/dsk_local/mr_y_smri__t1__gm50__dsk.tsv` | 33,795 × 73 | 19.9 MB |
| `processed/dsk_local/mr_y_smri__t1__wm1__dsk.tsv` | 33,795 × 73 | 22.0 MB |
| `processed/hcp/mr_y_smri__t2__gm50__hcp.tsv` | 28,689 × 365 | 98.6 MB |
| `processed/hcp/mr_y_smri__t2__wm1__hcp.tsv` | 28,689 × 365 | 98.0 MB |
| `processed/hcp/mr_y_smri__t1t2__gm50__hcp.tsv` | 28,689 × 365 | 73.6 MB |
| `processed/dsk_local/mr_y_smri__t2__gm50__dsk.tsv` | 28,689 × 73 | 19.8 MB |
| `processed/dsk_local/mr_y_smri__t2__wm1__dsk.tsv` | 28,689 × 73 | 19.7 MB |
| `processed/dsk_local/mr_y_smri__t1t2__gm50__dsk.tsv` | 28,689 × 73 | 14.9 MB |
| `processed/hcp/surf_metrics_session_qc.tsv` | 202,770 (session × measure; T2 rows add `bbr_mincost`) | 19.8 MB |
| `processed/hcp/PROVENANCE_surf_metrics.txt` | | <0.1 MB |
| `processed/hcp/surf_metrics_inventory.tsv` (discovery) | 33,795 | 6.1 MB |
| `processed/hcp/surf_metrics_pilot_sessions.tsv`, `surf_metrics_pilot/{gates.txt,report_stats.md,t2reg_pilot.tsv,t2reg_speed.tsv}` | | small |

Registrations and T2w bias fields (not release outputs; needed only to resample):
`legacy/hpc/work/d1_t2reg/<sub>/<ses>/{bbr.lta, bbr.dat.mincost, coreg.lta, t2_biasfield_8mm.nii.gz}`, with
per-session timing and cost in `legacy/hpc/work/d1_t2reg/status/` (first pass in `status_pass1/`).

Columns: `participant_id`, `session_id`, `<stem>__<region>__<hemi>_mean` (HCP region tokens exactly as in
`mr_y_smri__thk__hcp.tsv`, e.g. `V1`, `s6-8`; DK release codes, e.g. `bstmps`), `<stem>__lh_mean`, `<stem>__rh_mean`,
`<stem>_mean`. Not yet wired into `Release70Adapter` (its `_imaging_hcp` knows thk/area/vol only).

Rsync (as in the brief). It also picks up `surf_metrics_parts/*_qc.tsv` (169 small per-chunk QC files) and the pilot
TSVs. The parquet part files (the long tables, incl. unweighted means) stay on CSD3.
```bash
rsync -av --include='*/' --include='*.tsv' --include='*.txt' --exclude='*' \
  rajd2@login.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/abcd-data-release-7.0/processed/ \
  ~/Git/abcd_development/abcd-data-release-7.0/processed/
```

## Measures 4–5: T2w and T1w/T2w

### Method
- **Registration.** fMRIPrep's `desc-preproc_T2w` is already in its `desc-preproc_T1w` grid, so a single
  same-session T1w→T1w registration to the DAIRC `orig.mgz` places both images. It runs as `mri_coreg` (rigid) to
  initialise, then `bbregister --t1` (boundary-based, 6 dof). Output is `bbr.lta` per session under
  `legacy/hpc/work/d1_t2reg/<sub>/<ses>/`; `abcd.surf_metrics` reads it as a vox2vox LTA.
- **Faster `mri_coreg`.** On 19 pilot sessions, `--sep 4 --ftol 1e-4 --linmintol 1e-2` cut the time to 124 s per
  session (90 s coreg + 35 s bbregister), against ~270 s for the default. Its white-surface displacement from the
  default was 0.001 mm mean and 0.006 mm max. `--init-header` with no `mri_coreg` failed outright on 2 of 19 sessions
  (cost ≈ 0.9).
- **Sampling check.** The LTA sampling in Python vs `mri_vol2surf --reg bbr.lta` on the midpoint surface gives
  r = 1.00000, median |Δ| = 0.
- **The fMRIPrep T2w is not bias-corrected** (fMRIPrep applies N4 to the T1w only). A quadratic fit to log WM
  intensity on 9 pilot sessions showed a 25–111 % field spread in the T2w, against a few % in the fMRIPrep T1w.
  `t2reg.sbatch` therefore runs ANTs N4 on the T2w (fMRIPrep brain mask, `-s 4 -b [200]`, 4×50 iterations). It keeps
  only the field, at 8 mm, normalised to median 1 in the mask; `surf_metrics` divides the T2w samples by it. Residual
  spread after N4 was 6–13 %.
  - **N4 trap:** the T2w contains small negative values (resampling undershoot). N4 log-transforms its input, and
    given negatives it silently returns a perfectly flat field. A synthetic ±30 % field caught this. The fix is
    clipping to ≥ 1 (and `-r 1`), and the script now fails a session whose field is flat.
- **Ratio.** T1w/T2w = fMRIPrep `desc-preproc_T1w` / N4-corrected `desc-preproc_T2w`, computed per vertex at gm50
  through the same `bbr.lta`, then averaged like every other measure.
- **Field-of-view floor.** fMRIPrep resamples into its T1w grid with zero fill, so cortex outside the T2w field of view
  samples ≈ 0. In the first pass this let whole-cortex ratios reach 2,199 (median 3.9). Samples of an fMRIPrep image
  below 0.2 × that hemisphere's median gm50 sample are now set to NaN; this is far below any tissue. Median 0.02 % of
  cortex is excluded per session (99th percentile 0.4 %). Only 4 sessions lose more than 5 % of T2 vertices (max
  12.5 %). One session, `sub-EF2EVE42/ses-00A`, loses 48 % of its ratio vertices to the T1w floor; the counts are in
  the QC file (`n_nonfinite_*`).

### Registration outcome (28,765 sessions with fMRIPrep T1w + T2w)
- 28,763 registered. Two failed:
  - `sub-2ZLBBMTU/ses-04A`: the known incomplete recon (no `aparc+aseg.mgz`).
  - `sub-64PE3J4T/ses-06A`: failed on a FreeSurfer quirk. FreeSurfer picks the file format by finding `.mgz`
    anywhere in the path, and a random temp-dir suffix `.mgzt` made it misread `lh.thickness`. With dot-free temp
    names it reran fine.
- bbregister cost: median 0.024, 95 % 0.042, 99 % 0.064.
- **75 sessions are gross failures** (cost 0.84–0.98), with a clean gap between 0.3 and 0.8 below them. Tested on 8
  of them, `--init-header` and the default (two-scale) `mri_coreg` rescued only 1. `t2reg.sbatch` now retries every
  cost > 0.5 once with the default `mri_coreg` and keeps the better result; this was applied to all 77 high-cost
  sessions. `surf_metrics` treats cost > 0.5 as unregistered (`MAX_BBR_COST`); the QC status says "registration
  failed (bbr cost …)". The cause was not pinned down; the inputs (one T1w run each, same session) look ordinary.
  41 sessions with cost 0.1–0.3 are kept, and every T2 QC row carries `bbr_mincost` for filtering.
- Compute: ~1,590 core-hours for registration + N4 (median per session: coreg 93 s, bbregister 43 s, N4 32 s).
  Sampling took ~40 core-hours.

### Sessions (ok)
| measure | ses-00A | ses-02A | ses-04A | ses-06A | total |
|---|--:|--:|--:|--:|--:|
| t2__gm50, t2__wm1, t1t2__gm50 | 11,043 | 7,724 | 6,175 | 3,747 | 28,689 |

Not ok: 5,029 sessions have no fMRIPrep T2w (input gap, mostly ses-06A), 75 failed registration, and 2 failed as
above. That is 5,106 in all.

### Trends
Whole-cortex mean, median [IQR] by visit (HCP; DK near-identical):

| measure | ses-00A | ses-02A | ses-04A | ses-06A |
|---|---|---|---|---|
| t1t2__gm50 | 3.59 [3.10–4.07] | 3.81 [3.12–4.40] | 4.35 [3.20–5.09] | 4.89 [3.67–5.68] |
| t2__gm50 | 483 [402–872] | 507 [424–689] | 511 [425–681] | 488 [387–695] |
| t2__wm1 | 380 [318–695] | 398 [330–549] | 398 [328–545] | 376 [306–557] |

Within-subject change, consecutive visits (median per-child % change, fraction increasing):

| | 00A→02A (n 7,310) | 02A→04A (n 4,568) | 04A→06A (n 2,402) |
|---|---|---|---|
| t1t2__gm50 | +7.5 %, 0.66 | +9.5 %, 0.71 | +11.5 %, 0.77 |
| t2__gm50 | +3.8 %, 0.61 | +1.3 %, 0.54 | −2.9 %, 0.41 |

- **T1w/T2w rises, as expected with intracortical myelination**, consistently within subject. Cross-sectional
  per-visit medians also mix scanner composition: manufacturer 1 is flat until ses-06A.
- The raw **T2 level shows no consistent age trend.** It is in scanner units with a per-session receive gain:
  N4's field is normalised to median 1, so it removes spatial bias only, not the global scale. Do not read the T2
  level as a biological quantity.

### Scanner sensitivity
| measure | R² site | R² manufacturer | median by manufacturer 1 / 2 / 3 |
|---|--:|--:|---|
| t2__gm50 (hcp) | 0.862 | 0.840 | 298 / 179,200 / 506 |
| t2__wm1 (hcp) | 0.861 | 0.839 | 250 / 144,300 / 395 |
| t1t2__gm50 (hcp) | 0.477 | 0.413 | 3.40 / 0.43 / 4.23 |

- Manufacturer 2 writes its T2w on a scale ~400× larger, so its ratio is ~0.43 rather than ~4.
- The ratio cancels much of the scanner scale (R² 0.48 vs 0.86 for T2 alone, 0.83 for T1 gm50). It still needs
  manufacturer/site in the model, or a within-manufacturer or log-scale analysis. On a log scale the T2 manufacturer
  R² is 0.97.
- For the laptop LMMs, use `t1t2__gm50` (log, with manufacturer/site) and within-child change; treat the `t2__*`
  levels as QC/secondary.

## Decisions needed from the user
1. **Measures 6–7:** the same kind of registration for qsiprep T1w would make tensor FA/MD (and DKI/NODDI) sampleable
   without preprocessing raw diffusion. RSI is not available on CSD3.
