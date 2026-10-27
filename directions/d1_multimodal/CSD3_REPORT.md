# D1 on CSD3: report (measures 1–3 done; 4–7 held)

2026-10-27. Brief: `CSD3_INSTRUCTIONS.md`. Discovery and pilot: `CSD3_DISCOVERY.md`. Code:
`src/abcd/surf_metrics.py`, `tools/surf_metrics.sbatch`, `tools/surf_metrics_collect.py`. Helper scripts in this
directory: `csd3_discovery.py`, `csd3_pilot.sbatch`, `csd3_pilot_gates.py`, `csd3_report_stats.py`.
Slurm jobs: discovery 39326028, pilot 39327812, full array 39327815 (169 × 200 sessions, all COMPLETED), collect
39331883.

## What was produced
| # | measure | stem (HCP / DK twin) | status |
|---|---|---|---|
| 1 | GWC (`surf/?h.w-g.pct.mgh`) | `mr_y_smri__gwc__hcp` / `__dsk` | **done** |
| 2 | T1w (`nu.mgz`) at GM mid-depth | `mr_y_smri__t1__gm50__hcp` / `__dsk` | **done** |
| 3 | T1w (`nu.mgz`) 1 mm below white | `mr_y_smri__t1__wm1__hcp` / `__dsk` | **done** |
| 4–5 | T2w, T1w/T2w | — | **held for the user.** No T2 in the DAIRC FreeSurfer tree. fMRIPrep `desc-preproc_T2w` exists (28,765 sessions) but needs a fMRIPrep-T1w→DAIRC registration: ~4.5 core-min/session, ~2,150 core-hours, pilot OK (see discovery) |
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
| `processed/hcp/surf_metrics_session_qc.tsv` | 101,385 (session × measure) | 9.3 MB |
| `processed/hcp/PROVENANCE_surf_metrics.txt` | | <0.1 MB |
| `processed/hcp/surf_metrics_inventory.tsv` (discovery) | 33,795 | 6.1 MB |
| `processed/hcp/surf_metrics_pilot_sessions.tsv`, `surf_metrics_pilot/{gates.txt,report_stats.md,t2reg_pilot.tsv}` | | small |

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

## Decisions needed from the user
1. **Measures 4–5 (T2w, T1w/T2w):** run the fMRIPrep-T1w → DAIRC registration (~2,150 core-hours as piloted,
   probably much less with a faster `mri_coreg` setting) and sample fMRIPrep `desc-preproc_T2w`, with the ratio
   taken from fMRIPrep T1w/T2w? Given the scanner dominance above, the ratio may be the most useful of these
   measures.
2. **Measures 6–7:** the same kind of registration for qsiprep T1w would make tensor FA/MD (and DKI/NODDI) sampleable
   without preprocessing raw diffusion. RSI is not available on CSD3.
