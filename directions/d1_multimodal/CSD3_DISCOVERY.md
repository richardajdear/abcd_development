# D1 on CSD3: discovery (step 0)

Run on 2026-10-27 by `directions/d1_multimodal/csd3_discovery.py` (Slurm job 39326028). The inventory covers all
**33,795 sessions** of the `abcd.hcp_stats` discovery (24,921 primary + 8,874 backfill), the same set the thickness
table was built from. Per-session flags are in `abcd-data-release-7.0/processed/hcp/surf_metrics_inventory.tsv`.

## 1. FreeSurfer tree (`derivatives/freesurfer/<sub>/<ses>/`, DAIRC FS 7.1.1)

| file | sessions | missing |
|---|--:|---|
| `surf/{lh,rh}.{white,pial,area}` | 33,795 | — |
| `surf/{lh,rh}.w-g.pct.mgh`, `stats/{lh,rh}.w-g.pct.stats` | 33,794 | `sub-2ZLBBMTU/ses-04A` (it also has no `recon-all.done` CMDARGS line) |
| `mri/{orig,nu,T1,rawavg}.mgz` | 33,795 | — |
| `label/{lh,rh}.aparc.annot` | 33,795 | — |
| **any T2 file** (`mri/T2*.mgz`, `transforms/T2raw.*`, any name matching `t2`) | **0** | — |
| `recon-all` run with `-T2pial` | **0** of 33,794 with a CMDARGS line | — |

By visit (white/pial present): ses-00A 11,749; ses-02A 8,093; ses-04A 6,348; ses-06A 7,605.

## 2. Subject-space HCP-MMP annotation: already stored, no regeneration needed
- rr480's `parcellate.sh` wrote `label/?h.HCP.fsaverage.aparc.annot` **into the FreeSurfer session** (owner rr480,
  dated July 2026). The `.log` tables were computed from it. All **24,921 primary sessions** have both hemispheres.
- The backfill run (`tools/hcp_backfill_session.sh`) kept its `mri_surf2surf` output next to the `.log`:
  `parcellations_backfill/T1/<sub>/<ses>/HCP.fsaverage.aparc/?h.HCP.fsaverage.aparc.annot`. All **8,874
  backfill sessions** have both hemispheres.
- `abcd.surf_metrics` resolves the annotation with the same rule as `hcp_stats`: backfill wins, otherwise the
  FreeSurfer `label/` copy. Some primary-tree sessions that were backfilled also have an FS `label/` copy, which may
  be broken (the symlink race); it is never read.
- rr480 also left FreeSurfer's own per-parcel GWC for HCP (`<parc>/HCP.fsaverage.aparc/?h_HCP.fsaverage.aparc.w-g.pct.stats`,
  `mri_segstats` on `w-g.pct.mgh`), which gives a second GWC gate (A′).
- FreeSurfer modules: 6.0.0, 7.1.0, 7.3.2, 7.4.0, 7.4.1, 8.0.0, 8.2.0. `fsaverageSubP` is at
  `/rds/project/rds-CeXlNYOYMxw/userdata/rr480/fsaverageSubP` (the annotation source). FS 7.1.0 is used for the
  cross-checks.

## 3. Other imaging under `/rds/project/rds-CeXlNYOYMxw/`
Top level: `Code Data_Genetics Data_Imaging Data_Out Data_Phenotype Data_Raw derivatives logs recycle.bin Scratch
Temp Tools userdata workdir`. `Data_Raw/` is empty and `Data_Imaging/` holds one test subject. `find -iname '*T2w*'` over `Data_Imaging`,
`Data_Raw`, `Scratch`, `Temp` and `derivatives/scratch` finds nothing, so **no raw T2w** is on the project (only
fMRIPrep's processed one, below).
`derivatives/`: `abcc-fmriprep_v25.1.4 abcc-qsiprep_v0.21.4 abcc-qsirecon-{DIPYDKI,DSIStudioGQI,MSMTAutoTrack,wmNODDI}
abcc-xcp_d_v0.13.0 freesurfer parcellations scratch tabulated`.

| candidate | path pattern | per session | sessions | space | transform to the DAIRC FS subject? |
|---|---|---|--:|---|---|
| processed T2w | `abcc-fmriprep_v25.1.4/<sub>/<ses>/anat/*_desc-preproc_T2w.nii.gz` | T2w resampled to the fMRIPrep T1w (+ `desc-preproc_T1w`, brain mask, dseg, GIFTI surfaces) | 28,765 (00A 11,078; 02A 7,752; 04A 6,184; 06A 3,751) | fMRIPrep T1w | **no**. `from-T1w_to-fsnative_xfm.txt` is a near-identity to fMRIPrep's **own** FreeSurfer run, not to DAIRC's recon. |
| processed DWI | `abcc-qsiprep_v0.21.4/<sub>/<ses>/dwi/*_space-T1w_desc-preproc_dwi.nii.gz` | preprocessed DWI + qsiprep T1w, `from-T1wACPC_to-T1wNative` | 25,828 | qsiprep T1w (ACPC) | **no** |
| DTI/DKI maps | `abcc-qsirecon-DIPYDKI/.../dwi/*_space-T1w_model-{tensor,dki}_param-{fa,md,ad,rd,mk,...}_dwimap.nii.gz` | tensor FA; DKI MD/AD/RD/MK/...; MNI copies | 24,792 | qsiprep T1w | **no** |
| DTI/GQI maps | `abcc-qsirecon-DSIStudioGQI/.../dwi/*_space-T1w_model-tensor_param-{fa,md,ad,rd}_dwimap.nii.gz` | tensor FA/MD/AD/RD, GQI QA/GFA/ISO, RDI | 24,792 | qsiprep T1w | **no** |
| NODDI | `abcc-qsirecon-wmNODDI/.../dwi/*_space-T1w_model-noddi_param-{icvf,od,isovf}_dwimap.nii.gz` | NODDI | ~24.8k | qsiprep T1w | **no** |
| RSI (RNI/RND) | — | — | 0 | — | — |

## 4. Release tables on CSD3 (`$REPO/abcd-data-release-7.0/`)
`y/` holds only `mr_y_smri__thk__dsk`, **`mr_y_smri__t1__gm__dsk`**, **`mr_y_smri__t2__gm__dsk`**, `mr_y_qc__incl` and
`mr_y_qc__post__aut` (33,794 rows each). There is no T1 white-matter or contrast table, and no `mr_y_dti__*` or
`mr_y_rsi__*` table. Gate D can therefore check T1 only (and T2 later).

## Decisions
- **Measures 1–3 (GWC, T1 gm50, T1 wm1): go.** They need only the FreeSurfer tree and the annotations above.
- **Measures 4–5 (T2, T1w/T2w): stopped, waiting for the user.** No T2 is registered to the DAIRC subject. A
  *processed* T2w exists for 28,765 sessions, in fMRIPrep's T1w space. So the missing piece is one rigid
  same-session **T1w→T1w** registration (fMRIPrep T1w → DAIRC `orig.mgz`), not a T2→T1 bbregister from raw data.
  The pilot (below) times `mri_coreg` + `bbregister --t1` on the 20 pilot sessions. Two open questions for the user:
  the cost, and whether fMRIPrep's `desc-preproc_T2w` is bias-corrected. For the ratio, T1 and T2 should both come
  from fMRIPrep (same space, same treatment), and not from `nu.mgz`.
- **Measures 6–7 (DTI/RSI): skipped per the brief.** Processed tensor maps exist for 24,792 sessions, but there is
  no registration to the DAIRC subject, and RSI does not exist at all. Note: the same kind of T1w→T1w rigid
  registration (qsiprep T1w → DAIRC `orig.mgz`) would make FA/MD (and DKI/NODDI) sampleable without touching raw
  diffusion. If the user wants it, it is a second ~25k-session registration run.

## Pilot (step 2): Slurm job 39327812, `directions/d1_multimodal/csd3_pilot.sbatch`
20 sessions (`processed/hcp/surf_metrics_pilot_sessions.tsv`): 5 per visit, distinct sites within a visit, all
three manufacturers, 15 primary + 5 backfill. All 20 × 3 measures `ok`. Full output:
`processed/hcp/surf_metrics_pilot/gates.txt`.

| gate | test | result |
|---|---|---|
| **A** GWC reproduces FreeSurfer | own DK unweighted parcel mean vs `stats/?h.w-g.pct.stats` `Mean` (1,360 parcels) | r = 1.000000, max \|Δ\| = 0.0001 → **PASS** |
| A′ (extra) | own HCP GWC vs rr480's `?h_HCP.fsaverage.aparc.w-g.pct.stats` (30 hemispheres with that file, 5,400 parcels) | r = 1.000000, max \|Δ\| = 0.0001 |
| **B** annotation = thickness annotation | HCP per-parcel vertex counts vs `.log` `nverts` (7,200 parcel rows) | 0 unmatched, 0 mismatches → **PASS** |
| **C** sampler vs FreeSurfer (lh, 3 sessions) | `mri_vol2surf --mov nu.mgz --regheader --interp trilinear` | `--surf mid --projdist 0` (same point as gm50, the white–pial midpoint): r = 1.00000, median \|Δ\| 0.000 → **PASS**; `--projdist -1` vs wm1: r = 0.993–0.995, median \|Δ\| 0.02–0.04 → **PASS**; `--projfrac 0.5`: r = 0.936–0.962 (expected to differ: it steps 0.5 × thickness along the white normal, not to the midpoint) |
| **D** release comparison | own DK `t1__gm50` vs release `mr_y_smri__t1__gm__dsk` | per-region r across sessions: median 0.980 (IQR 0.968–0.987, min 0.896), all 68 > 0.8 → **PASS**; whole-cortex r = 0.994. The level differs (own median 67.5 vs release 81.4), as expected if the release uses another depth and/or intensity scaling (not verified here). No DTI/RSI tables on CSD3. |
| **E** timing | GWC + T1 gm50 + T1 wm1, HCP + DK, both hemispheres | 1.6 s/session (max 2.5) → ~15 core-hours for 33,795 sessions |

Other checks: the white-surface normal agrees with the white→pial direction at ≥ 99.86 % of vertices in every pilot
hemisphere (median 99.96 %). No non-finite samples.

### T2 registration pilot (for the measures 4–5 decision)
`mri_coreg --s <subj> --mov fMRIPrep desc-preproc_T1w` followed by `bbregister --t1 --init-reg` (FS 7.1.0, 1 thread):
- 18/20 registered; 2 have no fMRIPrep anat (`sub-8R9KBJ6M/ses-04A`, `sub-3ADKR0X9/ses-06A`).
- Time per session: `mri_coreg` median 231 s, `bbregister` 38 s; total median 272 s (max 324 s).
- bbregister minimum cost 0.016–0.061 (median 0.023). All are low, so no failed alignments.
- Projected: ~4.5 core-min × 28,765 sessions with fMRIPrep T2w ≈ **2,150 core-hours** (sl3 icelake).
  `mri_coreg` dominates. A coarser `--sep`, or bbregister alone with `--init-header`, should cut this a lot (the two T1s
  come from the same acquisition), but that is untested.
- **Measures 4–5 are not run. Waiting for the user** on the cost and the T2 source (fMRIPrep `desc-preproc_T2w`; the
  ratio would use fMRIPrep T1w/T2w, both in fMRIPrep T1w space, sampled on the DAIRC surfaces through this transform).
