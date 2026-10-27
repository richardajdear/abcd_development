# Task for the CSD3 agent: HCP-MMP (and DK) per-session cortical metrics beyond thickness (D1)

## Goal
D1 in `docs/DIRECTIONS.md` asks whether the risk-related part of adolescent thinning is intracortical
myelination (the grey/white boundary moving outward) or loss of neuropil. To answer it we need, for every ABCD
session, **per-parcel values of myelin- and microstructure-sensitive measures in the HCP-MMP1.0 parcellation**.
These must come from the same FreeSurfer surfaces and HCP annotation that gave the thickness table
`abcd-data-release-7.0/processed/hcp/mr_y_smri__thk__hcp.tsv`, which every current analysis uses.

The laptop pipeline then fits one single-LMM slope per whole-cortex measure (as Figure 1) and per parcel. It also
fits a cross-modal composite slope. **Nothing statistical is run on CSD3.** The deliverable is a set of
release-convention TSVs plus a report.

Measures, in priority order:

| # | measure | file stem (HCP; DK twin `__dsk`) | needs | expected age trend |
|---|---|---|---|---|
| 1 | grey/white contrast, FreeSurfer `pctsurfcon` | `mr_y_smri__gwc__hcp` | `surf/?h.w-g.pct.mgh` | falls |
| 2 | T1w intensity, GM mid-depth (projfrac 0.5) | `mr_y_smri__t1__gm50__hcp` | `mri/nu.mgz` + surfaces | rises |
| 3 | T1w intensity, superficial WM (1 mm below white) | `mr_y_smri__t1__wm1__hcp` | same | rises |
| 4 | T2w intensity, GM mid-depth and WM 1 mm | `mr_y_smri__t2__{gm50,wm1}__hcp` | T2w registered to the FreeSurfer subject | falls |
| 5 | T1w/T2w ratio, GM mid-depth | `mr_y_smri__t1t2__gm50__hcp` | 2 + 4 | rises |
| 6 | DTI FA and MD, GM mid-depth and superficial WM (1 mm and 2 mm) | `mr_y_dmri__{fa,md}__{gm50,wm1,wm2}__hcp` | processed DTI maps + registration to the FreeSurfer subject | FA rises in WM, MD falls |
| 7 | RSI restricted normalised isotropic / directional (RNI, RND), same depths | `mr_y_dmri__{rni,rnd}__{...}__hcp` | processed RSI maps | RND rises in WM |

Area and grey-matter volume already exist on CSD3 (`processed/hcp/mr_y_smri__{area,vol}__hcp.tsv` from
`abcd.hcp_stats`, and `processed/dsk_local/` from `abcd.dk_stats`). Do **not** recompute them; just confirm they
are there (see "Return").

## Environment (existing; do not reinstall)
- **Repo:** `REPO=/home/rajd2/rds/hpc-work/abcd_development`. Pull `main` first; this brief lives at
  `directions/d1_multimodal/CSD3_INSTRUCTIONS.md`.
- **Python env:** `legacy/hpc/work/envs/abcd`, as in `tools/hcp_extract.sbatch`. Activate it the same way, with
  `PYTHONPATH=src`. It needs `nibabel` and `scipy`. If `nibabel` is missing, `pip install nibabel` into that env
  and record the version.
- **FreeSurfer tree (DAIRC FS 7.1.1):** `/rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer/<sub>/<ses>/`.
- **HCP-MMP annotation tree:** `/rds/project/rds-CeXlNYOYMxw/derivatives/parcellations/T1/<sub>/<ses>/HCP.fsaverage.aparc/`,
  with backfill at `/rds/user/rajd2/hpc-work/abcd_development/legacy/hpc/work/parcellations_backfill/T1`. The backfill
  wins where both exist, exactly as `src/abcd/hcp_stats.py` resolves them. Reuse its session discovery; do not write a
  new one.
- **Slurm:** account `vertes-sl3-cpu`, partition `icelake`. Memory comes with cores, so give no `--mem` line (the old
  `tools/hcp_extract.sbatch` has one; do not copy that line). Logs go to `slurm/`.
- **Code goes in:** `src/abcd/surf_metrics.py` (module, CLI like `abcd.hcp_stats`), `tools/surf_metrics.sbatch`
  (array job over session chunks) and `tools/surf_metrics_collect.py`. Commit all three. Do not commit any data.

## Step 0 — discovery (report before any large job)
Write `directions/d1_multimodal/CSD3_DISCOVERY.md` and push it, then continue with step 1 without waiting.
Check on 20 sessions spread over sites and visits (`ses-00A` … `ses-06A`), then count over all sessions:
1. In the FreeSurfer tree, the presence of `surf/{lh,rh}.{white,pial,area,w-g.pct.mgh}`,
   `stats/{lh,rh}.w-g.pct.stats`, `mri/{orig,nu,T1,rawavg}.mgz`, `label/{lh,rh}.aparc.annot`, and **any T2**:
   `mri/T2.mgz`, `mri/T2.prenorm.mgz`, `mri/transforms/T2raw.*`, or a T2 file anywhere under the session. Note whether
   `recon-all` was run with `-T2pial` (check `scripts/recon-all.log` / `recon-all.done` command line).
2. In the HCP annotation folder, whether a **subject-space** `?h.HCP.fsaverage.aparc.annot` (or similar) is stored next
   to the `.log` tables. If it is not, find how `Code/parcellation_T1/parcellate.sh` creates it (`mri_surf2surf`
   from the fsaverage annotation). Step 1 then has to regenerate it per session, and FreeSurfer must be available
   (`module avail freesurfer`; record the version and where `fsaverage` lives).
3. `ls` the top two levels of `/rds/project/rds-CeXlNYOYMxw/` and of its `derivatives/`. List anything that looks
   like raw T2w (`*_T2w.nii*`), processed diffusion (`dmri`, `dwi`, `qsiprep`, `mproc`, `DTI`, `RSI`, `FA`, `MD`) or
   registration to the T1. For each candidate, give the path pattern, the files per session, a session count, the
   space (native dMRI, T1/FreeSurfer, MNI) and whether a transform to the FreeSurfer subject exists.
4. Whether the release tabulated tables on CSD3 (`$REPO/abcd-data-release-7.0/`) include
   `mr_y_smri__t1__gm__dsk`, any T1 white-matter or contrast table, and any `mr_y_dti__*` / `mr_y_rsi__*` table. These
   are the gate references below.

Decision rules:
- **Measures 1–3** need only the FreeSurfer tree (plus the HCP annotation). Always do them.
- **Measures 4–5:** do them if a T2 already registered to the subject exists (`mri/T2.mgz` or an equivalent with a
  transform). If only raw T2w exists, pilot `bbregister --t2` (or `mri_coreg` + `bbregister`) on 20 sessions.
  Report the per-session time and the projected total core-hours, then **stop measures 4–5 and wait for the user.**
- **Measures 6–7:** do them only if processed DTI/RSI maps **and** a registration to the FreeSurfer subject already
  exist. **Do not preprocess raw diffusion data.** If they are missing, report what exists and skip these measures;
  the laptop will use the DK release tables instead.

## Step 1 — sampling method (one module, all measures)
For each session and hemisphere, using nibabel and scipy (FreeSurfer binaries are only needed for the cross-check and
for annotation regeneration):
1. Read `surf/?h.white` and `surf/?h.pial` (vertex-matched meshes), `surf/?h.area` (vertex area, the weight), the
   subject-space HCP annotation and `label/?h.aparc.annot` (DK).
2. **Surface measures (GWC):** read `surf/?h.w-g.pct.mgh` directly; it already has one value per vertex.
3. **Volume measures:** sample with trilinear interpolation (`scipy.ndimage.map_coordinates`, order 1) at:
   - `gm50`: white + 0.5 × (pial − white);
   - `wm1`, `wm2`: white − d × n̂, with d = 1 and 2 mm and n̂ the outward unit normal of the white surface
     (vertex normals from the mesh; `nibabel.freesurfer` has no normal helper, so compute area-weighted face
     normals).
   
   Coordinates are FreeSurfer tkr-RAS. Convert them with the inverse of the volume's `vox2ras-tkr` (nibabel
   `MGHHeader.get_vox2ras_tkr()`). For volumes not in conformed space (registered T2, diffusion), compose with
   their registration.
4. **Parcel summary:** the vertex-area-weighted mean over the vertices of each parcel. Exclude the medial wall (`???`
   in HCP; `unknown`/`corpuscallosum` in DK), vertices with label −1, and non-finite samples. Also write hemisphere and
   whole-cortex means over all non-medial-wall vertices, weighted by vertex area (the release convention).
5. **T1w/T2w ratio:** compute it per vertex before averaging. Use `nu.mgz` for T1 in measures 2–3, and record the
   choice. For the ratio, take T1 from the same unnormalised source that the T2 registration uses (e.g. `rawavg`
   resampled to conformed space). The ratio exists to cancel the shared receive field, so do not divide a
   bias-corrected T1 by an uncorrected T2.

Output (release convention, so `abcd.io.Release70Adapter` can read it like the thickness table):
- `$REPO/abcd-data-release-7.0/processed/hcp/<stem>.tsv`, with columns `participant_id`, `session_id` and
  `<stem>__<region>__<hemi>_mean` for the 180 parcels × 2 hemispheres. Use the same region tokens as
  `mr_y_smri__thk__hcp.tsv` (copy the naming function from `abcd.hcp_stats`). Add `<stem>__lh_mean`,
  `<stem>__rh_mean` and `<stem>_mean`.
- The DK twin of each stem in `processed/dsk_local/`, using the release region codes (`abcd.dk_stats` naming).
- `processed/hcp/surf_metrics_session_qc.tsv`: per session and measure, the number of sampled vertices, the number of
  non-finite samples, and flags for a missing input or a regenerated annotation.
- `PROVENANCE_surf_metrics.txt`: inputs, sampling depths, intensity source per measure, software versions, date.

## Step 2 — pilot and gates (20 sessions, interactive or one short job)
Report each gate in `CSD3_DISCOVERY.md` (append a "Pilot" section). **Stop if gate A or B fails.**
- **A. GWC reproduces FreeSurfer:** own DK GWC parcel means vs `stats/?h.w-g.pct.stats` `Mean` column. Expect r > 0.999
  and max |Δ| < 0.5 (percent units). FreeSurfer's stats are unweighted vertex means, so compare the unweighted version.
- **B. Annotation agrees with thickness:** HCP per-parcel vertex counts equal the `nverts` in
  `processed/hcp/hcp_session_qc.tsv` / the `.log` tables for the same sessions (exact).
- **C. Sampler cross-check (if FreeSurfer is available):** `mri_vol2surf --mov nu.mgz --regheader <subj> --hemi lh
  --projfrac 0.5` vs own `gm50` on 3 sessions. Expect vertex-level r > 0.99.
- **D. Release comparison (if the tables exist on CSD3):** own DK `t1__gm50` vs release `mr_y_smri__t1__gm__dsk`, and
  own DK diffusion vs release `mr_y_dti__*` / `mr_y_rsi__*`, per region across the 20 sessions. Report r; the depths
  differ, so expect r > 0.8, not 1.
- **E. Timing:** seconds per session per measure group; projected core-hours for all sessions.

## Step 3 — full run
- Array job over all sessions found by the `hcp_stats` discovery (33,795 last time), in chunks of about 200 sessions
  per task, with 1 core per process and `--cpus-per-task` matched to a worker pool, as in `hcp_extract.sbatch`. Make
  it restartable: skip chunks whose part-file exists. Then run `surf_metrics_collect.py` to assemble the TSVs.
- Measures 1–3 (and 6–7 if allowed) go in the same pass, so each session's surfaces are read once.

## Checks to report (`directions/d1_multimodal/CSD3_REPORT.md`, pushed)
- Sessions processed per measure; per-visit counts; the sessions failing per measure, by reason.
- Gates A–E with numbers.
- Per measure: whole-cortex mean by visit (`ses-00A` … `ses-06A`; median and IQR). This is a sanity check that the
  trends go the way the table above expects (GWC falls; T1w and T1w/T2w rise).
- Per measure: the share of variance in the whole-cortex mean explained by site, and by scanner manufacturer from
  `ab_g_dyn`. This is the scanner-sensitivity figure we already know is a problem for the release T1w/T2w.
- Exact paths and sizes of every output file, and the rsync command (below) with its file list.

## Return
Commit code and reports only (never data). The user will rsync:
```bash
rsync -av --include='*/' --include='*.tsv' --include='*.txt' --exclude='*' \
  rajd2@login.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/abcd-data-release-7.0/processed/ \
  ~/Git/abcd_development/abcd-data-release-7.0/processed/
```
If the release directory on CSD3 is not at that path, put the correct path in the report.

## Not in scope
Preprocessing raw diffusion; VBM; any slope, GWAS, PRS or twin model; changing the existing thickness, area or volume
tables.
