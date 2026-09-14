# Cortical thickness in the HCP-MMP1.0 (Glasser) parcellation — investigation and plan

**Written 2026-09-12.** Investigation of what exists on `rds-abcd-CeXlNYOYMxw`
(CSD3, `/rds/project/rds-CeXlNYOYMxw`) and what it would take to produce
QC'd, release-7.0-vintage HCP-MMP thickness for the longitudinal model in
this repo. Everything below was checked against the filesystem on the date
above; counts are exact unless marked as a sample.

## 0. Summary

**Difficulty: low-to-moderate, and most of the compute is already done.**
Nothing needs to be downloaded and no `recon-all` needs to be run.

1. The official ABCD FreeSurfer 7.1.1 reconstructions (run by the ABCD Data
   Analysis, Informatics and Resource Center — DAIRC — at UCSD) are unpacked
   on rds for **11,823 subjects / 33,825 sessions** across ses-00A/02A/04A/06A.
   They are the very reconstructions the release tables were tabulated from:
   Desikan thickness in a local `aparc.stats` file is bit-identical to the
   release table value (checked on `sub-003RTV85/ses-00A`).
2. Rafael Romero-Garcia (`rr480`) has **already projected HCP-MMP1.0 to native
   space and run `mris_anatomical_stats`** (28–31 July 2026). Output
   directories exist for 30,360 sessions, but only **24,921 carry a complete
   180-parcel table per hemisphere**; 5,439 hold 90-byte stubs because the
   `mri_surf2surf` step failed (see §1.4 — a symlink race between concurrent
   sessions of the same subject). *(Corrected after the full parse: the
   earlier spot check tested only that the files were non-empty.)*
3. **3,465 sessions have no parcellation output at all, and 3,435 of them are
   parcellable** (complete surfaces, no HCP annotation — the array jobs never
   reached them; two ran into the time limit). Only 30 sessions are unusable
   reconstructions. Together with the 5,439 stubs that is **8,874 sessions to
   (re)run**, listed in `abcd-data-release-7.0/processed/hcp/sessions_to_reparcellate.txt`;
   roughly 750 CPU-hours, one array job (§3.2).
4. The `derivatives/tabulated/` copy on rds is **release 6.0, not 7.0**: it is
   byte-identical to `Data_Phenotype/ABCD60/`, its scan dates stop in
   February 2024, and it has 4,103 six-year sessions. The FreeSurfer derivatives
   carry **3,516 six-year sessions that are absent from those tables**, with
   scan dates running to July 2025. So the surfaces are 7.0-vintage but the
   QC/age/covariate tables on rds are not. **§1.3 raises the possibility that
   this repo's "7.0" imaging tables are also 6.0-vintage** — that needs to be
   settled before anything else.
5. Remaining work is therefore mostly *analysis engineering*: a census, a
   small gap-filling job, a parser, a QC policy that mirrors the DK one, and a
   `Release70Adapter._imaging_hcp` implementation. Roughly 2–4 working days,
   plus whatever it takes to obtain the 7.0 tabulated release.

## 1. What is on rds

### 1.1 Raw imaging: absent, and not needed

`Data_Imaging/` holds one subject. The per-session BIDS raw data were never
kept locally at scale (the 2023–24 pipeline of S. Orellana under `workdir/`
downloaded release-5 archives and processed them in place; that tree is gone).
Because the DAIRC reconstructions are on disk, raw T1s are not required for
cortical thickness.

### 1.2 FreeSurfer derivatives (the key asset)

```
derivatives/freesurfer/sub-<8char>/ses-{00A,02A,04A,06A}/{mri,surf,label,stats,scripts}
```

| | count |
|:--|--:|
| subject directories | 11,823 |
| sessions | 33,825 |
| ses-00A / 02A / 04A / 06A | 11,757 / 8,097 / 6,359 / 7,612 |
| sessions present in the 6.0 QC table | 30,301 |
| sessions **not** in the 6.0 QC table | 3,524 (3,516 of them ses-06A) |

Provenance: downloaded as 33,825 per-session zips (`Code/ziplist.txt`,
`Code/submit_unzip.sh`, run 28 May 2026 by rb643), then unpacked in place.
`scripts/build-stamp.txt` reads `freesurfer-linux-centos7_x86_64-7.1.1`
in both a baseline and a six-year session, run as user `abcdproc1` on
`mmil-compute-*` hosts — the DAIRC processing cluster. This is cross-sectional
`recon-all` per session (ABCD does not run the FreeSurfer longitudinal stream),
which is also what the release DK tables are.

Each `stats/aseg.stats` carries `lhSurfaceHoles`, `rhSurfaceHoles`,
`SurfaceHoles` — the Euler number is `2 − 2·holes` per hemisphere — so an
Euler-style QC variable is available locally for **every** session, including
the 3,516 that have no release QC row yet.

`label/` in each session already contains `{lh,rh}.HCP.fsaverage.aparc.annot`
(owner rr480, written 28 July 2026) alongside Schaefer 200/400, economo,
500.aparc, sjh, PALS lobes.

### 1.3 Tabulated release data: 6.0 on rds, and a question about this repo

`derivatives/tabulated/{e,g,l,p,t,y}` is `cmp`-identical to
`Data_Phenotype/ABCD60/`. Its `mr_y_adm__info` scan dates end 2024-02; it has
4,103 ses-06A rows in `mr_y_qc__incl` (4,083 with the T1 include flag set).

This repo's settled run (`out/thickness_dsk_70_*/manifest.json`, labelled
release 7.0) ends with **3,539** six-year scans after QC and the ≥2-visit
filter. That is consistent with a 4,083-row six-year input and **not** with
the ~7,600 six-year sessions the FreeSurfer derivatives show 7.0 to contain.
Either the tables vendored as `abcd-data-release-7.0/` are actually the 6.0 tabulation, or
7.0's imaging tables lag its imaging derivatives. The former is much more
likely. **Action:** count ses-06A rows in the `mr_y_smri__thk__dsk` table the
repo actually reads. If it is ~4,100, the "7.0" phenotype is 6.0-vintage and
the whole longitudinal model, not only the HCP work, stands to gain ~3,500
six-year scans from a proper 7.0 tabulated download.

Release QC variables available in the 6.0 tables (same names expected in 7.0):

| table | column | notes |
|:--|:--|:--|
| `mr_y_qc__incl` | `mr_y_qc__incl__smri__t1_indicator` | ABCD's recommended T1 inclusion; 0 for 789 of 30,391 sessions |
| `mr_y_qc__post__aut` | `mr_y_qc__post__aut__smri__topodfct_count` | topological defects before fixing; median 18, p95 46, p99 98, max 517 |
| `mr_y_qc__post__man__fsurf` | `..._fsurf_score` and sub-scores | manual review, populated for only 2,107 sessions (reviewed subset) |
| `mr_y_adm__info` | `..._dev_manufact`, `..._dev_model`, `..._dev__sftw_ver` | 18,838 Siemens, 7,930 GE, 3,615 Philips sessions |

### 1.4 The existing HCP-MMP parcellation pass (rr480, July 2026)

```
derivatives/parcellations/T1/sub-*/ses-*/HCP.fsaverage.aparc/
    lh.HCP.fsaverage.aparc.log      # mris_anatomical_stats table, 180 ROIs + ???
    rh.HCP.fsaverage.aparc.log
    {lh,rh}_HCP.fsaverage.aparc.w-g.pct.stats
    HCP.fsaverage.aparc.nii.gz, HCP.fsaverage.aparc_seq.nii.gz   # volumes, not needed
```

| | count |
|:--|--:|
| sessions with a parcellation directory | 30,360 |
| ses-00A / 02A / 04A / 06A | 10,564 / 7,292 / 5,684 / 6,820 |
| FreeSurfer sessions without one | 3,465 (1,193 / 805 / 675 / 792) |
| of which both hemisphere tables complete (180 parcels each) | **24,921** (8,705 / 5,880 / 4,633 / 5,703) |
| of which 90-byte stubs, no parcel rows (`mri_surf2surf` failed) | 5,439 |

Method (`Code/parcellation_T1/parcellate.sh`, the group's 2017 pipeline):
`mri_surf2surf --sval-annot` from `fsaverageSubP` (a copy of fsaverage whose
`label/` holds `{lh,rh}.HCP.fsaverage.aparc.annot`, dated 2016-10-28 — the
fsaverage projection of HCP-MMP1.0 used throughout this group's work and the
one the AHBA C1–C3 maps are defined on) through the subject's `sphere.reg`,
then `mris_anatomical_stats -a <annot> -b <subject> <hemi>`. That is the
standard way to get a Glasser parcellation out of FreeSurfer surfaces, and it
is the same code path that produced the 5.1 HCP files this repo already reads
(`Release51Adapter._imaging_hcp`). The `mris_anatomical_stats` per-parcel
thickness is the unweighted vertex mean, the same definition as the release
`aparc.stats` values.

The table has one row per parcel with: vertex count, surface area, grey-matter
volume, mean thickness, thickness SD, mean curvature, Gaussian curvature,
folding index, intrinsic curvature index. So SA, GMV and curvature in HCP-MMP
come free with the thickness.

Each per-hemisphere table has exactly 180 `*_ROI` rows plus one `???` row
(medial wall), as expected for HCP-MMP1.0.

Why 3,465 sessions were not parcellated: the submitted list
(`Code/subject_status_check/to_parcellate_t1.txt`) contained all 33,823
sessions, so nothing was deliberately skipped. A status check of all 3,465
(recon-all.done, `surf/{lh,rh}.{sphere.reg,thickness}`, `label/lh.aparc.annot`,
`label/lh.HCP.fsaverage.aparc.annot`):

| state | sessions |
|:--|--:|
| recon finished, surfaces + aparc present, **no HCP annotation** (parcellable) | 3,432 |
| recon-all.done present but no surfaces (only `mri/` + `scripts/`) | 28 |
| HCP annotation present, stats never written | 3 |
| no recon-all.done | 2 |

By wave the 3,432 are 1,185 / 801 / 664 / 785 (00A / 02A / 04A / 06A) — they
are simply the tail of the submission that the array jobs did not reach (two
task logs end with `CANCELLED ... DUE TO TIME LIMIT`).

**The stub failures.** In 5,439 sessions the `.log` files contain only the
three header lines: `mri_surf2surf` never wrote the HCP annotation, so
`mris_anatomical_stats` had nothing to summarise (these are the 2,597
subject-level `error.log` files, "could not read annot file"). The failure
rate is 6 % for subjects with one session in the tree, 12 % with two, 17 %
with three and 22 % with four. `parcellate.sh` does
`rm $SUBJECTS_DIR/fsaverageSubP; ln -s ...` inside the *subject* directory
at the start of every session and `rm` at the end; with sessions of one
subject spread over concurrent array tasks, one task deletes the link another
is reading. The fix when re-running is to link `fsaverageSubP` once, outside
the loop (or use a per-session scratch `SUBJECTS_DIR` as in §3.2), and to
delete the stubs first because the script's `-s` test treats them as done.

### 1.5 Legacy files to ignore

- `Data_Out/HCP/Struct_CT.csv` (May 2020, 11,289 rows, wide `lh_L_V1_ROI`
  columns) — release 2.0 baseline only, FreeSurfer 6.
- `workdir/Code_Extraction/FS_parameters/ses-*_HCP.fsaverage.aparc_CT.csv` —
  release 5.x, the files behind this repo's 5.1 HCP config.
- `Data_Phenotype/ABCD20`, `ABCD30`, `release4`, `abcd-data-release-5.0`.

None of these should feed the 7.0 model. The 5.1 files remain useful as a
cross-check on overlapping sessions (§4.3).

## 2. Permissions and tooling

- I am in `rds-CeXlNYOYMxw-users`, not `-managers`. `derivatives/` is
  group-`managers` with an ACL that grants **read** to users; I cannot write
  into `derivatives/freesurfer/*/label` or `derivatives/parcellations/`.
  Gap-filling therefore either runs as rr480/rb643, or writes to a location
  we own (see §3.2, which needs no write access to the FreeSurfer tree).
- CSD3 modules: `freesurfer/7.1.0`, `7.3.2`, `7.4.0`, `7.4.1`, `8.x`. Use
  **7.1.0** to stay as close as possible to the 7.1.1 reconstructions (the
  existing pass did the same).
- SLURM accounts available to me: `vertes-sl2-cpu`, `vertes-sl3-cpu`.
- Default `python3` lacks `nibabel`; the parser in §3.3 needs only the
  standard library + pandas (the stats logs are plain text), so no surface
  library is required. `data/hcp_centroids.csv` (360 rows) already exists for
  spin tests.
- Login-node rule: anything touching 30k session directories must be an
  `sbatch` job, not a foreground loop (an unbounded `find` over
  `parcellations/T1` was still running after five minutes here).

## 3. Plan

### 3.0 Settle the release vintage (blocking for the model, not for the parcellation)

1. Count ses-06A rows in the thickness table this repo reads under its
   `abcd-data-release-7.0/` directory. Check `mr_y_adm__info` max scan date there.
2. If it is 6.0-vintage, obtain the 7.0 tabulated release (imaging QC,
   `ab_g_dyn` ages/site/scanner, `ab_g_stc`) and place it as `abcd-data-release-7.0/`.
   Without it the 3,516 new six-year sessions have surfaces but no age,
   scanner, or release QC flag, and cannot enter the model.
3. Record the answer in `README.md` §Status and in `io.py`'s
   `Release70Adapter` docstring.

### 3.1 Census (one sbatch job, ~10 min on 16 cores)

Script: `tools/hcp_census.sbatch` (not yet written) → `docs/hcp_census/session_status.csv`,
one row per FreeSurfer session:

- `recon_done` (`scripts/recon-all.done` present), `recon_ok` (`finished
  without error` in `recon-all.log`), `fs_build_stamp`
- `has_surfaces` (`surf/{lh,rh}.{white,pial,thickness,sphere.reg}` non-empty)
- `has_hcp_annot` (both hemispheres)
- `has_hcp_stats`, `hcp_stats_rows_lh/rh` (must be 180), file sizes
- `lh_holes`, `rh_holes`, `total_holes` from `aseg.stats`; derived Euler
- `in_release_qc_table`, `t1_include`, `topodfct_count`, `manufacturer`
  (joined from the tabulated release, 6.0 now, 7.0 once available)

Deliverables: the CSV, a per-wave summary table, and the list
`to_parcellate.txt` = sessions with `has_surfaces & !has_hcp_stats`. A
first version of that list (3,432 sessions) already exists from this
investigation and is reproduced by the census; the census additionally
verifies the 30,360 existing stats tables have 180 rows per hemisphere, which
the 1,209-session sample suggests but has not been checked in full (an
unbounded `find` over the tree does not finish on the login node — this must
be an sbatch job).

### 3.2 Fill the gap (one array job; no write access to the FS tree needed)

For each session in `to_parcellate.txt`, in a scratch `SUBJECTS_DIR` we own:

```
SUBJECTS_DIR=$SCRATCH/fs_link ; mkdir -p $SUBJECTS_DIR
ln -s /rds/project/rds-CeXlNYOYMxw/userdata/rr480/fsaverageSubP $SUBJECTS_DIR/fsaverageSubP
ln -s /rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer/$sub/$ses $SUBJECTS_DIR/${sub}_${ses}
for hemi in lh rh; do
  mri_surf2surf --srcsubject fsaverageSubP --sval-annot $SUBJECTS_DIR/fsaverageSubP/label/$hemi.HCP.fsaverage.aparc.annot \
      --trgsubject ${sub}_${ses} --trgsurfval $OUT/$sub/$ses/$hemi.HCP.fsaverage.aparc.annot --hemi $hemi
  mris_anatomical_stats -a $OUT/$sub/$ses/$hemi.HCP.fsaverage.aparc.annot -b ${sub}_${ses} $hemi \
      > $OUT/$sub/$ses/$hemi.HCP.fsaverage.aparc.log
done
```

Only the HCP annotation and the stats table are produced; the volume
parcellation, MATLAB renumbering and the other atlases in rr480's script are
not needed for thickness. Expect ≈5 min per session, so the 8,874 sessions
in `sessions_to_reparcellate.txt` are ≈750 CPU-hours — one `--array` of
one-core tasks on `vertes-sl3-cpu` (SL2 is out of CPU-minutes), done within a
day or two of queueing. Output goes to
`Scratch/` or `hpc-work/abcd_development/processed/hcp_gapfill/`, and the
parser in §3.3 reads both trees.

Before running: tell rr480 and rb643. The tidy thing is for rr480 to re-run
their own array on `to_parcellate.txt` so everything stays in one tree; the
scratch route is the fallback if that is slow to arrange.

### 3.3 Extraction to a tidy table — DONE 2026-09-12 (`src/abcd/hcp_stats.py`, `tools/hcp_extract.sbatch`)

`src/abcd/hcp_extract.py` : parse every
`{lh,rh}.HCP.fsaverage.aparc.log`, drop the `???` row, emit

```
processed/hcp70_anatomical_stats.parquet
    subject (sub-NDARINV…), session (ses-00A…), hemi, region (V1, MST, …),
    label (lh_V1), nverts, area_mm2, gmv_mm3, thickness_mm, thickness_sd,
    mean_curv, gauss_curv, folding_index, intrinsic_curv, source (rr480|gapfill)
```

plus `processed/hcp70_session_qc.parquet` from the census (§3.1). Assert 360
parcels per session and fail loudly on any session with a different count.

### 3.4 QC policy for HCP thickness

Mirror the DK `coded` policy so the two parcellations are comparable, then add
the parcel-level checks that DK does not need:

1. **Session level, from the release** (once 7.0 tables exist):
   `t1_include == 1`; scanner policy identical to the DK configs (Philips are
   currently retained in `ct_70_*` with a site/scanner effect — keep that
   decision shared, do not fork it for HCP).
2. **Session level, local**: `recon_ok`; `total_holes` (or `topodfct_count`)
   below a threshold. Calibrate rather than assume: the 5.1 legacy Euler list
   is on disk and `Release70Adapter` already validated `topodfct_count`
   against it (AUC 0.73). Propose the same cut for DK and HCP; report attrition
   in the ledger as the `qc.py` predicates already do.
3. **Parcel level (new)**: flag parcels with `nverts < 30` or thickness outside
   0.5–5.0 mm as missing; sessions with any missing parcel fail
   `complete_regions` (existing predicate, `n_expected = 360`). Report how
   many parcels this affects per region — the small HCP parcels (V6, MST, PEF,
   A1 ≈ 250–400 vertices in a 10-year-old) are the ones at risk.
4. **Distributional**: per region, per wave, flag `|z| > 4` (robust z via MAD)
   as a diagnostic table, not an automatic exclusion.

### 3.5 Validation before any modelling

1. **DK consistency**: for a random 500 sessions, parse local `aparc.stats`
   and compare to the release DK table — must be bit-identical (it was for
   the one session checked). This proves the local surfaces *are* the release.
2. **Whole-cortex mean**: HCP area-weighted mean vs DK `_mean` column,
   expect r > 0.99 per session.
3. **Against the 5.1 HCP files**: for sessions present in both, per-region ρ of
   thickness; 7.0 re-uses 5.1 reconstructions for shared waves (99.94 %
   identical DK values), so ρ should be ≈1. Any region far from 1 indicates an
   annotation or parsing mismatch.
4. **Known anatomy**: V1 ≈ 1.7–1.9 mm, 3b ≈ 1.8–2.0, TE1m/STS > 2.9,
   left–right correlation of the group mean map > 0.9.
5. Group mean and age-slope maps on the fsaverage HCP surface for eyeballing.

### 3.6 Wire into the repo

- `Release70Adapter._imaging_hcp(metric)` reading
  `processed/hcp70_anatomical_stats.parquet`; expose `thickness`, `area`,
  `volume` at least. Return the same long format as the DK path, `is_global`
  row computed as the surface-area-weighted mean, the release's own
  convention (established 2026-09-14: the release ``__lh_mean`` equals the
  area-weighted mean over the 34 DK regions, not FreeSurfer's vertex-weighted
  cortex mean; the 5.1 HCP path used the unweighted parcel mean).
- `qc.py`: predicates `surface_holes_predicate(max_holes)` and
  `parcel_coverage_predicate(min_verts)`; the release-include and
  complete-regions predicates already exist.
- `data/region_labels.csv` has no `hcp` rows; either add 360 rows or make
  `region_labels("hcp")` derive them from `hcp_centroids.csv`. Spin tests
  already have HCP centroids.
- New config `configs/ct_70_hcp_noglobal_mv2_genetic.yaml` = the settled DK
  spec with `parcellation: hcp`; register it in `README.md`'s config table.
- Tests: parser round-trip on two committed example logs; 360-parcel
  invariant; adapter smoke test on a 20-session fixture.
- `make all` should then run unchanged for the HCP config.

## 4. Effort and risk

| step | wall time | who/what it depends on |
|:--|:--|:--|
| 3.0 vintage check | 1 h to check; days to weeks if a 7.0 download is needed | data access holder (rb643) |
| 3.1 census | 1 sbatch job, ~10 min | none |
| 3.2 gap fill | 8,874-task array, ≈750 CPU-h, 1–2 days incl. queue | coordination with rr480; scratch route otherwise |
| 3.3 extraction | half a day incl. tests | 3.1 |
| 3.4–3.5 QC + validation | 1–2 days | 3.3; 7.0 tables for the release flag on new sessions |
| 3.6 adapter/config/tests | 1 day | 3.3 |

Risks, in order of importance:

1. **Tabulated 7.0 is not on rds.** Until it is, the extra ~3,500 six-year
   sessions have no age or QC flag. Everything else can proceed and be
   validated on the 6.0-covered sessions in the meantime.
2. **Write access.** The FS and parcellation trees are managers-only; §3.2
   sidesteps this, but merging outputs back into `derivatives/parcellations`
   requires a manager.
3. **Small parcels.** HCP-MMP parcels are far smaller than DK regions and the
   slope reliability, already the binding constraint at DK resolution
   (median 0.15 at ≥2 visits), will be lower per parcel. Report per-region
   reliability alongside the maps; consider that the C3 test may be better
   powered on the group map than on subject-level parcel slopes.
4. **Annotation provenance.** The 2016 fsaverage HCP annot is the group's
   standard and matches the AHBA maps, which is the reason to use it; document
   its file hash in the census so nobody swaps it for another projection later.

## 5. Files referenced

| path (under `/rds/project/rds-CeXlNYOYMxw`) | what |
|:--|:--|
| `derivatives/freesurfer/` | DAIRC FreeSurfer 7.1.1, 33,825 sessions |
| `derivatives/parcellations/T1/` | rr480's native-space atlas outputs, 30,360 sessions |
| `derivatives/tabulated/` = `Data_Phenotype/ABCD60/` | release 6.0 tables |
| `Code/parcellation_T1/{parcellate.sh,parcellate_array_ind.sh,parcellate_array_run.sh}` | the parcellation pipeline |
| `Code/check_fs_parallel_parcellations.sh`, `Code/subject_status_check/` | rr480's completeness checker and lists |
| `Code/check_bids_derivatives.py`, `derivatives/test_qc.csv` | rb643's derivative checker (includes FS holes) |
| `userdata/rr480/fsaverageSubP/label/{lh,rh}.HCP.fsaverage.aparc.annot` | the fsaverage HCP-MMP1.0 annotation (2016) |
| `logs/parcellation_T1/` | 30,441 array-task logs, July 2026 |

## 6. Backfill run plan — 2026-09-14 (branch `hcp-backfill`; for review before the array is submitted)

rr480 has not replied, so the missing parcellations are run by us, read-only
against his tree and the FreeSurfer derivatives. Everything below runs on
`vertes-sl3-cpu` / `icelake` (SL2 has no CPU-minutes left; cclake is slow).

### 6.1 HCP-MMP backfill

**Input.** `abcd-data-release-7.0/processed/hcp/sessions_to_reparcellate.txt` — 8,874
sessions = 5,439 with empty stub stats + 3,435 never run. Sessions with
incomplete surfaces are skipped by the worker (exit 3), so the ~30 unusable
reconstructions cost nothing.

**Scripts** (committed on this branch):

| file | role |
|:--|:--|
| `tools/hcp_backfill_session.sh` | one session: private `SUBJECTS_DIR` of two symlinks (`fsaverageSubP`, the session), `mri_surf2surf` of the 2016 fsaverage HCP annotation through `sphere.reg`, `mris_anatomical_stats -a -b` per hemisphere; exit 0 only if both hemispheres have 180 `*_ROI` rows |
| `tools/hcp_backfill.sbatch` | array wrapper: `CHUNK=50` sessions per task, 1 core, 4 GB, 1 h 30 wall; skips sessions already complete, so a partial array is simply resubmitted |

Design decisions worth reviewing:

- **No writes to `derivatives/`.** Output mirrors rr480's layout under
  `legacy/hpc/work/parcellations_backfill/T1/<sub>/<ses>/HCP.fsaverage.aparc/`
  (gitignored). The annotation and stats live there, plus the surf2surf and
  stats stderr logs per hemisphere for post-mortems.
- **Same annotation, same FreeSurfer major version.** `fsaverageSubP`'s
  `sphere.reg` is byte-identical to the stock fsaverage, so it contributes
  only the HCP `.annot`; FreeSurfer 7.1.0 module against 7.1.1 surfaces, as in
  the July run. Backfilled and existing sessions are therefore the same
  measurement; §6.3 checks this on the sessions where both exist.
- **The race is designed out.** Each task gets its own `mktemp` `SUBJECTS_DIR`
  removed on exit; nothing is shared between tasks.
- **Idempotent and resumable.** Both the wrapper's completeness test and the
  worker's `-s` test on the annotation require real content (180 rows / a
  non-empty annot), so stubs can never be mistaken for done.

**Resources.** Measured 40 s per session in the test job (not the 5 min
assumed earlier), so 8,874 sessions ≈ 100 CPU-h: 178 array tasks of 50
sessions (`--array=0-177%100`), ≈35 min each, done within an hour or two of
scheduling. Disk: two annotations ≈ 2.4 MB per session → ≈21 GB under
`hpc/work/`; `/rds/user/rajd2` has 214 GB free of its 1 TB quota.

**Mechanism test — passed.** Job 35544957 ran `sub-00LH735Y/ses-00A` (a July
stub) with the identical scripts: 40 s, 180 `*_ROI` rows per hemisphere,
V1 = 1.974 mm and TE1m = 3.426 mm (plausible), stderr shows FreeSurfer
reading the session only through the private symlink, and the session's
`label/` directory in `derivatives/` keeps its 28 July mtime — nothing was
written there. Output at `legacy/hpc/work/parcellations_backfill_test/`.

**Then re-extract.** `abcd.hcp_stats` gains `--extra-parc-root`: sessions are
discovered in both trees, and for a session present in both the backfill tree
wins (it is only ever populated for sessions that were broken or absent in
rr480's). `hcp_session_qc.tsv` gains a `parc_source` column. Expected
coverage after the run: ≈33,790 of 33,825 sessions.

### 6.2 Desikan-Killiany — does the backfill apply, and can we check it against the release?

**No parcellation job is needed for DK.** `recon-all` writes the DK tables
itself (`stats/{lh,rh}.aparc.stats`), and every one of the 33,825 sessions
already has them; rr480's `parcellations/T1/<ses>/aparc/` output re-derives
the same thing and is redundant. What is missing is only the *extraction*
into release-shaped tables, and that is a parse, not a FreeSurfer run
(≈10 min, one 16-core job, same pattern as `hcp_extract.sbatch`).

Plan: `src/abcd/dk_stats.py` (+ `tools/dk_extract.sbatch`) parsing
`aparc.stats` — note the different layout from `mris_anatomical_stats -b`:
`StructName` is the *first* column, and the header carries
`# Measure Cortex, MeanThickness` and `NumVert`, which are the values the
release's `__lh_mean` columns should equal. Output to
`abcd-data-release-7.0/processed/dsk_local/`:

- `mr_y_smri__{thk,area,vol}__dsk.tsv` with the **exact release column names**
  (`mr_y_smri__thk__dsk__bstmps__lh_mean`, …) via the inverse of
  `Release70Adapter.REGION_CODES` and `data/region_labels.csv` (mapping
  completeness for all 34 regions is asserted, not assumed);
- a long parquet with every `aparc.stats` measure plus hemisphere `NumVert`,
  `MeanThickness`, and the `aseg.stats` hole counts.

The release table stays the canonical DK source for the pipeline. The local
table is (a) the comparison object, and (b) the only DK thickness that exists
for the ~3,500 six-year sessions newer than the tabulated release.

### 6.3 Validation report

`tools/validate_local_vs_release.py` → `docs/hcp_census/local_vs_release.md`
(+ TSVs), run once now against the 6.0 tables on rds and again when the 7.0
tables arrive:

1. **DK local vs release**, joined on participant × session: fraction of
   values equal at the release's 3-decimal precision, per-region max |Δ|,
   hemisphere mean vs `__lh_mean`/`__rh_mean`, and the list of sessions
   present on one side only. Expectation from the one-session check:
   identical. Disagreement would mean the release reprocessed or replaced
   sessions after these surfaces were produced, which matters for every DK
   result in this repo.
   Two known asymmetries to report rather than hide: 3,524 FreeSurfer sessions
   are absent from the 6.0 tables (newer scans), and **90 sessions are in the
   6.0 tables but have no FreeSurfer directory on rds** — worth knowing which.
2. **HCP backfill vs July run** on sessions with valid output in both trees
   (some stub sessions have a valid *other* session of the same subject; a
   direct same-session overlap exists only if we deliberately re-run a sample
   of 200 already-good sessions — proposed, cheap, and the only way to prove
   the two runs are interchangeable).
3. **HCP vs DK** on the same surfaces: area-weighted whole-cortex means
   close to within the medial-wall difference (the 0.9993 correlation in §3.5 becomes an
   equality check once DK comes from the same files).

### 6.4 Sequence

| step | what | gate |
|:--|:--|:--|
| 1 | one-session test job 35544957 | done: 180 rows per hemisphere, 40 s |
| 2 | check `hpc-work` quota headroom ≥ 30 GB | done: 214 GB free |
| 3 | submit the 178-task array (+ a 4-task array re-running 200 good sessions for §6.3.2) | **your go-ahead on this plan** |
| 4 | DK extraction job (independent of 3) | |
| 5 | re-run `hcp_extract.sbatch --extra-parc-root …` | array finished; FAILED lines reviewed and resubmitted once |
| 6 | validation report, README coverage numbers, commit, push `hcp-backfill`, PR | |

### 6.5 Failure modes considered

- `mri_surf2surf` cannot write into a read-only subject: it writes only to
  `--trgsurfval`, which is our path. `mris_anatomical_stats -b` writes only to
  stdout. Verified by the test job (nothing appears in `derivatives/`).
- A task dying mid-session leaves a partial annot → the worker's `-s` test
  would accept a non-empty partial file. Mitigation: the wrapper's 180-row
  test governs "done", and the worker regenerates stats every time; a corrupt
  annot fails the row count and is reported as FAILED for a manual re-run
  with the annot deleted.
- Sessions of the same subject in different tasks: independent
  `SUBJECTS_DIR`s, so no interaction.
