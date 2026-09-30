# abcd_development

**Which genes drive adolescent cortical development?**

Imaging genetics of *longitudinal change* in the ABCD study. The phenotype is a
per-subject **rate** of cortical thinning — a random slope of thickness on age
estimated from repeated scans — not thickness at one timepoint.

The hypothesis: genes driving this rate are (i) enriched for GWAS signal from
disorders with adolescent onset, schizophrenia and major depression, and
(ii) connected to this group's prior transcriptional results (AHBA components
C1–C3, and the leading component of an independent snRNA-seq analysis).

Why the design matters: the genetic architecture of a developmental *rate* need
not resemble that of a static measure. Almost all published brain-imaging GWAS
use cross-sectional phenotypes, which average over exactly the variation of
interest.

## Current state — 2026-09-28

ABCD release 7.0, 8,716 children with 2–4 QC-passing scans (26,946 scans), 8,596 genotyped
(European-ancestry arm 4,308). Primary specification: HCP-MMP parcellation, one mixed model
on the per-scan cortical mean (`global_slope`), 2025 SCZ GWAS. Every genetics number is a
row of [`genetic_analysis/current_results.tsv`](genetic_analysis/current_results.tsv) and is
explained in [`genetic_analysis/README_HPC.md`](genetic_analysis/README_HPC.md) §2; where this
summary and that file differ, that file wins. Figure 1: [`docs/figures/fig1.png`](docs/figures/fig1.png)
(caption `docs/figures/fig1_caption.md`).

## Key findings

**Phenotype**

- **Slope reliability, not sample size, is the limit.** Median regional slope reliability is
  0.21 at ≥ 2 visits (0.24 at ≥ 3, 0.26 at 4). The single LMM on the cortical mean is more reliable
  and less shrunk than averaging per-parcel BLUPs, so it is the primary trait.
- **Design rules for genetic phenotypes:** `min_visits: 2` (≥ 3 loses 25 % of children for no
  gain), no family random effect (it removes the between-family variance relatedness explains),
  no global covariate. Site explains 2.9 % of global-slope variance, scanner manufacturer 0.1 %.
- **The group thinning map is transcriptionally patterned**: thinning rate vs AHBA C3
  ρ = −0.55, p_spin = 0.002 (DK), reproduced across releases and tabulations.

**Genetics** (`genetic_analysis/`)

- **Heritable, but no locus at this n.** GREML h² = 0.18 ± 0.05 (DK 0.21). No genome-wide hit in
  any scan (λ_GC 0.99–1.03). LDSC h² z is ~0.4, so rg with any disorder is uninformative, not null.
- **Higher schizophrenia polygenic risk predicts faster thinning, not thinner cortex.** β ≈ −0.03
  to −0.04 SD/SD under C+T, PRS-CS and SBayesRC, in both arms and on both atlases, one DK cell borderline (e.g. pooled
  SBayesRC −0.034, p = 0.001). The same scores are null for baseline thickness.
- **It is not specific to schizophrenia.** Depression (−0.025 to −0.043), Alzheimer's with APOE
  (−0.03 to −0.045; gone without APOE) and, in the opposite direction, educational attainment
  (+0.03 to +0.04) have effects of the same size; autism, ADHD and intelligence are null, bipolar
  is significant only under C+T.
- **No robust gene-level signal.** The one enriched set (SCZ locus genes, MAGMA p = 0.005)
  does not survive swapping the 1000G LD panel for ABCD's own (p = 0.15). AHBA C1–C3 as gene
  properties are null. Imputed C4A expression is null (β −0.011, 95 % CI −0.039 to +0.017).
- **The transcriptomic link is group-level only.** No individual-level genetic readout connects
  to C3: the regional SCZ/MDD PRS maps do not resemble C3 or PLS2 (7 of 192 tests pass both
  nulls, about chance), and see *Between-child variation* below.

**Timing versus rate** (`directions/d2_prs_age_puberty/`, 2026-09-30)

- **Polygenic risk scales the thinning rate uniformly from 9 to 17; it does not shift its timing.**
  In a one-stage LMM on 26,597 scans the SCZ score changes the rate by −0.28 µm/yr per SD (p = 7 × 10⁻⁴;
  MDD −0.17, p = 0.04) with no thickness difference at 12.8. The normative rate peaks at 12–14, but the PRS
  effect is the same before and after (change after 13: p = 0.57 / 0.86), which rules out a phase advance
  beyond ~1 month per SD. Pubertal stage at the scan does not carry the effect; MDD risk does predict
  earlier puberty, and earlier puberty faster thinning, but adjusting for puberty leaves PRS → slope unchanged.

**Environment and direction of effect** (`directions/d5_adversity_direction/`, 2026-09-30)

- **Adversity and polygenic risk act on the rate largely separately.** SES and area deprivation relate to
  the rate only between sites (the slope keeps between-site rate differences); within site, SES predicts
  thickness *level* (+0.062) but not rate. Parent-reported negative life events predict faster thinning
  (−0.034, p = 0.003, robust to site and scan quality). SES + adversity explain 4 % of the SCZ polygenic
  effect, and 0 of 8 PRS × environment terms reach p < 0.05.
- **Thinning precedes symptoms, not the reverse, between children.** Faster thinning predicts later
  depressive symptoms (−0.033, p = 0.003, unchanged with environment adjusted); baseline symptoms predict
  the thinning rate for 0 of 4 scales. Within-child cross-lags are small in both directions and need a
  latent-curve check. The EA-with-SES test awaits the CSD3 EA score file.

**Symptoms** (exploratory; `ahba_pls/` scripts 23–25, 30; Figure 1g)

- Faster global thinning goes with more depressive symptoms by ages 15–17, given baseline
  (CBCL depressive problems +0.030 SD/SD, p = 0.006; parent-reported MDD OR 1.11 per SD). All effects
  are small (|β| ≤ 0.05).
- The spatial pattern of symptom-linked thinning follows the normative thinning map (with
  baseline CT adjusted, ρ = 0.25–0.41 for anxiety, internalising and p-factor), not PLS2 or C3.

**Imaging transcriptomics** (`ahba_pls/`)

- The NSPN-PLS2 / AHBA-C3 "signature of adolescent thinning" re-derives from the ABCD maps: PLS of
  AHBA expression on thinning rate and baseline thickness recovers a component matching both prior
  signatures in regional scores and gene weights, with the same neuronal-up / glial-down profile.
- Its gene weights carry SCZ and MDD GWAS signal at about half C3's effect and add nothing once C3
  is in the model: a developmental warrant for C3, not a better gene list. The HCP-MMP version
  carries more MDD signal than DK.

**Between-child variation in thinning** (exploratory, 2026-09-28; `genetic_analysis/c3axis/`,
README_HPC §8)

- Children's regional thinning rates covary as a global factor (PC1, 13 % of variance) plus the
  AHBA C1 axis (PC2; on DK it is also the T1w/T2w myelin axis) and the C2 axis (PC3). C3 appears only
  once each child's mean is removed, as a small left-hemisphere component (5th, ρ = 0.52 with C3,
  p_spin = 0.002; diffusion maps give the same) that is distinct from the normative dCT/PLS2 map.
- These components add no genetic or clinical signal: SCZ and MDD scores are null on the C3, C1
  and C2 axes (SCZ → C3 axis β −0.007 to −0.013); the SCZ effect is a uniform whole-cortex shift.
  Puberty relates to the global factor and the C1 axis (|β| ≤ 0.08), not to C3, and no component
  predicts CBCL beyond chance.
- The one open C3 test is at the gene level: whether the C3-gene partition of the SCZ score
  carries its global-slope effect (README_HPC §8.3, C3-D).

## Where to look

| you want | go to |
|:---|:---|
| **the cluster genetics: current results, what has been tried, next analyses, the steps to run** | [`genetic_analysis/README_HPC.md`](genetic_analysis/README_HPC.md); every number in [`genetic_analysis/current_results.tsv`](genetic_analysis/current_results.tsv) |
| what mechanism we are proposing, what would falsify it, and the ranked next directions | [`docs/DIRECTIONS.md`](docs/DIRECTIONS.md) |
| D2: one-stage PRS × age / PRS × puberty models (timing vs rate) | [`directions/d2_prs_age_puberty/`](directions/d2_prs_age_puberty/README.md) |
| D4: thinning rate vs cognitive gain, ages 10–16 (NIH Toolbox) | [`directions/d4_cognitive_gain/`](directions/d4_cognitive_gain/README.md) |
| D5: adversity/SES × polygenic risk on the thinning rate; direction of the thinning–symptom link | [`directions/d5_adversity_direction/`](directions/d5_adversity_direction/README.md) |
| D7: partitioned SBayesRC scores (SynGO, neuronal/glial, snRNA-seq, AHBA C3 gene sets) — CSD3 pipeline, pre-registered, not yet run | [`directions/d7_partitioned_prs/`](directions/d7_partitioned_prs/README.md) |
| Figure 1 (phenotype, PRS, symptoms) | [`docs/figures/fig1.png`](docs/figures/fig1.png), generator [`genetic_analysis/fig1.R`](genetic_analysis/fig1.R) |
| the PRS results on one slide (methods × parcellations, controls) | [`docs/figures/slide_prs_methods_1lmm.png`](docs/figures/slide_prs_methods_1lmm.png) |
| the imaging-transcriptomics PLS study and the symptom / PRS maps | [`ahba_pls/README.md`](ahba_pls/README.md), notebook [`ahba_pls/imaging_transcriptomics.qmd`](ahba_pls/imaging_transcriptomics.qmd) |
| slope components, the C3 axis, puberty and CBCL checks | [`genetic_analysis/c3axis/`](genetic_analysis/c3axis/), README_HPC §8 |
| the C4A (complement) test | [`c4_imputation/README.md`](c4_imputation/README.md) |
| how the mixed model works and why | [`notebooks/01_longitudinal_model.qmd`](notebooks/01_longitudinal_model.qmd) |
| spatial nulls and the map-to-gene tests | [`notebooks/02_maps_and_genes.qmd`](notebooks/02_maps_and_genes.qmd) |
| heritability and phenotype choice | [`notebooks/04_heritability.qmd`](notebooks/04_heritability.qmd) |
| the phenotype report (prose numbers are 6.0-vintage; tables and figures are current) | [`docs/REPORT_7.0.md`](docs/REPORT_7.0.md) |
| HCP-MMP thickness: where it comes from, coverage, how to regenerate | [`docs/PLAN_HCP_thickness.md`](docs/PLAN_HCP_thickness.md) and §HCP-MMP below |
| data-vintage history (the 6.0 tables once labelled 7.0, fixed 2026-09-14) | [`docs/RERUN_7.0_TABULATED.md`](docs/RERUN_7.0_TABULATED.md) |
| superseded pipelines, the 6.0-vintage genetics, the 5.1 draft | [`legacy/README.md`](legacy/README.md), [`docs/REPORT_5.1_legacy.md`](docs/REPORT_5.1_legacy.md) |

## HCP-MMP (Glasser) thickness

ABCD tabulates only Desikan (`dsk`). The HCP-MMP1.0 parcellation exists for 7.0
as a **derived** table at `abcd-data-release-7.0/processed/hcp/`, so
`parcellation: hcp` works (`configs/ct_70_hcp_noglobal_mv2.yaml`).

- **Source.** The release FreeSurfer 7.1.1 reconstructions on CSD3
  (`/rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer/`, 33,825 sessions).
  The fsaverage HCP-MMP1.0 annotation is carried to each session's sphere and
  summarised with `mris_anatomical_stats`: R. Romero-Garcia's July 2026 run
  covers 24,921 sessions (`derivatives/parcellations/T1/`), and
  `tools/hcp_backfill.sbatch` (2026-09-14) did the remaining 8,874 the same
  way into `legacy/hpc/work/parcellations_backfill/T1/`. A 200-session overlap
  re-run is bit-identical to the July output, so the two trees are one
  measurement. `src/abcd/hcp_stats.py` parses both (`--extra-parc-root`) into
  `mr_y_smri__{thk,area,vol}__hcp.tsv` in the DK column convention, plus
  `hcp_session_qc.tsv` (parcel/vertex counts, surface holes, `parc_source`).
- **Coverage is complete: 33,795 of 33,825 sessions**; the other 30 are
  reconstructions without surfaces. Regenerate on CSD3 with
  `sbatch tools/hcp_extract.sbatch --extra-parc-root legacy/hpc/work/parcellations_backfill/T1`
  (about 5 min).
- **Hemisphere and whole-cortex means are surface-area-weighted over parcels.**
  That is the release's own definition: the published DK `__lh_mean` equals
  the area-weighted mean of the 34 regions, not FreeSurfer's vertex-weighted
  `Cortex MeanThickness` (verified 2026-09-14; see `docs/hcp_census/`).
- **Parcel `H`** (hippocampus) lies on the medial wall and has thickness 0 in
  7,792 sessions (left; 8,668 in either hemisphere); it is excluded from the fits.
- **DK from the same surfaces, as a check on the release.** `src/abcd/dk_stats.py`
  (`sbatch tools/dk_extract.sbatch`) parses FreeSurfer's own `aparc.stats` for
  every session into `processed/dsk_local/` with the release column names.
  Against the 7.0 tables it is identical at 3 dp for **33,791 of 33,792 shared
  sessions**; the one exception (`sub-9RRGXBK5/ses-06A`) was re-tabulated in
  7.0 while the surfaces on rds still hold its 6.0-era reconstruction. So the
  DK phenotype and the HCP phenotype come from the same surfaces, and the
  release table stays the canonical DK source. Report:
  `docs/hcp_census/local_vs_release_7.0.md` (`tools/validate_local_vs_release.py`).
- Investigation, run plan and results: [`docs/PLAN_HCP_thickness.md`](docs/PLAN_HCP_thickness.md).

## Running it

```bash
export ABCD_CONFIG=ct_70_noglobal_mv2_genetic     # the settled specification
make all                                          # assemble -> fit -> phenotype -> gcta
```

`make help` prints the active config and the run directory it resolves to. Every
analysis choice is a `RunConfig` field that hashes into the `run_id`, so runs
cannot silently overwrite each other and `out/<run_id>/config.yaml` records what
produced the numbers. **The hash does not encode data vintage**: the manifest's
`data_vintage` key does, and assembly refuses a 7.0 run from 6.0-sized tables.

The same four steps run directly (needs `PYTHONPATH=src` or `pip install -e .`):

```bash
python -m abcd.assemble              # tidy long table
Rscript R/fit_lmm.R --cores 8        # per-region lme4 fits (reads $PY for the interpreter)
python -m abcd.phenotype             # BLUPs + reliability
python -m abcd.gcta_export           # GCTA/MAGMA inputs
```

To redo every specification the report needs (eleven fits, ~15 min):

```bash
tools/rerun_local.sh
```

Environments on this machine: `~/mambaforge/envs/abcd` has the Python
dependencies **and** an R with lme4/arrow/optparse/ggseg (the fitter and the
brain maps); `~/.claude-science/conda/envs/abcd-spatial` has pytest and the
spatial/statsmodels stack (tests, `ahba_pls/` analysis); `ahba-pls-r` holds the
R plotting stack for `ahba_pls/` figures. `tools/rerun_local.sh` and `make`
accept `PY`/`RSCRIPT` overrides.

`ABCD_ROOT` is **not** required: 7.0 is vendored at `abcd-data-release-7.0/`
and the 6.0 tables at `abcd-data-release-6.0/` (both gitignored —
access-controlled); 5.1 is found under `~/Git/ABCD`. A config with
`release: "6.0"` (`configs/ct_60_noglobal_mv2_genetic.yaml`) reproduces the
pre-2026-09-14 sample exactly, for comparison only.

## Layout

```
src/abcd/          # Python: assembly, QC, phenotypes, spatial stats, gene work
R/                 # model fitting (lme4)
configs/           # one YAML per specification; the run_id is a hash of it
docs/              # REPORT_7.0.md, RERUN_7.0_TABULATED.md + every table and figure they cite
tools/             # regenerators for every table and figure; rerun_local.sh; compare_vintage.py;
                   # CSD3 jobs: hcp_extract / hcp_backfill / dk_extract .sbatch; validate_local_vs_release.py;
                   # prs_assoc.R (used by genetic_analysis/)
genetic_analysis/  # the cluster genetics: README_HPC.md, config, GENESIS + PRS scripts;
                   # c3axis/ = slope components and the C3-axis phenotypes (exploratory)
c4_imputation/     # imputed C4A expression vs the thinning rate
ahba_pls/          # imaging transcriptomics: PLS of AHBA expression on the ABCD thinning maps,
                   # and its SCZ/MDD enrichment -- self-contained, see ahba_pls/README.md
directions/        # one sub-directory per direction in docs/DIRECTIONS.md (d2_prs_age_puberty, ...); code, tables, figure, README each
notebooks/         # explanatory documents, not analysis scripts
tests/             # pytest suite, incl. provenance and README checks
legacy/            # superseded: hpc/, hpc_v2/, hpc_v3/ (6.0-vintage genetics), handoff tables
out/               # run directories (gitignored); out/legacy_6.0_tabulated/ holds the old fits
```

The Python/R seam is Parquet in `out/<run_id>/`. Computation is kept separate
from plotting throughout.

| config | what it is for |
|:---|:---|
| `ct_70_noglobal_mv2_genetic.yaml` | **the settled specification** |
| `ct_70_noglobal_mv{2,3,4}.yaml`, `ct_70_noglobal_mv{3,4}_genetic.yaml` | the visit-filter and family-effect comparisons |
| `ct_70_global_mv3_genetic.yaml`, `ct_70_baseline.yaml`, `ct_70_genetic.yaml` | global-covariate contrasts |
| `ct_51_noglobal_mv2_matched.yaml` | 5.1 on the same specification |
| `ct_60_noglobal_mv2_genetic.yaml` | the settled spec on the 6.0 tables (comparison only) |
| `t1t2_70_noglobal_mv2_genetic.yaml` | T1w/T2w ratio, matched to the settled spec (maps for `ahba_pls/`) |
| `ct_70_hcp_noglobal_mv2.yaml` | the settled spec on HCP-MMP (Glasser); derived table, full coverage since 2026-09-14 |

## Reproducing the report

Every table and figure in `docs/` is produced by script from runs on disk and
committed data — never hand-entered:

```bash
export PYTHONPATH=src
python tools/regen_report_tables.py     # spatial/covariance/site tables
python tools/regen_h2_tables.py         # heritability tables (bootstrap; slow)
python tools/regen_report_figures.py    # table-based figures
python tools/regen_brain_maps.py        # DK surface maps
python tools/compare_vintage.py --old out/legacy_6.0_tabulated/thickness_dsk_70_139406217085 \
                                --new out/thickness_dsk_70_139406217085
python -m pytest tests/ -q              # 175 tests
```

The README test count check spawns pytest via `PYTEST_PY` (default `python`);
set it to an interpreter that has pytest if `python` is a bare shim here.

This is enforced, not trusted: `tests/test_docs_provenance.py` checks that every
cited figure and table exists and that every figure on disk has a generator. It
exists because 11 figures were once drawn in ad-hoc cells and silently survived
corrections to the numbers underneath them.

**If you work on this repo through Claude Science**, note that editing a file on
disk does not update its artifact, and neither does committing. Run
`make audit-artifacts` at the end of any session that edits deliverables.

## Related prior work

The transcriptional targets tested here are inputs to this project, not part of
it: **AHBA C1–C3** from this group's previous work, and the **snRNA-seq
maturation axis** in
[`transcriptional_maturation`](https://github.com/richardajdear/transcriptional_maturation).
