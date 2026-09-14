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

## Status — 2026-09-14: re-run on the true 7.0 tabulation

**What happened.** Every result in this repo dated before 2026-09-14 was
computed from the ABCD **6.0** tabulated tables, which sat in a directory
labelled 7.0. The two tabulations are column-identical; only the six-year row
count tells them apart (4,086 vs 7,607 thickness rows). The FreeSurfer
surfaces and the genotypes on the cluster were 7.0 throughout. The 7.0
tabulated tables were downloaded on 2026-09-14 and every **local** analysis
was re-run on them; the **cluster genetics has not yet been re-run** and is
specified, step by step, in
[`genetic_analysis/README_HPC.md`](genetic_analysis/README_HPC.md).

**What the 7.0 tables buy** (settled specification, `ct_70_noglobal_mv2_genetic`;
full table in [`docs/RERUN_7.0_TABULATED.md`](docs/RERUN_7.0_TABULATED.md)):

| | 6.0 tables | 7.0 tables |
|:--|--:|--:|
| six-year scans in the model | 3,539 | 6,510 |
| children with ≥2 QC-passing scans | 8,192 | **8,716** |
| … with all four scans | 1,830 | **3,005** |
| … phenotyped and genotyped | 8,082 | **8,596** |
| median slope reliability | 0.157 | **0.211** |
| effective N for the slope | 1,361 | **1,752** |
| group thinning map, old vs new (Spearman) | | 0.996 |

The map did not change; the per-child slopes did (r = 0.92 on shared children),
which is where heritability and polygenic-score analyses spend their power.

**What stands from the earlier work**, re-checked on the new tables:

- The design decisions — no global covariate, no family random effect, ≥2
  visits, the coded QC stack — hold (§Key findings; `docs/reliability_grid.csv`,
  `docs/h2_family_effect_contrast.csv`).
- **Absolute thinning rate vs AHBA C3: ρ = −0.546, p_spin = 0.002**, identical
  to before; the slope-component couplings (PC3–C2 +0.854, PC2–C1 −0.812) too
  (`docs/ahba_vs_maps_noglobal.csv`).
- The imaging-transcriptomics arm (`ahba_pls/`) re-derives the NSPN-PLS2 /
  AHBA-C3 signature from the new maps — see the dated section at the end of
  [`ahba_pls/README.md`](ahba_pls/README.md).

**What is superseded.** All cluster genetics numbers (h², GWAS, LDSC, MAGMA,
PRS) were computed on the 6.0-vintage phenotypes and now live under
[`legacy/`](legacy/README.md). The last canonical result there — SCZ polygenic
score → faster thinning, β −0.035 to −0.047 SD/SD, significant in 3 of 4 PRS
methods in both ancestry arms; MDD same direction, significant only pooled;
ASD and education null or opposite — is the benchmark the re-run must be read
against (`genetic_analysis/README_HPC.md` §5). Do not quote legacy genetics
numbers as current.

**Next step (needs the user's go-ahead):** rsync the tables to CSD3 and start
`genetic_analysis/README_HPC.md` step 0. HCP-MMP runs on the new tables as
well (`configs/ct_70_hcp_noglobal_mv2.yaml`); the DK arm remains primary until
the HCP re-parcellation on CSD3 is complete.

## Where to look

| you want | go to |
|:---|:---|
| **what changed with the 7.0 tables, old vs new** | [`docs/RERUN_7.0_TABULATED.md`](docs/RERUN_7.0_TABULATED.md), `docs/vintage_comparison.csv` |
| **the cluster genetics: state, what was learnt, the steps to run** | [`genetic_analysis/README_HPC.md`](genetic_analysis/README_HPC.md) |
| the findings, their caveats and the corrections (prose numbers are 6.0-vintage; tables and figures are current) | [`docs/REPORT_7.0.md`](docs/REPORT_7.0.md) |
| how the mixed model works and why | [`notebooks/01_longitudinal_model.qmd`](notebooks/01_longitudinal_model.qmd) |
| spatial nulls and the map-to-gene tests | [`notebooks/02_maps_and_genes.qmd`](notebooks/02_maps_and_genes.qmd) |
| heritability and phenotype choice | [`notebooks/04_heritability.qmd`](notebooks/04_heritability.qmd) |
| the imaging-transcriptomics PLS study (NSPN-PLS2 / AHBA-C3 re-derivation, SCZ & MDD enrichment) | [`ahba_pls/README.md`](ahba_pls/README.md), notebook [`ahba_pls/imaging_transcriptomics.qmd`](ahba_pls/imaging_transcriptomics.qmd) |
| HCP-MMP thickness: where it comes from, coverage, how to regenerate | [`docs/PLAN_HCP_thickness.md`](docs/PLAN_HCP_thickness.md) and §HCP-MMP below |
| superseded pipelines and the 5.1 draft | [`legacy/README.md`](legacy/README.md), [`docs/REPORT_5.1_legacy.md`](docs/REPORT_5.1_legacy.md) |

## Key findings so far

**Phenotype and modelling** (7.0 tables, 2026-09-14)

- **Slope reliability, not sample size, binds.** Median regional slope
  reliability is 0.211 at ≥2 visits, 0.242 at ≥3, 0.262 at 4. Effective N is
  1,752 at ≥2 visits — 2.1× release 5.1 at the same filter
  (`docs/handoff_release_comparison.csv`).
- **Use `min_visits: 2`, not 3.** The ≥3 filter discards 25 % of subjects and
  effective N falls (1,752 → 1,504); the group map is unchanged.
- **The family random effect must be omitted for genetic phenotypes.** It
  centres each family at zero, so the subject-level BLUPs lose the
  between-family variance that relatedness explains; the h² of the
  baseline-thickness control becomes impossible while the group map is
  untouched (`docs/h2_family_effect_contrast.csv`).
- **The developmental map is transcriptionally patterned.** Absolute thinning
  rate vs AHBA C3: ρ = −0.546, p_spin = 0.002, now reproduced across two
  releases, two tabulations and a pipeline rewrite. The strongest couplings are
  on the slope components (PC3–C2 ρ = +0.854, PC2–C1 ρ = −0.812).
- **Site explains 2.9 % of global-slope variance, scanner manufacturer 0.1 %**
  (`docs/site_scanner_icc.csv`); the model carries a site random effect.

**Imaging transcriptomics** (`ahba_pls/`; DK arm re-run on the 7.0 maps
2026-09-14, HCP-MMP arm likewise — numbers in its README)

- The NSPN-PLS2 / AHBA-C3 "signature of adolescent thinning" re-derives from
  the ABCD maps: PLS of AHBA expression on thinning rate plus baseline
  thickness recovers a component matching both prior signatures in regional
  scores and gene weights, with the same neuronal-up / glial-down cell-class
  profile.
- It carries SCZ and MDD GWAS signal at about half C3's effect size and adds
  nothing once C3 is in the model: a developmental warrant for C3, not a
  better gene list.
- Parcellation matters for MDD, not SCZ: the HCP-MMP version of the signature
  carries more MDD signal than the DK version. Whether the finer atlas' astrocyte
  sign flip is biology or parcel size is open.

**Genetics** — see the Status section: the legacy result is a benchmark, the
re-run is pending. The reasoning that governs it (PRS as the primary
disorder test, why rg is uninformative at this h², why controls are mandatory,
why region selection is not a lever) is in `genetic_analysis/README_HPC.md` §4.

## HCP-MMP (Glasser) thickness

ABCD tabulates only Desikan (`dsk`). The HCP-MMP1.0 parcellation exists for 7.0
as a **derived** table at `abcd-data-release-7.0/processed/hcp/`, so
`parcellation: hcp` works (`configs/ct_70_hcp_noglobal_mv2.yaml`).

- **Source.** The release FreeSurfer 7.1.1 reconstructions on CSD3
  (`/rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer/`, 33,825 sessions).
  R. Romero-Garcia projected the fsaverage HCP-MMP1.0 annotation onto them;
  `src/abcd/hcp_stats.py` parses the per-session tables into
  `mr_y_smri__{thk,area,vol}__hcp.tsv` in the DK column convention, plus
  `hcp_session_qc.tsv`. Regenerate with `sbatch tools/hcp_extract.sbatch` on
  CSD3 (about 10 min).
- **Coverage is incomplete and flagged.** 24,921 of 33,825 sessions parsed;
  3,435 were never reached by the July 2026 array job and 5,439 have empty
  stubs (a concurrency bug in the parcellation script). The re-run list for
  rr480 is `processed/hcp/sessions_to_reparcellate.txt`. Missing sessions have
  no row, so `assemble` reports fewer scans than DK.
- **Parcel `H`** (hippocampus) lies on the medial wall and has thickness 0 in
  ~6,350 sessions; exclude it or treat 0 as missing.
- With the 7.0 covariates the six-year HCP sessions that previously lacked an
  age row now enter the model; see `ahba_pls/README.md` for the resulting n.
- Details and the plan: [`docs/PLAN_HCP_thickness.md`](docs/PLAN_HCP_thickness.md).

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
                   # hcp_extract.sbatch (CSD3); prs_assoc.R (used by genetic_analysis/)
genetic_analysis/  # the cluster genetics re-run: README_HPC.md, config, GENESIS + PRS scripts
ahba_pls/          # imaging transcriptomics: PLS of AHBA expression on the ABCD thinning maps,
                   # and its SCZ/MDD enrichment -- self-contained, see ahba_pls/README.md
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
| `ct_70_hcp_noglobal_mv2.yaml` | the settled spec on HCP-MMP (Glasser); derived table, partial coverage |

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
python -m pytest tests/ -q              # 167 tests
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
