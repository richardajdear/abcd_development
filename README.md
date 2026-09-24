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

## Status — 2026-09-15: the 7.0 re-run is complete, on two parcellations

> **Genetics numbers below are dated 2026-09-15.** The current state (HCP-MMP,
> single-LMM slope, 2025 SCZ GWAS, full control panel) is
> [`genetic_analysis/README_HPC.md`](genetic_analysis/README_HPC.md) §2; where
> they differ, that file wins.

**What happened.** Every result in this repo dated before 2026-09-14 was
computed from the ABCD **6.0** tabulated tables, which sat in a directory
labelled 7.0 (column-identical; only the six-year row count tells them apart,
4,086 vs 7,607 thickness rows). The FreeSurfer surfaces and the genotypes were
7.0 throughout. On 2026-09-14 the 7.0 tables were installed and every local
analysis re-run; on 2026-09-14/15 the **cluster genetics were re-run end to
end** — every step of [`genetic_analysis/README_HPC.md`](genetic_analysis/README_HPC.md)
§3, on both the Desikan–Killiany (DK, 68 regions) and the HCP-MMP (Glasser, 358
parcels) parcellations of the same children. The dated run log in that file's
§8 (job ids, failures, diagnoses) is the authoritative record; its closing
statement (§7 item 5) is the one-paragraph verdict.

**What the 7.0 tables buy** (settled specification, `ct_70_noglobal_mv2_genetic`;
full table in [`docs/RERUN_7.0_TABULATED.md`](docs/RERUN_7.0_TABULATED.md)):

| | 6.0 tables | 7.0 tables |
|:--|--:|--:|
| six-year scans in the model | 3,539 | 6,510 |
| children with ≥2 QC-passing scans | 8,192 | **8,716** |
| … with all four scans | 1,830 | **3,005** |
| … phenotyped and genotyped | 8,082 | **8,596** (EUR arm 4,308) |
| median slope reliability | 0.157 | **0.211** |
| effective N for the slope | 1,361 | **1,752** |
| group thinning map, old vs new (Spearman) | | 0.996 |

The map did not change; the per-child slopes did (r = 0.92 on shared children),
and that is exactly where the re-run found its gains: subject-level precision.

**Genetics on the 7.0 phenotype — what changed against the 6.0 benchmark.**
Tables (6.0-vs-7.0 comparison, now legacy): `legacy/genetic_analysis/tables/summary_70tab_{dk,hcp}.tsv`;
side by side in `…/results_70tab_hcp/compare/table_dk_vs_hcp.tsv`.

- **Every matched SCZ/MDD polygenic-score SE tightened**, no sign flipped.
  **SCZ PRS → faster thinning is robust in the pooled arm under all four
  methods** (β −0.022 to −0.035 SD/SD, p_adj 0.004–0.037; min-p permutation
  0.0065), on **both** parcellations.
- **MDD PRS → faster thinning crossed into significance under three of four
  methods** in the pooled arm on DK (PRS-CS, SBayesR, SBayesRC), where 6.0 had
  one. On HCP it is one of four.
- **SCZ in the EUR arm attenuated to borderline on DK** (0/4 methods, p_adj
  0.065–0.074) and is 3/4 on HCP. A one-input-at-a-time decomposition showed
  the DK attenuation is the phenotype values themselves, not a pipeline
  difference. Effects of 0.02–0.03 SD/SD sit either side of 0.05 depending on
  the whole-cortex average taken; **only the atlas-invariant claims are
  licensed**, and the marginal cells are reported side by side.
- **SNP heritability of `global_slope` (GCTA GREML, dense imputed GRM, 6,011
  PC-AiR-unrelated children): 0.137 → 0.166 ± 0.045** with the SE unchanged —
  the signature of a de-attenuated phenotype. It is atlas-invariant (HCP
  0.156). Baseline thickness 0.223 (positive control holds).
- **No genome-wide hit on `global_slope`** in any of four scans (two atlases
  × two arms; λ_GC 0.99–1.05 everywhere). LDSC h² z is 0.57 (DK) / 1.78
  (HCP), far below the ≈4 needed, so **rg with SCZ/MDD is uninformative, not
  null**. The `SCZ_locus_pool` gene-set enrichment of the thinning GWAS
  strengthened (p 0.0074 → 0.0030 → 0.0005 on HCP) and is the most consistent
  gene-level result. AHBA C1–C3 gene-property and the ahba_pls H3 signature
  test are null in both directions on both atlases.
- The controls behave as documented: ASD null, ALZ associated at ~80 % of
  SCZ's magnitude, EA in the opposite direction.

**HCP-MMP.** The backfill is complete (33,795 of 33,825 sessions), and the
genetics pipeline is parcellation-agnostic: `PARC=hcp` in
`genetic_analysis/config.local.sh` runs every step unchanged on the HCP
export (`configs/ct_70_hcp_noglobal_mv2.yaml`, hippocampal parcel excluded
because it is thickness 0 in ~7,800 sessions). Whole-cortex readouts agree
across atlases; the slope PCs are atlas-specific decompositions except PC2
(`docs/figures/pc_loadings_dk_vs_hcp.html`), and structural-covariance PCs
are the same components by algebra. One candidate locus, chr7:35.5 Mb for
HCP's PC3 in the EUR arm, is recorded and not led with.

**What stands from the earlier work**, re-checked on the new tables:

- The design decisions — no global covariate, no family random effect, ≥2
  visits, the coded QC stack — hold (§Key findings; `docs/reliability_grid.csv`,
  `docs/h2_family_effect_contrast.csv`).
- **Absolute thinning rate vs AHBA C3: ρ = −0.546, p_spin = 0.002**, identical
  to before; the slope-component couplings (PC3–C2 +0.854, PC2–C1 −0.812)
  too, and both replicate on HCP-MMP (+0.72, −0.78).
- The imaging-transcriptomics arm (`ahba_pls/`) re-derives the NSPN-PLS2 /
  AHBA-C3 signature from the new maps — see the dated section at the end of
  [`ahba_pls/README.md`](ahba_pls/README.md).

**What is superseded.** All cluster genetics numbers computed before
2026-09-14 (h², GWAS, LDSC, MAGMA, PRS) were on 6.0-vintage phenotypes and
live under [`legacy/`](legacy/README.md). Do not quote them as current; the
7.0 tables above carry the 6.0 value alongside every row for the comparison.

## Where to look

| you want | go to |
|:---|:---|
| **what changed with the 7.0 tables, old vs new** | [`docs/RERUN_7.0_TABULATED.md`](docs/RERUN_7.0_TABULATED.md), `docs/vintage_comparison.csv` |
| **the cluster genetics: current results, what has been tried, next analyses, the steps to run** | [`genetic_analysis/README_HPC.md`](genetic_analysis/README_HPC.md); every number in [`genetic_analysis/current_results.tsv`](genetic_analysis/current_results.tsv) |
| the PRS results on one slide (4 methods × 2 parcellations, controls, all caveats) | [`docs/figures/slide_prs_methods.png`](docs/figures/slide_prs_methods.png), generator (legacy, PGC3-era) [`legacy/genetic_analysis/slide_prs_methods.py`](legacy/genetic_analysis/slide_prs_methods.py) |
| whether the slope construction matters (mean-of-BLUPs vs one LMM on the per-scan mean) | [`docs/figures/slide_prs_orderops.png`](docs/figures/slide_prs_orderops.png), generator (legacy) [`legacy/genetic_analysis/slide_prs_orderops.py`](legacy/genetic_analysis/slide_prs_orderops.py) |
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
why region selection is not a lever) is in `genetic_analysis/README_HPC.md` §3 and §6.

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
  ~6,350 sessions; exclude it or treat 0 as missing.
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
python -m pytest tests/ -q              # 170 tests
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
