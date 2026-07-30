# abcd_development

**Which genes drive adolescent cortical development?**

Imaging genetics of *longitudinal change* in the ABCD study. The phenotype is a
per-subject **rate** of cortical thinning — a random slope of thickness on age
estimated from repeated scans — not thickness at one timepoint. Release 7.0
supplies enough subjects with three or more imaging visits to make this
feasible.

The hypothesis is that genes driving this rate are (i) enriched for GWAS signal
from disorders with adolescent onset, schizophrenia and major depression, and
(ii) connected to prior transcriptional results from this group: AHBA component
**C3**, and the leading component of an independent snRNA-seq analysis.

Why the design matters: the genetic architecture of a developmental *rate* need
not resemble that of a static measure. Almost all published brain-imaging GWAS
use cross-sectional phenotypes, which average over exactly the variation of
interest.

> **Read [`docs/REPORT_7.0.md`](docs/REPORT_7.0.md) first.** It states what was
> built, what was found, what is blocked, and — in §11 — seven corrections to
> earlier versions of itself, two of which overturned substantive conclusions.
> Every number in it is sourced to a committed table.

## Status

Modelling framework settled and validated on release 7.0. Phenotype selection
for imaging genetics complete. Genetic analyses (GCTA heritability, GWAS, MAGMA)
specified and ready to dispatch to the cluster, not yet executed.

## Quick orientation

| you want | go to |
|:---|:---|
| the findings, their caveats and the corrections | [`docs/REPORT_7.0.md`](docs/REPORT_7.0.md) |
| how the mixed model works and why | [`notebooks/01_longitudinal_model.qmd`](notebooks/01_longitudinal_model.qmd) |
| spatial nulls and the map-to-gene tests | [`notebooks/02_maps_and_genes.qmd`](notebooks/02_maps_and_genes.qmd) |
| heritability and phenotype choice | [`notebooks/04_heritability.qmd`](notebooks/04_heritability.qmd) |
| running the genetic analysis on the cluster | [`hpc/README.md`](hpc/README.md) |
| the 5.1 draft this supersedes | [`docs/REPORT_5.1_legacy.md`](docs/REPORT_5.1_legacy.md) |

## Headline results (release 7.0)

- **Slope reliability is the binding constraint, and 7.0 relaxes it.** Median
  regional slope reliability is 0.149 at ≥2 visits, 0.229 at ≥3 and 0.255 at 4 —
  low enough that it, not sample size, sets the ceiling on any genetic analysis.
  Effective N rises from 833 on release 5.1 to 1,361 on 7.0 at the same ≥2-visit
  filter — a 1.63× gain from extra timepoints on nearly the same subjects.
- **Use `min_visits: 2`, not 3.** The ≥3 filter buys per-subject precision but
  discards 37% of subjects (8,192 → 5,195); effective N falls. The group map is
  unchanged either way (ρ = 0.998).
- **The family random effect must be omitted for genetic phenotypes.** With it
  enabled, the baseline-thickness control returns *h²* = 1.46 — impossible;
  omitting it gives 0.74, matching published twin estimates. The random-effects
  structure is part of the phenotype definition, not a fitting detail. It leaves
  the group-level map *exactly* unchanged (ρ = 1.000) while destroying the
  subject-level genetic signal, which is why it went undetected in the 5.1 draft.
- **The developmental map is transcriptionally patterned.** Absolute thinning
  rate versus AHBA C3: ρ = −0.546, p_spin = 0.0016, closely reproducing the
  thesis 5.1 result (−0.550, p = 0.004) across a release change and a full
  pipeline rewrite. C2 also tracks the rate (−0.330, p_spin = 0.0034); C3 is not
  the only component that does. This *overturns* the previous report's negative
  finding, which was an artefact of the global-thickness covariate.
- **The strongest associations in the project are on the slope components**:
  slope PC3 versus C2 at ρ = +0.854, slope PC2 versus C1 at ρ = −0.812 (both
  p_spin < 0.0002, at the 5,000-rotation floor).
- **Developmental slopes are globally coupled.** Of 2,278 region pairs only 13
  are negative (most negative −0.076); mean off-diagonal r = 0.205, homotopic
  pairs 0.455. Subjects who thin faster in one region thin faster nearly
  everywhere — the same fact that makes a global covariate so damaging.
- **The slope PCs *are* the structural-covariance PCs**, not merely similar to
  them: loading correlation r = 1.000 to numerical precision for PC1–PC3.
- **Site and scanner effects are small**: ICC 0.021 across 18 sites and 0.0013
  across scanner manufacturers. Manufacturer is the only scanner grouping either
  release ships — no device-serial column exists, so a serial-level test is not
  reproducible from these data.
- **The global mean slope should go to GWAS first**: held-out *h²* = 0.464 (SD
  0.10 over 40 family-level splits), statistically indistinguishable from the
  best component while being the only candidate whose definition does not consume
  the data. Ranking regions by heritability performs worse *and* is unstable.

The heritability estimates above are Falconer twin estimates, an upper bound
on SNP heritability; replacing them with GCTA GRM estimates is the first
cluster job.

## Layout

```
src/abcd/        # Python: assembly, QC, phenotypes, spatial stats, gene work
  io.py            release adapters (5.1, 7.0) — add a release by adding a class
  assemble.py      tidy long table from release files
  config.py        RunConfig; every analysis choice hashes into the run_id
  phenotype.py     BLUP extraction, reliability, slope PCA
  heritability.py  Falconer estimates, held-out splits, DZ/sibling handling
  spatial.py       spin nulls
  genemaps.py      AHBA components, map-to-gene tests
  covariance.py    structural covariance of slopes
  gcta_export.py   GCTA/MAGMA input files
  brainplot.py     DK surface rendering
R/               # model fitting: lme4, because it handles crossed random
                 # effects with correlated slopes and reports singularity honestly
configs/         # one YAML per specification; the run_id is a hash of it
docs/            # REPORT_7.0.md + every table and figure it cites
tools/           # regenerators for every table and figure in the report
hpc/             # SLURM pipeline for GRM, REML, GWAS, MAGMA
notebooks/       # explanatory documents, not analysis scripts
tests/           # 144 tests, incl. provenance and README checks
```

The Python/R seam is Parquet in `out/<run_id>/`. Computation is kept separate
from plotting throughout.

## Pipeline

One variable selects the run and every step picks it up:

```bash
export ABCD_CONFIG=ct_70_noglobal_mv2_genetic     # the settled specification
make all                                          # assemble -> fit -> phenotype -> gcta
```

`make help` prints the active config, the run directory it resolves to, and every
available config. The Makefile puts `src/` on the path itself, so it needs no
install.

The same four steps run directly, which is what the notebooks and HPC scripts
do. These need the package importable — either `pip install -e .` once, or
`PYTHONPATH=src` on each call:

```bash
python -m abcd.assemble              # tidy long table
Rscript R/fit_lmm.R --cores 8        # per-region lme4 fits
python -m abcd.phenotype             # BLUPs + reliability
python -m abcd.gcta_export           # GCTA/MAGMA inputs
```

`python -m abcd.run_dir` prints the run directory for the active config and
nothing else, which is how `R/fit_lmm.R` and the HPC scripts resolve the same
run without reimplementing the config hash.

Every analysis choice is a `RunConfig` field that hashes into the `run_id`, so
runs cannot silently overwrite each other and `out/<run_id>/config.yaml` records
what produced the numbers. The metric (`thickness`, `t1t2_ratio`, …) and
parcellation (`dsk`, `dst`, `fzy`, `hcp`) are config fields; adding a release
means one adapter class in `src/abcd/io.py`.

### Configs that matter

| config | what it is for |
|:---|:---|
| `ct_70_noglobal_mv2_genetic.yaml` | **the settled specification** — no global covariate, no family effect, ≥2 visits |
| `ct_70_noglobal_mv{2,3,4}.yaml` | the visit-filter comparison of §4.3 |
| `ct_70_global_mv3_genetic.yaml` | global-covariate contrast, for §4.1 only |
| `ct_51_noglobal_mv2_matched.yaml` | release 5.1 on the *same* specification, for the cross-release check |

`ABCD_ROOT` is **not** required: release 7.0 is vendored at
`abcd-data-release-7.0/` (gitignored — access-controlled), and 5.1 is found
under `~/Git/ABCD`. Releases are searched across all roots, so the two need not
share a parent; export `ABCD_ROOT` only to add a third location. Any step
accepts an explicit argument that overrides the export
(`python -m abcd.phenotype out/<run_id>`, `--config ct_70_baseline`).

## Reproducing the report

Every table and figure in `docs/` is produced by script from committed data —
never from a run directory — in this order:

```bash
python tools/regen_report_tables.py     # spatial/covariance/site tables
python tools/regen_h2_tables.py         # heritability tables (bootstrap; slow)
python tools/regen_report_figures.py    # table-based figures
python tools/regen_brain_maps.py        # DK surface maps
```

This is enforced rather than trusted: `tests/test_docs_provenance.py` checks
that every cited figure and table exists and that every figure on disk has a
generator. It exists because 11 figures were once drawn in ad-hoc cells and
silently survived corrections to the numbers underneath them.

```bash
python -m pytest tests/ -q                                 # 144 tests collected
python tools/check_notebook_chunks.py notebooks/01_*.qmd    # notebooks run
```

## What runs next, on the cluster

1. **GCTA GRM heritability** on the recommended phenotypes, replacing the
   Falconer estimates. `gcta_export.py` refuses to export from a run with the
   family effect enabled.
2. **GWAS** on the global mean slope first, then slope PC3 and PC2, with a mixed
   model to handle relatedness rather than pruning relatives.
3. **MAGMA gene-level enrichment**, then the SCZ and MDD tests — the primary
   hypothesis.
4. **Compare the gene sets** to AHBA C1–C3 and to the snRNA-seq leading
   component.

Step 3 should run on both the global mean *and* the components: they may
implicate different genes, and which one connects to the disorder GWAS is the
empirical question the project exists to answer. See
[`hpc/README.md`](hpc/README.md) for the SLURM scripts, resource estimates and
the input checks to run first.

## Related prior work

The transcriptional targets tested here come from two earlier projects, and are
inputs to this one rather than part of it:

- **AHBA C1–C3** — spatial transcriptomic components of the Allen Human Brain
  Atlas, from this group's previous work.
- **snRNA-seq maturation axis (PC1)** — a single-cell neuron⇄glia maturation
  contrast that recapitulates AHBA C3, in
  [`transcriptional_maturation`](https://github.com/richardajdear/transcriptional_maturation).
