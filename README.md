# transcriptional_maturation

A single-cell **PC1 neuronal-maturation axis** and its link to the AHBA **C3**
spatial transcriptomic programme and psychiatric GWAS risk.

This repo builds a six-panel headline figure for a chosen
(dataset, chemistry, brain region) subset of the Velmeshev developmental
snRNA-seq atlas, and quantifies whether the maturation axis carries GWAS
signal via MAGMA gene-property tests.

## The maturation axis (PC1)

For a subset of cells (neurons + macroglia; microglia excluded), we
log-normalise expression and take the first principal component. PC1 is
oriented so that **neurons load positive and glia negative** — it is a
neuron⇄glia maturation contrast. Each cell also gets an **AHBA C3 score** by
projecting its expression onto the C3 gene-weight vector (a spatial cortical
programme from the Allen Human Brain Atlas). The single-cell PC1 axis
recapitulates the AHBA C3 spatial programme (per-cell r ≈ 0.86–0.96 across
subsets), tying single-cell development to a cortical map that is enriched for
schizophrenia and MDD.

## Figure panels

- **a** Pseudotime vs age, faceted by major class (ExN / InN / Macroglia).
  An alternate layout (`assemble_figure(..., panel_a='umap')`) replaces this
  with a column of four UMAPs of all cells coloured by cell type, donor age,
  pseudotime and PC1. It reuses the precomputed global embedding
  (`obsm/X_umap`), which was computed across all regions/datasets — for a
  publication figure a subset-specific embedding should be recomputed.
- **b** Pseudotime vs PC1, coloured by cell type, with per-class fit lines.
- **c** AHBA C3 rendered on the HCP/Glasser cortical surface.
- **d** Per-cell PC1 vs AHBA C3 (overall Pearson r).
- **e** Top ± PC1 gene loadings.
- **f** MAGMA GWAS enrichment of the PC1 axis and AHBA C3 for SCZ, MDD, ASD, AZ
  (points annotated with significance stars).

## Layout

```
src/
  compute.py   # data loading, per-cell PCA + C3 projection, MAGMA runner/parser,
               # Ensembl->Entrez mapping, covariate builder, brain-map loader
  plotting.py  # palette (ExN yellow->red, InN pink->purple gradients),
               # per-panel drawers, compound six-panel assembler
  driver.py    # run(dataset, chemistry, region, pc1_col, out)  — uses precomputed
               # loadings/covariate; run_new_region(...) — computes region-specific
               # loadings + MAGMA covariate on the fly
data/          # gene weights, AHBA HCP cortical scores, Glasser annotation,
               # PC1 gene loadings, cached Ensembl->Entrez map
figures/       # rendered PNG + PDF headline figures
```

**Computation is kept separate from plotting**: `compute.py` produces
DataFrames/arrays, `plotting.py` only draws.

## Usage

```python
import sys; sys.path.insert(0, 'src')
import driver

# Region already covered by the precomputed loadings/covariate (Herring V3 FC, U01 V2 FC):
driver.run('Herring', 'V3', 'FC', 'PC1_herringV3', out='figures/velmeshev_C3_headline.png')
driver.run('U01', 'V2', 'FC', 'PC1_U01V2', out='figures/velmeshev_C3_headline_U01_FC.png')

# A new region — computes its own PC1 loadings + Entrez MAGMA covariate, runs MAGMA, renders:
driver.run_new_region('U01', 'V2', 'CC', out='figures/velmeshev_C3_headline_U01_CC.png',
                      cov_txt='.../input_PC1_U01CC.txt', magma_tag='PC1CC', pc1_col='PC1_U01CC',
                      cache_csv='data/ens2entrez_mygene.csv')
```

Or from the command line:

```
python src/driver.py --dataset U01 --chemistry V2 --region FC \
    --pc1-col PC1_U01V2 --out figures/velmeshev_C3_headline_U01_FC.png
```

## External dependencies (not vendored)

Paths are resolved at import time in `compute.py` and point at the local
analysis tree:

- **Velmeshev h5ad + metadata** — the developmental snRNA-seq atlas.
- **MAGMA** binary + `*.magma.genes.raw` gene-analysis files (SCZ, MDD2025,
  ASD, AZ). Educational-attainment sumstats were not available locally.

Set these via the constants at the top of `src/compute.py` if your paths differ.

## Results summary

| Subset | cells | PC1–C3 r | SCZ | MDD | ASD | AZ |
|---|---|---|---|---|---|---|
| Herring V3 FC | 132,947 | +0.96 | *** | *** | n.s. | n.s. |
| U01 V2 FC | 116,604 | +0.86 | *** | *** | *** | n.s. |
| U01 V2 CC | 70,664 | +0.90 | *** | *** | *** | * |

Stars are MAGMA gene-property p-values for the PC1 covariate
(\*\*\* p<0.001, \* p<0.05). The maturation axis carries schizophrenia and MDD
risk across all three subsets; the deeper U01 data additionally resolves ASD.

---

# ABCD: adolescent cortical development (branch `abcd-longitudinal`)

A second, larger strand of the same question: **which genes drive adolescent
cortical development?** Where the snRNA-seq work above establishes a maturation
axis and links it to AHBA C3 and psychiatric risk, this strand asks whether
*longitudinal imaging* in ABCD can identify the genetic drivers of developmental
change — as opposed to the baseline-only imaging genetics in the literature.

**Read [`REPORT.md`](REPORT.md) first.** It states what was built, what was
found, and what is blocked, with every number sourced.

## Quick orientation

| you want | go to |
|:---|:---|
| the findings and their caveats | [`REPORT.md`](REPORT.md) |
| how the mixed model works and why | [`notebooks/01_longitudinal_model.qmd`](notebooks/01_longitudinal_model.qmd) |
| spatial nulls and the gene tests | [`notebooks/02_maps_and_genes.qmd`](notebooks/02_maps_and_genes.qmd) |
| running the genetic analysis on CSD3 | [`hpc/README.md`](hpc/README.md) |

## Headline results (release 5.1)

- **Slope reliability is the binding constraint**: median 0.12–0.16 across DK
  regions with 2–3 visits, against 0.88–0.91 for the intercept. Release 5.1
  cannot answer the genetic question and more two-visit subjects would not help;
  7.0's extra timepoints attack the right term.
- **The family random effect must be omitted for genetic phenotypes.** Including
  it drives the sibling correlation of the intercept to −0.62 (impossible);
  omitting it recovers *h²* = 0.83, matching published twin estimates. The
  random-effects structure is part of the phenotype definition.
- **The group-mean developmental map is not transcriptionally patterned** along
  AHBA C1–C3, at map or gene level, with a reliable map (split-half ρ = 0.975)
  and a calibrated spin null. Positive control passes (baseline thickness vs C1,
  ρ = −0.58, p = 0.001).
- **Regional slope reliability *is* gene-linked** (p = 0.001, survives
  partialling noise, parcel size and baseline thickness) but fails a
  hemispheric-symmetry check — provisional pending 7.0.

## Pipeline

Install once, then one variable selects the run and every step picks it up:

```bash
pip install -e .                     # so `python -m abcd.*` works anywhere
export ABCD_CONFIG=ct_70_genetic     # the only thing you set

python -m abcd.assemble              # tidy long table
Rscript R/fit_lmm.R --cores 8        # per-region lme4 fits
python -m abcd.phenotype             # BLUPs + reliability
python -m abcd.gcta_export           # GCTA/MAGMA inputs
```

or `make all` for the same four steps. `make help` prints the active config,
the run directory it resolves to, and the available configs.

`ABCD_ROOT` is **not** required: release 7.0 is vendored at
`abcd-data-release-7.0/` in the repo (gitignored — 92 MB, access-controlled),
and 5.1 is found under `~/Git/ABCD`. Releases are searched across all roots, so
the two need not share a parent; export `ABCD_ROOT` only to add a third
location. Any step still accepts an explicit argument, which overrides the
export (`python -m abcd.phenotype out/<run_id>`, `--config ct_70_baseline`).

`python -m abcd.run_dir` prints the run directory for the active config and
nothing else, which is how `R/fit_lmm.R` and the HPC scripts resolve the same
run without reimplementing the config hash.

Python (`src/abcd/`) does assembly, QC, spatial statistics and gene work; R
(`R/`) does model fitting, because `lme4` handles crossed random effects with
correlated slopes and reports singularity honestly. The seam is Parquet.

Every analysis choice is a `RunConfig` field that hashes into the `run_id`, so
runs cannot silently overwrite each other and `out/<run_id>/config.yaml` records
what produced the numbers. Adding release 7.0 means one adapter class in
`src/abcd/io.py`; the metric (`thickness`, `t1t2_ratio`, …) and parcellation
(`dsk`, `dst`, `fzy`, `hcp`) are likewise config fields.

```bash
python -m pytest tests/ -q                                    # 79 tests
python tools/check_notebook_chunks.py notebooks/01_*.qmd       # notebooks run
```
