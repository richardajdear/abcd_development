# Notebooks

Explanatory documents, not analysis scripts. The analysis lives in `src/abcd/`
and `R/`; these render it with the reasoning attached.

| notebook | engine | what it explains |
|---|---|---|
| `01_longitudinal_model.qmd` | R (knitr) | The mixed model: every term, why it is there, what the phenotype means, and why the family random effect is a phenotype decision rather than a fitting detail. **Start here.** |
| `02_maps_and_genes.qmd` | Python (jupyter) | Spatial null calibration, the map-to-gene tests, and how the one positive result was interrogated. |

## Rendering

```bash
quarto render notebooks/01_longitudinal_model.qmd
quarto render notebooks/02_maps_and_genes.qmd
```

Both need a fitted run present in `out/`. `01` reads
`thickness_dsk_51_32bbb98e84e4` (DK, genetics config); `02` reads
`thickness_hcp_51_f408a620fde4` (HCP, absolute change) plus the summary CSVs in
`out/`. Regenerate a run with:

```bash
python -m abcd.assemble configs/ct_genetics.yaml
Rscript R/fit_lmm.R --run-dir out/<run_id> --cores 8
python -m abcd.phenotype out/<run_id>
```

`02` additionally needs the AHBA expression matrix; set `ABCD_AHBA_DIR` to the
directory holding `hcp_3d_ds5.csv` (see `src/abcd/paths.py`).

## Chunk verification without Quarto

Quarto was not installed in the environment these were written in, so the
chunks were verified by extraction and direct execution:

```bash
python tools/check_notebook_chunks.py notebooks/01_longitudinal_model.qmd
python tools/check_notebook_chunks.py notebooks/02_maps_and_genes.qmd
```

That runs every code chunk in order and fails on the first error. It does not
check that the *prose* renders, only that the code does — but code failure is
the usual reason a notebook does not render.
