# Task for the CSD3 agent: GWAS of the area, volume and thickness+area slopes (D1 step 9)

## Why
In ABCD (ages 9–17), cortical thinning is the primary phenotype, but its genetic signal is weak: EUR n 4,308, GREML h² 0.18 and LDSC h² z 0.37. The one gene-level result so far is that genes rising in adolescent oligodendrocytes are enriched in the thinning MAGMA output (gene-property p ≈ 3e-4, `transcriptional_maturation/gwas_cells/48*`).

D1 step 8 (`12_geometry_composite.R`) showed that per-scan composites of thickness with surface area (and GM volume) have more reliable slopes (0.40 vs 0.30). They also have slightly larger PRS effects (mean z over 8 cells 2.86 vs 2.60).

In twins, the area slope is only partly genetically correlated with thinning (rA −0.36 with faster thinning), and the volume slope much more (0.72). A composite GWAS therefore targets a mixed trait.

This task runs the **pre-specified decomposition**: the same single-LMM GWAS as `global_slope_1lmm` on four new slope phenotypes, followed by MAGMA. The oligodendrocyte enrichment test is then repeated locally on every output with identical settings. Roles are fixed now, before any result is seen:

| tag | phenotype (per-child BLUP slope) | role |
|---|---|---|
| `slope_thk` | whole-cortex mean thickness, recomputed on the Figure-1 scans | **reproduction control**: must match `global_slope_1lmm` |
| `slope_area` | total HCP-MMP surface area, z over scans | decomposition: area-change genetics, the contrast to thinning |
| `slope_vol` | total HCP-MMP GM volume, z over scans | conservative alternative (closest genetically to thinning) |
| `slope_thk_area` | per-scan z(thk) + z(area), re-standardised | secondary composite |
| `slope_thk_vol_area` | per-scan z(thk) + z(vol) + z(area), re-standardised | secondary composite |

Intercepts (`int_*`) are in the same file. Run them only as a second priority (positive controls, see below).

Expected outcomes:
- If the Oligo-rising enrichment is thinning-specific, it should shrink in the composites and be absent for area.
- If it strengthens in the composites, it reflects shared geometry.

Both outcomes are reported. No phenotype replaces thinning on the basis of its p-value.

## Input (individual-level; transferred by Richard, never committed)
- Laptop file: `directions/d1_multimodal/work/gwas_pheno_geom.tsv`, written by `13_export_gwas_pheno.R`.
  - It covers 8,716 children (Figure-1 scan set: 26,946 HCP scans).
  - Model: `y ~ age_c + sex + (1 + age_c | subject) + (1 | site)`, REML, age centre 12.797.
  - Columns: `subject` (`sub-NDARINV…`), `slope_<tag>` and `int_<tag>` for the five tags above. Values are raw BLUPs, not yet z-scored.
- Richard rsyncs it to `$REPO/directions/d1_multimodal/work/gwas_pheno_geom.tsv`. If it is absent, stop and ask.
- Laptop child-level r between the slopes:

  |  | thk | area | vol | thk_area |
  |---|---|---|---|---|
  | area | 0.18 | | | |
  | vol | 0.73 | 0.65 | | |
  | thk_area | 0.93 | 0.51 | 0.82 | |
  | thk_vol_area | 0.90 | 0.58 | 0.92 | 0.98 |

  Check that your merged file reproduces these to 2 dp.

## Environment (existing; do not reinstall)
- Repo: `REPO=/home/rajd2/rds/hpc-work/abcd_development`. Pull first with `git pull --ff-only origin main`.
- GENESIS pipeline: `genetic_analysis/run_all.sh 03 04`, `03_null_model.sbatch`, `04_assoc.sbatch`, `R/03–05`, `collect_1lmm.sbatch`, configured through `genetic_analysis/config.sh` (`PHENO_DIR`, `KIN_DIR`, `NULL_DIR`, `ASSOC_DIR`).
- Follow exactly how the existing single-LMM scans were run. Steps 11–13 of `legacy/genetic_analysis/README_HPC_2026-09-24.md` describe them: `PHENO_DIR=prs_final_1lmm/pheno`, `KIN_DIR=kinship_eur` for the EUR arm, outputs under `results_70tab_hcp/scan_1lmm/`. The pooled arm is the step-14 scan (35780487–8 / 35780507–8).
- Slurm: account `vertes-sl3-cpu`, partition `icelake`. Logs go to `slurm/`. Rule 10: `TMPDIR` on rds.
- MAGMA: as `step13_magma_panel.sbatch` (EUR: 1000G EUR, `ncbi37.window35-10.union3.genes.annot`) and `step14_magma_pooled_genes.sbatch` (pooled: in-sample ABCD reference).

## Work layout
- Root: `W=$REPO/genetic_analysis/work/results_70tab_hcp/scan_geom/`.
- Subfolders: `pheno/`, `nullmodel{,_eur}/`, `assoc{,_eur}/`, `reml_imp_pooled/`, `magma_eur/genes/`, `magma_pooled/genes/`, `README.md`.
- Commit the scripts to `genetic_analysis/geom_gwas/` (for example `00_build_pheno.py`, `run_geom.sh`, plus thin wrappers if needed).
- **Do not write anything into the existing `scan_1lmm/`, `magma_panel/` or `magma_pooled/` folders.**

## Steps
1. **Phenotype files** (`00_build_pheno.py`, modelled on `orderops/build_1lmm_pheno.py`).
   - Read the HCP `prs_final_1lmm/pheno/phenotypes_gcta.txt` used for `results_70tab_hcp/scan_1lmm`, i.e. the reference FID/IID set and order.
   - Match on the 8-character NDAR suffix.
   - Z-score each column over the matched set.
   - Write `pheno/phenotypes_gcta.txt` (FID IID + ten columns) and `pheno/phenotype_manifest.tsv` (priority 1 = the five slopes; 0 = intercepts).
   - Copy or symlink `covar_quant.txt` / `covar_categorical.txt` from the same `prs_final_1lmm/pheno` dir, with identical covariates (sex, site, baseline age, n_visits, 10 PC-AiR PCs).
   - Report: n matched of the reference set, n missing per column, and r(`slope_thk`, `global_slope_1lmm`) and r(`int_thk`, `baseline_thickness_1lmm`).
   - **Gate A:** r(slope) ≥ 0.99. The two fits use the same model, so it should be about 1; small differences can come from the scan set. If below 0.99, stop and report.
2. **EUR arm GWAS** (the primary arm; MAGMA's 1000G EUR panel matches it).
   - Run null model + 22-chromosome scan for the five slopes, with `KIN_DIR=kinship_eur` and the same sample set as `scan_1lmm/assoc_eur`. Then collect (`05_collect_assoc.R --maf 0.01 --mac 20`).
   - **Gate B:** for `slope_thk`, per-SNP Pearson r of BETA/SE z against `scan_1lmm/assoc_eur/global_slope_1lmm.sumstats.tsv.gz` ≥ 0.98, and λGC within 0.01 of 1.004. If not, stop and report.
3. **Pooled arm GWAS** for the same five slopes (n ≈ 8,596; the `scan_1lmm/assoc` configuration).
4. **GREML h²** for each slope (`step5_reml.sbatch` settings; same unrelated set as `scan_1lmm/reml_imp_pooled`). Also run bivariate rg of each with `slope_thk` if the existing script supports it; skip it otherwise.
5. **MAGMA gene analysis.**
   - EUR arm: 1000G EUR + 35/10 union annotation, as step 13 (default `snp-wise=mean`; P and N columns built the same way).
   - Pooled arm: in-sample reference + annotation, as step 14.
   - `slope_thk` acts as the reproduction check: gene-ZSTAT r against `magma_panel/genes/global_slope_1lmm.genes.out` (EUR) and `magma_pooled/genes/global_slope_1lmm.genes.out` (pooled). Expect ≥ 0.98.
6. **Second priority, only after 1–5 pass:** the five intercepts, EUR arm, GWAS + MAGMA.
   - These are positive controls: baseline area should be clearly heritable and show λGC > 1.
7. **Optional:** LDSC h² for the EUR slopes, as `step12_ldsc_panel.sbatch`. Report only; do not add to any panel.

Do **not** run the gene-set or gene-property tests here. They run locally with `gwas_cells/48*` so that every phenotype gets byte-identical covariates and sets.

## Report (`$W/README.md`, committed with the scripts)
- Gates A and B (values).
- Per phenotype × arm:
  - n, λGC, min P, and the number of SNPs at P < 5e-8 and < 1e-6;
  - GREML h² ± se;
  - rg with `slope_thk` if run;
  - MAGMA: genes tested, λ of gene P, top 10 genes, and gene-ZSTAT Spearman with `slope_thk`.
- Core-hours and wall time per stage.
- The exact rsync command for the deliverables.

## Deliverables (rsync back; gene-level summaries only)
- `magma_eur/genes/<pheno>.genes.{out,raw,log}` and `magma_pooled/genes/<pheno>.genes.{out,raw,log}`.
- `assoc{,_eur}/gwas_summary.tsv` and `reml_imp_pooled/reml_summary.tsv`.
- Local destination: `~/Git/abcd_development/genetic_analysis/work/results_70tab_hcp/scan_geom/` (same layout).
- The SNP-level sumstats stay on CSD3 unless asked.

## Do not
- Do not overwrite or rerun any existing scan, MAGMA, LDSC, GREML or PRS output.
- Do not change covariates, kinship, sample sets, MAF/MAC filters, LD references, annotations or the gene model.
- Do not commit anything under `work/` (it holds individual-level phenotypes and genotype-derived files; the repo is public). `git status` should show only the scripts and `README.md`.
- Do not rank-normalise the phenotypes (the existing scans don't).

## Local follow-up (laptop, not for the CSD3 agent)
`transcriptional_maturation/gwas_cells/48b` will be rerun on each `genes.raw`, with the same covariates and sets as for thinning.
- Primary readout: the Oligo 4–25 y excess-slope gene-property β and p for thk vs area vs composites, EUR arm, then pooled.
- Interpretation is fixed in advance:
  - thinning-specific if the area estimate is null and the composite estimate is ≤ the thinning estimate;
  - shared geometry if the composite is clearly larger and area is also positive.
