# legacy/ — superseded work, kept for the record

Nothing in here is needed to understand or run the project. Everything a new
agent needs is in the top-level [`README.md`](../README.md), in `docs/`, and in
[`genetic_analysis/README_HPC.md`](../genetic_analysis/README_HPC.md), which
states the current results and what has been tried.

**Why it is legacy.** Every genetics result in these trees (and every imaging
result dated before 2026-09-14) was computed from the ABCD **6.0** tabulated
tables, which sat in a directory labelled 7.0. The true 7.0 tabulation adds
3,520 six-year scans and about 520 subjects to the longitudinal sample, so the
phenotypes changed and the genetics has to be re-run
([`genetic_analysis/`](../genetic_analysis/)). The FreeSurfer surfaces and the
genotypes were 7.0 throughout; only the tabulated covariates and QC were stale.

| directory | what it was | still useful for |
|:--|:--|:--|
| `hpc/` | v1 cluster pipeline: GCTA GRM/REML, fastGWA, MAGMA, LDSC, C+T PRS (release-4.0 EUR genotypes, then 7.0 imputed). `README_HPC.md` §8 is the log of every pitfall. | the MAGMA / LDSC / PRS-association *scripts* that `genetic_analysis/` reuses; the disorder-side MAGMA gene results in `work/results/magma/` |
| `hpc_v2/` | v2 pipeline: GENESIS (KING → PC-AiR → PC-Relate) mixed-model GWAS keeping relatives, Zaitlen REML, within-family PRS, the four-method PRS grid (C+T, PRS-CS, SBayesR, SBayesRC). `README_HPC.md` §14 holds the last canonical PRS tables (`work/results_v2/prs_final/table_main.tsv`). | kinship files and PRS score profiles on CSD3 are phenotype-independent and are reused verbatim |
| `hpc_v3/` | Experiment A: region-subset phenotypes (top-ΔCT, top-C3, projections). Negative on every readout. `SETUP_CONTEXT.md` explains why MOSTest/JAGWAS were rejected and MTAG/genomic SEM were gated. | `align_export_v3.py` (the three silent export mismatches), `prs_paired_delta.py` |
| `genetic_analysis/` | superseded parts of the live `genetic_analysis/` arm, moved 2026-09-24: the full dated run log `README_HPC_2026-09-24.md` (every job id, failure and diagnosis of the 7.0 re-run, steps 1–14); the 6.0-vs-7.0 comparison tables and their builder (`tables/summary_70tab_{dk,hcp}.tsv`, `build_summary_70tab.py`); the PGC3-era PRS slides (`slide_prs_*.py`); v2/v3 setup scripts no longer called (`setup/`). Script paths inside assume the old location | the run log, for derivations and diagnoses |
| `handoff/` | cross-session hand-off tables from the 5.1 → 7.0(6.0) comparison | none; superseded by `docs/vintage_comparison.csv` |
| `hpc/tools/`, `hpc_v2/tools/` | local synthetic-fixture tests of the two pipelines (`local_test*.sh`, `make_test_genotypes*.py`); their relative paths assume the old layout | pattern for a fixture-based test of `genetic_analysis/` |

The `work/` trees are gitignored (per-subject data, genotypes, sumstats,
environments). On CSD3 they still sit at the **old** paths
(`~/rds/hpc-work/abcd_development/hpc*/work/`) because git does not move
untracked files; the migration snippet is in `legacy/genetic_analysis/README_HPC_2026-09-24.md` §2.2.
