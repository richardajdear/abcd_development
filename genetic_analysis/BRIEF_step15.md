# Brief for the CSD3 agent: run step 15 (EUR-arm MAGMA on ABCD LD, then set tests)

Repo on CSD3: `/home/rajd2/rds/hpc-work/abcd_development`, branch `main`. Pull
first: `git pull --ff-only origin main` (needs at least commit `f94965c`,
which adds `magma_gene_sets/genesets_wes.txt` and wires it into step 15b).

## Why this is needed
Panel g of Figure 1 (`genetic_analysis/fig1.R`) needs gene-set tests on four
cells per set: {thinning rate, baseline thickness} x {EUR arm, pooled arm}. The
EUR arm so far uses 1000 Genomes EUR (n = 503) as its LD reference, and the
pooled arm uses the ABCD analysis sample. Step 15a re-runs the EUR arm with
ABCD's own EUR children as the LD reference, so the two arms differ only in
sample and not in reference panel. That requires the controlled ABCD
genotypes, which is why it runs here.

## Run
```bash
cd /home/rajd2/rds/hpc-work/abcd_development
git pull --ff-only origin main
ls genetic_analysis/work/results_70tab_hcp/magma_pooled/genes/*_1lmm.genes.raw   # step 14 must be done
jid=$(sbatch --parsable genetic_analysis/step15_magma_eur_abcdld.sbatch)          # 15a, array 1-4, ~1-3 h
sbatch --dependency=afterok:$jid genetic_analysis/step15_magma_set_tests.sbatch     # 15b, ~10 min
```
- Array 15a: {dk, hcp} x {global_slope_1lmm, baseline_thickness_1lmm}. Task 1
  builds `work/magma_scz2025/ref_abcd_eur/abcd_eur.{bed,bim,fam}` (step-14
  `abcd_analysis` restricted to `legacy/hpc/work/results/ancestry/eur_anchor.keep`);
  the other tasks wait on its lock.
- Output 15a: `work/results_70tab{,_hcp}/magma_panel_abcdld/genes/<pheno>.genes.{raw,out}`
- Output 15b: `work/results_70tab{,_hcp}/magma_set_tests/table_magma_set_tests.tsv`,
  with columns version (EUR_1000G / EUR_ABCD / pooled_ABCD), phenotype,
  variable, ngenes, beta, se, p.

## Checks before committing
1. `abcd_eur.fam` has ~4,300 rows (the EUR arm is n = 4,308 in the PRS grid).
2. Each `.genes.out` has ~18-19k genes; any task log with `FATAL` means the job failed.
3. Sanity: the EUR_1000G rows of `table_magma_set_tests.tsv` for
   `SCZ_locus_pool` x phenotype `global_slope` (hcp; the table drops the `_1lmm` suffix) must reproduce p = 0.005 (the
   value already in Figure 1g). The local EUR_1000G run of the exome sets on
   the hcp thinning rate gave SCZ_WES p = 0.026 and MDD_WES p = 0.97; 15b
   should match these exactly.
4. Report EUR_ABCD vs EUR_1000G side by side for SCZ_locus_pool, MDD_highconf,
   SCZ_WES and MDD_WES (hcp, both phenotypes).

## Commit
Commit ONLY the two `table_magma_set_tests.tsv` files plus a short Step 15
result paragraph in `README_HPC.md` (replace the "scripted, not yet run" note).
Never commit `.genes.raw`/`.genes.out` or anything under `ref_abcd_eur/`
(the reference is individual-level ABCD genotype data; the repo is public).
All of these are already gitignored, so `git status` should show only the tables.
Once the tables are pushed, `fig1.R` switches panel g from the local table to
`table_magma_set_tests.tsv` automatically; no code change is needed.
