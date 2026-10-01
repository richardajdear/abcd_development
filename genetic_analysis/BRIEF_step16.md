# Brief for the CSD3 agent: step 16, multivariate (MOSTest) discovery on regional thinning

Repo on CSD3: `/home/rajd2/rds/hpc-work/abcd_development`, branch `main`.
Pull first: `git pull --ff-only origin main`. Everything is in
`genetic_analysis/mostest/`; README_HPC.md §9 summarises it. Rules 1–19 of
README_HPC §6 apply. Rule 16 matters most here: `pheno/`, `geno/`, `regenie/`,
`zmat/`, `sumstats/` and `magma/` are per-subject or SNP-level and gitignored.
Commit only `table_mostest_*.tsv`.

## Why this is run

Every ABCD GWAS so far scans one scalar trait per run: global slope, baseline
CT, or the slope PCs. The 68 DK regional slopes have an effective
dimensionality of about 50, and global slope, their mean, discards most of
it. MOSTest (van der Meer et al. 2020, Nat Commun 11:3512) combines the 68
per-region z-scores per SNP into z′R⁻¹z. It detects variants with distributed,
mixed-sign regional effects, which cancel in the mean and never clear the bar
in any single region. Caseras et al. 2026 (Nat Commun, s41467-026-77476-2)
used MOSTest followed by MAGMA on resting-state networks in 30,820 UK Biobank
adults (114 independent SNPs, 315 genes).

**This is a discovery arm with written-down expectations, not a heritability
optimiser.** This was settled in the 2026-09-08 assessment (`legacy/hpc_v3/SETUP_CONTEXT.md`):
- **Unsigned output.** MOSTest gives one unsigned p per SNP. No LDSC rg, no
  LAVA, no h² "of the MOSTest phenotype". The downstream is MAGMA, which needs
  only unsigned SNP p, and gene-property tests.
- **Expected power is low.** Effective N per regional slope is about
  8,596 × reliability 0.2 ≈ 1,700, versus about 31k in Caseras et al. Few or
  no genome-wide loci are expected for the slopes. The deliverable that
  matters is the MAGMA gene-level z (`genes.raw`), which feeds laptop
  gene-property tests against developmental snRNA-seq.
- **What will dominate.** The most heritable directions of slope variation are
  the global factor and the C1/myelin axis, so they are likely to dominate
  any signal.

## Design (why each choice)

| choice | value | why |
|:--|:--|:--|
| measures | DK 68 regions × 3 families: `slope` (per-region LMM slope BLUP), `ct` (per-region intercept BLUP, regional CT at age 12.8), `slopeols` (per-child OLS slope) | `ct` is the **positive control** (cross-sectional, reliable, whole-cortex LDSC h² z ≈ 4.4). `slopeols` is the unshrunk sensitivity for the region-dependent BLUP shrinkage caveat. HCP-MMP (358) comes only after the gates (see below). |
| sample | the 8,596-child pooled analysis set, relatives kept | same set as every pooled result; ~1,300 sibling pairs would be lost by pruning to unrelated |
| univariate engine | **REGENIE** (LOCO ridge; handles relatedness), `--apply-rint` | GENESIS (steps 3–4) needs one null model and 22 tasks per phenotype, i.e. about 10k tasks here. REGENIE scans all 68 measures of a run in one pass. RINT follows MOSTest. |
| covariates | sex, site (categorical), baseline_age, n_visits, PC1–10 | identical to the GENESIS null model |
| MOSTest null | one permutation (phenotype + sex, site, age, n_visits rows shuffled jointly; PC1–10 stay with the genotypes; seed 16), full REGENIE re-run | this is MOSTest's genotype permutation. R = corr(permuted z). The p for the MOSTest statistic comes from a gamma fitted to the permuted statistic. MinP is calibrated on the permuted minP. Code: `mostest_core.py`, unit-tested in `tests/test_mostest_core.py`. |
| regularisation | eigenvalue floor 1e-4 × max | DK R is well conditioned: the laptop phenotype-correlation spectrum has condition number 120 (slope) and 179 (ct), so the floor should not bite (`n_floored` = 0) |
| step-2 genotypes | `genetic_analysis/work/inputs/geno/abcd_imp_prs` (7.07 M SNPs), MAF ≥ 1 %, MAC ≥ 20 | the same fileset as the step-14 MAGMA LD reference, so SNP IDs join |
| loci / genes | PLINK clump (5e-8 and 1e-6, r² 0.1, 500 kb) and MAGMA 35/10 kb on `ref_abcd/abcd_analysis` | rule 19: the LD reference is the GWAS sample itself |

## Setup (once)

```bash
cd /home/rajd2/rds/hpc-work/abcd_development && git pull --ff-only origin main
conda env create -p ~/rds/hpc-work/envs/regenie -f genetic_analysis/envs/regenie_env.yml
~/rds/hpc-work/envs/regenie/bin/regenie --version
```
If conda cannot solve, use the static binary (instructions in the yml) and
set `REGENIE=` in `genetic_analysis/config.local.sh`. Check that the inputs
named in `mostest/paths.sh` exist. In particular:
- `REF_1LMM` = `work/results_70tab/prs_final_1lmm/pheno/` (DK single-LMM
  export; FID = family id). If it is not there, point it at the DK
  single-LMM export the step-14 pooled scan was built from.
- `GENESIS_1LMM` = `work/results_70tab/scan_1lmm/assoc/` (pooled DK
  `global_slope_1lmm.sumstats.tsv.gz`, which step 14 read).
- `DK_RUN` = `out/thickness_dsk_70_139406217085/` (`fits/blups.parquet`,
  `model_table.parquet`).
- `magma_scz2025/ref_abcd/abcd_analysis.{bed,bim,fam}` and its `.genes.annot` (step 14).

Also confirm that the Python in `paths.sh` imports numpy, scipy, pandas and
pyarrow. Before submitting, check that `pd.read_parquet` works in that env.

## Run

```bash
bash genetic_analysis/mostest/run_step16.sh          # submits stages 1-6 with afterok
FROM=4 bash genetic_analysis/mostest/run_step16.sh   # resume at a stage after a fix
```

| stage | script | array | rough cost | output (under `work/results_70tab/mostest/`) |
|:--|:--|:--|:--|:--|
| 1 build | `01_build.sbatch` → `01_build_pheno.py` + PLINK array QC | 1 | 30 min | `pheno/`, `geno/step1_array.*` |
| 2 REGENIE step 1 | `02_regenie_step1.sbatch` | 7 runs | 1–6 h each, 16 cpu | `regenie/step1_<run>_pred.list` |
| 3 REGENIE step 2 | `03_regenie_step2.sbatch` → `02_to_zmat.py` | 154 = 7 × 22 | 0.5–2 h each, 8 cpu | `zmat/<run>_chr<c>.npz` |
| 4 MOSTest + G1 | `04_mostest.sbatch` → `03_mostest.py`, `04_engine_check.py` | 4 | < 1 h, 96 G | `sumstats/mostest_<fam>{,_perm}.sumstats.tsv.gz` |
| 5 loci + genes | `05_clump_magma.sbatch` | 6 (3 real + 3 perm) | 2–12 h each | `clump/`, `magma/genes/mostest_<scan>.genes.{raw,out}`, `magma/sets/` |
| 6 collect | `06_collect.sbatch` → `06_collect.py` | 1 | minutes | `table_mostest_*.tsv` (tracked) |

The 7 runs are `slope slope_perm ct ct_perm slopeols slopeols_perm global`.
`global` puts `global_slope_1lmm` and `baseline_thickness_1lmm` through
REGENIE only for the engine check. The total is roughly 3–5k CPU-hours on
SL3; these are estimates, so check `sacct` after stage 2 and re-size if
needed. Disk: about 10 GB of REGENIE text per real run (kept as the
per-region univariate sumstats) and about 2 GB of z per run.

## Gates (in `table_mostest_gates.tsv`; read NO result until G0–G3 pass)

| gate | test | threshold | if it fails |
|:--|:--|:--|:--|
| G0 | analysis set in both genotype sets | n ≥ 8,500 (expect 8,596); regional-mean slope vs `global_slope_1lmm` r ≥ 0.7 (laptop: 0.91 against HCP single LMM) | ID join (rule 1) or wrong run dir |
| G1 | REGENIE `global` vs the GENESIS pooled scan | \|r(z)\| ≥ 0.90 and \|Δλ\| ≤ 0.05, both traits | covariates, IDs or RINT mismatch; do not proceed |
| G2a | MOSTest λ on the **permuted** scan | 0.95–1.05 per family | null fit or R estimate broken |
| G2b | genome-wide loci in the permuted scan | ≤ 1 per family | tail miscalibrated |
| G2c | MAGMA gene λ on the permuted scan | 0.90–1.10 | gene-level p not usable |
| G3 | R from real vs permuted z, max \|Δ\| | ≤ 0.05 | population structure or strong polygenic inflation; inspect `univariate_<fam>.tsv` λ |
| G4 | **positive control**: `ct` genome-wide loci | ≥ 1 | see decision rule 1 |

## Decision rules (fixed before seeing results)

1. **If G4 fails** (`ct` has no genome-wide locus), MOSTest at this n cannot
   find even cross-sectional structure. Report every slope result as "not
   assessed at this power", not as absence. Do not extend to HCP-MMP. Still
   pull the `genes.raw` files (below), because a gene-level enrichment can
   exist without loci.
2. **If G4 passes and `slope` has loci**, report them only together with
   (a) their `slopeols` p, (b) their `ct` p, and (c) the univariate minP.
   A slope locus that is also a `ct` locus is a thickness locus. Also report
   `table_mostest_gene_z_corr.tsv`: gene-level slope vs ct Spearman near the
   slope vs slope_perm floor means the slope signal is its own.
3. **HCP-MMP extension** (358 measures; the eigenvalue floor will then bite):
   only if `slope` shows loci or a gene-level enrichment on the laptop test.
   It needs an HCP copy of the pipeline with `DK_RUN` and `REF_1LMM` changed;
   do not start it without asking.
4. Compare MOSTest against minP in `table_mostest_summary.tsv`. MOSTest
   should beat minP on `ct`. If it does not, the distributed-effect premise
   does not hold for these measures, and that is a finding.

## Report back (in the commit message and README_HPC §9)

- Job ids, n, SNP count, and the gates table.
- `table_mostest_summary.tsv` rows for the three real families: MOSTest and
  minP λ, loci at 5e-8 and 1e-6, gene λ, Bonferroni genes.
- `table_mostest_loci.tsv` and `table_mostest_gene_z_corr.tsv` in full.
- The set tests (`table_mostest_set_tests.tsv`) for `SCZ_locus_pool`,
  `MDD_highconf`, `SCZ_WES` and `MDD_WES` on slope and ct.

## Commit

Commit ONLY `genetic_analysis/work/results_70tab/mostest/table_mostest_*.tsv`
plus a short result paragraph in README_HPC §9. `git status` should show
nothing else under `work/`. Run `git diff --cached --name-only` before
committing, because other sessions leave files staged.

## After the cluster (laptop, not for the CSD3 agent)

The user pulls the gene results himself. No individual-level data is
involved, but they stay gitignored:
```bash
cd ~/Git/abcd_development
rsync -avh login-q-1.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/results_70tab/mostest/magma/genes/ \
  genetic_analysis/work/results_70tab/mostest/magma/genes/
```
The laptop agent then runs MAGMA gene-property tests (`--gene-covar`, seconds
per test) of `mostest_{slope,ct,slopeols}.genes.raw` against developmental
snRNA-seq properties: cell-type specificity from Velmeshev and Herring, and
the maturation PC1 loadings. `ct` and the permuted scans are the comparators.
This is the step that links ABCD genetics to single-cell transcriptomics.
