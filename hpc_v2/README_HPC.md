# hpc_v2 — sibling-aware genetics pipeline

**Written 2026-08-23, for an agent starting fresh on CSD3 with no access to the
conversation that produced it.** Everything you need is in this file; read
`../hpc/README_HPC.md` only where this one points you at it.

`hpc_v2` does not replace `hpc/`. v1's results stay reproducible and nothing
here writes into v1's output tree (`$OUT_V2` is separate). v2 exists to fix one
design flaw and add one analysis v1 could not do.

---

## 1. Why this pipeline exists, in one paragraph

**ABCD is enriched for siblings, and v1 threw that away.** v1's heritability
step (`hpc/02_reml.sbatch`) ran GCTA REML on a relatedness-pruned subset —
5,649 of 8,082 phenotyped-and-genotyped children — because a single-GRM REML
with relatives in it loads shared environment onto the additive term. The
association step (`hpc/03_gwas.sbatch`, fastGWA) *did* keep relatives, so the
GWAS was not losing N; but its sparse GRM and its ancestry PCs both came from a
**pooled multi-ancestry GCTA GRM**, whose off-diagonals confound ancestry with
kinship. That confounding is not hypothetical: `--grm-cutoff 0.05` on that GRM
returned a **94 %-European "unrelated" set** while reporting nothing wrong
(`../hpc/README_HPC.md` §8.6 — read it, it is the single most important finding
in v1). v2 replaces the relatedness machinery with estimators that separate
kinship from ancestry, keeps every sibling in every step, and adds the one test
the sibling design uniquely enables.

## 2. What v2 does that v1 did not

| | v1 | v2 |
|---|---|---|
| relatedness estimator | GCTA pooled GRM | **KING-robust → PC-AiR → PC-Relate** (GENESIS) |
| ancestry PCs | GCTA `--pca` on the pooled GRM (family structure leaks into the axes) | **PC-AiR** (PCs from unrelateds, relatives projected in) |
| unrelated set | `--grm-cutoff 0.05` → 94 % EUR | **not needed** — no step requires one |
| association | fastGWA-MLM, sparse GCTA GRM | **GENESIS `assocTestSingle`**, sparse PC-Relate kinship |
| heritability | single-GRM REML, relatives **dropped** | **Zaitlen two-GRM REML**, relatives **kept**; reports SNP h² *and* pedigree h² |
| PRS ↔ disorder | population association only | **+ within-family (Fulker) decomposition** |

**What v2 does NOT change:** the phenotypes (five standardised BLUPs from the
same `python -m abcd.gcta_export`), the genotypes (7.0 array for kinship, the
v1-converted TOPMed imputed set for association), and the PRS scores
themselves (v2 consumes v1's `score_*.profile`; it does not re-score).

**Three honest expectations, so nobody over-sells this:**

1. **Keeping relatives in REML improves h²'s SE less than raw N suggests** —
   relatives carry partly redundant genotype information. The genuinely new
   quantity is **pedigree h²** from the same fit, which bounds how much of the
   family-design signal the SNPs tag.
2. **The pooled multi-ancestry GWAS is still a compromise**, just a smaller
   one. PC-Relate removes ancestry from the *kinship*; it does not give the
   score test ancestry-specific allele frequencies. Everything in
   `../hpc/README_HPC.md` §8.16.2 about **European LD references** (LDSC, MAGMA,
   PRS clumping) still applies downstream and is not fixed here.
3. **The within-family PRS test is underpowered for confirmation** at ABCD's
   pair count. Its interpretable output is `beta_diff` (between minus within),
   which *bounds the confounding* in the population estimate. A null
   `beta_within` is **not** evidence against the effect. The script prints this
   itself; do not remove that text.

## 3. Layout

```
hpc_v2/
  config.sh                  every path/parameter; sourced by every script
  config.local.sh.example    copy to config.local.sh (gitignored) and edit
  00_check_inputs.sh         preflight — run it first, it is seconds
  01_gds.sbatch              PLINK → GDS  (array + 22 imputed chromosomes)
  02_kinship.sbatch          KING → PC-AiR → PC-Relate
  03_null_model.sbatch       GENESIS fitNullModel, one per phenotype
  04_assoc.sbatch            assocTestSingle, one (phenotype × chromosome)
  05_reml_zaitlen.sbatch     GCTA two-GRM REML (independent of 01-04)
  06_prs_family.sbatch       within-family PRS (independent of 01-04)
  run_all.sh                 preflight + submit/run with dependencies
  R/                         the statistical work: 01_make_gds, 02_kinship,
                             03_null_model, 04_assoc, 05_collect_assoc,
                             06_prs_family
  envs/genesis_env.yml       conda spec for the R stack
```

Outputs land under `$OUT_V2` (default `$ABCD_HPC_ROOT/results_v2`):
`gds/ kinship/ nullmodel/ assoc/ reml_zaitlen/ prs_family/ logs/`.

## 4. Setup on CSD3 — do these three things first

**(a) Build the R environment.** GENESIS is not on `PATH` and CSD3's system R
lacks it (v1 lesson: a bare `python` on CSD3 is 3.7.4 and nothing useful is on
`PATH`; the same holds for R).

```bash
mamba env create -p ~/rds/hpc-work/envs/genesis -f hpc_v2/envs/genesis_env.yml
```

**Note:** `bioconductor-genesis` failed to solve on conda-forge/bioconda in
August 2026 (`bioconductor-snprelate` was absent from the channel). If the yml
fails, install from Bioconductor inside a plain R:

```r
install.packages(c("gdsfmt","SNPRelate","GWASTools","GENESIS"),
  repos = c(BioCsoft = "https://bioconductor.org/packages/release/bioc",
            CRAN = "https://cloud.r-project.org"), type = "source")
```

**(b) Write `hpc_v2/config.local.sh`** from the example. It is gitignored —
never commit it. The variables that matter:

| variable | what it must point at |
|---|---|
| `GENO_ARRAY` | the **merged** 7.0 Smokescreen array fileset (11,670 × 515,228, hg19) |
| `GENO_IMP_DIR` / `GENO_IMP_TPL` | v1's converted per-chromosome imputed filesets (rsIDs restored) |
| `GRM_FULL` | v1's dense imputed GRM over the **full** sample — **not** a `.unrel` one |
| `PRS_PROFILE_DIR` | v1's `results/prs_imp/` (`score_*.profile`) |
| `EUR_IDS` | the EUR keep list behind v1's EUR-primary PRS arm |
| `RSCRIPT` | `~/rds/hpc-work/envs/genesis/bin/Rscript` |

**(c) Run the preflight.**

```bash
bash hpc_v2/00_check_inputs.sh
```

It checks `.bed`/`.bim`/`.fam` size consistency (the 4.0 defect detector), that
phenotype IDs intersect the genotype `.fam` on the 8-char NDAR token, that the
phenotype export still carries **real family IDs** (step 06 is impossible
without them), that the manifest's `--mpheno` indices are in range, and that
`RSCRIPT` can load GENESIS. It exits non-zero on anything fatal, and
`run_all.sh` refuses to submit if it fails.

## 5. Running it

```bash
bash hpc_v2/run_all.sh                # everything, with dependencies
bash hpc_v2/run_all.sh 01 02 03 04    # the GENESIS chain only
bash hpc_v2/run_all.sh 05             # Zaitlen REML only
DRY_RUN=1 bash hpc_v2/run_all.sh      # print every command, execute nothing
```

Dependency shape — **05 and 06 do not wait for the GENESIS chain**, they depend
only on the GRM/scores and the export, so `run_all.sh` submits them
immediately and they run in parallel:

```
01 gds ──► 02 kinship ──► 03 null models ──► 04 assoc ──► (final collect)
05 reml_zaitlen        (independent)
06 prs_family          (independent)
```

Array sizing comes from the manifest, never from a hand-edited `--array`
range: `03` is `1-n_phenotypes`, `04` is `1-(n_phenotypes × 22)` with
`task = (pheno_idx − 1) × 22 + chr`, `05` is `1-n_phenotypes`.

**Expected cost** (11,670 subjects; extrapolated from v1's timings, **not yet
measured for v2** — record the real numbers when you run it):

| step | shape | guess |
|---|---|---|
| 01 gds | 23 tasks | minutes each; I/O-bound |
| 02 kinship | 1 job, 8 CPU, 64 G | a few hours (KING + PC-Relate over ~100k pruned SNPs) |
| 03 null model | 5 tasks | minutes each |
| 04 assoc | 110 tasks | tens of minutes each |
| 05 reml | 5 tasks | as v1's REML (~minutes at this N) |
| 06 prs_family | 1 job | seconds |

## 6. Reading the outputs

| file | what to look at |
|---|---|
| `kinship/kinship_summary.tsv` | `n_pcair_unrelated` should be in the same range as ABCD's shipped `genesis/unrelateds_individuals.txt` (**8,181** of 11,670). First-degree pair count should be in the thousands — ABCD is sibling-enriched. |
| `nullmodel/<pheno>_null_summary.tsv` | **`n` must exceed v1's 5,649.** If it equals it, relatives were dropped and the whole point was lost. |
| `assoc/gwas_summary.tsv` | same columns as v1's (`phenotype n_snps n_mean lambda_gc n_p5e8 n_p1e5 min_p`) so the two GWAS are directly comparable. λ_GC near 1.00; above ~1.10 means residual structure — check the PCs came from `pcair_pcs.tsv`, not the release PCs. |
| `assoc/<pheno>.sumstats.tsv.gz` | `SNP CHR POS A1 A2 AF1 N BETA SE P`, rsID-keyed, ready for MAGMA/LDSC/PRS. |
| `reml_zaitlen/reml_zaitlen_summary.tsv` | `h2_snp` (compare with v1's `reml_summary.tsv`) and `h2_ped` (new). `n` must be the **full** sample. |
| `prs_family/prs_withinfamily.tsv` | `beta_within` / `beta_between` / `beta_diff` with `n_pairs`. Read §2.3 above before interpreting. |

**Verify outputs, not `sacct`.** A v1 run had 13 jobs reporting `COMPLETED`
while producing nothing. Every step here writes a `*_summary.tsv` and the
sbatch wrappers assert their own outputs exist and are non-empty before
exiting 0 — but check the row counts yourself.

## 7. Gotchas — all of these have already bitten this project once

- **The "sparse" kinship matrix can silently come out dense — check it.**
  `pcrelateToMatrix(thresh=)` does **not** zero individual entries; it
  **clusters samples by transitive closure** (any pair above the threshold
  joins a cluster, all within-cluster pairs are then kept, and only
  *between*-cluster entries become zero). A chain of noise-level pairs can
  therefore merge the whole cohort into one block and hand `fitNullModel` a
  dense matrix — **~1.1 GB at n = 11,670**, with no error anywhere. Measured on
  the local fixture (1,500 subjects, only 4,977 pruned SNPs): 783 true sibling
  pairs, but 2,545 pairs cleared the threshold and their closure pulled 1,472
  of 1,500 subjects into a single block, giving **96 % density**. With ~100k+
  pruned SNPs on real data the noise is far smaller and this should not happen,
  but do not assume it: `kinship_summary.tsv` reports `sparse_density` and
  `sparse_largest_block`, and `02_kinship.R` prints a loud warning above
  `SPARSE_KIN_MAX_DENSITY` (default 0.10) explaining the fix. If it fires,
  raise `SPARSE_KIN_THRESH` (0.0884 keeps only 2nd-degree-and-closer pairs) and
  check `n_pruned_snps` first — too few pruned SNPs is the usual cause.
  **Note the factor of 2:** `thresh` is compared against `2 × kinship`, so the
  0.0442 default corresponds to raw kinship ~0.0221.
- **IDs.** ABCD spells subject IDs three ways (`NDAR_INVxxxxxxxx`,
  `sub-NDARINVxxxxxxxx`, bare `sub-xxxxxxxx`). Join on the 8-char NDAR token;
  `R/03_null_model.R` and `R/06_prs_family.R` both normalise with the same
  `norm_id()`. A mismatched ID convention gives a **zero-subject analysis with
  no error**.
- **`FID = IID` in every genotype `.fam` on CSD3.** Family structure is
  invisible to anything reading the `.fam`; it lives in the phenotype export.
  Step 06 refuses to run when all `FID == IID` rather than silently returning
  `b_W = 0`, which would look like a clean null.
- **`plink2 --vcf` writes `FID = 0`** unless `--double-id` is given
  (`../hpc/README_HPC.md` §8.14). If you regenerate any fileset, re-check the
  ID convention against a keep-list *before* trusting a step with no `--keep`
  to fail for you.
- **Name the partition.** Five v1 scripts had no `#SBATCH --partition`, took
  the cluster default `cclake` while it was in maintenance, and sat PENDING
  overnight showing only `(Priority)`. Every script here declares
  `--partition=icelake`; `run_all.sh` passes `$SLURM_PARTITION` too. Run
  `squeue -o '%P'` before concluding a job is merely queued.
- **Size jobs for backfill.** Dropping v1's GRM array from 16 CPU × 3 h to
  8 CPU × 1 h 45 moved its start from 3.8 hours away to minutes. Queue time
  dominates compute time on this pipeline.
- **`/home` is nearly full.** Everything goes under `/rds/user/<user>`;
  `ABCD_HPC_ROOT` relocates the whole tree.
- **Bash brace trap.** `${GENO_IMP_TPL:-abcd_imp_chr{CHR}}` does **not** work —
  bash ends the expansion at the `}` closing `{CHR}`, appending `}` as literal
  text. `config.sh` uses an `if [[ -z ... ]]` assignment; keep `{CHR}`-style
  placeholders out of `${...:-...}` defaults.
- **bash 3.2.** No `mapfile`, no `declare -A` (macOS ships 3.2 and the local
  test runs there). Use `while read` and `case`.
- **No absolute paths in scripts.** Every path is a `config.sh` variable. If
  you add a literal, you have broken portability for the next account.

## 8. Local validation before you submit anything

```bash
python tools/make_test_genotypes.py       # v1 fixture (source genotypes)
python tools/make_test_genotypes_v2.py    # v2 fixture
cp hpc_v2/config.local.sh.example hpc_v2/config.local.sh
bash tools/local_test_v2.sh
```

Unlike v1's harness — where GCTA's MKL abort on Apple Silicon meant only option
parsing could be checked — **most of v2 executes for real** locally, because it
is R: GDS conversion, KING, PC-AiR, PC-Relate, `fitNullModel`,
`assocTestSingle`, the collector and the between/within decomposition all
produce genuine numbers on a 1,500-subject / 5,000-SNP fixture with 739
multi-member families.

The strongest check is a **planted positive control**: the fixture writes a
`global_slope` phenotype with a known within-family PRS effect
(`b_W = −0.20`) plus family-level confounding (family-mean coefficient
`−0.35`). Step 06 must recover `beta_within ≈ −0.20`, must return a
`beta_between` more negative than it, and must flag `beta_diff` as significant.
A script that confuses the two coefficients, loses the family IDs, or
standardises at the wrong point cannot pass. The harness also asserts the
FID=IID guard fires, that the integrity gate rejects a `.bed`/`.bim` mismatch,
that the collector is idempotent, and that the null model's `n` **exceeds** the
unrelated count — the v2-specific claim.

GCTA (step 05) remains the exception: on Apple Silicon it aborts inside MKL, so
the harness validates that its options parse and that it opens the GRM, plus
the two-GRM `.hsq` collector against a synthetic file in GCTA's exact format
(including a negative case: a single-GRM `.hsq` must yield no row rather than
be mislabelled as Zaitlen).

**Result as committed: 73 checks, 0 failures**
([`../docs/hpc_v2_local_test.log`](../docs/hpc_v2_local_test.log)). The
load-bearing ones, with their measured values on the fixture:

| check | measured |
|---|---|
| PC-Relate recovers the true pair count | 783 first-degree pairs vs **783** implied by the `.fam` FIDs |
| null model keeps relatives | fitted on **n = 1,500** vs 642 PC-AiR unrelated |
| association uses relatives | per-variant N up to **1,500** |
| model is calibrated | λ_GC **0.994** (`global_slope`), 1.011 (`baseline_thickness`) |
| seeded-phenotype control | `sim_h2_50` (20 causal SNPs, h² = 0.5) → **8 hits, min p = 2.2e-23**; unseeded `global_slope` → 0 hits |
| within-family PRS positive control | β_W = **−0.163** (planted −0.20, p = 1.2e-03); β_B = **−0.329** (planted −0.35); p_diff = 4.4e-03; MDD null score p = 0.76 |
| sparsity guard | density **0.963** > 0.10 → warning **fired** (expected on a 5k-SNP fixture; see §7) |

The seeded-phenotype row is worth understanding: without it, a pipeline that
shuffled genotypes against phenotypes — an ID-join bug, the most expensive class
of error in this project — would still produce a clean λ_GC and pass every other
check.

**Step 05's numerics are not validated anywhere yet.** GCTA does not run on the
development machine, so the two-GRM REML has only been checked for option
validity and output parsing. **The first cluster run is its first real
execution** — treat `h2_ped ≥ h2_snp` and an `n` equal to the full
phenotyped-genotyped sample as the two things to verify before believing it.

## 9. What is NOT in v2, and why

- **No MAGMA / LDSC / PRS scoring.** They are unchanged from v1 and consume
  sumstats; `assoc/<pheno>.sumstats.tsv.gz` is written in the layout v1's
  `04_magma` / `05_ldsc_rg` / `06_prs` already read. Run them from `hpc/`
  against v2 sumstats rather than duplicating them here. The European
  LD-reference problem (`../hpc/README_HPC.md` §8.16.2, §8.16.6) is untouched
  by v2 and still governs how those results may be read.
- **No within-family GWAS.** Sib-difference association is not powered at ABCD
  scale — the effective N is the pair count, and the within-pair variance is
  roughly half the total. It is a robustness check for specific hits, not a
  discovery arm. If you want it for the top loci from step 04, add it then.
- **No snipar / direct-vs-indirect effect decomposition.** Realistic at PRS
  level, not SNP level, at this N.
- **No ancestry-stratified GWAS meta-analysis.** The non-EUR strata are
  individually small (v1 §8.9d: SEs of 0.23–0.45 span the parameter space), and
  the pooled GENESIS design is the reason v2 exists. `STRATA_FILE` optionally
  enables heteroscedastic residual variances by stratum inside the *pooled*
  model, which is the middle ground.

## 10. How to document what you do

Same contract as `../hpc/README_HPC.md` §7, which earned its place:

**Update this file in place**, with a dated section for what you ran. Do not
start a third README; if this grows unwieldy, move superseded parts to
`hpc_v2/legacy/` with a datestamp.

**For every result record:** the job ID and script, the **N** it ran on (this
project's history is one long argument about N, and a result without its N
cannot be compared to anything), the file it landed in, and whether it
**replaces or supplements** an existing number — if it replaces one, say which
and keep the old value visible.

**Record failures and dead ends with the diagnosis.** Several v1 defects were
*silent*, producing plausible numbers rather than errors. If something looked
right and was not, that is the highest-value thing you can write down.

**Distinguish verified from assumed.** State which claims you checked against a
file. The most expensive class of error here has been a confident statement
nobody re-derived.

**Flag anything that changes a conclusion** at the top of your section — a
material h² disagreement with v1, or the within-family PRS estimate diverging
from the population one, is a headline, not a table row.
