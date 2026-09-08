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

---

## 11. First cluster run — 2026-08-23 (CSD3, account VERTES-SL3-CPU)

**Status: setup complete, preflight passed, all six steps submitted.** This
section records the setup decisions (two of them resolve ambiguities §4 left
open), one **defect in `run_all.sh` found and fixed here**, and the results as
each step lands. Everything below was checked against a file unless it says
otherwise.

### 11.1 Headline

1. **Steps 01–04 and 06 ran clean and deliver what v2 promised.** The
   association scan keeps every relative (**n_mean 7,893** per variant,
   9,413,281 SNPs) and is better calibrated than v1's: λ_GC for
   `baseline_thickness` falls from **1.107 to 1.036**. Null models fit on
   **n = 8,082**, comfortably above v1's 5,649 — §6's central check.
2. **v1's two genome-wide loci survive but drop just below the line.** chr2
   ~26.93 Mb (rs7599286, **p = 8.2e-08**) and chr12 ~69.38 Mb (rs2603090,
   **p = 9.5e-08**) are still the top two signals under PC-AiR PCs and a
   PC-Relate kinship. They did not evaporate — which is what an ancestry
   artefact would have done — but neither arm should now be called
   genome-wide significant.
3. **Step 05 as shipped is invalid, and not for a numerical reason.** It
   defines its "close relative" component by thresholding the *pooled* GRM,
   which is v1 §8.6's confound. It kept **3,837,371 pairs where PC-Relate finds
   29,371**. Every phenotype degenerated. Details and the fix in 11.8–11.9.
4. **`h2_ped` is not obtainable from this design at all.** The pooled GRM's
   **diagonal correlates −0.990 with PC-AiR PC1** (R² = 0.98), both Zaitlen
   components inherit it, and their sum is therefore an ancestry-heteroscedasticity
   parameter wearing a pedigree-h² label. This is a design conclusion, not a
   re-run: see 11.9.
5. **The within-family PRS test ran and reproduces v1's population estimate as
   its between-family term** — the sharpest available check that the two
   pipelines agree where they should. 11.10.
6. **`run_all.sh` had a defect that made `DRY_RUN=1` submit the pipeline** and
   silently hollow out every job. Fixed here; 11.5.

### 11.2 The environment (§4a)

`hpc_v2/envs/genesis_env.yml` **solved without modification** — the August 2026
note in §4 about `bioconductor-snprelate` being absent no longer applies.

| | |
|---|---|
| built by | `hpc_v2/work/setup/build_genesis_env.sbatch` (job 34276944, 10 min, icelake) |
| prefix | `/rds/user/rajd2/hpc-work/envs/genesis` (274 packages, 407 MB download) |
| versions | R 4.4.3, GENESIS 2.36.0; lme4/lmerTest/data.table/optparse present |
| package cache | `hpc/work/mamba/pkgs` via `MAMBA_ROOT_PREFIX` — **not** `/home`, which is at 48.8 of 52.4 GB |

Built under SLURM rather than on the login node: the solve+download is well
past what a login node tolerates. v1's conda R (`hpc/work/envs/abcdR`) was not
reusable — it has lme4 but no GENESIS/SNPRelate/GWASTools.

### 11.3 Setup decision 1 — which phenotype export, and the GCTA FID clash

§4b says `$PHENO` must carry real family IDs. On CSD3 **two exports exist and
they differ exactly in that column**, because v1's two consumers need different
things:

| export | FID | n | used by v1 for |
|---|---|---|---|
| `hpc/work/pheno_allanc/` | = IID | 8,082 | GCTA (01–03): GCTA matches an individual on the **FID+IID pair**, and every `.fam`/`.grm.id` on CSD3 spells FID = IID |
| `hpc/work/pheno_allanc_prs/` | family_id | 8,082 | `prs_assoc.R`: family structure comes from FID |

v2 has **one** `$PHENO` serving both GENESIS (joins on the 8-char NDAR token,
FID irrelevant) and step 06 (reads family_id straight from FID and refuses to
run when FID == IID). So `$PHENO` is the **family-ID copy** —
and step 05's GCTA is made to agree rather than being left to fail:

**`GRM_FULL` points at an id view, not at v1's GRM.**
`hpc_v2/work/setup/make_grm_famid.sh` writes `hpc_v2/work/grm/abcd_imp_famid`
whose `.grm.bin`/`.grm.N.bin` are **symlinks** to
`hpc/work/results/grm_imp_pooled/abcd_imp` and whose `.grm.id` is the same file
with column 1 replaced by family_id **in the same row order** — the row order is
the matrix order, and reordering it would silently permute the GRM. The script
asserts the row count is unchanged and that the IID column is byte-identical to
the source. 11,670 rows, 8,082 FIDs replaced (the 3,588 genotyped-but-not-
phenotyped subjects keep FID = IID; GCTA drops them either way). **Nothing in
v1's tree is modified.**

**Verified, not assumed:** GCTA read the id view and reported
`8082 individuals are in common in these files` — the full phenotyped-genotyped
sample, not zero and not v1's 5,649 unrelated subset.

The phenotype sample: 8,082 subjects, 6,789 families, **1,250 multi-member
families, 1,339 sibling pairs** (all 8,082 present in the array `.fam`).

### 11.4 Setup decision 2 — which covariates, and which imputed filesets

- `COVAR_QUANT` = `pheno_allanc_prs/covar_quant.txt`, which carries the
  **array-GRM PCs** — the covariate set behind v1's *primary* REML
  (`results/reml_imp_pooled/`), so v2's `h2_snp` is directly comparable with it.
  The imputed-PC sensitivity arm sits beside it as `covar_quant_imp10.txt`.
  The GENESIS steps drop these PC columns and substitute PC-AiR PCs anyway.
- `GENO_IMP_TPL` = `imp_union_chr{CHR}` (**not** the `abcd_imp_chr{CHR}` guess
  in `config.sh`'s default). 22 filesets, 11,670 subjects each, rsID-keyed,
  **25,952,942 variants total**. `genotype_imputed/` also holds a stray
  `imp_pooled_chr22` from an earlier v1 attempt — not used.
- All 23 filesets (array + 22 imputed) pass the `3 + ceil(n/4)*m` bed-size and
  `6c1b01` magic checks.

### 11.5 Defect found and fixed: `DRY_RUN=1` submitted the pipeline

**This is the failure mode §6 warns about, and it was in `run_all.sh` itself.**

`DRY_RUN` was wired only into `config.sh`'s `run()`, which the SLURM branch
never calls — it invoked `sbatch` directly. So `DRY_RUN=1 bash run_all.sh`
**submitted all seven jobs**. Worse, `sbatch` exports the submitting
environment, so every job inherited `DRY_RUN=1` and echoed its commands
instead of running them: the 23-task GDS array reported **COMPLETED in one
second per task having produced nothing**. Jobs 34277587–34277593 were
cancelled; no output had been written.

Fixed in `run_all.sh` with a `submit()` wrapper: under `DRY_RUN=1` it prints
the `sbatch` command to stderr and returns the sentinel job id `DRYRUN`; a real
submission passes `--export=ALL,DRY_RUN=0`, so a stale `DRY_RUN` export in the
submitting shell can never hollow out a real run. Re-tested: `DRY_RUN=1` now
leaves `squeue` empty.

### 11.6 Submission (jobs 34277745–34277818)

| step | job | shape | measured cost |
|---|---|---|---|
| 01 gds | 34277745 | array 0–22 | **25 s – 1 m 34 s per task** (23/23 COMPLETED; 73 GB of GDS, 26,468,170 variants) |
| 02 kinship | 34277746 | 1 job, 8 CPU, 64 G | *(running)* — 98,528 LD-pruned SNPs (maf≥0.05, r²≤0.1) |
| 03 null model | 34277747 | array 1–5 | *(queued, afterok 02)* |
| 04 assoc | 34277748 | array 1–110 | *(queued, afterok 03)* |
| 04c collect | 34277749 | 1 job | *(queued, afterok 04)* |
| 05 reml | 34277750 | array 1–5 | *(running)* |
| 06 prs_family | 34277818 | 1 job | **1 m 11 s, COMPLETED** |

98,528 pruned SNPs is in the range §7 says makes the sparse-kinship transitive-
closure blow-up unlikely (the local fixture's 96 % density came from 4,977).
The density is reported in `kinship_summary.tsv` — check it, do not assume it.

### 11.7 Steps 02–03 — the relatedness machinery works, and the v2 claim holds

**Step 02, kinship (job 34277746, 55 m 53 s, 8 CPU / 64 G).** 98,528 LD-pruned
SNPs (MAF ≥ 0.05, r² ≤ 0.1) from the 515,228-SNP array set.

| quantity | value | read against |
|---|---|---|
| PC-AiR unrelated / related | **8,536 / 3,134** | ABCD's shipped `genesis/unrelateds_individuals.txt` is 8,181 of 11,670 — same range, gate passed |
| dup/MZ pairs (κ ≥ 0.354) | 421 | |
| 1st-degree pairs | 1,504 | sibling-enriched, as §6 expects |
| 2nd-degree pairs | 982 | |
| 3rd-degree pairs | 26,464 | |
| `sparse_density` | **0.00407** | §7's blow-up threshold is 0.10 — **no warning fired** |
| `sparse_largest_block` | 726 (6.2 % of sample) | |

So §7's transitive-closure hazard did not materialise at 98.5k pruned SNPs,
exactly as it predicted. **Verified from `kinship_summary.tsv`, not assumed.**

**Step 03, null models (job 34277747, 1 m 15 s per phenotype).** All five fit
on **n = 8,082** with 14 covariates. §6's test — *"`n` must exceed v1's
5,649"* — **passes**: every relative is in the model.

| phenotype | n | varcomp_kin | varcomp_resid | prop_kin |
|---|---|---|---|---|
| baseline_thickness | 8,082 | 0.8416 | 0.1155 | 0.879 |
| global_slope | 8,082 | 0.3400 | 0.5415 | 0.386 |
| slope_PC1 | 8,082 | 0.3131 | 0.6165 | 0.337 |
| slope_PC2 | 8,082 | 0.3446 | 0.6244 | 0.356 |
| slope_PC3 | 8,082 | 0.3598 | 0.5688 | 0.387 |

`prop_kin` is **not** an h². The covariance matrix here is the *sparse*
PC-Relate kinship (density 0.004, diagonal ≈ 1), so this component is largely
absorbing per-individual variance and is not comparable with a REML h². It is
reported because `fitNullModel` returns it, not because it answers anything.

### 11.8 HEADLINE — step 05 as specified is invalid on this GRM, and why

**This changes a conclusion, so it is here rather than in a table.** §8 warned
that step 05's numerics had never been executed. They have now, and the step as
shipped **does not estimate what it claims to**.

`05_reml_zaitlen.sbatch` builds its second component with GCTA
`--make-bK 0.05` on `$GRM_FULL` — i.e. it decides who is a close relative from
the **pooled multi-ancestry GRM's off-diagonals**. Those off-diagonals confound
kinship with ancestry. That is v1 §8.6, the finding §1 of this README calls
"the single most important finding in v1", and step 05 walks straight back into
it.

**Measured on the bK GRM this run produced** (counted directly out of the
binary, 11,670 × 11,670 lower triangle):

| | pairs retained as "close relatives" | share of all 68,088,615 pairs |
|---|---|---|
| `--make-bK 0.05` on the pooled GRM | **3,837,371** | 5.64 % |
| PC-Relate, κ ≥ 0.0442 (3rd degree+) | 29,371 | 0.043 % |
| PC-Relate, κ ≥ 0.025 (= the same 0.05 GRM-value threshold) | 90,365 | 0.13 % |

**42–130× too many pairs.** G2 is therefore an ancestry-similarity component,
not a pedigree one, and it swallows the additive variance. Every phenotype
shows it (jobs 34277750_1–5):

| phenotype | what AI-REML did | v1's h² for comparison |
|---|---|---|
| `global_slope` | **FAILED** at iteration 8: `V(G1)`→0, `V(e)`→0, `V(G2)`→1.043, *"more than half of the variance components are constrained"* | 0.137 ± 0.046 |
| `slope_PC3` | `V(G1)` constrained to **0.000**, G2 takes 0.255 | 0.041 ± 0.045 |
| `slope_PC2` | `V(G1)` dropped out entirely — fit reduced to G2 + residual | 0.161 ± 0.046 |
| `slope_PC1` | `V(G1)` → 0.037, G2 holds 0.193 | 0.075 ± 0.046 |
| `baseline_thickness` | ran away the other way: `V(G1)` = 0.705 → h²_snp ≈ **0.67** | 0.246 ± 0.049 |

Both directions are the same artefact. A component containing 3.8 M
ancestry-driven pairs either absorbs the additive variance (four phenotypes) or
leaves the dense G1 free to absorb residual population structure that the ten
PCs did not remove (`baseline_thickness`). **Do not report any number from
`results_v2/reml_zaitlen/`.**

**The fix, and it is the one v2's own logic dictates:** define the second
component from **PC-Relate** — kinship estimated *conditional on the ancestry
PCs* — instead of from the pooled GRM. The model, the GRM values, the diagonal
and the threshold semantics are unchanged; only the *set of pairs called close*
changes, and it changes to the estimator this pipeline exists to use.

Implemented in three files under `hpc_v2/work/setup/`, all config-driven:

| file | what it does |
|---|---|
| `export_pcrelate_pairs.R` | reads `kinship/pcrelate.rds`, writes the pairs above `BK_THRESH/2` (GCTA GRM values are ≈ 2 × kinship — the same factor-of-2 §7 flags for `pcrelateToMatrix`) |
| `make_bk_pcrelate.py` | streams the dense GRM's lower triangle and zeroes every off-diagonal PC-Relate does not call related; asserts the output is byte-for-byte the same size |
| `reml_zaitlen_pcrelate.sbatch` | the same GCTA `--reml --mgrm` call, against `reml_zaitlen_pcrel/mgrm.txt`, writing to a separate directory so both arms stay visible |

Build: job 34279445 (1 m 24 s) — 90,365 off-diagonals kept, 67,998,250 zeroed.
Re-run: job 34279502. Results in 11.9.

### 11.9 The PC-Relate second component fixes the pair set — and reveals a deeper problem

Two corrected arms were run, differing only in which PC-Relate pairs define G2:

| arm | κ floor | pairs in G2 | job | outcome |
|---|---|---|---|---|
| `reml_zaitlen_pcrel/` | 0.025 (= `BK_THRESH`/2, 3rd degree+) | 90,365 | 34279502 | still unstable — `global_slope` drove `V(G1)` to 0 again |
| `reml_zaitlen_pcrel2nd/` | 0.0884 (2nd degree and closer) | **2,907** | 34280243 | **all five converged, ~1 m 15 s each** |

The 2nd-degree arm converges and passes both of §6's checks — `n` = 8,082 (the
full sample) and `h2_ped ≥ h2_snp` everywhere:

| phenotype | v2 `h2_snp` | v1 h² (n = 5,649) | v2 `h2_ped` | agreement |
|---|---|---|---|---|
| `baseline_thickness` | 0.548 ± 0.052 | 0.246 ± 0.049 | 0.846 ± 0.015 | **disagrees, ~6 SE** |
| **`global_slope`** | **0.115 ± 0.056** | **0.137 ± 0.046** | 0.273 ± 0.033 | consistent |
| `slope_PC2` | 0.193 ± 0.058 | 0.161 ± 0.046 | 0.286 ± 0.032 | consistent |
| `slope_PC1` | 0.139 ± 0.056 | 0.075 ± 0.046 | 0.172 ± 0.033 | within ~1 SE |
| `slope_PC3` | 0.099 ± 0.056 | 0.041 ± 0.045 | 0.165 ± 0.033 | within ~1 SE |

§2's honest expectation 1 is confirmed and then some: the SEs are **larger**
than v1's (0.056 vs 0.046 for `global_slope`) despite n rising from 5,649 to
8,082. Relatives carry redundant genotype information, and the second component
costs precision on the first.

**But `h2_ped` from this fit must not be reported as a pedigree heritability,
and here is the evidence.** Two facts:

1. `--make-bK` and its PC-Relate replacement both keep the source GRM's
   **diagonal**, so G1 and G2 have the *same* diagonal and differ only in
   off-diagonals. `h2_ped = (V(G1)+V(G2))/Vp` is therefore identified largely
   by the diagonal-vs-identity contrast — which is why its SE (0.015 for
   `baseline_thickness`) is three times *smaller* than the SE of either
   component it sums (0.052 and 0.053). A quantity supposedly estimated from
   2,907 relative pairs cannot be that precise.
2. **That diagonal is an ancestry axis.** Measured on this GRM (11,670
   subjects): diagonal mean 1.0346, **SD 0.4546**, range 0.681–2.419, and

   > **corr(GRM diagonal, PC-AiR PC1) = −0.990  (R² = 0.980;  0.988 on PC1+PC2)**

So the pooled multi-ancestry GRM's diagonal is essentially PC1 rescaled, both
Zaitlen components inherit it, and `h2_ped` is measuring ancestry-driven
heteroscedasticity rather than pedigree resemblance. The same shared diagonal
makes the `h2_snp`/`h2_ped` split ill-conditioned, which is why the answer moves
with the κ floor.

**Conclusion on step 05.** The pair-set fix is necessary but not sufficient.
The Zaitlen decomposition cannot be estimated cleanly from *any* thresholding of
this pooled GRM, because the confound is in its diagonal as well as its
off-diagonals. What can be taken from the run:

- `h2_snp` for `global_slope`, `slope_PC1`, `slope_PC2`, `slope_PC3` **supplements**
  v1's REML numbers — same direction, overlapping intervals, now with every
  relative in the sample. It does not replace them, and it is not more precise.
- `baseline_thickness`'s 0.548 **conflicts** with v1's 0.246 and with the
  published range for cortical thickness. Do not report it. v1's 0.246 ± 0.049
  stands as this project's estimate.
- **`h2_ped` is not usable from this design.** Getting it would need a second
  component whose diagonal is not an ancestry variable — a within-ancestry GRM,
  or GENESIS with the PC-Relate kinship plus an explicit family matrix (the
  sparse kinship from step 02 is already built and ancestry-adjusted). That is
  new work, not a re-run.

### 11.10 Step 04 — the association scan, and step 06 — within-family PRS

**Step 04 (job 34277748, 110/110 tasks COMPLETED, zero failures).** Measured
cost, replacing §5's guess of "tens of minutes each": **3 m (chr22) to 31 m
(chr1)** per (phenotype × chromosome), 4 CPU / 24 G — the 6 h budget is roughly
12× more than the largest task needed. Collection (job 34277749) took 14 m 58 s
for all five phenotypes.

`assoc/gwas_summary.tsv`, against v1's pooled scan (§8.16.1, n = 7,932):

| phenotype | v2 λ_GC | v1 λ_GC | v2 p<5e-8 | v1 p<5e-8 | v2 min p |
|---|---|---|---|---|---|
| `baseline_thickness` | **1.0363** | 1.107 | 0 | **20** | 8.23e-08 |
| `global_slope` | 1.0139 | 1.013 | 0 | 0 | 5.11e-07 |
| `slope_PC1` | 1.0142 | 0.999 | 0 | 0 | 1.63e-06 |
| `slope_PC2` | 1.0358 | 1.030 | 1 | 0 | 2.90e-08 |
| `slope_PC3` | 1.0048 | 0.993 | 2 | 3 | 2.03e-08 |

9,413,281 SNPs after MAF ≥ 0.01 / MAC ≥ 20 (v1 tested 9,077,609), n_mean 7,893.
Every λ_GC is below §6's 1.10 alarm line, and `baseline_thickness`'s drop from
1.107 to 1.036 is the clearest single sign that PC-AiR PCs plus a PC-Relate
kinship control structure better than a pooled GCTA GRM did.

**What happened to v1's 20 hits — checked per-SNP, not inferred from λ.** Both
v1 loci are still v2's top two signals, at just above threshold:

| locus | v2 top SNP | v2 p | v1 |
|---|---|---|---|
| chr2 ~26.93 Mb | rs7599286 | **8.23e-08** | 8 SNPs at p < 5e-8 |
| chr12 ~69.38 Mb | rs2603090 | **9.47e-08** | 12 SNPs at p < 5e-8 |

An ancestry artefact would have collapsed. These moved by a factor of ~2 in p
and stayed at the top of the genome — consistent with v1 §8.16.1's per-SNP
argument, while no longer clearing 5e-8. **This supplements v1's claim; it does
not replace it. Neither run is a replication of the other** (shared subjects,
different model).

**Step 06 (job 34277818, 1 m 11 s).** 160 rows in
`prs_family/prs_withinfamily.tsv`; **1,339 informative sibling pairs** in the
pooled arm (686 in the EUR arm), against the ~688 §2 anticipated — the 7.0
release's sibling structure is bigger than v1's EUR-only sample had.

The load-bearing consistency check: for SCZ at the threshold carrying v1's
population signal (`0p5`, pooled, `global_slope`),

| quantity | value |
|---|---|
| v1 population estimate (`prs_imp/prs_association.tsv`) | −0.0413 ± 0.0145 (p = 0.0043) |
| v2 `beta_between` | **−0.0456 ± 0.0151** (p = 0.0025) |
| v2 `beta_within` | +0.0070 ± 0.0497 (p = 0.89) |
| v2 `beta_diff` (between − within) | −0.0526 ± 0.0518 (p = 0.31) |

`beta_between` reproducing v1's population estimate is exactly what should
happen and is the best available evidence that the two pipelines agree where
they are meant to. **`beta_diff` is not significant, so this run gives no
evidence that the population estimate is inflated by confounding** — and,
per §2.3, the null `beta_within` is *not* evidence against the effect: the
script itself reports ~18 % power at 1,339 pairs.

One EUR-arm caveat, verified against the files: v2's EUR stratum is
`results/ancestry/eur_anchor.keep` ∩ phenotyped = **4,116**, whereas v1's PRS
EUR arm used `prs_assoc.R`'s 3-SD PC-distance cut and had **5,361**. The two EUR
arms are not the same subjects; only the pooled arm (8,082 in both) is directly
comparable.

### 11.11 Measured costs, for the next person sizing these jobs

| step | shape | §5's guess | **measured** |
|---|---|---|---|
| 01 gds | 23 tasks, 4 CPU / 16 G | "minutes each" | 25 s – 1 m 34 s; 73 GB of GDS |
| 02 kinship | 1 job, 8 CPU / 64 G | "a few hours" | **55 m 53 s** |
| 03 null model | 5 tasks, 4 CPU / 32 G | "minutes each" | 1 m 15 s each |
| 04 assoc | 110 tasks, 4 CPU / 24 G | "tens of minutes" | 3–31 m (chr-dependent) |
| 04c collect | 1 job | — | 14 m 58 s |
| 05 reml | 5 tasks, 8 CPU / 64 G | "~minutes" | **wrong: 22 m to failure**; the converging PC-Relate 2nd-degree arm took 1 m 15 s |
| 06 prs_family | 1 job | "seconds" | 1 m 11 s |
| env build | 1 job, 4 CPU / 16 G | — | 10 m |

Total wall clock from first submission to the last GENESIS output: **~1 h 45 m.**
Queue time was negligible on icelake at this size, except that the 110-task
scan hit `QOSMaxCpuPerUserLimit` and held up the REML re-runs behind it — worth
knowing before submitting a wide array and a separate job in the same window.

### 11.12 The decisive test: `h2_ped` does not depend on the relatives at all

Both corrected arms eventually converged (the κ ≥ 0.025 arm needed 1 h 10 m –
1 h 35 m per phenotype against the 2nd-degree arm's 75 s — slow convergence is
itself a symptom). Comparing them is the cleanest possible test of 11.9's
claim, because **the only difference between the two fits is the number of
relative pairs in G2: 90,365 versus 2,907, a 31-fold cut.**

| phenotype | v1 h² | `h2_snp` 90,365 pairs | `h2_snp` 2,907 pairs | Δ | `h2_ped` 90,365 | `h2_ped` 2,907 | Δ |
|---|---|---|---|---|---|---|---|
| `baseline_thickness` | 0.246 | 0.539 ± 0.055 | 0.548 ± 0.052 | +0.009 | 0.846 ± 0.0153 | 0.846 ± 0.0153 | **+0.0000** |
| `global_slope` | 0.137 | 0.183 ± 0.049 | 0.115 ± 0.056 | −0.068 | 0.252 ± 0.032 | 0.273 ± 0.033 | +0.021 |
| `slope_PC1` | 0.075 | 0.178 ± 0.051 | 0.139 ± 0.056 | −0.039 | 0.178 ± 0.032 | 0.172 ± 0.033 | −0.006 |
| `slope_PC2` | 0.161 | 0.230 ± 0.050 | 0.193 ± 0.058 | −0.037 | 0.279 ± 0.031 | 0.286 ± 0.032 | +0.007 |
| `slope_PC3` | 0.041 | 0.096 ± 0.052 | 0.099 ± 0.056 | +0.003 | 0.165 ± 0.033 | 0.165 ± 0.033 | −0.0002 |

**Delete 87,458 of the 90,365 relative pairs and `h2_ped` does not move** —
0.846 ± 0.0153 both times for `baseline_thickness`, to four decimal places and
with the same SE; 0.165 both times for `slope_PC3`. A pedigree heritability
estimated *from relatives* cannot be invariant to removing 97 % of them.
`h2_ped` here is the shared GRM diagonal — which is PC-AiR PC1 at R² = 0.98 —
and nothing else. 11.9's diagnosis is confirmed by direct experiment rather
than by argument.

`h2_snp` does respond to the pair set (`global_slope` 0.183 → 0.115, a 1.2 SE
swing on the primary phenotype), which is the other half of the same problem:
with G1 and G2 sharing a diagonal, only the off-diagonals separate them, and
the answer follows whatever you put there.

**Bottom line for step 05.** Two things are robust across both arms and worth
keeping: `n` = 8,082, and `h2_snp` for `baseline_thickness` sitting at
**0.54 in both arms against v1's 0.246** — the two-GRM fit inflates the
positive control by more than a factor of two regardless of the threshold, so
the inflation is not a threshold artefact. Everything else in step 05 is
threshold-dependent or diagonal-driven. **v1's `reml_imp_pooled/reml_summary.tsv`
remains this project's heritability result**; v2 contributes the demonstration
that the Zaitlen route is closed on a pooled GRM, and a specification for what
would open it (11.9's last paragraph).

---

## 12. EUR-stratified arm, MAGMA and LDSC — 2026-09-07

**Question this section answers: did v2 add signal, or only rigour?**
Answer: **rigour in the pooled arm (a real confound removed), precision in the
EUR arm, one strengthened MAGMA result, and no new discoveries.** Details below;
every number is from a summary file named in the row.

Jobs: EUR kinship 34985236 (15 m 38 s), null models 34985244 (37 s each),
association 34985245 (110/110, 1 m 36 s – 15 m 40 s), collect 34985246
(15 m 49 s), MAGMA+LDSC 34985266 (1 h 30 m), prioritised gene sets 34985572
(5 m 14 s), ancestry-matched rg 34985917 (2 m 5 s). Plus the pooled LDSC
diagnostic 34983075 (5 m 8 s).

### 12.1 How the EUR arm was built

Same 4,116 subjects as v1's EUR GWAS (`eur_anchor.keep` ∩ phenotyped), so the
arms are directly comparable. PC-AiR and PC-Relate were **re-run within the
European anchor set** (5,656 genotyped) rather than reusing the full-sample
axes: inside one continental group the leading pooled PCs are near-constant and
carry no information about within-EUR structure. v1's EUR arm used EUR-specific
PCs for the same reason.

MAGMA and LDSC were run by invoking **v1's own `04_magma.sbatch` and
`05_ldsc_rg.sbatch`** with `GWAS_DIR`/`MAGMA_DIR`/`LDSC_DIR` redirected, so any
difference between the arms is the GWAS and not the downstream code.

**Within-EUR relatedness, against the pooled run** — an independent
confirmation of §11.9's diagnosis:

| pair class | pooled (n = 11,670) | within EUR (n = 5,656) | expected within EUR if the pooled counts were real |
|---|---|---|---|
| 1st degree | 1,504 | 748 | ~353 |
| 2nd degree | 982 | **14** | ~231 |
| 3rd degree | 26,464 | **106** | ~6,218 |

First-degree pairs *exceed* the naive scaling (siblings are real and
concentrated), while the 2nd/3rd-degree bands collapse by 16–60×. **The pooled
scan's large distant-relative bands were cross-ancestry residue**, which is
exactly why the κ ≥ 0.025 Zaitlen arm was unstable and the 2nd-degree arm was
not.

### 12.2 GWAS

| | | v1 | v2 | | v1 | v2 |
|---|---|---|---|---|---|---|
| arm | phenotype | λ_GC | λ_GC | | p<5e-8 | p<5e-8 |
| **pooled** | `baseline_thickness` | 1.107 | **1.036** | | 20 | 0 |
| | `global_slope` | 1.012 | 1.014 | | 0 | 0 |
| | `slope_PC1` | 0.999 | 1.014 | | 0 | 0 |
| | `slope_PC2` | 1.030 | 1.036 | | 0 | 1 |
| | `slope_PC3` | 0.993 | 1.005 | | 3 | 2 |
| **EUR** | `baseline_thickness` | 1.024 | 1.027 | | 0 | 0 |
| | `global_slope` | 1.000 | 1.000 | | 0 | 0 |
| | `slope_PC1` | 1.007 | 1.005 | | 0 | 0 |
| | `slope_PC2` | 1.017 | 1.019 | | 0 | 0 |
| | `slope_PC3` | 0.999 | 0.995 | | 1 | 0 |

n: pooled 7,932 (v1) vs 7,893 (v2); EUR 4,039 vs 4,017. SNPs: pooled 9.08M vs
9.41M; EUR 7.18M vs 7.47M.

**The whole GWAS difference is in the pooled arm.** In EUR the two pipelines
agree to within noise on every λ — as they should, because within one
continental group there is little structure for PC-AiR/PC-Relate to fix. That
is itself the control that makes the pooled difference interpretable.

### 12.3 HEADLINE — v1's LDSC intercept was confounding, not misspecification

v1 §8.16.2 recorded an LDSC intercept of **1.0995 ± 0.0085** for
`baseline_thickness` on the pooled sumstats with EUR LD scores, and attributed
it to *"the EUR LD scores applied to a multi-ancestry sample … not population
structure in the GWAS"*. That attribution does not survive v2:

| pooled sumstats, EUR LD scores | v1 | v2 |
|---|---|---|
| intercept, `baseline_thickness` | **1.0995 ± 0.0085** | **1.0036 ± 0.0060** |
| h² | 0.424 ± 0.080 | 0.298 ± 0.056 |
| **in-sample** LD scores, h² | **−0.062 to −0.495** (z to −8.04) | **+0.313, +0.107, +0.059, +0.034, −0.015** |

v2's scan is **equally multi-ancestry** and used the **same EUR LD scores**. If
the intercept were LD-score misspecification it would be unchanged; it is 1.00.
**It was residual population structure, and PC-AiR + PC-Relate removed it.**

The in-sample LD panel says the same thing from the other side. Run against v1's
sumstats it produced strongly *negative* heritability (`slope_PC3` −0.495,
z = −8.04) — recorded here for the first time, from the previously undocumented
`results/ldsc_imp_insample/`. Against v2's sumstats the same panel behaves. The
mechanism (interpretation, not measured here): in-sample LD scores are inflated
at SNPs with large between-ancestry frequency differences; v1's χ² was inflated
at those SNPs too, but not proportionally, so the regression slope went
negative. Remove the structure and the pathology goes.

**Consequence for v1's results.** §8.16.6 recommends in-sample LD scores as
"the fix for LDSC's multi-ancestry arm". That recommendation was written before
the result arrived and is **contradicted by it on v1 sumstats** — though it is
sound advice for v2 sumstats. Update it rather than repeating it.

### 12.4 LDSC genetic correlation (EUR arm)

| disorder | phenotype | rg v1 | rg v2 | h² z v1 | h² z v2 |
|---|---|---|---|---|---|
| SCZ | `baseline_thickness` | 0.019 ± 0.057 | 0.026 ± 0.052 | 3.89 | **4.49** |
| SCZ | **`global_slope`** | −0.199 ± **0.218** | −0.149 ± **0.109** | 0.57 | 1.16 |
| SCZ | `slope_PC1` | 0.163 ± 0.369 | 0.094 ± 0.150 | 0.32 | 0.66 |
| SCZ | `slope_PC2` | −0.066 ± 0.062 | −0.052 ± 0.060 | 3.18 | 3.46 |
| MDD | `baseline_thickness` | −0.027 ± 0.058 | −0.019 ± 0.053 | 3.62 | **4.29** |
| MDD | `global_slope` | −0.043 ± 0.151 | −0.051 ± 0.107 | 0.58 | 1.04 |
| MDD | `slope_PC1` | −0.018 ± 0.466 | 0.013 ± 0.168 | 0.26 | 0.43 |
| MDD | `slope_PC2` | −0.019 ± 0.065 | −0.008 ± 0.061 | 3.26 | 3.77 |

**Every rg is null in both arms, and no conclusion changes.** What changes is
precision: rg SEs fall by 2–2.8× on the slope phenotypes, and
`baseline_thickness` crosses LDSC's h² z ≥ 4 usability line for the first time
(LDSC's own `underpowered` flag flips from yes to no).

**Read the mechanism honestly.** The SE fall is not extra precision on h² — the
h² SEs are nearly unchanged (0.108 → 0.104). It is that the **h² point
estimates rose**, and rg's SE scales as 1/√(h²₁h²₂). For `global_slope`, LDSC h²
goes 0.062 → 0.124, which now sits on top of the GREML estimate of 0.137 where
before it was half of it. EUR-arm intercepts are ~1.00 in *both* versions
(0.993/1.005/1.020 v1, 0.993/1.000/1.009 v2), so this is not inflation. Whether
it is less attenuation or chance is not settled by one comparison.

`slope_PC3` returns negative h² in both arms (−0.31 v1, −0.23 v2); its rg is
not estimable either way.

### 12.5 A mismatch neither version had noticed: the disorder GWAS are themselves multi-ancestry

v1 §8.16.3 labels its EUR arm "properly specified". That was true only of the
ABCD side. The disorder sumstats in `work/sumstats/` are:

- **SCZ** = `PGC3_SCZ_wave3.**primary**` — the multi-ancestry release
  (74,776 cases; European **and** East Asian).
- **MDD** = `pgc-mdd2025_no23andMe_**div**` — the trans-ancestry release
  (Adams et al. 2025).

The European-only SCZ release, `PGC3_SCZ_wave3.**core**` (55,085 cases /
78,957 controls), was already on disk and unused. **So no arm of this project,
v1 or v2, has ever been ancestry-matched on both sides.**

Re-running the EUR rg against the EUR-only SCZ (`ldsc_eur_matched/`, N taken as
2 × NEFF to match `prepare_sumstats.sh`'s SCZ convention):

| phenotype | rg, v1 EUR × SCZ-core | rg, v2 EUR × SCZ-core |
|---|---|---|
| `baseline_thickness` | 0.018 ± 0.059 | 0.022 ± 0.054 |
| `global_slope` | −0.153 ± 0.168 | −0.121 ± 0.099 |
| `slope_PC2` | −0.052 ± 0.061 | −0.042 ± 0.059 |
| `slope_PC1` | 0.132 ± 0.250 | 0.077 ± 0.132 |

**Matching the disorder side changes nothing** — every estimate moves less than
a third of an SE. The `SCZ × global_slope` sign stays negative in all four
specifications now on record (v1 pooled −0.140, v1 EUR −0.199, v2 EUR −0.149,
v2 EUR matched −0.121), consistent and never significant.

An equivalent EUR-only MDD release is **not** on disk; PGC ships one and it
would cost only a download.

### 12.6 MAGMA — AHBA C1–C3

| phenotype | set | β v1 | p v1 | β v2 | p v2 |
|---|---|---|---|---|---|
| `slope_PC1` | **C1−** | −0.042 | 0.204 | −0.035 | 0.283 |
| `global_slope` | C1− | −0.046 | 0.163 | −0.030 | 0.359 |
| `slope_PC3` | C1− | −0.043 | 0.196 | −0.043 | 0.193 |
| `baseline_thickness` | C3− | −0.091 | 0.216 | −0.136 | 0.060 |

**§8.16.4's ruling stands and v2 does not overturn it.** The C1− →
`slope_PC1` result (pooled: β = −0.107, p = 5.1e-04) remains **null in the
ancestry-matched arm** under v2 as it was under v1. The sign is consistent
across slope phenotypes in both, and that is all. Nothing in C1/C2/C3 reaches
significance anywhere in the EUR arm; the strongest cell is
`baseline_thickness` × C3− at p = 0.060, one of 30 tests. **Do not report the
C1− result.**

### 12.7 MAGMA — SCZ and MDD gene sets (EUR arm)

| set → phenotype | β v1 | p v1 | β v2 | p v2 |
|---|---|---|---|---|
| **`SCZ_locus_pool` → `baseline_thickness`** | 0.153 | 1.4e-03 | **0.187** | **1.1e-04** |
| `SCZ_locus_pool` → `global_slope` | 0.127 | 7.0e-03 | 0.124 | 7.4e-03 |
| `SCZ_prioritised` → `baseline_thickness` | 0.110 | 0.116 | 0.154 | **0.048** |
| `SCZ_prio_finemap` → `baseline_thickness` | 0.119 | 0.150 | 0.116 | 0.158 |
| `MDD_hc_finemap` → `baseline_thickness` | 0.174 | 8.6e-03 | 0.147 | 0.023 |
| `MDD_highconf` → `baseline_thickness` | 0.115 | 0.034 | 0.075 | 0.117 |
| `MDD_pool` → `slope_PC2` | 0.046 | 0.042 | 0.034 | 0.098 |

**The one result that genuinely strengthens is the project's most robust one.**
`SCZ_locus_pool → baseline_thickness` goes from p = 1.4e-03 to **1.1e-04**, and
now clears Bonferroni over all 35 set × phenotype tests (threshold 1.4e-03) with
an order of magnitude to spare, where v1 sat exactly on the line.
`SCZ_locus_pool → global_slope` is unchanged at p ≈ 7e-03 — the §8.16.4 finding
that the SCZ locus pool reaches the *slope* replicates under a different
association model.

**MDD moves the other way**, uniformly and modestly: every MDD row weakens by
roughly a factor of 2–4 in p. §8.16.4's "MDD is no longer null" claim rests on
`MDD_hc_finemap → baseline_thickness`, which survives at p = 0.023 but no longer
at 8.6e-03. Read that claim as weaker than v1 stated.

### 12.8 PRS (unchanged by design, plus the new within-family arm)

PRS scoring does not use our GWAS — it applies PGC weights to our genotypes — so
**v2 changes nothing about the population PRS estimates**. What v2 adds is the
sibling decomposition (§11.10):

| disorder | phenotype | thr | v1 population β | v2 β_between | v2 β_within | p_diff |
|---|---|---|---|---|---|---|
| SCZ | `global_slope` | 0p5 | −0.041 (p 4.3e-03) | −0.046 (p 2.5e-03) | +0.007 (p 0.89) | 0.31 |
| SCZ | `slope_PC1` | 0p5 | 0.044 (p 3.3e-03) | 0.045 (p 3.8e-03) | 0.031 (p 0.55) | 0.80 |
| SCZ | `slope_PC2` | 0p5 | 0.033 (p 0.030) | 0.036 (p 0.023) | −0.001 (p 0.99) | 0.50 |
| MDD | `global_slope` | 0p01 | −0.046 (p 0.016) | −0.040 (p 0.047) | **−0.115 (p 0.071)** | 0.26 |
| MDD | `slope_PC2` | 0p01 | 0.051 (p 0.013) | 0.055 (p 9.3e-03) | 0.004 (p 0.95) | 0.46 |
| MDD | `baseline_thickness` | 1 | −0.059 (p 0.032) | −0.074 (p 0.011) | 0.039 (p 0.57) | 0.12 |

`β_between` reproduces every v1 population estimate, which is the check that the
two pipelines agree where they must. **No `p_diff` is significant**, so this run
gives no evidence that the population PRS estimates are inflated by
between-family confounding — and, at 1,339 pairs and ~18 % power, no evidence
against it either. The only within-family estimate approaching nominal
significance is MDD → `global_slope` (β_W = −0.115, p = 0.071), *larger* than
its between-family counterpart and in the same direction; at this power that is
a curiosity to re-test in a larger release, not a result.

### 12.9 Verdict: what v2 added

**Added, and worth keeping:**
1. A demonstrated confound removal in the pooled arm — intercept 1.0995 → 1.0036,
   λ_GC 1.107 → 1.036, and in-sample LD scores that work instead of returning
   negative h². This reclassifies v1's pooled inflation from "misspecification"
   to "structure", which is a correction to the record, not a new finding.
2. `SCZ_locus_pool → baseline_thickness` strengthened from p = 1.4e-03 to
   1.1e-04, clearing multiple-testing correction with room.
3. 2–2.8× tighter rg standard errors and `baseline_thickness` crossing LDSC's
   h² z ≥ 4 line for the first time.
4. The within-family PRS arm, which has no v1 counterpart.
5. Two documentation corrections: the in-sample-LD recommendation (§8.16.6) and
   the "properly specified EUR arm" claim (§8.16.3), neither of which held.

**Not added:**
1. **No new loci** — the pooled arm loses v1's 20 hits (they attenuate to
   p = 8.2e-08 and 9.5e-08, still the top two signals; §11.10).
2. **No rg becomes significant**, in any of four specifications.
3. **No AHBA component becomes significant**; C1− stays null in the matched arm.
4. **MDD gene-set results weaken** by 2–4× in p across the board.

**Bottom line: v2 bought correctness and precision, not discovery.** The
project's scientific claims are where v1 left them, with one strengthened
(SCZ locus pool → baseline thickness), one weakened (MDD high-confidence sets),
and one retired (v1's 20 genome-wide hits).

---

## 13. Control disorders, local rg and PRS-CS — 2026-09-08

Three additions aimed at the same question: **is the SCZ association specific to
adolescent-onset illness, and is it real?**

### 13.1 HEADLINE — the onset-age control panel does not behave as hoped, but the reason matters

Two control disorders were run through the identical pipeline (rg, MAGMA, P+T
polygenic scores): **ASD** (childhood onset; SPARK+iPSYCH+PGC, N = 58,794) and
**Alzheimer's** (late onset; PGC-ALZ2, Wightman 2021 excl. 23andMe).  Both are
already on disk; neither needed downloading.

Naively the panel fails: ASD's polygenic score predicts `baseline_thickness` at
**p = 8.4e-05**, *stronger* than any SCZ result, and `global_slope` at
p = 1.9e-03.  Alzheimer's (APOE excluded) predicts `global_slope` at p = 0.021.
Taken at face value that would say the association is not developmental at all.

**But the threshold profile discriminates them, and it is decisive.**  A genuine
polygenic overlap strengthens as more SNPs enter the score.  An artefact does
not.  Pooled sample, β (p) by p-value threshold, `global_slope`:

| disorder | 5e-8 | 1e-5 | 0.001 | 0.01 | 0.05 | 0.1 | 0.5 | 1 |
|---|---|---|---|---|---|---|---|---|
| **SCZ** | +0.006 | +0.003 | −0.005 | −0.020 | **−0.033** | **−0.033** | **−0.041** | **−0.039** |
| **MDD** | −0.018 | −0.021 | −0.031 | **−0.046** | **−0.048** | **−0.057** | **−0.055** | **−0.056** |
| ASD | +0.006 | **+0.036** | +0.007 | −0.012 | −0.010 | −0.011 | −0.009 | −0.009 |
| ALZ (no APOE) | +0.008 | −0.010 | −0.020 | −0.014 | **−0.034** | −0.025 | −0.029 | −0.029 |

(bold = p < 0.05; `baseline_thickness` shows the same pattern — ASD significant
only at 5e-8 and 1e-5, sign reversing at looser thresholds.)

- **SCZ and MDD are monotone**: null at stringent thresholds, strengthening
  steadily as SNPs accumulate.  That is the signature of distributed polygenic
  overlap.
- **ASD is the opposite**: significant at 1e-5 on **69 SNPs**, dead null at
  31,740 / 183,700 / 271,508 SNPs, and the sign flips.  That is a handful of
  loci or a threshold-selection artefact, not polygenic sharing.  It is present
  in EUR too (p = 3.2e-03), so it is not simply stratification.
- **Alzheimer's is the genuinely awkward one**: weakly monotone in the SAME
  direction as SCZ (p = 0.02 at 0.05, p ≈ 0.06–0.07 at the loosest thresholds),
  never surviving correction.  It is not clean evidence against specificity, but
  it is not the flat null the hypothesis wants either.

**Reading.** The panel supports SCZ/MDD being polygenic in a way ASD is not, and
that is a real discriminator the single-threshold summary would have missed.  It
does **not** establish onset-age specificity, because Alzheimer's does not
cleanly fail.  Both controls should be reported.

**Caveat on the rg arm of the controls:** Alzheimer's LDSC h² on these sumstats
is 0.005 ± 0.004, so its *genetic-correlation* null (all p > 0.6) is
uninformative — the test could not have detected anything.  Only the polygenic
arm of the ALZ control carries weight.

### 13.2 Local genetic correlation (LAVA) — the MAGMA finding becomes testable

Genome-wide rg averages over the whole genome and was null.  MAGMA said the
signal concentrates in the SCZ **locus pool**.  LAVA tests that directly by
estimating rg within each of 2,495 approximately independent LD blocks.

**223 of 2,495 blocks contain a SCZ locus-pool gene.**  Restricting to those is
pre-specified by MAGMA — computed before any local rg existed, so it is not
selection on the outcome — and cuts the testing burden 11-fold.  1,053 bivariate
tests ran (LAVA requires significant univariate local h² in both traits first;
100–124 blocks passed per phenotype, 213–222 for the disorders).

| phen1 | phen2 | LOC | rho [95% CI] | p | q (BH) | |
|---|---|---|---|---|---|---|
| baseline_thickness | MDD | 864 | −1.00 [−1.00,−0.54] | 2.5e-05 | 0.014 | **boundary — unstable** |
| **slope_PC1** | **SCZ** | **335** | **+0.62 [+0.35,+0.91]** | **2.6e-05** | **0.014** | |
| global_slope | SCZ | 335 | −0.63 [−1.00,−0.33] | 1.3e-04 | 0.046 | |
| slope_PC3 | SCZ | 335 | −0.90 [−1.00,−0.47] | 1.9e-04 | 0.050 | |
| slope_PC3 | MDD | 2483 | −1.00 [−1.00,−0.49] | 2.5e-04 | 0.052 | boundary |
| slope_PC2 | SCZ | 1719 | −0.50 [−0.81,−0.24] | 3.5e-04 | 0.062 | |

2 pass Bonferroni (p < 4.75e-05), 4 pass FDR q < 0.05.  What the loci contain:

| LOC | region (GRCh37) | SCZ locus-pool genes |
|---|---|---|
| **335** | chr2:144.5–146.0 Mb | **ZEB2** |
| **1719** | chr11:112.8–113.9 Mb | **DRD2** |
| 2483 | chr22:38.7–40.4 Mb | ATF4, **CACNA1I**, MGAT3, MIEF1, RPS19BP1 |
| 864 | chr5:106.4–107.3 Mb | EFNA5 |

**Locus 335 is the result.**  Three phenotypes independently show local rg with
SCZ in the same block, two of them past Bonferroni, and the block's only
locus-pool gene is **ZEB2** — a transcription factor required for cortical
interneuron migration and the gene behind Mowat-Wilson syndrome.  A
developmental transcription factor is exactly the kind of locus this hypothesis
predicts.  DRD2 at locus 1719 is the canonical SCZ gene and reaches p = 3.5e-04.

**Three cautions, none of which the tables should hide.**
1. The three locus-335 rows are **not independent** — the slope PCs and
   `global_slope` come from the same measurements.  They are one finding seen
   three ways, which is reassuring for consistency but is not triplicate
   evidence.
2. `slope_PC1` is +0.62 while `global_slope` is −0.63 at the same locus.  That
   is coherent only if PC1 loads opposite to the global slope; **this has not
   been checked and must be before the sign is interpreted.**
3. Estimates of exactly ±1 with a CI touching the bound are boundary solutions,
   not perfect correlations.  Both MDD rows above are of that kind.

The genome-wide LAVA run is still going; it is needed to establish whether 4
hits in 1,053 targeted tests exceeds the genome-wide background rate.  **Until
that lands, the enrichment claim is not established** — only the individual loci.

### 13.3 PRS-CS — set up, and four input defects found

PRS-CS replaces clumping+thresholding with a Bayesian continuous-shrinkage prior
over all HapMap3 SNPs, learning the global shrinkage from the data
(PRS-CS-auto).  It removes the eight-threshold multiple-testing burden entirely
and typically buys 1.5–2× the variance explained.  All four disorders are run
through it, deliberately: a specificity claim resting on a better-powered SCZ
score against weaker control scores would be an artefact of method.

Four defects were found and fixed getting there, all of the silent kind:

1. **One all-tab blank row** in `SCZ_core.tsv` made PRS-CS's parser raise
   `IndexError` on every chromosome (it splits on whitespace and indexes
   `ll[1]`).
2. **`rs283815` in PGC-ALZ2 has z = inf.**  PLINK loaded it happily
   ("50,440 valid predictors") and returned `SCORESUM = 0` for all 11,670
   subjects.  The APOE-inclusive Alzheimer's score in `control_ALZ/` is
   therefore void; the APOE-excluded score, which drops that variant, is
   correct and is the one reported in §13.1.
3. **A race across array tasks.**  All 22 chromosome tasks per disorder built
   the same input file concurrently; each read a partial write.  SCZ came out
   as **391 SNPs instead of 7,585,076** and would have produced a meaningless
   score without ever failing.  Input preparation is now its own job
   (`prscs_prep.sbatch`).
4. **A "skip if it exists" guard** around that preparation also enclosed the
   `n_gwas` assignment, so a stale file left the sample size unset and it fell
   back to 100,000 — for MDD, whose real effective N is **927,713**.  PRS-CS
   scales its shrinkage by n, so this would have quietly mis-weighted the score.

Corrected inputs: SCZ 7,585,076 SNPs (n = 130,000), MDD 7,363,302 (927,713),
ASD 2,910,895 (58,794), ALZ 1,121,059 (74,004).

**Alzheimer's needed one further correction to be a fair control.**  PGC-ALZ2
supplies Z, not beta, and no allele frequency.  Feeding Z to PRS-CS is worse
than to P+T, because PRS-CS models `beta_hat ~ N(beta, sigma²/n)` and Z rescales
every SNP by 1/SE, which varies with MAF.  Beta is reconstructed as
`z / sqrt(2p(1-p)(n + z²))` using MAF from PRS-CS's own 1000G EUR reference, so
the control is weighted on the same scale as the other three.

### 13.4 CORRECTION to 13.2 — the ZEB2 locus is not enriched, and the local-rg nulls are not calibrated

The genome-wide LAVA run (job 35061400, 3 h 21 m) has now landed, and it
**retracts the enrichment reading of 13.2**. 10,412 bivariate tests across
2,028 blocks, of which the 223 SCZ locus-pool blocks contribute 1,053.

**The targeted blocks show no excess of local genetic correlation whatever:**

| threshold | in SCZ locus-pool blocks | elsewhere | ratio | Fisher p |
|---|---|---|---|---|
| p < 0.05 | 118/1,053 (11.21 %) | 1,040/9,359 (11.11 %) | 1.01× | 0.92 |
| p < 1e-3 | 7/1,053 (0.66 %) | 55/9,359 (0.59 %) | 1.13× | 0.67 |
| p < 1e-4 | 2/1,053 (0.19 %) | 13/9,359 (0.14 %) | 1.37× | 0.66 |
| Bonferroni (4.8e-06) | 0/1,053 | 1/9,359 | — | 1.00 |

So the MAGMA-motivated hypothesis — that SCZ/thinning sharing concentrates in
the SCZ locus pool — **is not supported**. The four FDR hits reported in 13.2
are what 1,053 tests drawn from anywhere in this genome produce.

**The nulls are also anti-conservative, which inflates every p-value in 13.2:**

| p-value bin | observed | expected under a calibrated null |
|---|---|---|
| < 0.01 | 3.48 % | 1 % |
| 0.01–0.05 | 7.65 % | 4 % |
| < 0.05 (total) | **11.1 %** | **5 %** |

A ~3.5× excess in the tail, *genome-wide* — not confined to the targeted
blocks. The likely cause is LAVA's own univariate pre-filter: conditioning on
local h² being significant in both traits selects blocks where the statistics
are already large, so the bivariate null is no longer uniform. Whatever the
cause, **the genome-wide empirical distribution, not the nominal p-value, is
the right reference**, and against it the targeted set is unremarkable.

**Where that leaves ZEB2.** `slope_PC1 × SCZ` at locus 335 survives *genome-wide*
FDR (q = 0.044) — but as the **7th of 7** such hits, the other six sitting in
blocks with no SCZ locus-pool gene, and two of those seven being boundary
solutions (|rho| = 1). It is a plausible locus with a developmentally sensible
gene, and it is not distinguishable from background. **Do not report ZEB2 as a
finding.** The internal consistency noted in 13.2 (all four slope phenotypes
agreeing in sign once their −0.583 correlation with `global_slope` is accounted
for) remains true and remains the reason the locus is worth revisiting in a
better-powered sample — but consistency across correlated measures of one
phenotype is not independent evidence.

**What 13.2 got right and what it got wrong.** Right: the method, the
pre-specification, the multiple-testing arithmetic within the targeted set, and
the sign-consistency check. Wrong: treating "passes FDR within the targeted
set" as evidence of enrichment without the genome-wide comparison that was still
running. That comparison was flagged as outstanding in 13.2 and it has now
falsified the reading. The lesson is the one this project keeps relearning — a
result computed on a selected subset means nothing until the background rate is
measured on the same pipeline.

**One thing genuinely survives**: `baseline_thickness` failed the univariate
local-h² gate at locus 335 (p = 0.096) while all four slope phenotypes passed
(p = 9.3e-07 to 0.034). Local heritability there is specific to the rate of
change rather than the cross-sectional measure. That is a statement about our
phenotypes' own local h², independent of the disorder correlation and of the
enrichment failure.

### 13.5 PRS-CS — the headline SCZ association does not survive it

PRS-CS-auto completed for SCZ, MDD and ASD (ALZ still running at time of
writing). Posterior weights: SCZ 1,085,902 SNPs, MDD 1,029,904, ASD 369,623 —
the HapMap3 intersection, against 7M for clumping+thresholding.

**The comparison is not like-for-like in the direction that flatters C+T**, and
that is the point. v1's reported p-values are the *best of eight thresholds*,
already selected on the outcome and needing a Bonferroni factor of 8. PRS-CS-auto
is a **single pre-specified test** with no threshold to tune.

| disorder | phenotype | C+T best-of-8 β (p) | PRS-CS β (p) | PRS-CS r² |
|---|---|---|---|---|
| SCZ | **`global_slope`** | −0.041 (**4.3e-03**) | −0.026 (**0.195**) | 0.00066 |
| SCZ | `slope_PC1` | +0.044 (3.3e-03) | **+0.052 (0.011)** | **0.00273** |
| SCZ | `slope_PC2` | +0.033 (0.030) | +0.028 (0.186) | 0.00076 |
| SCZ | `slope_PC3` | −0.029 (0.035) | −0.033 (0.105) | 0.00111 |
| SCZ | `baseline_thickness` | +0.016 (0.282) | −0.021 (0.302) | 0.00044 |
| MDD | `global_slope` | −0.046 (0.016) | −0.029 (0.111) | 0.00085 |
| MDD | `slope_PC2` | +0.051 (0.013) | +0.030 (0.128) | 0.00087 |
| MDD | `baseline_thickness` | −0.059 (0.032) | −0.032 (0.088) | 0.00104 |
| ASD | `baseline_thickness` | −0.047 (**8.4e-05**) | −0.015 (**0.250**) | 0.00024 |
| ASD | `slope_PC3` | −0.034 (5.4e-03) | −0.034 (**0.012**) | 0.00117 |

**Two readings, and they point opposite ways.**

1. **The project's headline result weakens badly.** SCZ → `global_slope` goes
   from p = 4.3e-03 (p_adj = 0.034 after the eight-threshold correction) to
   **p = 0.195** under a method that needs no correction at all. In the EUR arm
   it is p = 0.284. A result that survives one scoring method and not the other,
   at this effect size, is a lead and not a finding.
2. **The ASD control result evaporates**, exactly as 13.1's threshold-profile
   argument predicted: p = 8.4e-05 → 0.250. That is a genuine confirmation that
   the ASD signal was a 69-SNP threshold artefact, and it is reassuring about
   the diagnostic logic in 13.1.

**What survives PRS-CS**: SCZ → `slope_PC1` (p = 0.011, r² = 0.0027 — the
largest r² in the table, and the only cell where PRS-CS *beats* C+T), and
ASD → `slope_PC3` (p = 0.012). Neither survives correction across the 15 cells
above (threshold 0.0033).

**Do not read PRS-CS's weakness as proof of absence.** Three caveats:
- PRS-CS is restricted to HapMap3 (~1.1M SNPs) where C+T used the full 7M set.
- Its LD reference is 1000G EUR while the pooled target sample is 32 %
  non-European; the reference is mismatched for those subjects. (C+T's clumping
  reference was our own sample, so this cuts the other way.)
- PRS-CS-auto learns one global shrinkage parameter; if it over-shrinks at this
  discovery power, real signal is attenuated.

The EUR arm removes the ancestry mismatch and gives SCZ → `global_slope`
β = −0.017, p = 0.284 — same direction, no significance, n = 4,116.

**Net effect on the project's claims.** Combined with 13.4, two of the three
strands supporting SCZ → faster thinning are now weaker than reported: the local
genetic correlation is not enriched, and the polygenic association does not
survive a method without threshold selection. The MAGMA locus-pool enrichment
(p = 7.4e-03, reproduced across v1 and v2) is the strand still standing.

### 13.6 The Alzheimer's arm completes the panel — and the panel does not discriminate

All 88 PRS-CS tasks completed (0 failures). The Alzheimer's score is the last
piece, and it needed one extra step to interpret.

**Raw result:** ALZ → `global_slope`, β = −0.0248 ± 0.0126, **p = 0.0486** —
nominally significant, same direction as SCZ, and with a *smaller* p-value than
SCZ's own PRS-CS result. Taken at face value that is fatal to the specificity
claim.

**But 77.4 % of the score's squared posterior weight sits on 786 APOE-region
variants** (chr19:44.4–46.5 Mb, of 1,085,921 total). PRS-CS concentrates weight
where the evidence is, and in Alzheimer's that is overwhelmingly APOE. So the
score is an APOE score wearing a polygenic label, and "does APOE predict
adolescent cortical thinning" is a different question from "does polygenic
Alzheimer's risk predict it".

Re-scoring with the APOE region dropped (1,085,135 variants, same posteriors):

| disorder (PRS-CS, pooled n = 8,082) | β | SE | p |
|---|---|---|---|
| SCZ | −0.0256 | 0.0198 | 0.195 |
| MDD | −0.0292 | 0.0184 | 0.111 |
| ASD | −0.0026 | 0.0131 | 0.841 |
| **ALZ, APOE included** | **−0.0248** | 0.0126 | **0.0486** |
| **ALZ, APOE excluded** | **−0.0201** | 0.0141 | **0.155** |

**Two conclusions, and the second is the important one.**

1. **The nominal Alzheimer's association is APOE-driven** — it halves in
   significance and loses nominal significance once APOE is removed. Worth
   recording as its own observation (APOE effects on adolescent cortex are a
   real and separate literature), not as evidence about polygenic AD risk.
2. **Under the untuned method, schizophrenia and non-APOE Alzheimer's are
   statistically indistinguishable**: β = −0.026 (p = 0.195) against β = −0.020
   (p = 0.155), same direction, overlapping intervals. The onset-age panel does
   not discriminate — **not because the disorders differ, but because none of
   them is significant.** A control panel can only demonstrate specificity when
   the target effect is itself detectable, and here it is not.

**Where the project's central claim now stands.** Of the three strands:

| strand | status after this session |
|---|---|
| Local genetic correlation in SCZ loci | **retracted** — no enrichment over genome-wide background (§13.4) |
| Polygenic score, SCZ → `global_slope` | **does not survive** an untuned method (p = 0.195), and is matched by a late-onset control (§13.5, §13.6) |
| MAGMA SCZ locus pool → `global_slope` | **stands** — p = 7.0e-03 (v1) and 7.4e-03 (v2), reproduced under different association models |

One strand, not three. The MAGMA result is also the one whose statistic depends
least on our phenotype's heritability, which is consistent with h² being the
binding constraint rather than the effect being absent. Nothing here refutes the
hypothesis; what it removes is the impression that several independent lines of
evidence supported it. They were not independent, and two of them do not hold up
under their own follow-up tests.

**What this argues for next.** Every failure in §13.4–§13.6 traces to the same
root: `global_slope` has h² ≈ 0.137 with an LDSC h² z of 1.16, so no
cross-trait method has the power to separate a real small effect from zero.
Raising that number — MOSTest across parcels, a better-estimated slope, more
waves — is worth more than any further re-analysis of the current phenotype.
