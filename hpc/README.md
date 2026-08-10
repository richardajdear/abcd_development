# HPC pipeline: heritability, GWAS, and gene-set enrichment

**Status: locally validated, not yet run on real genotypes.** CSD3 refused key
authentication during development (`Permission denied
(publickey,keyboard-interactive,hostbased)`), so nothing here has touched the
ABCD genotypes. Every script is complete, parameterised, and exercised against
a synthetic fixture — see *Local validation* below for exactly what that does
and does not establish. Treat the first cluster run as a dry run
(`DRY_RUN=1 bash hpc/01_grm.sbatch` prints commands without executing) and
check `00_check_inputs.sh` before launching anything expensive.

These scripts assume the CSD3 conventions from the previous analysis:
account `VERTES-SL2-CPU`, partition `cclake`, working tree
`/home/rajd2/rds/hpc-work/ABCD/`. Override via `hpc/config.sh`.

## Configuration: no hard-coded paths

`hpc/config.sh` is sourced by every script and defines every path as
`${VAR:-default}`. Two consequences worth knowing:

- **Exporting `ABCD_HPC_ROOT` alone relocates the whole tree.** The defaults
  derive from each other (`GRM_DIR` from `OUT`, `OUT` from `ABCD_HPC_ROOT`), and
  `hpc/config.local.sh` is sourced *before* them so that a root set there
  propagates. Overriding one path individually is also fine.
- **`require_paths` validates by variable name** and knows which variables are
  file-set *prefixes* (`GENO`, `MAGMA_REF`, `GRM`, `GRM_UNREL`, `GRM_SPARSE`)
  rather than files, so a missing input is reported as
  `ERROR: GRM=... -- neither ....grm.bin nor ....grm.gz found` up front rather
  than as a tool error an hour in.

Copy `hpc/config.local.sh.example` to `hpc/config.local.sh` (gitignored) to
point the pipeline at a different tree.

## Local validation

```bash
python tools/make_test_genotypes.py     # synthetic fixture: 1500 subjects, 5000 SNPs
bash tools/local_test.sh                # 44 checks
```

The fixture carries real ABCD FIDs/IIDs and real family structure, so sibling
and twin pairs are present and relatedness pruning has something to do. One
phenotype (`sim_h2_50`) is simulated from 20 causal SNPs with h²=0.5 as a
positive control for the cluster run.

**What this establishes.** GCTA ships statically linked against Intel MKL, which
aborts on any CPU without AVX — under Rosetta on Apple Silicon that is every
invocation, so the three numerical steps (GRM element computation, REML,
fastGWA) cannot execute here. But GCTA parses every option and opens every input
file *before* entering the linear algebra, so an invocation that dies in MKL has
already proved its options are valid, its files exist and parse, its IDs
intersect, and `--mpheno` is in range. `tools/local_test.sh` distinguishes those
two failure points across all 15 GCTA invocations. This is what caught
`--geno`/`--hwe` being passed to GCTA — they are PLINK flags, and GCTA
implements neither.

Two GCTA steps are *not* MKL-bound and do run for real: `--grm-cutoff` and
`--make-bK-sparse` read an existing GRM and threshold it. So the relatedness
logic — the part most likely to be silently wrong — is genuinely verified:
1500 → 726 unrelated at 0.05, sparse GRM of 2,486 pairs with mean off-diagonal
0.41. The GRM they read is built by PLINK 1.9, which writes GCTA's own format
and runs natively; note PLINK's `--rel-cutoff` and GCTA's `--grm-cutoff` use
different pruning algorithms, so the retained subset is not bit-identical to
what the cluster will produce.

`04_magma.sbatch` is exercised against `scratch/bin/magma`, a stub that records
its arguments and emits MAGMA-format outputs. That tests invocation construction
and output parsing only — every number it produces is the stub's.

**LDSC and PRS run for real.** Both are pure Python / PLINK, so stages 5 and 6
execute the actual tools. `05_ldsc_rg` munges all eight sumstats files and
produces the rg summary; the fixture's 5,000 SNPs and synthetic LD scores are
far too little information for a *defined* rg, so rg comes back `NA` by design
and the checks assert the machinery (sign check, log parsing, the
`underpowered` flag) rather than any value. `06_prs` clumps, scores, and fits
every association model — 72 models over 6 phenotypes × 2 disorders × 3
thresholds × 2 strata — and includes a null-calibration check: with random
weights against a simulated phenotype, nothing should survive correction, and a
pile of hits there would mean the family random effect is not being fitted.

Stage 6's R step needs `Rscript` with `lme4`/`lmerTest` on `PATH`; it skips with
a message rather than failing when they are absent, so run the harness from an
environment that has both if you want those four checks to execute.

**What it does not establish.** No heritability estimate, no association
statistic, and no gene-set enrichment here is real. On the cluster this is a
pre-flight for the code, not evidence about the data.

## What runs where, and why

The local pipeline (`src/abcd/`, `R/fit_lmm.R`) produces subject-level
phenotypes. Everything downstream of that is genotype work: it needs the
imputed genotypes, several GB of LD reference, and GCTA/MAGMA binaries, all of
which live on the cluster. Nothing here needs a GPU; all of it is
CPU-and-IO-bound, and the GRM step is the only one that needs real memory.

## Order of operations

1. `00_check_inputs.sh` — verify genotypes, phenotype file, and binaries exist and that phenotype IDs intersect the `.fam` file. Cheap, run interactively.
2. `01_grm.sbatch` — build the genetic relatedness matrix (GCTA), then prune to unrelated at 0.05. ~2 h, 32 GB.
3. `02_reml.sbatch` — SNP heritability of each phenotype (GCTA REML). Array job over phenotypes.
4. `03_gwas.sbatch` — per-phenotype GWAS (GCTA fastGWA mixed model). Array job.
5. `04_magma.sbatch` — SNP→gene aggregation, then competitive gene-set tests against AHBA components and the snRNA-seq PC1 gene sets.
6. `05_ldsc_rg.sbatch` — genetic correlation (LDSC) with SCZ and MDD. ~1 h.
7. `06_prs.sbatch` — polygenic score association with SCZ and MDD. ~2 h.

Steps 04, 05 and 06 are **siblings, not a chain**: each consumes 03's summary
statistics and none reads another's output, so `run_all.sh` gives all three the
same dependency on 03 rather than serialising them.

Or just `bash hpc/run_all.sh`, which runs the preflight and then submits 01–06
with the right dependencies (or runs them sequentially where `sbatch` is absent).

## Testing association with SCZ and MDD

Three different questions, and they are not interchangeable. The power argument
in *Expected power* below decides which one to lead with.

| step | method | asks | power here |
|---|---|---|---|
| `06_prs` | polygenic score → phenotype regression | does disorder genetic risk predict the developmental phenotype in these subjects? | **good** — discovery N is large (SCZ ~160k, MDD ~674k); only the target regression uses our N |
| `05_ldsc_rg` | LDSC genetic correlation | do the two traits' genome-wide signals covary? | poor — needs adequate h² on *both* traits |
| `04_magma` | competitive gene-set enrichment | are disorder risk genes enriched in this phenotype's gene-level signal? | moderate — aggregates SNPs into genes |

**`06_prs` is the primary test.** rg requires the ABCD side to be
well-estimated, and it is not: at N≈8,000 with slope reliability ≈0.21, the
slope h² z-score is expected around 2, where LDSC's own guidance (z > 4) says rg
is unreliable. A polygenic score moves the statistical burden onto the discovery
GWAS, so the score enters the model as a single well-estimated regressor and a
null becomes interpretable rather than merely inconclusive. `05_ldsc_rg` is
still worth running — `baseline_thickness` (reliability ≈0.9) *is* adequately
powered and acts as the positive control, and the h²/intercept LDSC reports are
a useful check on the GWAS itself. The summary table flags each row
`underpowered=yes/no` from its own h² z-score so a null rg cannot be read as
evidence of no overlap.

**Ancestry.** The discovery GWAS are predominantly European and polygenic scores
transfer poorly across ancestry, so the **EUR subset is primary and the full
multi-ancestry sample is a sensitivity analysis**. PC adjustment controls
stratification but does not repair transferability: in the full sample a null is
ambiguous between "no effect" and "the score does not transfer". `prs_assoc.R`
fits the strata separately and labels them; it never pools. The EUR definition
is recorded in the output — by default a PC-space cut, but pass `--eur-ids` to
use a reference projection instead, which is preferable when one is available.

**Two things `06_prs` does deliberately.** It clumps in the *reference panel*
rather than in ABCD, because using the target sample's own LD would let target
noise choose the index SNPs. And it fits `(1 | family_id)`, because ABCD
contains twins and siblings and treating them as independent understates the
standard errors. Scores are computed at several p-value thresholds since the
best one is not known a priori; the threshold is a nuisance choice rather than a
hypothesis, so `p_adj` is Bonferroni across thresholds within each
phenotype × disorder × stratum. A result significant at one threshold only, and
not after correction, is a threshold-selection artefact.

### Reference data these two steps need

| what | config variable | where to get it |
|---|---|---|
| HapMap3 SNP list | `HM3_SNPLIST` | `w_hm3.snplist`, from the LDSC repo |
| LD scores (EUR) | `LD_REF`, `LD_WEIGHTS` | `eur_w_ld_chr/`, from the LDSC repo |
| LDSC executables | `LDSC`, `LDSC_MUNGE` | `ldsc.py`, `munge_sumstats.py` |

LDSC is Python 2 in its original form; the maintained Python 3 port is on PyPI
as `ldsc`. Two packaging notes if you install it that way rather than from the
Broad repo: the 2.0.1 wheel installs `ldsc.py` only as a console script, so
`munge_sumstats.py`'s `from ldsc import ...` fails until a copy is importable on
`PYTHONPATH`; and it calls `pd.read_csv(delim_whitespace=...)`, removed in
pandas 3.0, so it needs pandas < 3. Both are handled locally by a dedicated
environment — see `hpc/config.local.sh.example`.


## The phenotype decision that matters most

Use the `ct_genetics.yaml` config, **not** `ct_baseline.yaml`, for anything
here. The two differ in whether the mixed model carries a `family_id` random
intercept. With it, family variance is partialled out of the subject BLUPs and
the DZ/sibling correlation of the intercept phenotype goes to −0.62, which is
biologically impossible and would destroy any heritability estimate. Twin
correlations on the whole-cortex mean (207 MZ / 617 DZ pairs):

| phenotype | family RE | MZ *r* | DZ *r* | Falconer *h²* |
|---|---|---|---|---|
| intercept | yes | 0.14 | −0.62 | — |
| intercept | no | 0.88 | 0.47 | 0.83 |
| slope | yes | — | — | 0.54 |
| slope | no | — | — | 0.58 |

The random-effects structure *is* part of the phenotype definition, not a
fitting detail.

## Expected power, stated up front

On release 5.1 the subject-level slope phenotype has reliability 0.15 (2 visits)
to 0.19 (3 visits) per region, and 0.12 for the whole-cortex mean. Averaging
regions does not help because the dominant noise is per-occasion, not
per-region. A slope GWAS on 5.1 is therefore attenuated roughly 5–17× relative
to a perfectly measured phenotype and is **not** expected to yield genome-wide
hits; run it to establish the pipeline, not to discover loci. The intercept
phenotype (reliability 0.86–0.90 without the family term) is the sanity check —
it should recover the known thickness heritability of ~0.8.

This is the main argument for release 7.0: the added third and fourth timepoints
raise slope reliability, and it is reliability, not sample size, that currently
binds.

## Files this pipeline needs that are not in the repo

Each is addressed by a config variable, so none of these paths appears inside a
script — set the variable if your copy lives elsewhere.

| what | config variable | notes |
|---|---|---|
| imputed genotypes | `GENO` | PLINK prefix; expects `$GENO.bed/.bim/.fam` |
| SCZ sumstats | `SCZ_SUMSTATS` | with `SCZ_SNP_COL=ID`, `SCZ_P_COL=PVAL` |
| MDD sumstats (2025) | `MDD_SUMSTATS` | with `MDD_SNP_COL=rsid`, `MDD_P_COL=p_value` |
| LD reference (1000G EUR) | `MAGMA_REF` | PLINK prefix, not a directory |
| gene locations | `MAGMA_GENE_LOC` | e.g. `NCBI37.3.gene.loc` |
| MAGMA SNP→gene annot | `MAGMA_ANNOT` | built by `04_magma` if absent |
| gene-set covariates | `GENESET_DIR` | written by `python -m abcd.magma_export` |

The 2018 MDD sumstats (`PGC_UKB_23andMe_depression_10000.txt`, SNP col
`MarkerName`, p col `P`) were used in the previous analysis and are superseded by
the 2025 release; add a third case in `04_magma.sbatch` if you want both.

## Interpreting `04_magma.sbatch` output

The gene-set test is **competitive**, not self-contained: it asks whether the
AHBA/snRNA-seq gene sets are more associated with the disorder than other genes
of similar size, density, and MAC — those covariates are conditioned internally
(see the `CONDITIONED_INTERNAL` header line in any `.gsa.out`). A
self-contained test would be confounded by gene length alone.

Note that this is a route to the SCZ/MDD hypothesis that is *independent* of the
map-level result — independent, not a fallback. The map-level spatial analysis is
**positive**: 8 of 18 map × component pairs survive the spin test
(`docs/ahba_vs_maps_noglobal.csv`, §7 of the report). The group-mean thinning
rate correlates with C3 (ρ = −0.55, p_spin = 0.001) and C2 (−0.33, p = 0.006),
and the slope PCs are coupled more strongly still: PC3–C2 ρ = +0.85, PC2–C1
ρ = −0.81, both p_spin = 0.001, plus PC1–C1 (−0.51), PC3–C3 (+0.44) and
PC1–C3 (−0.41). Only C1-versus-group-mean is null (−0.09, p = 0.49).

An earlier version of this file claimed the opposite ("no correlation between
the group-mean developmental change map and AHBA C1–C3"), citing a spatial-nulls
figure that is no longer in the repo. That null came from
correlating the **global-adjusted** age coefficient, and the project has since
settled on the no-global specification. The difference is not subtle:
`docs/ahba_by_global_spec.csv` puts adjusted-versus-C3 at ρ = −0.14
(p_spin = 0.45) against −0.55 (p = 0.002) for the total rate. Adjusting for the
cortex-wide mean removes the shared thinning gradient that carries most of the
transcriptional signal, so it answers "which regions thin faster *than average*"
rather than "where does cortex thin", and the AHBA components speak to the
latter.

The distinction still worth keeping is the one the old text was reaching for:
a map-level correlation is about *where* cortex changes on average, while the
MAGMA test is about which genes carry disorder risk in a phenotype's own
association signal. They are separate claims and separate evidence — but the
map-level one is now supporting evidence, not a dead end.

## Ancestry PCs (added automatically for 7.0)

`python -m abcd.gcta_export <run-dir>` now writes the first 10 ancestry PCs into
`covar_quant.txt` from the 7.0 static table (`ab_g_stc__gen_pc__01..32`); 5.1's
local copy has none, and the CLI prints a loud warning in that case. Note the
7.0 `participant_id` is `sub-0A4P0LWM` — it carries **no** `NDARINV` token,
unlike the imaging tables, so the join goes through the adapter rather than a
string transform. Verified exact against the source table (max abs diff 0.0).

Two properties to be aware of:

- The PCs are standardised, so they are **not** ordered by decreasing variance
  in the released table. Do not use variance ordering as a sanity check.
- They are orthogonal in the full genotyped cohort (max |r| = 0.09) but not in
  the imaging subset (max |r| = 0.26), which is expected when subsetting. This
  is harmless collinearity for covariate adjustment; if you want strictly
  orthogonal covariates, recompute PCs within the analysis sample.
- 115 of 8,192 subjects have no PCs and appear as `NA`. GCTA drops them.
