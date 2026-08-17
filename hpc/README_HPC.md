# HPC genetics pipeline — state and plan

**Last updated: 2026-08-17.** Written for an agent starting fresh on CSD3 with
no access to the conversations that produced it.

> **You should not need to read `legacy/`.** Those two files are the full
> working record through 2026-08-17, kept because they contain the diagnostic
> detail behind everything asserted here. Read them **only** if the repo state
> contradicts this file, or if you need the derivation of a specific claim:
> - [`legacy/README_HPC_2026-08-17.md`](legacy/README_HPC_2026-08-17.md) — the
>   cluster run log: 17 sections, every defect found and fixed, full result tables.
> - [`legacy/README_pipeline_2026-08-10.md`](legacy/README_pipeline_2026-08-10.md)
>   — the pre-cluster design document, written before anything had run.
>
> Bugs recorded there are **fixed**. Do not re-litigate them; the only ones that
> can still bite you are reproduced in §6 below.

---

## 1. The project, in one paragraph

**Which genes drive adolescent cortical development?** The phenotype is a
per-subject *rate* of cortical thinning — a random slope of thickness on age
from repeated scans — not thickness at one timepoint. The hypothesis is that
genes driving this rate are enriched for schizophrenia and depression GWAS
signal, and connect to this group's prior transcriptional work (AHBA components
C1–C3). Local pipeline (`src/abcd/`, `R/fit_lmm.R`) produces subject-level
phenotypes; everything genotype-related runs here.

Five phenotypes go to genetics, all standardised BLUPs so every β is per SD:
`baseline_thickness` (model-predicted thickness at age 12.44, the high-reliability
positive control), `global_slope` (cortex-wide thinning rate, **the primary
phenotype**), and `slope_PC1/PC2/PC3`.

## 2. What has been run, and what it found

The full pipeline (steps 01–06) **has executed end to end on CSD3**. Do not
treat this as unrun code. Genotypes were `abcd_eur` (13,697,177 SNPs × 5,678
subjects), **N = 4,119** after intersecting phenotypes, covariates and genotypes.

**Heritability** — GCTA REML on 3,329 unrelated:

| phenotype | h² | SE | p |
|---|---|---|---|
| `baseline_thickness` | **0.575** | 0.147 | 4.3e-05 |
| `global_slope` | 0.228 | 0.141 | 0.051 |
| `slope_PC2` | 0.113 | 0.141 | 0.214 |
| `slope_PC1`, `slope_PC3` | 0.000 | ~0.14 | 0.500 |

LDSC independently gives `baseline_thickness` h² = 0.584 ± 0.131 against REML's
0.575 — two estimators, different samples, agreeing to 0.01. That jointly
validates the GRM, ID alignment, covariates and sumstats munging.

**GWAS** — fastGWA-MLM, 8.78M SNPs: **no genome-wide-significant loci for any
phenotype**, which is expected at N = 4,119, not a surprise. λ_GC 1.00–1.03, so
the mixed model is properly calibrated.

**Genetic correlation** — null and, for the slope, *uninformative*: only
`baseline_thickness` has h² z > 4 (z = 4.45). SCZ rg = 0.036 ± 0.047. Every
slope row is flagged `underpowered` in `results/ldsc/ldsc_rg_summary.tsv`.

**Polygenic scores** — 160 models; nothing survives correction. The one coherent
pattern is SCZ → `global_slope`, monotone across thresholds, reaching
β = −0.038 (p = 0.015 EUR, 0.009 full) at p<0.5. Negative β means **higher SCZ
risk → faster thinning**, the predicted direction. `p_adj` = 0.119 Bonferroni
across 8 thresholds. **A lead to power up, not a result.**

**Gene sets** — the SCZ locus pool is enriched in `baseline_thickness`
(β = 0.183, p = 1.9e-04) and nominally in `global_slope` (p = 0.017); the 104
*prioritised* genes carry none. MDD high-confidence sets are null throughout.
Caveat that still stands: SCZ loci are brain-expressed and cortical thickness is
a brain phenotype, so a non-brain negative-control set is needed to call this
SCZ-specific. It has not been built.

### Follow-up analyses done locally since (not on the cluster)

Worth knowing so you do not redo them; all committed:

- **Regional SCZ PRS maps**: 63/68 regions negative, 10 at p<0.05, and the map
  correlates with thickness slope PC1 (ρ = +0.28, spin p = 0.014). MDD's map is
  uncorrelated with SCZ's — disorder-specific, not a generic PRS artefact.
- **Age × PRS interaction**: SCZ × age p = 0.035, MDD null, both main effects
  null. **The effect is on the rate of thinning, not the level attained.**
- **Within-family test**: no detectable confounding (between-minus-within
  ≈ 0), but only 688 pairs, so 11 % power — it bounds confounding, it cannot
  confirm the effect.
- **Multiple-comparison correction was wrong**: Bonferroni across 8 *nested*
  score thresholds treats them as independent. Three effective-test estimators
  give m_eff = 2.7–4.6, putting the corrected p at 0.024–0.041.
- **T1w/T2w ratio** was fitted as an alternative metric: slope h² = 0.595 (vs
  0.444 for thickness) but site ICC 0.153 (vs 0.045) — 3.4× more site-confounded,
  and **no** PRS association on either its slope or its intercept. Not a
  replication; read it as a different phenotype.
- **Selecting the top-10 PRS-associated regions is a dead end**: h² is *lower*
  than the whole-cortex mean and at the median of random 10-region sets, and
  split-half cross-validation shows the selected phenotype performs worse out of
  sample. Do not build a GWAS on a selected region set.

## 3. The one thing that matters now: N

**Every result above is limited by N, not by method.** The analysis sample is
4,119: those both phenotyped (8,192) and in the EUR-only prebuilt GRM (5,678).

| scenario | N | h² z | PRS power |
|---|---|---|---|
| current (EUR GRM, 4.0 genotypes) | 4,119 | 1.45 | 76 % |
| cross-ancestry GRM only | 5,678 | 2.0 | 88 % |
| 7.0 genotypes only (EUR) | 5,719 | 2.0 | 88 % |
| **both** | **~7,278** | **2.6** | **94 %** |

**Set expectations honestly.** An interpretable LDSC rg needs h² z ≈ 4, which at
h² = 0.165 requires N ≈ 11,300 — above ABCD's *phenotyped* ceiling of 8,192. So
rg on the slope stays out of reach in ABCD alone even after this succeeds. What
the extra N buys is a well-powered PRS test and a usable h² point estimate.
Do not sell this as unlocking rg.

## 4. Where the missing subjects went — checked, and the answer is clean

A real question was whether subjects absent from the 4.0 genetics were **dropped
by genotyping QC** or **never genotyped**. From ABCD's published documentation
(`docs.abcdstudy.org` genetics page for 7.0; 6.0 release notes):

| | subjects |
|---|---|
| ABCD enrolled cohort | 11,868 |
| **7.0 curated Smokescreen PLINK** (~515k variants) | **11,670** |
| 7.0 TOPMed r3 imputed | 11,670 |
| 7.0 WGS (separate product) | 8,710 |
| never genotyped (6.0 notes: "still not been genotyped") | ~200 |
| removed — consent withdrawal / relatedness inconsistency | 4 |

**11,868 − 11,670 = 198, against ~200 never genotyped plus 4 removed. The
arithmetic closes.** The curated 7.0 genotype file is essentially the whole
enrolled cohort, so what we have been working with is *not* QC attrition:

| what we had | subjects | shortfall vs 11,670 | why |
|---|---|---|---|
| `abcd_eur` (every result above) | 5,678 | 5,992 | **ancestry restriction** |
| 4.0 all-ancestry per-chromosome | 10,072 | 1,598 | **release vintage** |

Both are recoverable by moving to 7.0. Neither represents subjects who failed
genotyping.

**What this means for you:** if the 7.0 fileset you find has ~11,670 subjects,
it is the full curated set and nothing is missing. **If it has appreciably
fewer, someone has pre-filtered it — find out on what basis before building a
GRM from it.** An undocumented subset is exactly how the 4.0 attempt failed.

*Provenance: these counts are from ABCD's public documentation, not from files
on CSD3. Confirm against the actual `.fam`; the table tells you what "correct"
looks like.*

## 5. YOUR TASK: the cross-ancestry GRM on the 7.0 genetics

The 7.0 genetics has been acquired. A previous attempt on the **4.0** data
failed — that failure was a property of those files, not of the approach.

### 5.1 Do this first — one second, and it would have saved 22 jobs

```bash
bash hpc/work/check_bfile_integrity.sh "$GENO_ALLANC_DIR"/ABCD_chr{1..22}_hg19
```

All 22 array tasks in the 4.0 attempt died in seconds on `Unexpected PLINK 1
.bed file size`. That is diagnosable from file sizes alone without submitting
anything, because a PLINK 1 `.bed` is exactly `3 + ceil(n/4) * m` bytes
(n = `.fam` rows, m = `.bim` rows). The script inverts that and reports **which
file is wrong**, because the three causes need different fixes:

- **`.bed` holds more variants than the `.bim` lists** → the `.bim` was filtered
  without regenerating the `.bed`. This was the 4.0 defect (`.bim` at ~78 % of
  the `.bed`). **Not recoverable by filtering further** — nothing records which
  `.bed` columns the surviving `.bim` rows correspond to. You need the original
  `.bim` or a re-derived fileset.
- **`(size-3)` does not divide by `ceil(n/4)`** → the `.fam` does not match;
  usually the wrong `.fam` was copied and the right one still exists.
- **magic bytes ≠ `6c1b01`** → not variant-major PLINK 1. `6c1b00` is
  sample-major; fix with `plink --make-bed`.

It is wired into `grm_allanc.sbatch` as a gate, so a bad fileset aborts that
chromosome in about a second with a diagnosis rather than a GCTA error. Tested
against a valid fixture and a deliberately reconstructed copy of the 4.0 defect.

**If 7.0 shows the same mismatch, stop and report it.** Do not work around it by
filtering: misaligned variants produce plausible-looking, meaningless h². GCTA
refused loudly here; that was the good outcome.

### 5.2 Running it

Point the scripts at the 7.0 tree in `hpc/config.local.sh` (gitignored — never
commit it; `config.local.sh.example` is the template):

```bash
GENO_ALLANC_DIR="/path/to/7.0/per-chromosome/filesets"
GENO_ALLANC_TPL='ABCD_chr{CHR}_hg19'   # {CHR} is substituted; 7.0 naming may differ
```

```bash
bash hpc/work/check_bfile_integrity.sh "$GENO_ALLANC_DIR"/<prefixes>   # gate
sbatch hpc/work/grm_allanc.sbatch                                      # array 1-22
sbatch --dependency=afterok:<jobid> hpc/work/grm_allanc_merge.sbatch   # merge + PCA
```

The merge job produces everything downstream consumes, so switching over is a
config edit rather than another round of jobs:

| output | config variable | consumed by |
|---|---|---|
| `abcd_all` | `GRM_ALLANC` | `02_reml` |
| `abcd_all_sp` | `GRM_ALLANC_SPARSE` | `03_gwas` (fastGWA) |
| `abcd_all.unrel` | `GRM_ALLANC_UNREL` | `02_reml` unrelated subset |
| `abcd_all_pca.eigenvec` | `GRM_ALLANC_PCA` | ancestry covariates (20 PCs) |

To switch the pipeline over, set `GRM`, `GRM_SPARSE` and `GRM_UNREL` to the
`GRM_ALLANC*` values in `config.local.sh`. They are deliberately separate
variables so the EUR results stay reproducible instead of being overwritten.

The merge refuses to run on an incomplete set of 22 rather than quietly building
from fewer — that quiet-wrong failure has already happened twice here.

### 5.3 Two things that are NOT optional on the all-ancestry set

1. **Use the in-sample PCs** from `abcd_all_pca.eigenvec`. An audit found that
   swapping the release's multi-ancestry PCs for in-sample EUR PCs moved h² by
   only 0.009 — but **that reassurance does not transfer here**. Within EUR
   there was no structure left for a PC to correct; across ancestries it is real
   and large. Replace `PC1..PC10` in `covar_quant.txt` with the first 10–20
   columns of the eigenvec.

2. **A single pooled GRM across ancestries is a real methodological compromise.**
   It assumes common allele frequencies and LD structure, which is exactly what
   does not hold across these groups; in-sample PCs mitigate but do not repair
   it. Compute the pooled estimate, but **also** run ancestry-stratified REML and
   compare. If pooled h² differs materially from the EUR estimate, the
   defensible design is stratified REML meta-analysed across groups, and the
   pooled number should not be the headline.

### 5.4 Definition of done

1. `check_bfile_integrity.sh` passes on all 22 chromosomes.
2. `abcd_all.grm.id` has appreciably more than 5,678 subjects.
3. **The phenotyped ∩ genotyped join exceeds 4,126** — this is the number that
   matters. `hpc/work/align_ids.py` already does the token join.
4. `02_reml` re-run on `global_slope` and `baseline_thickness` with the new GRM
   and in-sample PCs, reported **alongside** the EUR estimates, not replacing them.
5. Ancestry-stratified REML run and compared, per §5.3.2.
6. If N materially improves, re-run `06_prs` — the SCZ → `global_slope` lead is
   the result most likely to change.

## 6. Gotchas that can still bite you

Everything else in `legacy/` is fixed. These are live:

- **IDs.** ABCD spells subject IDs three ways (`NDAR_INVxxxxxxxx`,
  `sub-NDARINVxxxxxxxx`, and 7.0's bare `sub-xxxxxxxx`). GCTA identifies an
  individual by the **FID+IID pair**, so a mismatched FID gives a
  **zero-subject analysis that reports no error**. Always join on the
  8-character NDAR token; `work/align_ids.py` does this and derives the target
  spelling *from the `.fam`* rather than assuming one.
- **Verify outputs, not `sacct`.** A previous run had 13 jobs reporting
  `COMPLETED` while producing nothing. Every step writes a `*_summary.tsv`;
  check it exists with the expected row count before believing a green state.
- **A bare `python` on CSD3 is 3.7.4.** Use the environments under `work/`
  (`work/envs/abcd`, `work/bin/Rscript`). Nothing usable is on `PATH` — no
  `gcta64`, `magma`, `plink`, no Python ≥ 3.10, and system R lacks `lme4`.
- **`/home` is nearly full** (48.6 of 52.4 GB). Everything goes to
  `/rds/user/<user>`; `ABCD_HPC_ROOT` relocates the whole tree.
- **`06_prs` needs PLINK 1.9 specifically**, for `--clump`/`.clumped` and
  `--score f 1 2 3 sum`. PLINK 2's `--clump` has a different output format.
- **Bash brace trap.** `${GENO_ALLANC_TPL:-ABCD_chr{CHR}_hg19}` does **not**
  work: bash ends the expansion at the `}` closing `{CHR}`, so `_hg19}` is
  appended as literal text and any override silently gains that suffix.
  `${VAR:=...}` has the same problem. `config.sh` uses an `if [[ -z ... ]]`
  assignment. Keep `{CHR}`-style placeholders out of `${...:-...}` defaults.
- **All paths are config variables.** `hpc/config.sh` defines every path as
  `${VAR:-default}` and is sourced by every script; `config.local.sh` is sourced
  *first* so a root set there propagates through the derived defaults. No script
  should contain an absolute path — if you add one, you have broken portability
  for the next person on a different account.

## 7. How to document what you do — please read this

Your work gets pulled back and read by an agent that **was not present for it**
and cannot see your terminal, your job output, or your reasoning. What follows
is what makes that handover work.

**Update this file, in place.** Add a dated section for what you ran. Do not
start a third README; if this file grows unwieldy, move the superseded parts to
`legacy/` with a datestamp and keep this one current — that is what happened to
produce the two files in there now.

**Commit messages should carry reasoning, not a file list.** The diff already
shows which files changed. What it cannot show is *why*, what you ruled out, and
what you would have done differently — and that is what the next agent needs.

**For every result, record:**

- the **job ID** and the script that produced it, so it can be traced in `sacct`
  and the logs;
- the **N** it ran on — this pipeline's history is one long argument about N,
  and a result without its N cannot be compared to anything;
- the **file** it landed in (`results/<step>/*_summary.tsv`), committed;
- whether it **replaces or supplements** an existing number. If it replaces one,
  say which and keep the old value visible rather than silently overwriting.

**Record failures and dead ends too, with the diagnosis.** The most useful parts
of the legacy file are the defects: seven were found and fixed, and several were
*silent* — producing plausible numbers rather than errors. If something looked
right and was not, that is the highest-value thing you can write down. If you
abandoned an approach, say why, so nobody re-attempts it.

**Distinguish what you verified from what you assumed.** State which claims you
checked against a file and which you inferred. The single most expensive class
of error in this project has been a confident statement nobody re-derived — a
stale table with no generator, prose asserting a result no test covered, and a
figure caption overstating its own panel. If you did not check it, say so.

**Flag anything that changes a conclusion**, prominently, at the top of your
section. If the cross-ancestry h² materially disagrees with the EUR estimate, or
the PRS lead strengthens or evaporates, that is the headline — not a detail in a
results table.
