# HPC genetics pipeline — state and plan

**Last updated: 2026-08-17.** Written for an agent starting fresh on CSD3 with
no access to the conversations that produced it.

> **§5 has been executed and the pipeline extended to the imputed genotypes.**
> Start at **§8.1** (findings), then **§8.15** (heritability) and **§8.16**
> (GWAS / rg / MAGMA / PRS). Headlines:
>
> - **N nearly doubled**: 4,126 → 8,082 phenotyped ∩ genotyped.
> - **`global_slope` is heritable**: h² = 0.137 ± 0.046, p = 0.0012 — the primary
>   phenotype, significant for the first time. It needed *both* the extra N and
>   the imputed variant set; neither alone was enough (§8.15).
> - **First genome-wide-significant loci**: 2 loci for `baseline_thickness`,
>   verified not to be ancestry artefacts (§8.16.1).
> - **The SCZ → faster-thinning PRS lead survives correction** in the full
>   sample (p_adj = 0.034) but *not* in the EUR arm the design makes primary
>   (§8.16.5). Strengthened, not settled.
> - **Two things §5 asked for are unsafe as specified.** `--grm-cutoff 0.05` on a
>   pooled multi-ancestry GRM returns a **94 % European** "unrelated" set (§8.6),
>   and the European LD reference behind LDSC/MAGMA/PRS produces an LDSC
>   intercept that reads as confounding but is misspecification (§8.16.2, §8.16.6).
> - **Run every ancestry-sensitive analysis twice** — pooled and EUR-stratified.
>   Where they disagree, the ancestry-matched arm is the one to believe.
>
> §5 is kept as written so the two can be compared.

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
β = −0.038, p = 0.015 (EUR, n = 3,725) and β = −0.040, p = 0.009 (full,
n = 4,126) at p<0.5. Negative β means **higher SCZ
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

---

## 8. 2026-08-17 — the cross-ancestry GRM on the 7.0 genetics

Written by the agent that executed §5, for one that was not present. Every N is
stated, because this pipeline's history is one long argument about N and a
result without its N cannot be compared to anything.

### 8.1 Headline

1. **The 4.0 defect does not recur.** The 7.0 curated fileset passes
   `check_bfile_integrity.sh` exactly, on the source and on all 22 splits.
2. **N nearly doubles: the phenotyped ∩ genotyped join is 8,082**, against 4,126
   on `abcd_eur`. That is *above* the 7,278 §3 projected, because §3 treated the
   7.0 gain and the cross-ancestry gain as more overlapping than they are.
3. **`--grm-cutoff 0.05` is not safe on a pooled cross-ancestry GRM** (§8.6).
   It returns a 94 %-European "unrelated" subset while reporting nothing wrong.
   This is the single most important thing in this section.
4. **The 7.0 array data has a genotyping-batch artefact** §5 did not anticipate
   (§8.5). It does not invalidate the GRM; it is why two GRMs were built.
5. **The extra N did NOT make `global_slope` detectable, and §3's forecast was
   optimistic for a reason worth understanding.** Moving from imputed to array
   SNPs halves every h², so the 2.6x precision gain from N was cancelled by a
   matching fall in the estimate: z went 1.62 -> 1.42. §3 assumed h² would hold
   while N grew. See §8.9(a,b).
6. **New: `slope_PC2` is significantly heritable** (h² = 0.179 ± 0.055,
   p = 4.4e-04, n = 5,649) — the first phenotype here besides the positive
   control to survive correction. It holds under 20 ancestry PCs and reproduces
   in the QC track, but it is a **lead, not a result**: §8.9 says exactly why.

### 8.2 What the 7.0 tree holds — verified against the files, not the docs

| product | contents | build | used |
|---|---|---|---|
| `smokescreen/merged_chroms` | 11,670 × 515,228, **one merged fileset** | hg19 | **yes** |
| `imputed/chr*_dose.vcf.gz` | TOPMed r3 dosages, ~890 GB | **GRCh38** | no |
| `genesis/pcrelate_*`, `unrelateds_individuals.txt` | PC-Relate kinship, PC-AiR unrelated set | — | the unrelated set, yes |

11,670 is exactly §4's documented curated-genotype count, so **nothing has been
pre-filtered** and the §4 warning did not trigger. Checked directly: 0 duplicate
variant IDs, 0 duplicate `chr:pos`, 0 duplicate IIDs, FID == IID on every row,
502,528 autosomal variants (the remainder are chr23/25/26, dropped — a GRM is
conventionally autosomal).

The imputed set's build was read off its own header — `##mis_panel=topmed-r3`,
`##contig=<ID=chr22>`, first chr22 record at position 10,557,776, where hg19
chr22 begins near 16.05 Mb. **It is GRCh38 while everything else in this
pipeline is hg19**, so it is not a drop-in for any existing step: it needs
liftover on top of the conversion. Worth knowing before anyone reaches for it.

### 8.3 Why the array data and not the imputed data

Converting ~890 GB of VCF to PLINK would cost days of wall-clock and hundreds of
GB before a single GRM existed. GREML does not need imputed density — SNP-h² is
estimated from the LD-tagging of ~500k common array SNPs, which is the standard
input for it.

**The cost of that choice, stated plainly:** the multi-ancestry GRM rests on a
different variant set from the EUR GRM (13.7M imputed at MAF 0.001, against
456,015 array SNPs surviving MAF 0.01 here). A pooled-vs-EUR h² difference
therefore confounds ancestry with SNP set. That is what the EUR-only GRM in §8.8
is for: it re-estimates EUR h² on **this** variant set, so the two confounds can
be separated instead of being argued about.

`genesis/pcrelate_*` was not used as the GRM: PC-Relate is a PC-adjusted kinship
estimator meant for association-model relatedness control, not a GCTA additive
GRM, and substituting it would make every h² here incomparable with §2. Its
companion PC-AiR unrelated set *is* used, for the reason in §8.6.

### 8.4 Compute is not the constraint; the queue is

A per-chromosome GRM over 11,670 subjects takes **17 seconds**. The whole
22-chromosome array is under 10 minutes of CPU. Every delay in this run was SLURM
queue time on SL3, and `mybalance` shows SL2-CPU exhausted (0 hours available)
against ~154,000 on SL3 — so SL3 is not a preference, it is the only option.

**Size jobs for backfill.** Dropping the GRM array from 16 CPUs × 3 h to 8 CPUs ×
1 h 45 moved its estimated start from 3.8 hours away to minutes. Ask for what the
job needs, which on this pipeline is far less than the headers suggest.

**And name the partition explicitly — this cost fourteen hours.** `02_reml`,
`03_gwas`, `04_magma`, `05_ldsc_rg` and `06_prs` carry no `#SBATCH --partition`
line, so they take the **cluster default, `cclake`** — *not* `config.sh`'s
`SLURM_PARTITION`, which is only read by scripts that pass it to sbatch
themselves. Four REML arrays submitted this way sat PENDING overnight while 641
cclake nodes were in maintenance, showing nothing but `(Priority)`. Resubmitted
with `--partition=icelake` they ran at once: 60 tasks, all COMPLETED, under two
minutes. **`squeue -o '%P'` before concluding a job is merely queued behind
others** — a job pending on a drained partition looks identical to a job pending
on a busy one.

### 8.5 A defect §5 did not anticipate: missingness *is* genotyping batch

`work/qc_missing_allanc.sbatch`, `work/qc_keeplist_allanc.sbatch`.

Genome-wide autosomal call rate averages 99.08 %, which looks unremarkable. The
distribution is not:

| BATCH | n | mean F_MISS | frac > 5 % |
|---|---|---|---|
| BATCH_5_Saliva | 203 | 0.099 | 100 % |
| BATCH_3_WB | 83 | 0.099 | 100 % |
| BATCH_2_WB | 192 | 0.082 | 100 % |
| BATCH_6 | 33 | 0.074 | 52 % |
| *other six batches* | 11,159 | 0.004–0.017 | ~0 % |

**87.6 % of the variance in per-sample missingness is explained by the ten
Smokescreen batches** (`smokescreen/batch.info`), against **23.5 %** by the ten
release ancestry PCs jointly. So the ancestry correlation that exists
(PC5 r = −0.36) is mostly batch–ancestry confounding — but it does mean the
artefact runs **along the ancestry axes**, which is where a pooled GRM can least
afford it.

**Correcting the obvious-but-wrong mechanism, because it is what a reader will
assume.** GCTA does *not* mean-impute missing genotypes. Verified against
`chr1.grm.N.bin`, whose per-pair SNP count **varies** (29,324–37,800) rather than
being constant: GCTA estimates each pair from the SNPs non-missing in both
members. So the damage is not a shared imputation bias but ~6 % fewer SNPs behind
every relatedness estimate touching a bad-batch subject (35,256 against 37,683 on
chr1) — batch-structured noise in the GRM, which attenuates h². Milder than
imputation would have been; not nothing.

Note what does **not** fix this: genotyping batch as a *covariate*. The
distortion is in the GRM, and a covariate on the outcome cannot repair a
relatedness matrix. Filtering can, so two tracks were built:

| track | filter | subjects | phenotyped ∩ |
|---|---|---|---|
| primary | `--maf 0.01` only | 11,670 | 8,082 |
| QC sensitivity | `--mind 0.05`, `--geno 0.02`, `--maf 0.01` | 11,167 | 7,756 |

`--mind` would drop 2,550 (at 0.01), 956 (0.02), 503 (0.05) or 31 (0.10)
subjects, so 0.05 is a choice and its cost is visible rather than implied.
`--mind` is computed **once, genome-wide**: applied per chromosome inside the
split it would drop a different sample set on each, leaving 22 filesets with 22
different `.fam` files for `--mgrm` to combine anyway.

### 8.6 **`--grm-cutoff` does not mean "unrelated" across ancestries**

**This is the finding. Read it before reusing any part of §5.**

§5.3.2 warned that a pooled GRM assumes common allele frequencies and LD. The
concrete consequence is worse than a bias in h²: it breaks the *relatedness
pruning step that every REML in this pipeline depends on*, and does so silently.

On the pooled GRM the off-diagonal averages **+0.0402 within stratum** and
**−0.0348 between**. A 0.05 threshold therefore sits inside the within-stratum
bulk rather than out in the tail of genuine relatives. Retention through
`--grm-cutoff 0.05`:

| stratum | in GRM | kept | retained |
|---|---|---|---|
| EURlike | 7,448 | 5,633 | **75.6 %** |
| cluster1 | 2,309 | 61 | **2.6 %** |
| cluster2 | 1,390 | 308 | **22.2 %** |
| cluster3 | 523 | 9 | **1.7 %** |
| total | 11,670 | 6,011 | 51.5 % |

**The resulting "unrelated" subset is 94 % EUR-like.** A pooled REML on
`abcd_all.unrel` is an essentially European analysis with a larger N, a plausible
h², and no error anywhere — reported as cross-ancestry. It would have passed
review, because `--grm-cutoff` is correct within one ancestry and nobody thinks
to question it.

Corroborating the mechanism, not just the symptom: genome-wide over 456,015 SNPs
the sampling noise in an off-diagonal is ~0.0015, yet the observed off-diagonal
SD is **0.0643, 43× that**. And while 1,659,570 pairs exceed 0.2, only 21,895
exceed 0.4. The 0.2–0.4 mass is not relatives; it is same-ancestry pairs scored
against pooled allele frequencies.

**The replacement.** ABCD 7.0 ships a GENESIS **PC-AiR** unrelated set
(`genesis/unrelateds_individuals.txt`, 8,181 subjects), built with
ancestry-adjusted kinship, which is exactly the estimator this situation calls
for. It prunes evenly:

| stratum | `--grm-cutoff` retained | PC-AiR retained |
|---|---|---|
| EURlike | 75.6 % | 73.9 % |
| cluster1 | 2.6 % | 65.2 % |
| cluster2 | 22.2 % | 68.6 % |
| cluster3 | 1.7 % | 39.5 % |
| **non-EUR share of the unrelated set** | **5.7 %** | **32.1 %** |

PC-AiR retains **5,649 phenotyped unrelated subjects** against 3,329 in the EUR
analysis. Both definitions are run and both are reported: `--grm-cutoff` because
it is what the published EUR h² used and comparability matters, PC-AiR because it
is the one that is correct here. **A gap between them is the finding, not a
nuisance to resolve.**

Tooling: `work/grm_diagnostics.py` reports all of the above for any GRM;
`work/grm_pcair_subset.sbatch` cuts the PC-AiR unrelated GRM.

### 8.7 Ancestry strata — derived, and here is exactly how

§5.3.2 requires a stratified REML, which needs a stratum label. **ABCD 7.0 ships
none.** Checked, not assumed: the only genetics columns under
`derivatives/tabulated` are `ab_g_stc__gen_pc__01..32` and `gn_y_genrel_*`
(pihat/zygosity) — no ancestry-proportion or ancestry-group variable anywhere.
Nor is a multi-ancestry 1000 Genomes panel on this account (`hpc-work` holds
`g1000_eur` and `g1000_eas` only), so "project onto 1000G and take the nearest
reference centroid" was not available.

So labels are derived, two ways, so neither is trusted alone
(`work/assign_ancestry.py`):

- **`EUR_anchor`** — membership of the `abcd_eur` fileset behind every published
  result here (5,656 of its 5,678 are in the 7.0 GRM). The one label not inferred
  by us.
- **k-means on the in-sample PC scores**, reported with each cluster's overlap
  against that anchor.

Two independent checks that the split is real rather than merely self-consistent:
all **5,656** anchor-EUR subjects land in one cluster at every k in 3–5 and every
non-EUR cluster is **0.0 %** anchor; and an independent k-means on the *release's*
own 32 PCs reproduces the sizes to about 1 % (7,448/2,309/1,390/523 against
7,515/2,286/1,353/509).

**A bug worth recording because it produced a plausible answer.** GCTA writes
unit-norm eigen*vectors*, so all columns have similar SD (0.0091, 0.0088, 0.0088,
0.0093 for PC1–4) even though their eigenvalues are 727/175/49/15 (6.04 %,
1.45 %, 0.41 %, 0.13 % of variance). Clustering the raw columns weighted a
0.13 %-of-variance PC as heavily as a 6.04 % one, and the noise PCs drove the
split: k=4 put 8,385 subjects in one cluster at 54.8 % anchor-EUR and named a
separate 1,477-subject cluster "EURlike" at 72 %. **Cluster PC scores
(eigenvector × √eigenvalue), not eigenvectors.**

**This is not a continental-ancestry classifier and must not be reported as
one.** ABCD contains a large admixed group and k-means gives no partial
membership: an admixed subject lands at whichever centroid is nearest. These
strata are regions of PC space. That supports "is h² stable across strata?" and
does **not** support "h² in African-ancestry participants is X".

### 8.8 In-sample PCs — §5.3.1 confirmed, and it matters here

`insample_pcs.py`'s docstring records that inside EUR the swap moved h² by 0.009,
because there was no structure left for a PC to correct. Across ancestries there
is: **in-sample PC1 explains 6.04 % of variance here against 0.13 % within EUR**,
PC2 1.45 %, PC3 0.41 %. §5.3.1's "real and large" is now measured, not asserted.
All 8,082 subjects kept their PCs, and 11 who had no release PCs at all rejoined
the analysis sample.

### 8.9 Results

Jobs 33823211–14 (`02_reml`), 33823242 (20-PC sensitivity), 33793546
(`strat_reml`). Summaries committed at `results/reml_*/reml_summary.tsv`.
**These supplement §2; they are not comparable with it** — see (a).

**(a) Array SNPs halve h². Isolated, not inferred.** Same ancestry definition,
N differing by 3, allele frequencies estimated within EUR in both — only the
variant set differs:

| | n | SNPs | `baseline_thickness` | `global_slope` |
|---|---|---|---|---|
| §2 published EUR | 3,329 | 13.7M imputed | 0.575 ± 0.147 | 0.228 ± 0.141 |
| EUR-only, this GRM | 3,332 | 456k array | 0.246 ± 0.091 | 0.074 ± 0.087 |

Expected (array SNPs tag less causal variation), but it governs everything below.

**CORRECTION (2026-08-20): that table understates the difference by one term.**
The published EUR GRM's build log was overwritten by a later `--pca` run, so its
MAF was unrecorded and the row above implicitly assumed it matched ours. It does
not. Recovered from `abcd_eur.grm.N.bin`, whose per-pair SNP count is
**13,695,250** against the fileset's 13,697,177: **no MAF 0.01 filter was
applied** — the published GRM used essentially every variant down to the
fileset's own MAF 0.001 floor. Ours used MAF 0.01.

So 0.575 vs 0.246 confounds **two** differences, not one: array-vs-imputed
tagging *and* MAF 0.01 vs 0.001. Including 0.001–0.01 variants generally raises
GREML h², because they add tagging and because the 1/(2p(1−p)) standardisation
upweights them. **Phase 4 must therefore run EUR-imputed at MAF 0.001 as well as
0.01**; only the 0.001 arm is like-for-like with §2, and only it can test the
array-vs-imputed explanation cleanly. The `MAC ≥ 10` union fileset (pooled MAF
≥ 0.00043) supports it: a variant at MAF 0.001 within EUR carries ~11 copies, so
it survives the pooled MAC floor.

**(b) §3's forecast fails: precision improved, z did not.**

| phenotype | arm | n | h² | SE | z |
|---|---|---|---|---|---|
| `baseline_thickness` | §2 EUR imputed | 3,329 | 0.575 | 0.147 | 3.91 |
| | **pooled PC-AiR** | **5,649** | **0.269** | **0.056** | **4.77** |
| | pooled PC-AiR + QC | 5,422 | 0.257 | 0.058 | 4.43 |
| | pooled `--grm-cutoff` (94 % EUR) | 4,299 | 0.346 | 0.077 | 4.47 |
| `global_slope` | §2 EUR imputed | 3,329 | 0.228 | 0.141 | 1.62 |
| | **pooled PC-AiR** | **5,649** | **0.076** | **0.054** | **1.42** |
| | pooled PC-AiR + QC | 5,422 | 0.081 | 0.056 | 1.46 |
| | pooled `--grm-cutoff` (94 % EUR) | 4,299 | 0.042 | 0.078 | 0.53 |

SE on `global_slope` fell 0.141 → 0.054 (2.6×, as §3 forecast) but h² fell by the
same factor, so **z went 1.62 → 1.42**. §3's arithmetic held h² fixed while N
grew; on array data it does not. **More N will not make `global_slope`
detectable — the variant set has to go back to imputed.** Update expectations
accordingly; do not quote §3's "z = 2.6".

**(c) New: `slope_PC2` is heritable — a lead, not a result.**

| | §2 EUR (n=3,329) | pooled PC-AiR (n=5,649) | + QC (n=5,422) | 20 PCs (n=5,649) |
|---|---|---|---|---|
| `slope_PC1` | 0.000 ± ~0.14 | 0.137 ± 0.054 | 0.144 ± 0.056 | 0.128 ± 0.055 |
| `slope_PC2` | 0.113 ± 0.141 | **0.179 ± 0.055** (p=4.4e-04) | 0.184 ± 0.057 | 0.186 ± 0.055 |
| `slope_PC3` | 0.000 ± ~0.14 | 0.111 ± 0.054 | 0.138 ± 0.056 | 0.101 ± 0.055 |

`slope_PC2` survives Bonferroni across the five phenotypes — the first phenotype
here besides the positive control to do so. It is significant pooled and **not**
in the EUR-only arm (0.090 ± 0.089), which is also what residual population
structure looks like, so it was tested rather than announced:

- **20 ancestry PCs instead of 10 → 0.186.** Structure-driven inflation shrinks
  under more PCs; this does not (`results/reml_allanc_pc20/`, job 33823242).
- **QC track reproduces it** (0.184) on a different 5,422-subject sample, so it
  is not the §8.5 batch artefact.
- **Caveat that remains:** pooled and EUR-only differ by ≈1 SE, so "the EUR arm
  is underpowered" explains it equally well. **EUR at imputed density settles
  it** — the §8.11 project.

**(d) Stratified (§5.3.2), job 33793546.**

| stratum | n | `baseline_thickness` | `global_slope` |
|---|---|---|---|
| EURlike | 3,968 | 0.365 ± 0.085 | 0.112 ± 0.084 |
| EUR_anchor | 3,159 | 0.310 ± 0.107 | 0.101 ± 0.105 |
| cluster1 | 966 | 0.000 ± 0.231 | 0.000 ± 0.235 |
| cluster2 | 571 | 0.063 ± 0.427 | 0.000 ± 0.453 |
| cluster3 | skipped — 221 unrelated, below `STRAT_MIN_N=300` | | |

**The non-European strata are uninformative, not null** — SEs of 0.23–0.45 span
the parameter space, and a 0.000 there must not be read as evidence of absence.
Pooled (0.269) sits inside the EUR estimates' error bars, so **pooled and
stratified do not materially disagree** and the pooled number is defensible as
the headline *provided* it uses the PC-AiR unrelated set (§8.6). EUR_anchor's
0.310 vs EUR-only's 0.246 is the pooled-allele-frequency effect, measured:
restricting a pooled GRM to one group ≠ building that group's own GRM.

**(e) The batch artefact does not matter.** QC track moves h² by <1 SE
everywhere (0.269 → 0.257) at a cost of 227 subjects. **Use the unfiltered track
as primary; the QC track is the recorded sensitivity.** A defect measured and
found negligible is a result.

**(f) What §8.6 would have cost.** The `--grm-cutoff` arm gives 0.346 vs 0.269
and 0.042 vs 0.076 — about 1 SE apart, in a sample 94 % European instead of
68 %. Not obviously wrong; just quietly a different analysis from the one it
claims to be.

### 8.10 Verified vs assumed

**Verified against a file:** every count in §8.2; the batch/missingness tables in
§8.5; GCTA's pairwise-complete handling (`chr1.grm.N.bin`); every retention and
off-diagonal figure in §8.6; the cluster concordances in §8.7; the PC variance in
§8.8; that all 22 per-chromosome `.grm.id` files are byte-identical; that the
split conserved exactly 502,528 autosomal variants and 11,670 subjects per
chromosome.

**Assumed / inferred, and not checked:** that the four high-missingness batches
are an assay or array-version difference rather than sample quality — the
direction of the effect on the GRM is the same either way, so it was not chased;
that ABCD's "curated" designation implies upstream QC beyond call rate;
that the k-means clusters correspond to conventional continental groupings —
deliberately *not* claimed (§8.7).

### 8.11 Not done, and why

**`06_prs` was not re-run, and should not be re-run on these genotypes.** The
scores would be computed over ~456k array SNPs instead of the 8.78M imputed ones
the published PRS used. The SCZ index SNPs are clumped in `g1000_eur` at imputed
density, so most would simply be absent from the array, and the score would be
built from a small and non-random subset. A weakened or vanished SCZ →
`global_slope` signal would then be uninterpretable: it would be consistent with
"the lead was noise" and equally with "the score lost its SNPs", and those are
not distinguishable after the fact. **That is a worse outcome than not running
it**, because the §2 lead is the result most likely to be over-read.

Powering up the PRS properly needs the imputed 7.0 genotypes, which are GRCh38
(§8.2) against this pipeline's hg19 — so it needs a liftover and a ~890 GB
conversion, or a targeted extraction of the clumped index SNPs by position after
liftover. That is a project, not a step, and it is the obvious next one.

### 8.14 `plink2 --vcf` writes FID = 0 — and the loud failure was luck

Jobs 34058416/34058417 (EUR-only imputed GRMs) died on all 22 tasks with

```
Get 5656 samples from list [.../eur_anchor.keep].
After keeping individuals, 0 subjects remain.
```

`plink2 --vcf` sets **FID = 0** for every sample unless `--double-id` is given,
while the array filesets and every keep list in this project use **FID = IID**.
GCTA identifies an individual by the FID+IID **pair** — §6's ID gotcha, arriving
through a door §6 does not mention.

**The important part is which job did *not* fail.** Those two used `--keep`, so
they errored. The pooled GRM (34058415) has no `--keep`: it was running happily
and would have completed, producing a GRM whose `.grm.id` carried FID = 0 — after
which *every* later keep-based step, including the PC-AiR unrelated subset that
§8.6 exists to enforce, would have silently matched zero subjects. **A loud
failure in the arm that happened to use `--keep` is the only reason this was
caught before the pooled result existed.**

Fixed two ways: `--double-id` added to `work/convert_imputed.sbatch` so a rerun
cannot reintroduce it, and the 22 existing `.fam` files rewritten in place
(FID := IID; the `.fam` does not affect `.bed` layout, and the integrity gate
still passes). Verified after the fix, before resubmitting 66 array tasks:
`eur_anchor.keep` matches 5,656/5,656 FID+IID pairs and `pcair_unrelated.keep`
8,178/8,181 (the 3 are the known non-genotyped).

**Generalisable lesson:** when a fileset changes provenance, check the ID
convention against a keep-list *before* trusting any step that has no `--keep`
to fail for you.

### 8.15 IMPUTED GENOTYPES — `global_slope` is heritable (jobs 34089870–3, 34092188–9)

**Headline: the primary phenotype reaches significance for the first time.**
`global_slope` h² = **0.137 ± 0.046, p = 0.0012** (n = 5,649 unrelated), which
survives Bonferroni across the five phenotypes. It took **both** changes — the
extra N halved the SE, the imputed variants restored the point estimate — and
neither alone was enough, which is why the array run (§8.9b) came back flat.

**The pipeline validates against a known number.** EUR-only at imputed density
with MAF 0.001 — matched to §2 on ancestry, N and frequency floor — gives
`baseline_thickness` **0.474 ± 0.142** against the published **0.575 ± 0.147**,
within ~0.7 SE, on a *different* imputation panel (TOPMed r3) and genome build.

**§8.9(a) is now properly decomposed, and its attribution was half wrong.** Same
children throughout:

| arm | n | `baseline_thickness` | attributable to |
|---|---|---|---|
| array, MAF 0.01 | 3,332 | 0.246 ± 0.091 | — |
| imputed, MAF 0.01 | 3,336 | 0.378 ± 0.107 | **variant set: +0.132** |
| imputed, MAF 0.001 | 3,312 | 0.474 ± 0.142 | **frequency floor: +0.096** |
| §2 published | 3,329 | 0.575 ± 0.147 | different panel/build |

**Both terms are real and roughly equal.** §8.9(a) attributed the whole halving
to the variant set; that was about half right.

**Primary results — imputed, pooled multi-ancestry, PC-AiR unrelated set:**

| phenotype | h² | SE | p | vs array (§8.9b) |
|---|---|---|---|---|
| `baseline_thickness` | 0.246 | 0.049 | 1.5e-07 | 0.269 |
| **`global_slope`** | **0.137** | **0.046** | **0.0012** | 0.076 (p=0.075) |
| `slope_PC2` | 0.161 | 0.046 | 1.6e-04 | 0.179 |
| `slope_PC1` | 0.075 | 0.046 | 0.046 | 0.137 |
| `slope_PC3` | 0.041 | 0.045 | 0.18 | 0.111 |

**Robustness — three checks, all passed:**

| check | `global_slope` | `slope_PC2` |
|---|---|---|
| array-derived PCs, 10 | 0.137 ± 0.046 | 0.161 ± 0.046 |
| imputed-derived PCs, 10 | 0.135 ± 0.047 | 0.165 ± 0.047 |
| imputed-derived PCs, 20 | 0.133 ± 0.047 | 0.158 ± 0.047 |
| drop cluster3 (n = 5,505) | 0.160 ± 0.047 | 0.163 ± 0.047 |

Residual population structure inflates under too few PCs and shrinks under more.
Neither result moves. (The PC-source mismatch was a genuine setup error — the
covariates carried PCs from the *array* GRM while the analysis used the
*imputed* one — caught and corrected here rather than left as an assumption.)

**§8.9(c)'s open caveat is resolved.** On array data `slope_PC2` was significant
pooled (0.179) and not in the EUR arm (0.090 ± 0.089), which is also what a
structure artefact looks like. At imputed density the EUR arm gives
**0.205 ± 0.108**, consistent with the pooled 0.161. So *"the EUR arm was
underpowered"* was the correct explanation, not residual structure. Likewise
`global_slope`: EUR-only gives 0.157 ± 0.107 (MAF 0.01) and 0.151 ± 0.138
(MAF 0.001), both bracketing the pooled 0.137 — **pooling is not inflating it.**

**Answer to "should we drop the poorly-imputed cluster3?" — no.** Dropping its
523 children moves `global_slope` by half an SE (0.137 → 0.160), *strengthening*
it. The group adds a little noise and nothing qualitative depends on it. Keep,
with the sensitivity on record (`results/reml_imp_pooled_noc3/`).

**What still limits this.** The pooled GRM uses pooled allele frequencies — the
compromise `docs/multianc/README.md` §7 flags as unresolved, and the reason no
single GCTA GRM can have ancestry-specific frequencies *and* cross-ancestry
entries. And these numbers are **not** comparable with §2: different panel,
build, variant set and sample. The EUR-MAF-0.001 arm is the only bridge.

### 8.16 GWAS, rg and MAGMA on the imputed genotypes (jobs 34108582, 34108774, 34108739–40, 34109942–3, 34109963, 34109989, 34110020–2)

**Every analysis is run twice: pooled multi-ancestry and EUR-stratified.** That
is not thoroughness for its own sake. LDSC's LD scores, MAGMA's LD panel and the
PRS clumping reference are all `g1000_eur` / `eur_w_ld_chr` — **European**. On a
32 %-non-European sample each carries an untested LD model, and §8.16.2 shows
that mismatch producing a result that looks exactly like a finding.

#### 8.16.1 GWAS — the first genome-wide-significant loci in this project

| phenotype | multi-ancestry (n = 7,932) | | EUR (n = 4,039) | |
|---|---|---|---|---|
| | λ_GC | p<5e-8 | λ_GC | p<5e-8 |
| `baseline_thickness` | 1.107 | **20** | 1.024 | 0 |
| `slope_PC3` | 0.993 | **3** | 0.999 | 1 |
| `global_slope` | 1.013 | 0 | 1.000 | 0 |
| `slope_PC2` | 1.030 | 0 | 1.017 | 0 |
| `slope_PC1` | 0.999 | 0 | 1.007 | 0 |

9,077,609 SNPs tested. The published run (N = 4,119, 8.78M SNPs) found **none**.
The 20 hits are **2 loci** — chr2 ~26.9 Mb (8 SNPs) and chr12 ~69.4 Mb (12) —
the SNP count being LD within loci.

**These are not ancestry artefacts, and the check that settles it is not λ_GC.**
λ_GC rises with N under genuine polygenicity too, so 1.107 → 1.024 at half the N
is consistent with either story. The decisive test is per-SNP: **all 20 hits
appear in the EUR-only scan at p = 4.4e-05 to 7.9e-04 — 20 of 20 below 0.05,
where ~1 is expected.** Had they been driven by between-ancestry structure,
restricting to Europeans would have destroyed them. It does not; they sit at
exactly the strength a true effect shows at half the sample size.

#### 8.16.2 The LDSC intercept was misspecification, not confounding

Worth recording in full, because the intermediate state looked alarming and a
reader stopping there would draw the wrong conclusion.

| phenotype | intercept, multi-ancestry sumstats + **EUR** LD scores | intercept, EUR sumstats + EUR LD scores |
|---|---|---|
| `baseline_thickness` | **1.0995 ± 0.0085** | **0.9927 ± 0.0062** |
| `global_slope` | 1.0283 | 1.0050 |
| `slope_PC2` | 1.0350 | 0.9949 |
| `slope_PC1` | 1.0223 | 1.0018 |
| `slope_PC3` | 1.0284 | 1.0196 |

An LDSC intercept above 1 conventionally indicates confounding, and 1.0995
accounts for essentially all of λ_GC = 1.107. **But when the LD scores match the
sample's ancestry, every intercept sits at 1.00.** The inflation was the EUR LD
scores applied to a multi-ancestry sample — LDSC's intercept absorbs exactly
that kind of model misspecification — not population structure in the GWAS.

**So the multi-ancestry LDSC h² and rg estimates should not be used.** Their
model is wrong in a way the intercept makes visible.

#### 8.16.3 Genetic correlation — still out of reach for the slope, as §3 predicted

| | SCZ × `baseline_thickness` | SCZ × `global_slope` |
|---|---|---|
| §2 published | 0.036 ± 0.047 | not estimable |
| multi-ancestry *(mis-specified, §8.16.2)* | 0.003 ± 0.045, h² z 5.30 | −0.140 ± 0.083, p = 0.092, h² z 2.26 |
| **EUR (properly specified)** | 0.019 ± 0.057, h² z 3.89 | −0.199 ± 0.218, p = 0.36, h² z 0.57 |

**Every slope row is flagged underpowered in both arms**, and §3's arithmetic is
why: an interpretable rg needs h² z ≈ 4, and LDSC estimates h² less precisely
than GREML on a HapMap3 subset. The GREML z of 2.96 does not carry over.
`baseline_thickness` reaches h² z 3.89 in EUR — borderline — and its rg with both
disorders is a clean null. MDD is null throughout.

The one thing worth noting: SCZ × `global_slope` is **negative in both arms**
(−0.140, −0.199), the direction the PRS and the age×PRS interaction also give —
higher SCZ risk, faster thinning. It is not significant in either and should not
be reported as support; it is consistent with the lead, nothing more.

#### 8.16.4 MAGMA — the SCZ locus pool moves onto the slope, and MDD stops being null

Gene sets tested against our phenotypes' gene-level signal:

| set → phenotype | §2 published (4.0 EUR) | **EUR imputed** | multi-ancestry |
|---|---|---|---|
| `SCZ_locus_pool` → `baseline_thickness` | β 0.183, **p 1.9e-04** | β 0.153, **p 1.4e-03** | β 0.079, p 0.051 |
| `SCZ_locus_pool` → `global_slope` | p 0.017 | β 0.127, **p 7.0e-03** | β 0.119, **p 6.5e-03** |
| `SCZ_prioritised` (101 genes) | none | none | none |
| `MDD_hc_finemap` → `baseline_thickness` | null | β 0.174, **p 8.6e-03** | — |
| MDD gene covariate → `baseline_thickness` | null | **p 5.1e-05** | p 1.7e-03 |

Three readings:

1. **The published `baseline_thickness` result replicates** in the ancestry-matched
   arm (0.183 → 0.153) and **attenuates in the multi-ancestry arm** (0.079,
   p = 0.051) — the LD mismatch degrading MAGMA's gene Z-scores, same mechanism
   as §8.16.2.
2. **`global_slope` strengthened in both arms**, from nominal (p = 0.017) to
   p ≈ 7e-03. That is the change the extra N and denser variants bought.
3. **MDD is no longer null** — new, and present in the properly-specified arm.
   §2's "MDD high-confidence sets are null throughout" no longer holds.

The prioritised genes still carry nothing, reproducing §2: the signal is in the
broad locus pool, not the fine-mapped genes.

**AHBA C1–C3: do NOT report the C1− result.** The multi-ancestry arm gives
C1− → `slope_PC1` β = −0.107, p = 5.1e-04, with the same sign across all four
slope phenotypes — which reads as a coherent finding. In the EUR arm the same
test gives β = −0.042, p = 0.20. MAGMA's gene-covariate SE depends on the gene
count (6,940, identical in both arms), not on sample size, so **this is a real
attenuation of the effect, not a loss of power**. C2 and C3 are null everywhere.
The ancestry-matched arm is the one to believe.

#### 8.16.5 PRS — the SCZ lead survives correction, in the arm that is not primary

SCZ → `global_slope` at p<0.5, the threshold where it has always peaked:

| arm | n | β | p | p_adj (Bonferroni × 8) |
|---|---|---|---|---|
| §2 published, full | 4,126 | −0.0397 | 0.0090 | 0.072 |
| **imputed, full** | **8,082** | **−0.0413** | **0.0043** | **0.034** |
| §2 published, EUR | 3,725 | −0.0385 | 0.0149 | 0.119 |
| imputed, EUR | 5,361 | −0.0300 | 0.030 | 0.240 |

**The full-sample result survives Bonferroni across all 8 thresholds for the
first time.** Effect size is unchanged (−0.040 → −0.041), monotone across
thresholds, same direction: higher SCZ risk → faster thinning.

**But `06_prs`'s design makes the EUR arm primary**, because European-discovery
scores transfer poorly across ancestry — in the full sample a result is
ambiguous between a real effect and a transferability artefact. The EUR arm gives
p_adj = 0.24 and is *weaker* than published despite 44 % more children. The two
β's (−0.030 vs −0.041) differ by ~0.6 SE and are compatible, but the full-sample
estimate being the larger one is also what residual ancestry confounding
produces. **Strengthened, not settled.**

**New:** SCZ → `slope_PC1`, β = **+0.0440**, p = 0.0033, p_adj = 0.026 —
*positive*, i.e. opposite in sign to `global_slope`. This is consistent with the
local regional-PRS finding (§2 follow-ups: the SCZ PRS map correlates with slope
PC1, ρ = +0.28, spin p = 0.014). **It appears only in the full sample**, so it
carries the same caveat. All effects are small: r² ≈ 0.17 % of variance.

#### 8.16.6 What to do about the European LD reference

Three analyses depend on it and none of them should. Options, in order of what
was actually done:

- **In-sample LD scores (implemented, `work/insample_ldscores.sbatch`).** We hold
  11,670 genotypes; LD scores computed from the analysis sample describe its LD
  by construction, and 11,670 is >20× the 503 EUR individuals behind
  `eur_w_ld_chr`. This is the fix for LDSC's multi-ancestry arm.
- **In-sample MAGMA panel.** Same logic — MAGMA's `--bfile` only needs a panel
  representing the sample's LD, and the 7.07M-SNP PRS fileset can serve.
- **What neither fixes:** in an admixed sample, admixture creates long-range LD
  that a 1 cM window misses, so in-sample LD scores are themselves biased.
  **cov-LDSC** (Luo et al. 2021 — 20 cM windows plus PC covariates) is the method
  for that and is not in this LDSC build.
- **For rg specifically the problem is deeper still:** SCZ and MDD discovery GWAS
  are themselves European, so a cross-ancestry rg needs a method modelling *two*
  LD structures — **Popcorn** (Brown et al.). Standard LDSC assumes one.
  **The EUR-stratified rg therefore remains the only arm with no ancestry
  mismatch anywhere, and is the one to report.**

#### 8.16.7 Two defects found, both mine, both worth the warning

**`--set-all-var-ids` destroyed every rsID.** `convert_imputed.sbatch` set all
variant IDs to `chr:pos:ref:alt`; the intent was `--set-MISSING-var-ids`. MAGMA,
LDSC and the PRS all match on rsID, so one flag blocked all three — while the GRM
work, which needs no IDs, proceeded looking healthy. Recovered without
re-converting by `work/extract_rsids.sbatch`, which re-reads the VCFs taking only
columns 1–5 and never parses a genotype: 25,950,182 of 25,952,942 mapped
(99.99 %), 0 duplicate IDs, 97.0 % of HapMap3 present. Original `.bim` files kept
in `genotype_imputed/bim_backup_chrposID/`.

**This also removed the liftover requirement.** An rsID names the same variant in
GRCh38 and hg19, and all three tools match on rsID — so the build mismatch §8.2
flagged never bites. That was luck: any position-matched step would have needed a
chain file.

**A race condition in `insample_ldscores.sbatch`.** All 22 array tasks wrote the
same `hm3.snps` path; five read it mid-truncation and died on "0 variants
remaining". They failed loudly, but **the same race could have delivered a
partial list and produced LD scores over a silently truncated SNP set**, which
nothing downstream would flag. Verified after the fact that the surviving
chromosomes' SNP counts scale correctly with chromosome size, so none were
truncated. Now writes a per-chromosome list.

### 8.13 Multi-ancestry diagnostics — see the dedicated report

**[`docs/multianc/README.md`](../docs/multianc/README.md)** is the full
methodological record behind §8.5–§8.7 and §8.12: why each pre-GRM check exists,
what it found, every threshold we chose and on what basis, a k-sensitivity sweep
with figures, and a comparison against published ABCD genetics methodology.
**Read it before repeating or extending this analysis.** In brief:

- **The `--grm-cutoff` finding does not depend on k.** At k=3/4/5 it retains
  75.8/75.6/75.4 % of the European-like cluster and 0–33 % of every other,
  leaving an unrelated set that is 93.6/93.7/90.8 % European. Expected rather
  than lucky: the failure is a property of the GRM's pooled allele frequencies,
  not of any partition. The strata are a lens for seeing it.
- **The clusters are validated, not asserted.** Blind to self-report, k=4 yields
  clusters that are 91 % White, 80 % Black, 95 % Hispanic and 45 % Asian.
- **k=4 is a judgement, and k=5 is defensible.** k=3 merges the Hispanic and
  Asian groups; k=5 isolates an admixed cluster and thereby raises `cluster1`
  from 79.6 % to 94.1 % Black — better purity, but its 221 unrelated children sit
  below `STRAT_MIN_N`. `strata_k2..k6.tsv` are all written; **prefer k=5 if you
  need cleaner group definitions rather than more fittable strata.**
- **We are not the first here.** PC-AiR is what the ABCD release itself uses, and
  `FastSparseGRM` exists because multi-ancestry GRMs are confounded with
  population structure. Two published choices are **better** than ours and are
  named as such: 1000 Genomes projection for ancestry labels (no multi-ancestry
  panel on this account) and PC-Relate as the GRM (we use GCTA for comparability
  with §2, which is a trade, not a win).
- **The unresolved compromise:** the pooled GRM still uses pooled allele
  frequencies. No single GCTA GRM can have ancestry-specific frequencies *and*
  cross-ancestry entries. State this wherever the pooled h² is reported.

### 8.12 Phase 0 of the imputed project — sizing it (job 33838952)

Before committing to the §8.11 project, `work/probe_imputed.sbatch` counted what
actually survives filtering on chr22. The raw variant count is not the relevant
number; the surviving one is.

| chr22 | variants | of total |
|---|---|---|
| total imputed | 5,789,636 | — |
| R² ≥ 0.8 | 1,019,601 | 17.6 % |
| **MAF ≥ 0.01** | **144,524** | **2.5 %** |
| R² ≥ 0.8 **and** MAF ≥ 0.01 | 135,050 | 2.3 % |

Scaling by chr22's share of the autosome (1.642 %, from the array variant
counts) projects **~8.2M genome-wide** — against the 8.78M the published EUR
GWAS used, which is a useful independent check that the filter is sane. As
PLINK that is **~24 GB**, ~31 GB including GRMs, against 184 GB free. **The
project fits; it does not need the ~890 GB the raw VCFs suggest**, provided the
filters are applied in the same pass as the conversion so the unfiltered form
never lands on disk.

**The finding that redirects the risk.** MAF, not R², does nearly all the
filtering: 97.5 % of TOPMed variants are too rare to use here, and of those that
clear MAF 0.01, **93.4 % already clear R² 0.8**. So ancestry-differential
*imputation quality* — the thing §8.11 was most worried about — is barely
binding among common variants.

What is binding is the MAF filter, and the `MAF` in the INFO field is the
**pooled whole-cohort frequency**. A variant at 4 % in one ancestry group and
absent elsewhere is ~0.8 % pooled and is discarded, despite being common in the
children it is informative for. **That is the §8.6 problem for the third time:
a single pooled number standing in for a quantity that differs by ancestry.**
`work/impqual_by_ancestry.sbatch` measures per-group MAF and per-group dosage
ambiguity so the threshold is chosen on evidence rather than convention.

Cost note: the chr22 census took 57 minutes single-threaded (12.6 GB of gzip);
chr1 at 69 GB scales to ~5 h. Phase 2 should use
`plink2 --extract-if-info` rather than an awk pass.

**`03_gwas` / `05_ldsc` were not re-run** either: the sparse GRM for fastGWA
exists (`abcd_all_sp`) but a GWAS on 456k array SNPs is not comparable with the
published 8.78M-SNP scan, and §3 already establishes that rg stays out of reach
in ABCD regardless of this N.

