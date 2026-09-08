# hpc_v3 — setup context for a multivariate GWAS of adolescent cortical thinning

**Written 2026-09-08 for an agent starting fresh.** This is not a pipeline; it is
the context you need before writing one. `hpc/` (v1) and `hpc_v2/` are large and
contain a lot of superseded back-and-forth — you do not need to read them. The
findings from them that bear on *this* task are below.

Everything numeric here was recomputed from committed result files, not quoted
from prose. Where a number is uncertain or unverified, it says so.

---

## 1. Read this first: the premise needs correcting

The suggestion that prompted `hpc_v3` was to use **MOSTest** to obtain a
"higher-heritability multivariate phenotype". Half of that is right and half is
a category error, and getting it wrong will waste weeks.

**MOSTest does not produce a phenotype.** It is a multivariate *association
test*. For each SNP it computes one omnibus p-value against the null "this SNP
is associated with none of the input measures", by combining the per-measure
z-scores through a regularised inverse of their null correlation matrix, with
the tail calibrated by permutation on individual-level data. What comes out is:

- one **unsigned p-value per SNP**;
- **no** scalar per-individual phenotype;
- **no** signed effect size, and therefore **no h²** for "the MOSTest phenotype".

That matters because it determines which downstream analyses survive:

| downstream | works on MOSTest output? | why |
|---|---|---|
| Locus discovery | **yes** — this is the point | omnibus p per SNP is all it needs |
| MAGMA (gene, gene-set) | **yes** | MAGMA consumes unsigned SNP p-values |
| PRS (disorder → ABCD) | **unaffected** | runs the other direction; MOSTest is irrelevant to it |
| LDSC genetic correlation | **no** | needs signed z for one scalar trait |
| LAVA local rg | **no** | needs signed sumstats + per-trait local h² |

So MOSTest can strengthen the discovery and gene-set arms and **cannot** feed
the genetic-correlation arms that currently carry the disorder hypothesis.

## 2. Why it is still worth doing — the honest case for

**The global mean is a genuinely poor summary of the regional slopes.** Computed
from `out/thickness_dsk_70_139406217085/fits/blups.parquet` (8,192 subjects ×
68 Desikan regions, no missing cells):

- median pairwise correlation between regional slopes **r = 0.186** (5–95%: 0.036–0.449)
- **PC1 explains only 24.1%** of the variance; PC1–5 40.6%; PC1–10 52.4%
- effective dimensionality **≈ 47 (Li & Ji) to 52 (Galwey)** of 68

`global_slope` is the unweighted mean of those 68 slopes, so it captures roughly
the PC1 direction and discards most of the multivariate structure. That is
precisely the regime MOSTest was designed for — distributed effects across many
imperfectly correlated measures — and it is a stronger argument for a
multivariate approach than anything in the disorder analyses.

Supporting evidence already in hand: **`slope_PC2` has LDSC h² z = 3.46**
(h² = 0.377 ± 0.109) against `global_slope`'s **z = 1.16** (0.124 ± 0.107). A
plain PCA of the same regional slopes already yields a 3× more heritable
phenotype. The premise that a multivariate combination beats the mean on
heritability is therefore **already demonstrated** — no new method required to
establish it.

And discovery is currently **zero**: 0 genome-wide significant variants for
`global_slope` in both v1 (fastGWA, min p 3.16e-07) and v2 (GENESIS, 5.11e-07).
MOSTest is the most credible route to any locus at all.

## 3. Why it does not solve the problem as stated — the honest case against

**Heritability and disorder-relevance are anti-ranked across the phenotypes we
already have.** From `hpc_v2/work/results_v2/ldsc_eur/ldsc_rg_summary.tsv` and
the v2 PRS tables (EUR arm, C+T p<0.5):

| phenotype | LDSC h² z | rg with SCZ | SCZ PRS p |
|---|---|---|---|
| `baseline_thickness` | **4.49** | +0.026 | 0.576 |
| `slope_PC2` | **3.46** | −0.052 | 0.233 |
| `global_slope` | 1.16 | **−0.149** | **0.0039** |
| `slope_PC1` | 0.66 | +0.094 | 0.086 |
| `slope_PC3` | −2.55 | n/a | 0.334 |

Spearman(h² z, |rg|) = **−0.80** and Spearman(h² z, PRS −log₁₀p) = **−0.80**
over the four positively-heritable phenotypes. With n = 4 that is not a
hypothesis test, and it must not be reported as one — but the direction is the
whole point: **the two most heritable phenotypes have essentially no genetic
overlap with schizophrenia, and the phenotype carrying the SCZ signal is the
least heritable of the four.**

MOSTest maximises multivariate association with genotype. On this evidence that
optimises a direction the existing data suggest is close to orthogonal to the
disorder hypothesis. Expect a more heritable, better-powered readout of
*something*, with no guarantee — and some contrary evidence — that it is the
something you care about.

**Relatedness is a second, concrete problem.** MOSTest as published fits
per-measure linear regressions on individual-level data and calibrates by
permutation; it has **no relatedness model**. ABCD is sibling-enriched
(1,339 sib pairs in the pooled genotyped-phenotyped sample, 686 in the EUR
arm), and the whole of `hpc_v2` exists to exploit that rather than discard it.
Running MOSTest naively would either reintroduce inflation or force you back
onto an unrelated subset. Two mitigations, in order of preference:

1. **Permute whole families, not individuals** (this is what the v2 threshold
   work did, and it worked: see §5). This fixes the *null calibration* but not
   the inflation of the per-measure z-scores themselves.
2. **Use a method that accepts summary statistics from a relatedness-aware
   GWAS.** `JAGWAS` (2024) is explicitly built for many phenotypes on
   *possibly related* individuals and works from per-phenotype summary
   statistics — meaning it could consume `hpc_v2`'s GENESIS per-region output
   directly and inherit the sibling-aware design for free. **This is probably
   the better first thing to try, and it has not been evaluated. Verify its
   assumptions and current implementation before committing to it.**

## 4. Alternatives that target the stated goal more directly

If the aim remains *the disorder link* rather than *locus discovery*, these
preserve the signed scalar trait that LDSC and LAVA require:

- **MTAG** — boosts power for one primary trait (`global_slope`) by conditioning
  on genetically correlated secondary traits (the regional slopes or the PCs).
  Outputs **signed summary statistics**, so LDSC, LAVA and everything else works
  unchanged. If you want a better-powered version of the *same* trait, this is
  the tool, not MOSTest.
- **Genomic SEM** — model the genetic covariance among the regional slopes *and*
  SCZ jointly, and define a genetic factor by its relationship to SCZ rather
  than by its own heritability. This directly addresses the anti-ranking in §3:
  it optimises for the covariance you care about instead of for h².
- **CPC** (shipped in the MOSTest package) — closed-form null, no permutation.
  Cheaper than MOSTest and equivalent to it when the measure covariance is the
  identity. Worth knowing about if permutation cost becomes the bottleneck.

**Recommended framing for `hpc_v3`:** run MOSTest (or JAGWAS) as a
**discovery arm answering a different question** — *which loci shape adolescent
cortical thinning?* — with MAGMA downstream, and report it as such. Do **not**
present it as a fix for the schizophrenia analysis, and do not attempt LDSC rg
or LAVA on its output. If the disorder link is the priority, do MTAG or genomic
SEM first.

## 5. Inputs, and what already exists

**Regional slopes (the MOSTest phenotype matrix)** —
`out/thickness_dsk_70_139406217085/fits/blups.parquet`: long format, columns
`subject`, `label`, `re_slope`, `se_re_slope`; 8,192 subjects × 68 regions,
complete. These are per-region random-effect slopes (BLUPs) from
`lmer(value ~ age_c + sex + (1 + age_c | subject) + (1 | site))` fitted per
region. `model_table.parquet` in the same directory has the underlying
observations (23,409 subject-visits) if you need to refit.

**Note the BLUP caveat.** These are shrunken estimates, not raw slopes. Shrinkage
differs by region (regions with noisier data shrink more), which induces
heteroscedasticity across the 68 measures. MOSTest rank-inverse-normalises each
measure, which helps, but **the differential shrinkage has not been assessed as
a source of artefact in a multivariate test.** Consider running on raw
per-region OLS slopes as a sensitivity check. The design is 100% balanced (every
subject-visit has all 68 regions), which makes raw slopes well defined.

**Genotypes** — hpc_v2's `config.local.sh.example` documents the CSD3 layout:
array fileset `~11,670 × 515,228` (hg19) and per-chromosome imputed filesets.
Kinship and PC-AiR PCs from hpc_v2 step 02 are in
`hpc_v2/work/results_v2/kinship/` (`kinship_sparse.rds`, `pcair_pcs.tsv`,
`pcair_unrelated.txt`). **Reuse these rather than recomputing** — v1's pooled
GCTA GRM with a relatedness cutoff returned a 94%-European "unrelated" set while
reporting nothing wrong, which is why v2 replaced it with KING → PC-AiR →
PC-Relate.

**Per-region GWAS summary statistics do not exist yet.** hpc_v2 ran GENESIS on
5 phenotypes, not 68. If you take the JAGWAS route you will need to run the v2
association step over the 68 regional slopes first — that is 68 × 22 tasks, and
the v2 `04_assoc.sbatch` array pattern extends to it directly.

**Software** — MOSTest is MATLAB (original release) with a permutation scheme on
individual-level data; there are third-party reimplementations of varying
fidelity. **Confirm which implementation you are using and validate it on a
subset before trusting genome-wide output.** Do not assume a Python port is
equivalent to the published method.

## 6. Findings from v1/v2 you must not re-derive or contradict

- **Heritability.** `global_slope` SNP h² is **0.183 ± 0.049** (Zaitlen two-GRM,
  relatives kept, n=8,082) but **unstable to the kinship sparsity threshold**
  (0.115 ± 0.056 at the 2nd-degree threshold — a 60% swing). Pedigree h² is
  stable at **0.252–0.273 ± 0.033**. `baseline_thickness` is stable throughout
  (SNP h² ≈ 0.54, pedigree 0.846). Keeping relatives raised the point estimate
  but **did not improve the standard error** (0.046 → 0.049).
- **Calibration.** v1 fastGWA reported 20 genome-wide hits for
  `baseline_thickness` at λ_GC = 1.107; v2 GENESIS gives λ = 1.036 and **zero**.
  They were stratification. Any new method must report λ_GC, and a hit count
  without it is not interpretable.
- **`slope_PC3` has negative LDSC heritability** (−0.240 ± 0.094, z = −2.55). It
  should not be treated as a heritable trait, and any method that reports signal
  for it is telling you about its own calibration, not about biology.
- **LAVA local rg is not usable as it stands.** Its bivariate test only runs
  where local h² clears an *uncorrected* nominal p<0.05 in both traits, and that
  gate passes 53–58% of blocks where 5% is expected — including 54% for
  `slope_PC3`, which has no genome-wide heritability at all. Every genome-wide
  FDR survivor has |ρ| ≥ 0.62 with 6 of 7 confidence intervals running into
  ρ = ±1, and within one block the sign flips across correlated phenotypes. Treat
  all existing LAVA output as provisional pending permutation calibration.
- **Threshold multiplicity.** For nested C+T PRS thresholds, the correct
  correction is **family-level permutation of the min-p statistic**, not
  Bonferroni. Measured on the array run: 8 nested scores carry ~4.1 effective
  tests, not 8, and the verdict is correction-dependent (survives permutation at
  0.035, fails Bonferroni at 0.068). Use permutation, or avoid the multiplicity
  entirely with a continuous-shrinkage score.
- **A negative control is not clean.** The Alzheimer's PRS is itself associated
  with `global_slope` (C+T p<0.5, EUR: β = −0.0366, p = 0.020), same sign as SCZ
  (−0.0460, p = 0.0039) and ~80% its magnitude. It *does* attenuate within
  families, which is the confounding signature; SCZ and MDD do not. So a
  non-specific confound is present in the population arm, and **any new
  phenotype should be tested against ASD and Alzheimer's controls from the
  start**, not at the end.
- **ID convention.** Subject IDs are spelled differently across the genotype
  filesets and the phenotype exports; join on the normalised 8-char NDAR token.
  The genotype `.fam` has FID = IID, so **real family IDs must come from the
  phenotype export**, never from the `.fam`.

## 7. Data-use rule (non-negotiable)

This repository is public. **Per-subject files derived from genotypes —
`.profile`, `.sscore`, `pcair_pcs.tsv`, and any per-individual score or
component — are never committed.** Copy them to `scratch/` (gitignored) for
local work; only summary-level tables go into the repo. The regional slope
BLUPs are phenotype-derived and already tracked under `out/`, so they are not
affected by this rule.

## 8. Definition of done for a first `hpc_v3` pass

1. A decision recorded on **MOSTest vs JAGWAS vs MTAG**, with the relatedness
   handling stated explicitly and justified.
2. Multivariate association run over the 68 regional slopes, reporting
   **λ_GC or the permutation-calibrated equivalent**, the locus count, and the
   number of loci that replicate the univariate `global_slope` signal.
3. MAGMA gene and gene-set analysis on the output, including the ASD and
   Alzheimer's comparators.
4. An explicit statement of what the result does **not** license — specifically
   that no genetic correlation or local-rg claim can be made from an omnibus
   unsigned statistic.
5. A sensitivity check on **raw per-region slopes vs BLUPs** (§5).
