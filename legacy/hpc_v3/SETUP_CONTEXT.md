# hpc_v3 — setup context for borrowing power across the regional thinning phenotypes

**Written 2026-09-08. For an agent starting fresh.** This is not a pipeline; it
is the context you need before writing one. `hpc/` (v1) and `hpc_v2/` are large
and contain a lot of superseded back-and-forth — you do not need to read them.
The findings that bear on *this* task are below.

Every number here was recomputed from committed result files or from the
phenotype BLUPs. Where something is unverified or assumed, it says so.

---

## 1. The problem this directory exists to solve

`global_slope` — the rate of cortical thinning across adolescence — is the
project's target phenotype, and it is **underpowered rather than uninformative**.
Its SNP heritability is real but barely distinguishable from zero:
LDSC h² = 0.124 ± 0.107, **z = 1.16**. Every downstream analysis inherits that
ceiling. Zero variants reach genome-wide significance (best p = 5.1e-07 in v2),
and the schizophrenia association, while nominally significant
(β = −0.046, p = 0.004), sits on a phenotype whose heritability we cannot
confidently bound above zero.

More children would fix this. Absent that, the remaining lever is to **borrow
statistical strength from the 68 individual regional slopes**, which were
measured in the same children and are currently collapsed into a single
unweighted average.

## 2. Scope decision: MOSTest and JAGWAS are OUT

This directory was opened to evaluate **MOSTest**, and an earlier draft of this
document was written around MOSTest and **JAGWAS**. Both are now ruled out, for
the same reason. Read this before reinstating either.

**The goal is the disorder link, not locus discovery.** The question is whether
genetic risk for schizophrenia and depression relates to thinning *rate*. That
needs **signed, trait-specific** effect estimates: you must know a variant
pushes thinning *faster* rather than merely that it does *something*, because
the analysis correlates a direction in the brain phenotype with a direction in
disorder risk.

**MOSTest and JAGWAS are the same statistic and cannot supply that.** Both
compute, per variant, `T = z' R⁻¹ z` — the vector of per-phenotype z-scores
squared up through the inverse phenotype correlation matrix — tested against a
chi-square null with k degrees of freedom (JAGWAS analytically; MOSTest by
permutation). **Squaring destroys the sign.** The output is one *unsigned*
p-value per variant against the null "associated with none of the k phenotypes":
no scalar per-individual phenotype, no signed effect, no h².

Which downstream analyses survive that:

| arm | omnibus test (MOSTest/JAGWAS) | MTAG / genomic SEM |
|---|---|---|
| Locus discovery | **yes** — the point | yes |
| MAGMA gene / gene-set | **yes** (unsigned p is enough) | yes |
| PRS (disorder → ABCD) | unaffected — runs the other direction | unaffected |
| LDSC genetic correlation | **no** — needs signed z for one trait | **yes** |
| LAVA local rg | **no** — needs signed sumstats | **yes** |

The JAGWAS authors draw this contrast themselves: they group MOSTest, MultiABEL
and metaUSAT as computing multivariate p-values via MANOVA-F or chi-square
tests, and note that MTAG instead does *not* test the multivariate null but
leverages evidence from related traits to boost a primary trait.

For completeness, because it is genuinely good at what it does: JAGWAS reports
6× more replicated loci than single-phenotype GWAS on UK Biobank imaging
phenotypes, runs 128 phenotypes × 8.1M SNPs in ~20 minutes, and — unlike
MOSTest — **inherits relatedness handling from whatever mixed-model GWAS
produced its input summary statistics**, which would have let it consume
`hpc_v2`'s GENESIS per-region output and keep ABCD's sibling design. That
property is real and was why JAGWAS was originally preferred over MOSTest. It
does not rescue either method here, because the limitation is the *output type*.
**If the priority ever shifts to locus discovery, JAGWAS is the right tool.**

## 3. The two methods that do fit

### MTAG (Multi-Trait Analysis of GWAS; Turley et al. 2018)

**The idea, in plain terms.** You have one trait you care about and several
others measured in the same people that are genetically correlated with it. If a
variant affects thinning in the frontal lobe it probably also affects the
parietal lobe, so evidence about the neighbours is partial evidence about your
trait. MTAG makes that precise: it takes ordinary per-trait GWAS summary
statistics plus the genetic covariance matrix between traits (estimated by
LDSC), and returns a **revised, signed estimate for each trait** in which each
variant's effect has been shrunk toward what the correlated traits imply. It is
an empirical-Bayes/shrinkage step, not a new test.

The analogy: estimating a student's maths ability from their maths exam alone,
versus also using their physics and chemistry scores. You still report *maths
ability* — same units, same meaning — you just estimate it more precisely.

**Why it fits this project.** The output is a set of summary statistics for
`global_slope` itself — signed betas and standard errors, same interpretation as
before, just better powered. LDSC, LAVA, and the PRS comparison all keep working
unchanged. This is the tool if you want **a better-powered version of the same
scalar trait**.

**Preconditions and risks — read these before running it:**

1. **MTAG needs well-estimated genetic covariances, and ours are marginal.**
   It obtains them from LDSC, whose reliability depends on per-trait h² z. The
   MTAG documentation advises against traits with low h² z (guidance in the
   region of z > 4; **verify the current recommendation**). `global_slope` is at
   **1.16**, and individual regional slopes will be *noisier still* than their
   own average. This is the single biggest threat to the whole approach and
   should be tested before any production run — see §5, step 1.
2. **Homogeneity assumption.** MTAG assumes the variant-level genetic
   covariance between traits is the same for all variants. Where that fails, the
   per-trait estimates can be biased. MTAG ships a **`maxFDR`** diagnostic that
   bounds the resulting false-discovery inflation; report it, do not skip it.
3. **The "effective N" is not a real sample size.** MTAG output has inflated
   apparent power by construction. Do not describe an MTAG result as a GWAS of
   N children, and do not feed MTAG's effective N into anything that assumes an
   independent sample.
4. **Circularity risk in the rg step.** If you MTAG `global_slope` using the
   regional slopes and then compute rg between the boosted `global_slope` and
   schizophrenia, that is legitimate — schizophrenia was not used in the boost.
   But do **not** boost using any trait that itself involved the disorder
   statistics, and state explicitly in any write-up which traits informed the boost.

### Genomic SEM (Grotzinger et al. 2019)

**The idea, in plain terms.** Instead of boosting one trait, model the whole
genetic covariance structure at once. LDSC gives you a matrix of genetic
covariances among all your traits; genomic SEM fits a **factor model** to that
matrix, exactly as classical factor analysis fits one to a phenotypic
correlation matrix. You get latent genetic factors, and you can then run a GWAS
*on a factor* and obtain signed summary statistics for it.

**Why it may fit this project better than MTAG.** You can include schizophrenia
in the model and define a factor **by its relationship to the disorder rather
than by its own heritability**. That directly attacks the central obstacle in
§4: the most heritable phenotypes here are not the disorder-relevant ones, so
any method that simply maximises heritability optimises the wrong direction.
Genomic SEM lets you specify the direction you care about.

**Costs.** More modelling choices to defend (factor structure, fit indices,
identification), and it inherits the same dependence on well-estimated LDSC
covariances as MTAG — so §3's precondition 1 applies here too.

### How to choose

| | MTAG | genomic SEM |
|---|---|---|
| What you get | signed sumstats for **your existing trait** | signed sumstats for a **latent factor you define** |
| Optimises for | precision of that trait | whatever covariance you specify, incl. with SCZ |
| Modelling burden | low — nearly turnkey | substantial |
| Best first move if… | you trust `global_slope` and just want power | you suspect `global_slope` is the wrong summary |

**Recommendation: MTAG first**, because it is cheap, it changes nothing about
interpretation, and its failure mode is informative — if MTAG cannot improve
`global_slope`'s precision, that tells you the regional slopes carry no
independent genetic signal, which is itself a finding worth having. Then genomic
SEM if MTAG works and you want to reshape the phenotype.

## 4. Feasibility evidence from these data — read before committing

Three facts from this dataset bear directly on whether power-borrowing can work.
They do not all point the same way.

**(a) In favour: the global mean discards a lot.** From
`out/thickness_dsk_70_139406217085/fits/blups.parquet` (8,192 subjects × 68
Desikan regions, no missing cells): median pairwise correlation between regional
slopes **r = 0.186** (5–95%: 0.036–0.449), **PC1 explains only 24.1%** of the
variance, PC1–10 just 52.4%, and effective dimensionality is **≈ 47 (Li & Ji) to
52 (Galwey)** of 68. There is substantial independent regional signal that a
single average cannot represent. Eigenspectrum: `regional_slope_eigen.csv`.

**(b) Against, and important: `global_slope` and `slope_PC1` are the same
trait.** Recomputed from the real BLUPs, their correlation is **−0.991**, while
`global_slope`'s correlation with PC2 and PC3 is only −0.077 and +0.083
(`phenotypic_corr_real.csv`). This is expected — all 68 inter-regional
correlations are positive, so PC1 has same-sign loadings and is nearly
collinear with the equally-weighted mean — but it has two consequences:

- The existing five phenotypes contain **essentially one informative axis**
  (`global_slope` ≈ `slope_PC1`) plus two near-orthogonal ones. So do **not**
  run MTAG on the five existing phenotypes: PC1 has nothing to lend
  (it is the same trait) and PC2/PC3 are nearly orthogonal to the target.
  **The 68 individual regional slopes are the right MTAG input.**
- Any v1/v2 result that treats `global_slope` and `slope_PC1` as separate
  findings is double-counting one trait. Their apparently different SCZ
  associations (`global_slope` β = −0.046, p = 0.004; `slope_PC1` β = +0.028,
  p = 0.086) are one association with a sign flip from PC1's negative loading.
  *Caveat:* the PC1 above was recomputed here with its own standardisation and
  may not be identical to the pipeline's `slope_PC1`; the magnitude difference
  suggests it is not. **Confirm against the pipeline's own PC scores before
  relying on this.**

**(c) Against: heritability and disorder-relevance are anti-ranked.** From
`hpc_v2/work/results_v2/ldsc_eur/ldsc_rg_summary.tsv` and the v2 PRS tables
(EUR arm, C+T p<0.5) — see `h2_vs_scz_relevance.csv`:

| phenotype | LDSC h² z | rg with SCZ | SCZ PRS p |
|---|---|---|---|
| `baseline_thickness` | **4.49** | +0.026 | 0.576 |
| `slope_PC2` | **3.46** | −0.052 | 0.233 |
| `global_slope` | 1.16 | **−0.149** | **0.0039** |
| `slope_PC1` | 0.66 | +0.094 | 0.086 |
| `slope_PC3` | −2.55 | n/a | 0.334 |

Spearman(h² z, |rg|) = **−0.80** and Spearman(h² z, PRS −log₁₀p) = **−0.80**
over the four positively-heritable phenotypes. **With n = 4 this is not a
hypothesis test and must never be reported as one** — but the direction is the
point: the two most heritable phenotypes have essentially no genetic overlap
with schizophrenia, while the phenotype carrying the SCZ signal is the least
heritable of the four. Any method that maximises heritability alone optimises a
direction these data suggest is close to orthogonal to the hypothesis. This is
the argument for genomic SEM over a purely power-maximising approach — and the
reason to keep the disorder association, not h², as the success criterion.

## 5. Required first steps, in order

1. **Compute the genetic correlation matrix among the ABCD phenotypes.**
   This does not exist — `ldsc_rg_summary.tsv` is phenotype × *disorder* only.
   MTAG and genomic SEM both require phenotype × phenotype genetic covariances.
   Run LDSC rg across `global_slope`, `baseline_thickness`, the PCs, and a pilot
   set of regional slopes. **This is also the go/no-go test:** if per-trait h² z
   is too low for LDSC to return stable rg estimates, both methods are blocked
   and you should say so rather than proceeding on unstable covariances.
2. **Run per-region GWAS.** No per-region summary statistics exist; `hpc_v2` ran
   5 phenotypes, not 68. The v2 `04_assoc.sbatch` array pattern extends directly
   — 68 × 22 tasks. Reuse the v2 kinship and PC-AiR files (§6) so the sibling
   design is preserved and MTAG's inputs are relatedness-aware from the start.
3. **Pilot MTAG on a subset** — e.g. `global_slope` primary with 10 regional
   slopes — and check `maxFDR` and whether the standard error on `global_slope`
   actually falls. Only scale to 68 if it does.
4. **Sensitivity check: raw per-region OLS slopes vs BLUPs** (§6).

## 6. Inputs, and what already exists

**Regional slopes (the phenotype matrix)** —
`out/thickness_dsk_70_139406217085/fits/blups.parquet`: long format, columns
`label`, `subject`, `re_intercept`, `re_slope`, `se_re_intercept`,
`se_re_slope`; 8,192 subjects × 68 regions, complete. These are per-region
random-effect slopes from
`lmer(value ~ age_c + sex + (1 + age_c | subject) + (1 | site))` fitted per
region. `model_table.parquet` alongside has the 23,409 subject-visit
observations if you need to refit.

**BLUP caveat.** These are shrunken estimates, not raw slopes, and shrinkage
differs by region (noisier regions shrink more), inducing heteroscedasticity
across the 68 measures. **This has not been assessed as an artefact source in a
multi-trait analysis.** The design is 100% balanced (every subject-visit has all
68 regions), so raw per-region OLS slopes are well defined and make a clean
sensitivity check.

**Genotypes and relatedness** — `hpc_v2/config.local.sh.example` documents the
CSD3 layout: array fileset ~11,670 × 515,228 (hg19) plus per-chromosome imputed
filesets. Kinship and ancestry axes from `hpc_v2` step 02 are in
`hpc_v2/work/results_v2/kinship/` (`kinship_sparse.rds`, `pcair_pcs.tsv`,
`pcair_unrelated.txt`). **Reuse these; do not recompute.** v1's pooled GCTA GRM
with a relatedness cutoff returned a 94%-European "unrelated" set while
reporting nothing wrong, which is why v2 replaced it with KING → PC-AiR →
PC-Relate.

**Software** — MTAG is a Python tool from the Turley et al. release; genomic SEM
is an R package. Both are actively maintained but **confirm current versions and
argument syntax rather than assuming; validate on a subset before trusting
genome-wide output.**

## 7. Findings from v1/v2 you must not re-derive or contradict

- **Heritability.** `global_slope` SNP h² = **0.183 ± 0.049** (Zaitlen two-GRM,
  relatives kept, n=8,082) but **unstable to the kinship sparsity threshold**
  (0.115 ± 0.056 at the 2nd-degree threshold — a 60% swing). Pedigree h² is
  stable at **0.252–0.273 ± 0.033**. `baseline_thickness` is stable throughout
  (SNP h² ≈ 0.54, pedigree 0.846). Keeping relatives raised the point estimate
  but **did not improve the standard error** (0.046 → 0.049).
- **Calibration.** v1 fastGWA reported 20 genome-wide hits for
  `baseline_thickness` at λ_GC = 1.107; v2 GENESIS gives λ = 1.036 and **zero**.
  They were stratification. Report λ_GC for anything new — a hit count without
  it is not interpretable.
- **`slope_PC3` has negative LDSC heritability** (−0.240 ± 0.094, z = −2.55). It
  is not a heritable trait, and a method reporting signal for it is telling you
  about its own calibration. It is a useful negative control for exactly that.
- **LAVA local rg is not usable as it stands.** Its bivariate test only runs
  where local h² clears an *uncorrected* nominal p<0.05 in both traits, and that
  gate passes 53–58% of blocks where 5% is expected — including 54% for
  `slope_PC3`. Every genome-wide FDR survivor has |ρ| ≥ 0.62 with 6 of 7
  confidence intervals running into ρ = ±1, and within one block the sign flips
  across correlated phenotypes. Treat all existing LAVA output as provisional
  pending permutation calibration.
- **Threshold multiplicity.** For nested C+T PRS thresholds the correct
  correction is **family-level permutation of the min-p statistic**, not
  Bonferroni: 8 nested scores carry ~4.1 effective tests, and the verdict is
  correction-dependent (survives permutation at 0.035, fails Bonferroni at
  0.068). Or avoid the multiplicity with a continuous-shrinkage score.
- **A negative control is not clean.** The Alzheimer's PRS is itself associated
  with `global_slope` (C+T p<0.5, EUR: β = −0.0366, p = 0.020), same sign as SCZ
  (−0.0460, p = 0.0039) and ~80% its magnitude. It *does* attenuate within
  families, which is the confounding signature; SCZ and MDD do not. A
  non-specific confound is present in the population arm, so **test any new
  phenotype against the ASD and Alzheimer's controls from the start.**
- **ID convention.** Subject IDs are spelled differently across the genotype
  filesets and the phenotype exports; join on the normalised 8-char NDAR token.
  The genotype `.fam` has FID = IID, so **real family IDs must come from the
  phenotype export**, never from the `.fam`.
- **Do not read `scratch/hpc_test/` or `scratch/hpc_v2_test/` as data.** They
  are synthetic test fixtures with deliberately planted effects, used by the
  local smoke tests. A phenotype file there will load cleanly and give
  plausible-looking wrong answers.

## 8. Data-use rule (non-negotiable)

This repository is public. **Per-subject files derived from genotypes —
`.profile`, `.sscore`, `pcair_pcs.tsv`, and any per-individual score or
component — are never committed.** Copy them to `scratch/` (gitignored) for
local work; only summary-level tables go into the repo. The regional slope BLUPs
are phenotype-derived and already tracked under `out/`, so are unaffected.

## 9. Definition of done for a first `hpc_v3` pass

1. The phenotype × phenotype genetic correlation matrix, with an explicit
   **go/no-go verdict** on whether per-trait h² z supports MTAG at all (§5.1).
2. Per-region GWAS summary statistics for the 68 regional slopes, produced with
   the v2 relatedness-aware pipeline, with λ_GC reported per region.
3. An MTAG run on `global_slope`, reporting: the standard error before and after
   (the whole point is whether it falls), `maxFDR`, and which traits informed
   the boost.
4. The SCZ/MDD association recomputed on the MTAG-boosted trait, **with the ASD
   and Alzheimer's controls alongside** (§7).
5. A raw-slope vs BLUP sensitivity check (§6).
6. An explicit statement of what the result does not license — in particular
   that MTAG's effective N is not a real sample size, and that a precision gain
   is not evidence of a larger effect.
