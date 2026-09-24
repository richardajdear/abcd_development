# genetic_analysis/ — genetics of the adolescent cortical thinning rate

State of knowledge as of **2026-09-24**, ABCD release 7.0. This file states
what is currently known from the best data and methods, what has been tried,
and what to do next. It is not a run log. The dated log that produced these
results (job ids, failures, diagnoses, the 6.0 → 7.0 re-run) is retired to
[`legacy/genetic_analysis/README_HPC_2026-09-24.md`](../legacy/genetic_analysis/README_HPC_2026-09-24.md);
read it only if the repo contradicts this file or you need a specific derivation.

**Every number below is a row of [`current_results.tsv`](current_results.tsv)**
(1,104 rows, built by `python genetic_analysis/build_current_results.py` from
the committed cluster tables; each row names its source table). Update the
table first, then this file.

---

## 1. Question, trait, primary specification

**Question.** Which genes drive the rate of cortical thinning across
adolescence, and do they overlap the genetic risk for adolescent-onset
disorders (schizophrenia, depression)?

**Trait.** `global_slope`: each child's rate of change in mean cortical
thickness, the random age slope from a linear mixed model on 2–4 scans per
child (ages ~9–16). `baseline_thickness` (the intercept) is carried alongside
as the positive control for every genetic readout.

**Sample.** 8,716 children with ≥ 2 QC-passing scans (26,946 scans); 8,596 of
them genotyped = the **pooled arm**; 4,308 of those in the European-ancestry
anchor set = the **EUR arm**. GREML uses the 6,011 PC-AiR-unrelated children.

**Primary specification** (the one `fig1.R` uses; everything else is a
sensitivity analysis):

| choice | primary | kept as sensitivity | why |
|:--|:--|:--|:--|
| parcellation | HCP-MMP1.0, 358 parcels (hippocampal parcel dropped) | Desikan–Killiany, 68 | finer; whole-cortex results agree across atlases |
| slope construction | **single LMM on the per-scan cortical mean** (`*_1lmm`) | mean of 358 per-parcel slope BLUPs (`perregion`) | averaging separately-shrunk BLUPs over-shrinks (SD 0.6×, r 0.82 with the single-LMM slope on HCP); the single LMM is less shrunk and more heritable |
| SCZ discovery GWAS | 2025 multi-ancestry SCZ GWAS (Nature 2025; EUR Neff 174 k) | PGC3 (2022) | larger; EUR-arm estimate no longer borderline |
| PRS methods | C+T (min-p permutation over 8 thresholds), PRS-CS, SBayesRC | SBayesR, PRS-CSx | three methods of different construction; agreement across them is the robustness criterion |
| ancestry matching | EUR GWAS → EUR arm; multi-ancestry GWAS → pooled arm, score z-scored within ancestry cluster | EUR weights → pooled, z within cluster | a European GWAS scored in a pooled sample is confounded (rule 4) |
| model | `scale(trait) ~ scale(PRS) + sex + age + 10 PCs + (1 \| family)` (`tools/prs_assoc.R`) | Fulker within/between family | β in SD of trait per SD of score |

## 2. Current state of knowledge

### 2.1 Summary

1. **The thinning rate is heritable, but no single variant is detectable at
   this n.** GREML h² = 0.18 ± 0.05 (HCP, single LMM; DK 0.21). No
   genome-wide-significant locus for `global_slope` in any of six scans (two
   atlases × EUR and pooled per-region, EUR single-LMM), λ_GC 0.99–1.03.
2. **Higher polygenic risk for schizophrenia predicts faster thinning, not
   thinner baseline cortex.** Significant under every method in both arms on
   both atlases (one exception: DK pooled C+T, p_adj 0.077), β ≈ −0.03 to
   −0.04 SD/SD; the same scores give null β for
   baseline thickness. This is the most robust genetic result in the project.
3. **The association is not specific to schizophrenia.** On the primary
   trait, APOE-carrying Alzheimer's scores (β −0.03 to −0.045), depression
   (−0.025 to −0.043) and, in the opposite direction, educational attainment
   (+0.03 to +0.04) have effects of the same size; autism is null. The
   Alzheimer's signal disappears when the APOE region is removed. SCZ is one of
   several polygenic influences on the rate, and whether it is independent of
   the others has not been tested on the current data (§4, item 1).
4. **At the gene level the only signal is the curated schizophrenia locus
   genes.** The Trubetskoy 2022 ST12 locus-gene set is enriched in the
   thinning-rate GWAS (MAGMA β 0.135, p 0.0046 HCP; 0.138, p 0.0039 DK), while
   gene sets built simply as "genes under SCZ GWS peaks" are weaker (p 0.04
   for the 2025 peaks, 0.14 for PGC3) and the ST12 set stays enriched
   conditional on the PGC3 peak set, so most of the enrichment belongs to the
   curated gene mapping rather than to proximity. Every
   other gene-level test on the slope is null or nominal on one atlas.
5. **Genetic correlation is not measurable here.** LDSC h² z for the slope is
   0.4 (HCP) / 0.15 (DK) against the ≳ 4 required; every slope rg is
   uninformative, not null.
6. **The transcriptomic link is group-level only.** The group thinning map
   tracks AHBA C3 (ρ −0.55, p_spin 0.002; see `ahba_pls/`), but no
   individual-level genetic readout connects to it: AHBA C1–C3 and the
   ahba_pls signature are null as MAGMA gene properties, and the per-parcel
   SCZ-PRS β map does not resemble C3 or the PLS lead (`prs_beta_map/`).

### 2.2 Heritability and GWAS

| | HCP single LMM | HCP per-region | DK single LMM | DK per-region |
|:--|--:|--:|--:|--:|
| GREML h² `global_slope` (n 6,011) | **0.183 ± 0.045** | 0.156 ± 0.045 | 0.207 ± 0.045 | 0.166 ± 0.045 |
| GREML h² `baseline_thickness` | 0.228 ± 0.047 | 0.225 ± 0.047 | 0.225 ± 0.046 | 0.223 ± 0.046 |
| LDSC h² z, slope / baseline (EUR) | 0.37 / 4.34 | 1.8 / 4.4 | 0.15 / 4.4 | 0.6 / 4.5 |
| EUR GWAS λ_GC, hits (slope) | 1.004, 0 | 0.994, 0 | 1.004, 0 | 0.995, 0 |
| pooled GWAS λ_GC, hits (slope) | not run | 1.028, 0 | not run | 1.027, 0 |

GREML (dense imputed GRM, unrelated set) is the estimate to quote. LDSC on
4,300 EUR children has SE ≈ 0.10 on h², so it cannot resolve anything below
~0.2; the disagreement between the estimators is expected.

One candidate locus exists, for an atlas-specific component rather than the
primary trait: chr7:35.55–35.63 Mb (rs62452241 and two SNPs in LD, MAF 5–8 %)
for HCP-EUR `slope_PC3`, p 1.2e-9, weaker in the pooled arm of twice the n
(1.5e-5), which is the pattern of winner's curse. It needs an independent
cohort. The slope PCs are atlas-specific decompositions (except PC2) and are
not primary phenotypes.

### 2.3 Polygenic scores: 2025 SCZ → `global_slope` (primary)

β (SD/SD), SE, p (threshold-adjusted for C+T). EUR arm n 4,308; pooled n 8,596.

| cell | method | HCP single LMM | DK single LMM | HCP per-region |
|:--|:--|:--|:--|:--|
| EUR GWAS → EUR | C+T | −0.044 (0.015), 0.020 | −0.041, 0.043 | −0.051, 0.0045 |
| | PRS-CS | −0.036 (0.015), 0.013 | −0.034, 0.022 | −0.041, 0.0057 |
| | SBayesRC | −0.035 (0.015), 0.016 | −0.030, 0.038 | −0.042, 0.0043 |
| multi GWAS → pooled (z within cluster) | C+T | −0.033 (0.010), 0.010 | −0.027, 0.077 | −0.033, 0.016 |
| | PRS-CS | −0.029 (0.010), 0.0044 | −0.027, 0.0098 | −0.023, 0.027 |
| | SBayesRC | **−0.034 (0.010), 0.0010** | −0.030, 0.0036 | −0.028, 0.0079 |

- **Baseline thickness:** every SCZ 2025 cell null (|β| ≤ 0.021, p_adj ≥ 0.17
  on HCP). The association is with the rate.
- **Min-p permutation (C+T, 2,000 family-block permutations, per-region
  construction, HCP):** EUR 0.0040, pooled 0.0080, ancestry-matched composite
  0.011.
- **Ancestry.** PRS × ancestry-cluster LRT p 0.19 (META) / 0.92 (EUR weights);
  the stratum-fixed-effects pooled β equals the headline β. Within clusters
  the SE is 0.03–0.05, so this is absence of evidence of heterogeneity, not
  evidence of homogeneity. Scores built from the AFR (Neff 35 k) and EAS
  (29 k) GWAS predict nothing in any cluster, including their own; EUR weights
  transfer a same-sign signal into the African-American-like and
  Hispanic-like clusters. The 2025 meta has no Latino cohort.
- **Within family (Fulker; per-region construction, HCP).** EUR arm, 726
  sibling pairs: β_within −0.015 to −0.095 (SE ≈ 0.047), β_between −0.030 to
  −0.046, p_diff 0.33–1.0; pooled arm (META weights), 1,449 pairs: β_within +0.013 to
  +0.026 (SE 0.03–0.06). Power to detect β = 0.04 within families is ~15 %, so this
  neither confirms nor contradicts a direct genetic effect.

### 2.4 Specificity panel (PGC-era scores, primary trait, HCP single LMM)

Count of methods (of 4: C+T, PRS-CS, SBayesR, SBayesRC) with p_adj < 0.05,
matched stratum only.

| score | `global_slope` | β range | `baseline_thickness` |
|:--|:--|:--|:--|
| SCZ pooled (PGC3) | 4/4 | −0.025 to −0.036 | 0/4 |
| SCZ EUR (PGC3) | 0/4 (p_adj ≥ 0.061) | −0.026 to −0.039 | 0/4 |
| MDD pooled | 3/4 | −0.024 to −0.031 | 0/4 |
| MDD EUR | 3/4 | −0.029 to −0.043 | 0/4 |
| ALZ (Wightman, with APOE) | 4/4 | −0.033 to −0.045 | 0/4 |
| ALZ without APOE | 0/4 | −0.015 to −0.027 | 0/4 |
| ALZ_IGAP (Kunkle) with / without APOE | 3/4 / 0/4 | | 0/4 |
| EA (Okbay) | 3/4, **positive** | +0.032 to +0.040 | 0/4 |
| ASD | 0/4 | −0.010 to +0.005 | 1/4 (PRS-CS −0.039) |

On DK the pattern is the same (ALZ 4/4, EA 3/4, MDD 3/4 pooled / 2/4 EUR). On
the per-region construction the non-SCZ signals are weaker (ALZ 1–2/4, EA 0–3/4),
which is why the single-LMM slope, not SCZ specifically, is what gains.

### 2.5 Gene level (MAGMA, EUR arm, 1000G EUR LD; single LMM, HCP)

| test on `global_slope` | β | p |
|:--|--:|--:|
| ST12 `SCZ_locus_pool` (455 curated locus genes) | 0.135 | 0.0046 |
| … minus the 101 prioritised genes (354) | 0.185 | 0.00074 |
| ST12 prioritised genes (101) | −0.035 | 0.65 |
| `SCZ25_locus_pool` (1,103 genes under 2025 GWS peaks) | 0.074 | 0.039 |
| `PGC3_locus_pool` (858 genes under PGC3 peaks) | 0.051 | 0.14 |
| SCZ gene-level-significant genes (2025 / PGC3) | −0.008 / −0.040 | 0.57 / 0.79 |
| MDD pool / high-confidence genes | 0.011 / −0.014 | 0.34 / 0.59 |
| reverse gene-property, SCZ25_EUR (disorder gene Z ~ thinning gene Z) | 0.022 | 0.045 (DK 0.11) |
| MDD, ASD, ALZ, EA gene-property | | all ≥ 0.20 |
| AHBA C1–C3, ahba_pls signature (H3), marginal and C3-conditional | | all ≥ 0.066 |

`baseline_thickness` is enriched in the SCZ ST12 pool (p 1e-4) and associated
with SCZ, MDD and EA gene signal (p 3e-4 to 0.04): at the gene level the
disorders track thickness, not its change. The multi-ancestry SCZ gene
analysis is per-ancestry with matched 1000G panels, combined with
`magma --meta`; running the meta file against the EUR panel alone changes ~16 %
of the significant gene list but no phenotype-level conclusion.

**Step 14 (pooled-arm MAGMA, ABCD sample as its own LD reference; 8,596
children, both constructions, both atlases; tables committed in
`work/results_70tab*/magma_pooled/table_magma_pooled.tsv`).** With the sample
doubled and LD matched, SCZ is the one trait with a gene-level slope signal on
both atlases and both constructions (SCZ25_META p 0.011–0.029, PGC3_EUR
0.010–0.046; SCZ25_EUR weaker, 0.13–0.38), and the 2025 peak-window pool
`SCZ25_locus_pool` is enriched on the pooled slope (p 0.004–0.017), joining the
ST12 pool (0.03–0.13). None survives correction over the 11 × 4 panel; it is
consistent in direction with the PRS result. MDD, ASD, ALZ and EA stay null on
the slope. The EUR-arm baseline gene-property signals (SCZ, MDD, EA) vanish in
the pooled arm (every p ≥ 0.07): present them as EUR-arm-only. Single-LMM and
per-region agree within noise. Run detail in the archived log (§8 step 14).

**Step 15 (run 2026-09-24) — EUR-arm MAGMA on ABCD's own LD, and gene-set
tests on three gene-result versions (tables committed in
`work/results_70tab*/magma_set_tests/table_magma_set_tests.tsv`).**
15a re-ran the EUR-arm single-LMM gene analysis with the ABCD EUR children
(step-14 reference ∩ `eur_anchor.keep`, n = 4,308; 17,901 genes) as LD
reference, so the EUR and pooled arms now differ only in sample. 15b tested 31
sets (step-13 SCZ/MDD sets, `genesets_topgenes.txt`, `genesets_wes.txt`) on
EUR_1000G, EUR_ABCD and pooled_ABCD gene results. EUR_1000G rows reproduce
step 13 exactly (SCZ_locus_pool, HCP slope, p 0.0046; SCZ_WES 0.026, MDD_WES
0.97). **Swapping the LD panel removes the headline ST12-pool result:**
SCZ_locus_pool on the HCP slope goes from p 0.0046 (1000G) to 0.15 (ABCD EUR;
pooled 0.13), and on baseline thickness from 1e-4 to 0.43. DK is the same (0.004
→ 0.23; 4e-4 → 0.61). The §2.5 ST12-pool rows therefore reflect the
503-person 1000G panel, not the sample: treat them as not robust. HCP, EUR_ABCD
vs EUR_1000G, p (slope / baseline): SCZ_locus_pool 0.15 / 0.43 vs 0.0046 /
1e-4; MDD_highconf 0.31 / 0.045 vs 0.59 / 0.075; SCZ_WES (24 genes) 0.050 / 0.061
vs 0.026 / 0.021 (pooled 0.18 / 0.055; DK EUR_ABCD slope 0.027); MDD_WES 0.55 /
0.20 vs 0.97 / 0.49. What does strengthen on the ABCD panel is the
gene-level-significant MDD set on the slope (MDD_genesig p 0.0044 HCP, 0.0021 DK;
top-250 0.003 / 0.006; pooled 0.013), same direction in all three versions.
Nothing survives Bonferroni over 31 sets (0.0016). Run order: step14 prep → 15a
(4 × 1.5 h) → 15b (9 min).

## 3. What has been tried

Closed means do not repeat without a new reason; the reason is given.

| approach | outcome | status |
|:--|:--|:--|
| GENESIS mixed-model GWAS keeping relatives (sparse PC-Relate kinship, PC-AiR PCs) | 0 hits on the slope, λ in band | **current** |
| GCTA GREML on dense imputed GRM, unrelated set | h² 0.16–0.21 | **current** (canonical h²) |
| fastGWA / GCTA on release-4.0 EUR genotypes (v1) | 20 "hits" at λ 1.107 were stratification | closed: superseded by GENESIS |
| Zaitlen two-GRM "pedigree h²" | GRM diagonal ≈ PC1 (r −0.99); insensitive to deleting 97 % of relative pairs | closed: not identifiable in this design |
| four-method PRS grid, ancestry-matched cells, within-cluster standardisation | §2.3–2.4 | **current** |
| PRS-CSx and per-child ancestry-matched weights | no gain over META SBayesRC or EUR weights; AFR/EAS GWAS too small | sensitivity; revisit when non-EUR GWAS grow |
| within-family (Fulker) | uninformative, ~15 % power | current, reported with its power |
| slope construction: mean of parcel BLUPs vs single LMM | single LMM less shrunk, more heritable, all PRS signals stronger | **single LMM adopted** |
| parcellation DK vs HCP-MMP | whole-cortex results agree; marginal cells move across p = 0.05 | HCP primary, DK sensitivity |
| region-subset phenotypes (top-ΔCT, top-C3 regions, projections; v3) | share 80–90 % of variance with the global slope, less heritable, C3 set weaker than 80 % of random sets | closed |
| slope PCs as GWAS phenotypes | atlas-specific; one candidate locus (chr7, HCP PC3) | secondary only |
| LDSC h² and rg (EUR GWAS files only, full disorder panel) | slope h² z < 2; rg uninformative | closed until n grows |
| LAVA local rg | anti-conservative (11 % of tests p < 0.05 genome-wide) | closed unless permutation-calibrated |
| MOSTest / JAGWAS multivariate GWAS | output is unsigned, so no disorder direction; only useful for locus discovery | closed for the disorder question (`legacy/hpc_v3/SETUP_CONTEXT.md`) |
| MTAG / genomic SEM across regional slopes | needs LDSC h² z ≳ 4 per trait; slope is at 0.4 | gated |
| MAGMA gene-property and gene sets, SCZ/MDD/ASD/ALZ/EA, both directions | §2.5 | **current** |
| AHBA C1–C3 and ahba_pls signature as MAGMA gene properties (H3) | null, both atlases | closed |
| per-parcel SCZ-PRS β map vs C3 / PLS lead, spin test | null; the map is a uniform shift plus noise (lh–rh ρ 0.3) | closed at parcel level; the projected-phenotype version (H4) is open |
| imputed C4A expression → slope (`c4_imputation/`, Sekar panel, EUR arm) | null: β −0.011 (95 % CI −0.039, +0.017), excludes the pre-specified MDE 0.048; imputation QC fine | closed; only the pooled arm (Kamitaki panel, dbGaP) is untried |
| PRS conditional on EA (`R/07_prs_conditional.R`) | run only on 6.0-vintage data | **open**: re-run on the current trait (§4, item 1) |
| PRS × age interaction in a one-stage LMM (`tools/age_prs_interaction.R`) | run only on 6.0-vintage data | **open** (§4, item 2) |

## 4. Next analyses, from first principles

The question is "which genes drive the thinning rate". The GWAS cannot answer
it directly at n ≈ 8,600 (h² ≈ 0.18 spread over many variants; a
single-variant hit needs something like 10× the sample). The polygenic score
is the only individual-level signal that works, so the most informative next
steps either (a) establish what that signal is, or (b) decompose it into
genes. Ranked by information per unit of effort:

1. **Joint multi-score model (specificity).** Fit SCZ 2025 + MDD + EA +
   ALZ-without-APOE + APOE ε4 dosage (rs429358/rs7412) in one model on the
   primary trait. If SCZ survives adjustment, the claim "SCZ risk predicts
   thinning" stands on its own; if the scores share one component, the finding
   is a general polygenic influence and should be written that way. Code
   exists (`R/07_prs_conditional.R`); scores exist; minutes of compute. Also
   report APOE ε4 on its own: the ALZ result is APOE, so the direct genotype
   test is one interpretable number.
2. **One-stage longitudinal model.** Fit `thickness ~ PRS × age + … + (1 + age
   | child) + (1 | site) + (1 | family)` on the scans directly. This removes the
   BLUP-shrinkage question entirely (the reason the two constructions give
   different answers) and gives correct SEs. `tools/age_prs_interaction.R` is
   the template. It is also the natural model for §4 items 3–5.
3. **Measurement confounds of change.** A PRS can predict apparent thinning
   through scan quality: motion and FreeSurfer surface-reconstruction quality
   (Euler number) reduce estimated thickness, and if they change with age
   differently by genotype they produce a spurious slope. EA and ADHD
   polygenic scores are the usual carriers. Test (i) PRS → per-scan Euler
   number and its age slope; (ii) the primary PRS association with time-varying
   Euler number as a covariate in the one-stage model. Add negative-control
   scores that should not act on the cortex (e.g. height) and an ADHD score.
   The EA effect in the opposite direction makes this more, not less, worth
   checking. Also test PRS → number of usable scans: how much a child's slope
   is shrunk depends on how many scans they have, so a score that predicts
   attrition or scan failure can shift the phenotype without acting on the
   cortex; an inverse-probability-weighted sensitivity closes this.
4. **Partitioned polygenic scores: which genes carry the SCZ signal.** Split
   the SCZ 2025 score into gene-set partitions and test each against the
   thinning rate, with size- and LD-matched random sets as the null (PRSet, or
   SBayesRC annotation-partitioned weights). Candidate partitions:
   the ST12 locus genes (the one enriched set), cell-type-specific genes from
   the group's developmental snRNA-seq analysis (Velmeshev), AHBA C3
   positive/negative poles, SynGO synaptic genes. This is the most direct route
   from the one robust signal to genes, and it is powered by the discovery GWAS
   rather than by ABCD.
5. **The C4 hypothesis.** Complement component 4A, imputed from SNPs (Sekar
   2016 reference haplotypes), is the best-characterised SCZ mechanism tied to
   synaptic pruning in adolescence. Imputed C4A expression → thinning rate is
   one pre-specified test and would name a gene. The SCZ scores do carry MHC
   SNPs (SBayesRC 2,532, C+T 980; `c4_imputation/results/mhc_in_scz_scores.tsv`),
   so §2.3 partly contains it; C4A is null with or without the score.
   **Done, null** (§3; [`c4_imputation/`](../c4_imputation/README.md)).
6. **Individual-level projected phenotypes (ahba_pls H4).** Project each
   child's 358 slopes onto C3 and the PLS lead component (one phenotype per
   child) and regress on the SCZ score. Averages out the parcel noise that
   sank the β-map test; `ahba_pls/FOLLOWUP_GENETICS.md` H4 has the design.
7. **External positive controls and replication.** (i) A score for adult
   cortical thickness (e.g. Warrier et al. 2023 or ENIGMA) should predict
   `baseline_thickness`; whether it predicts the rate tells us how much of the
   rate's genetics is shared with adult structure. (ii) The ENIGMA-Plasticity
   longitudinal change GWAS (Brouwer et al. 2022, 15,640 people aged 4–99) is
   the only published GWAS of cortical change, but **ABCD release 3.0 is one of
   its cohorts**, so its summary statistics cannot be scored in, or
   correlated with, ABCD without a leave-ABCD-out version from ENIGMA. Its
   thickness-change GWAS is also weak (h² z < 4 for every phenotype, no
   thickness-change locus), so it is a poor positive control. What transfers
   is its APOE rs429358 result (age-dependent effects on change rates,
   including in development) and its finding that the genetics of change
   overlap little with the genetics of level. (iii) Replicate the SCZ-PRS → thinning association in
   an independent longitudinal adolescent cohort (IMAGEN, Generation R). Of all
   items here, replication is what a reviewer will ask for first.
8. **Other measures of the same process.** Thinning partly reflects
   intracortical myelination rather than grey-matter loss. Slopes of surface
   area and of grey/white contrast from the same scans would test whether the
   SCZ association is specific to thickness. Lower priority: new phenotype
   pipelines, not new genetics.

Not recommended at this n: rg/MTAG/genomic SEM on the slope (h² z 0.4),
TWAS/colocalisation on the slope GWAS (no locus to colocalise), further
region-selection phenotypes (closed, §3).

## 5. Pipeline

### 5.1 On CSD3

Repo `~/rds/hpc-work/abcd_development`; account `VERTES-SL3-CPU`, partition
`icelake` (the default `cclake` once cost 14 h); never write to `/home`.
`genetic_analysis/config.local.sh` (gitignored; template
`config.local.sh.example`) sets the run root to `genetic_analysis/work` and
names each reused asset at its real path. It is the authority for paths; do
not re-derive them. `PARC=dsk|hcp` is the only switch between atlases
(outputs `work/results_70tab/` and `work/results_70tab_hcp/`).

Reused genotype-only products (do not rebuild): GDS and kinship
(`legacy/hpc_v2/work/results_v2/{gds,kinship,kinship_eur}`), the dense imputed
GRM, the 22 imputed `imp_union_chr*` filesets, every PGC-era PRS score and
posterior-weight file, normalised sumstats, LD references, the disorder-side
MAGMA gene results. Environments: `~/rds/hpc-work/envs/genesis` (R 4.4.3,
GENESIS 2.36; spec `envs/genesis_env.yml`), `legacy/hpc/work/envs/abcd`
(Python), `~/.conda/envs/ldsc` (Python 2.7 LDSC). Binaries: `gcta64`, `magma`,
PLINK 1.9 (required for `--clump`/`--score`), `plink2`, GCTB, PRScs/PRScsx.

### 5.2 Steps → scripts → outputs

Run from the repo root on CSD3. Each step writes a summary table under
`work/results_70tab*/`; those tables are tracked and feed
`build_current_results.py`.

| step | what | scripts | output (per atlas root) |
|:--|:--|:--|:--|
| 1 | phenotype export (per-region BLUP construction), align IDs, z-score over the genotyped set, verify | `step1_export_pheno.sbatch`, `step1_align_and_preflight.sbatch`, `step1_verify_export.py`, `setup/align_export.py`, `setup/make_grm_famid.sh` | `pheno_70tab*/` (per-subject, never committed) |
| 1b | single-LMM phenotypes (primary) | `orderops/fit_globalmean.sbatch`, `orderops/check_blups.py`, `orderops/build_1lmm_pheno.py` | `prs_final_1lmm/pheno/` |
| 2 | preflight | `00_check_inputs.sh` | — |
| 3–4 | GENESIS null model + 22-chromosome scan, pooled and EUR | `run_all.sh 03 04`, `03_null_model.sbatch`, `04_assoc.sbatch`, `R/03–05`, `setup/04_collect.sbatch`, `collect_1lmm.sbatch` | `assoc{,_eur}/gwas_summary.tsv`, `scan_1lmm/assoc_eur/` |
| 5 | GREML h² | `step5_reml.sbatch` | `reml_imp_pooled/`, `scan_1lmm/reml_imp_pooled/` |
| 6 | PGC-era PRS association grid, controls, Fulker, min-p | `step6_prs_assoc.sbatch` → `setup/prs_final.sbatch`, `setup/collect_final.py`, `setup/standardise_within_ancestry.py`, `06_prs_family.sbatch`, `step6_minp_permutation.{py,sbatch}` (uses `setup/prs_paired_delta.py`), `step6_prs_1lmm_controls.sbatch`, `step6_orderops_collect.py`, `orderops/prs_1lmm.sbatch`, `orderops/compare_1lmm.py` | `prs_final/`, `prs_final_1lmm/table_order_of_operations_all.tsv` |
| 7 | LDSC + MAGMA on the ABCD GWAS (EUR), prioritised gene sets | `step7_ldsc_magma.sbatch`, `setup/02_ldsc_magma.sbatch`, `setup/05_prio_gsa.sbatch`, `setup/eur_{magma_ldsc,prio_gsa}.sbatch` | `ldsc_eur/`, `magma_eur/`, `magma_prio_eur/` |
| 8 | ahba_pls H3 gene-property | `step8_ahba_pls_h3.sh` | `magma_ahba_pls_h3/table_h3.tsv` |
| 9 | 2025 SCZ GWAS: normalise (Neff fix), score all methods, gather, associate, strata, min-p | `run_scz2025.sh`, `step9_scz2025_*.{sbatch,py}`, `setup/normalise_gwas.py`, `setup/build_matched_scores.py`, `R/09_prs_ancestry_strata.R`, `step9_scz2025_assoc_1lmm.sbatch`, `step9_scz2025_orderops.py` | `prs_scz2025/table_*.tsv` |
| 10 | MAGMA on the 2025 SCZ GWAS, per ancestry + `--meta` | `run_scz2025_magma.sh`, `step10_scz2025_magma_*`, `setup/build_scz2025_genesets.py` | `magma_scz2025/` |
| 11–13 | disorder-side MAGMA, LDSC panel, MAGMA panel on both constructions | `step11_magma_disorders.sbatch`, `step12_ldsc_panel.{sbatch,_collect.py}`, `step13_magma_panel.{sbatch,_collect.py}` | `ldsc_1lmm/`, `magma_panel/` |
| 14 | pooled-arm MAGMA with in-sample LD | `step14_magma_pooled_{prep,genes,tests}.sbatch` | `magma_pooled/` |
| — | tables, comparison, figures (laptop) | `build_current_results.py`, `compare_parcellations.py`, `fig1.R` (R, ggplot2 + patchwork; `LC_ALL=en_US.UTF-8 Rscript`, env ahba-pls-r; reads `fig1_inputs/`), `fig_genetics_panel_1lmm.py` | `current_results.tsv`, `results_70tab_hcp/compare/`, `docs/figures/` |

Unused R entry points kept for §4: `R/07_prs_conditional.R` (joint model with
EA). `R/01_make_gds.R` and `R/02_kinship.R` are provenance only; do not re-run.

## 6. Rules

Numbering is kept from the retired README so that script comments citing
"rule N" still resolve.

1. **IDs.** ABCD spells a child three ways (`NDAR_INVxxxxxxxx`,
   `sub-NDARINVxxxxxxxx`, `sub-xxxxxxxx`); CSD3 genotype files use
   `sub-xxxxxxxx`. GCTA and PLINK match FID+IID literally and return zero
   subjects without an error. Join on the 8-character token.
2. **Two FID conventions.** Genotype `.fam` has FID = IID. GCTA wants FID
   matching the GRM (`make_grm_famid.sh` writes a family-id view in the GRM's
   row order; never reorder it); `prs_assoc.R` and Fulker want FID = family id.
3. **Scale.** Every β is per SD. Z-score over the phenotyped-and-genotyped set
   *after* subsetting.
4. **Ancestry matching decides the stratum.** Multi-ancestry GWAS → pooled
   (score z within cluster); European GWAS → EUR. European GWAS → raw pooled
   is confounded and once manufactured a false ASD hit. Read only `matched`
   cells.
5. **Allele frequencies from the target sample** (`maf_union.tsv`), never a
   reference panel.
6. **C+T multiplicity is a family-block min-p permutation**, not Bonferroni.
7. **Controls are mandatory** (ASD, ALZ ± APOE, Kunkle, EA) and are reported
   with every disorder result.
8. **Report λ_GC with every scan.**
9. **Fresh output directories for LDSC/MAGMA**: their collectors overwrite
   `*_summary.tsv` without warning.
10. **SLURM:** `--export` splits on commas; `DRY_RUN=1` once submitted a
    pipeline of no-ops; verify outputs (row counts, `n_chr == 22`), never
    `sacct` state; export `TMPDIR` onto rds and avoid here-documents (node-local
    `/tmp` fills and kills bash); `afterok` on an array dies if any task fails.
11. **Zaitlen two-GRM h² is not obtainable here.** Quote GREML only.
12. **LAVA is anti-conservative here.** Do not use uncalibrated.
13. **rg needs h² z ≳ 4 on both traits.** Report the z; below it rg is
    uninformative, not null.
14. **Age is a covariate in every PRS model** (the export calls it
    `baseline_age`; the R scripts derive `age_c`).
15. **Region selection is not a lever** (closed, §3).
16. **Never commit per-subject files** (`.profile`, `.sscore`, `pcair_pcs.tsv`,
    phenotype exports, per-child scores). The repo is public; summary tables
    only.
17. **`scratch/hpc_test*` are synthetic fixtures with planted effects.** They
    load cleanly and give plausible wrong answers.
18. **Effective N in meta-analysed sumstats.** bcftools +metal sums each
    study's N column; PGC files carry Neff/2. Check max N against the largest
    single contributing study before scoring (`step9_scz2025_normalise.sbatch`).
19. **MAGMA LD reference must match the GWAS sample.** For multi-ancestry
    discovery GWAS, run per ancestry with the matched panel and combine with
    `--meta` (both `raw=` and `genes=` calls are needed).

## 7. How to record new work

A new result goes into its step's tracked table, then into
`current_results.tsv` (extend `build_current_results.py`), then into §2 here,
replacing the number it supersedes. Keep the job id, n and output path in the
commit message, not in this file. A result that changes a conclusion in §2.1
goes at the top of that section. Record dead ends in §3 with their diagnosis
in one line: several past defects produced plausible numbers rather than
errors, and the diagnosis is what stops the next agent repeating them.
