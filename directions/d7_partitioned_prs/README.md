# D7 — partitioned SBayesRC polygenic scores from cell-type and synaptic gene programmes

*Implements `docs/DIRECTIONS.md` D7. Written 2026-09-30; **run on ABCD
2026-10-01 (CSD3; result below).** It absorbs README_HPC §8.3 **C3-D** (the
C3-partitioned SCZ score): the C3, C1, C2 and PLS2 deciles are its secondary
sets.*

## Result (2026-10-01)

**Disorder side, yes; thinning side, no detectable enrichment.** Neuronal and
SynGO genes carry 1.25–1.35× their size-expected share of every score, but on
the HCP single-LMM global slope no set carries more of the PRS–thinning
association than its share of the score predicts. The gene-level link stays
group-level (DIRECTIONS §1), with the power caveat below: only ER_rel ≳ 3
was detectable.

- **Gates pass.** Step 2: row sums reproduce every Figure-1 `.profile` at
  r = 1.000000 (mean-imputed; 1,129,989 of 1,154,522 weighted SNPs matched,
  n = 11,670). Step 3: |Δβ| vs the step-9 lmer β ≤ 0.0016 on `global_slope_1lmm`
  for all four primary arms (worst cell anywhere: EA baseline, −0.0047).
- **K1c positive control (non-brain-expressed genes) depletes** on the
  disorder side in all four SCZ and MDD arms (f_enrich 0.75–0.86, p_f
  0.96–0.99 for enrichment), but not for EA (1.07, p_f 0.23). On the thinning side ER_rel is 0.0–0.56, the predicted
  direction, but not significant (p_ER 0.27–0.66).
- **Disorder side (`p_f`, null A).** K2 neuronal f_enrich 1.25–1.35 and K3
  SynGO 1.26–1.31, p_f ≤ 0.003 in all five arms, EA included: a property of
  brain-trait scores, not of SCZ. Glial and oligodendrocyte genes are not
  enriched (0.86–1.18). The S0 mechanical check holds: the ST12 SCZ locus pool
  carries 3.8× (SCZ25_META) and 3.4× (SCZ25_EUR) its share.
- **Thinning side (primary, 7 tests per arm, Bonferroni within arm).** Nothing
  survives in any arm. SCZ25_META (full β −0.035, p 4.7e-4, n 8,596):
  SynGO ER 0.98 (95 % CI −0.66 to 2.60; ER_rel 1.19), neuronal 0.31, glial
  0.14, oligodendrocyte 0.51; SynGO − oligodendrocyte ΔER 0.47 (−2.7 to 4.0).
  SCZ25_EUR (n 4,308): SynGO ER 1.47, every p ≥ 0.19. The nearest cell is
  MDD_eur oligodendrocyte genes, ER −3.8 (−27.9 to 0.83), p 0.011,
  p_bonf 0.078. No K2 or K3 contrast is nominal in any arm (p ≥ 0.11).
- **Secondary.** AHBA C3/C1/C2 and PLS2 deciles: no ER signal for SCZ. On the
  disorder side C3-bottom is depleted (f_enrich 0.67–0.71) and C1-top is
  nominally enriched for SCZ25_META only (1.33, p_f 0.013). The scattered
  nominal p_ER on symptom and baseline outcomes (e.g. MDD_eur glial on
  internalising, p 4e-4) sit in cells where ER is unstable (|ER| up to 57);
  they are not read.

Tables: `results/table_d7_{primary,partition,pairs,score_gate}.tsv`,
`results/SUMMARY.md`; rows in `genetic_analysis/current_results.tsv`
(section `d7_partition`). Step 4 (annotated SBayesRC refits) and K4 (D6 sets)
are not run.

**Run notes.** Steps 1–2 took 24 min, not 3–6 h. `config.sh` now reads the
MDD/EA weights from `legacy/hpc_v2/work/results_v2/prs_final/SBayesRC`, the
scores Figure 1 used. The old `$RES/prs_final` holds only association tables,
so those three arms would have been skipped with only a warning. `env_d7.yml`
pins `pandas<3`: copy-on-write and the new str dtype in pandas 3 are untested
against this code.

## Question

The SCZ and MDD polygenic scores predict a faster per-child rate of cortical
thinning (Figure 1; pooled SBayesRC SCZ 2025 β −0.034, SE 0.010). Is that
association carried by the part of each score that sits in a named biological
programme — synaptic genes, neuronal or glial cell-type genes, and later the
adolescent-window maturation genes of D6 — **more than that part's share of the
score predicts**? If so, the single-cell programme links to the individual-level
thinning effect. If not, the gene-level link stays group-level (DIRECTIONS §1).

## Method

**Partitioning with SBayesRC.** Each Figure-1 SBayesRC score is fitted once,
genome-wide and jointly with LD, with the baseline 2.2 annotation. Its per-SNP
posterior weights are subset by gene, so a set's partial score is the Figure-1
score restricted to that set's SNPs. This is the valid partitioned score that
README_HPC §8.3 describes. No refit per set is needed, and the parts add up
exactly to the full score. Every SNP goes to at most one gene: MAGMA 35/10 kb
windows (NCBI37.3, GRCh37), with SNPs in several windows going to the gene whose
body is nearest. The extended MHC (chr6:25–34 Mb) is its own bucket and is
excluded from every set.

**Per-gene partial scores.** One pass over the genotypes
(`work/inputs/geno/abcd_imp_prs`, the fileset Figure 1 scored) writes an
n × genes matrix per arm. Every set score, and every random-set null, is then a
column sum. A hard gate requires the row sums to reproduce the Figure-1
`.profile` (r > 0.9999).

**Statistics** (`code/partition_core.py`). y and the scores are residualised on
the `tools/prs_assoc.R` covariates (sex, age, PC1–10). The symptom arm adds the
Figure-1 panel-h terms: baseline score, age at follow-up, site. Two quantities
are computed for each set T:

- `f` = share of the score's variance carried by T;
- `s` = share of the score's association with y carried by T.

Both sum to 1 over any partition. **ER = s / f** is the enrichment ratio, and
β_T ≈ β_full × s / √f (the DIRECTIONS D7 table). Nulls are 5,000 random gene sets
from T's universe:

- **null A** is matched on the per-gene SNP-count decile and gives `p_f`. This is
  the disorder side: does T carry more of the disorder score than its size
  predicts?
- **null B** is additionally matched on each gene's covariance with the full
  score and gives `p_ER`. This is the thinning side: given its variance share,
  does T carry more or less of the association with y?

CIs are family-resampled bootstraps. OLS with family-clustered SEs replaces
lmer; a gate checks the full-score β against the step-9 lmer value
(|Δ| < 0.005).

**Read ER against its matched null, not against 1.** Intergenic SNPs are part of
the full score, so any genic set has ER > 1 if genic SNPs carry more than their
share. `ER_rel = ER / median(null B)` is the enrichment beyond matched random
genes.

## Pre-registration

These are fixed before any ABCD result is seen (2026-09-30). The roles follow
DIRECTIONS.md D7.

| id | set(s) | n genes | role |
|:--|:--|--:|:--|
| K1 | brain-expressed genes (17.3k that pass the developmental-cortex snRNA-seq expression filter) vs the rest of the genome | 17,293 | **positive control; must enrich.** K1 covers about 90 % of the genome, so the informative side is its complement, `K1c_not_brain_expressed`. The prediction is ER_rel < 1 and f_enrich < 1 for K1c. |
| K2 | neuronal markers vs glial markers (Seidlitz 2020; made disjoint) | 1,484 / 2,033 | primary |
| K3 | SynGO 1.3 synaptic genes vs oligodendrocyte markers (made disjoint) | 1,696 / 694 | primary; the pruning-vs-myelination discrimination of DIRECTIONS §2.2 |
| K4 | adolescent-window maturation genes from D6 | slot | primary once D6 exists; append rows to `gene_sets/d7_gene_sets.tsv` and rerun steps 3 and 5 |
| S0 | SCZ locus pool (Trubetskoy 2022 ST12) | 398 | secondary; also a mechanical check that f_enrich ≫ 1 for the SCZ score |
| S1 | snRNA-seq maturation PC1 top vs bottom decile (top = neuronal/synaptic, bottom = oligodendrocyte), snRNA-seq universe | 1,729 each | secondary |
| S2–S4 | AHBA C3, C1, C2 and ABCD PLS2 (HCP base) top vs bottom deciles, AHBA universe | 663 / 650 | secondary (the former C3-D). C1 and C2 are comparison axes, **not negative controls**: all three AHBA components are hypothesised mediators, and several poles are themselves SCZ-enriched. The null is the matched random sets. |
| S5 | union of Seidlitz cell-type markers | 4,591 | secondary; the background of K2 |

- **Primary readout.** HCP-MMP `global_slope_1lmm`, matched cells (rule 4):
  SCZ25_META and MDD_pooled → pooled (`_zanc`); SCZ25_EUR and MDD_eur → EUR (the
  replication). The tests are K1c, K2 ×2, K3 ×2 and the K2 and K3 contrasts
  (ΔER), seven per arm, with Bonferroni correction within arm. The direction is
  two-sided, because the tempo account does not say whether synaptic or myelin
  genes should carry the effect.
- **Interpretability rule.** ER and ΔER are read only where the full score
  predicts the outcome (`er_interpretable`, full p < 0.05). SCZ and MDD on the
  global slope qualify from Figure 1. Baseline thickness does not, and there
  only β_T is reported.
- **Secondary outcomes.** DK global slope, baseline thickness, the five
  Figure-1 CBCL change outcomes (p-factor, internalising, externalising,
  depressive, thought), and `c3axis_rc` if C3-A has been run (`D7_C3AXIS=1`).
  EA (EUR) is run as the opposite-direction trait.

**Expected power; read this before interpreting.** The β_T SE is about 0.011
for any partition. For a SynGO-sized set (f ≈ 0.08) and β_full ≈ 0.034, the SE
of ER is about 1.1. So only a set carrying roughly a quarter or more of the
signal in under a tenth of the variance, ER_rel ≳ 3, is detectable. A null
result means "not assessed at this power", **not** "the programme is not
involved". The disorder-side test (`p_f`, and step 4) is powered by the
discovery GWAS and is the better-powered half.

## Run order (CSD3)

Rules 1–19 of README_HPC §6 apply. Rule 16 matters most: `work/` holds
SNP-level and per-subject files and is gitignored. Commit only `results/`.

```bash
cd /home/rajd2/rds/hpc-work/abcd_development && git pull --ff-only origin main
conda env create -p ~/rds/hpc-work/envs/d7 -f directions/d7_partitioned_prs/env_d7.yml   # or legacy/hpc/work/bin/micromamba create -y -p ... (system conda is 4.7)
~/rds/hpc-work/envs/d7/bin/python -m pytest -q directions/d7_partitioned_prs/code/test_partition_core.py   # 5 pass
# laptop step L1 (the CBCL outcomes), then copy the file up: see below
bash directions/d7_partitioned_prs/run_d7.sh        # steps 1-2 -> 3 (15 tasks) -> 5
```

| step | script | where / time | output |
|:--|:--|:--|:--|
| 0 | `code/00_build_gene_sets.py` | laptop, done | `gene_sets/*.tsv` (committed) |
| L1 | `code/L1_export_symptom_outcomes.py` | laptop, done 2026-09-30 (8,716 children) | `work/d7_symptom_outcomes.tsv`: **scp it to the same path on CSD3** |
| 1–2 | `d7_scores.sbatch` → `code/01_assign_snps.py`, `code/02_gene_scores.py` | CSD3, 16 cores, ≈ 3–6 h | `work/gene_scores/<ARM>.npy`; `results/table_d7_score_gate.tsv` |
| 3 | `d7_partition.sbatch` (array; arms × {imaging_hcp, imaging_dk, symptoms}) → `code/03_partition_test.py` | CSD3, ≈ 15–30 min per task | `work/partition/*.tsv` |
| 5 | `code/05_collect.py` | anywhere | `results/table_d7_{partition,pairs,primary}.tsv`, `results/SUMMARY.md` |
| 4 (secondary, after 5) | `code/04_build_annot.py` → `d7_annot.sbatch --array=1-4` → rerun 2–3 with the `<ARM>_annot` arm strings the log prints (`d7_scores.sbatch '<arm string>'`, `D7_EXTRA_ARMS=...`) → `code/04_annot_readout.py` | CSD3, ≈ 6 h per fit | `results/table_d7_annot_enrichment.tsv` |

Laptop → CSD3 copy (run on the laptop, from `~/Git/abcd_development`):

```bash
scp directions/d7_partitioned_prs/work/d7_symptom_outcomes.tsv \
    rajd2@login.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/directions/d7_partitioned_prs/work/
```

**Inputs assumed on CSD3** (`config.sh`; check they exist before submitting):
the five SBayesRC weight and `.profile` pairs (SCZ 2025 in
`work/scores_scz2025/SBayesRC/`; MDD and EA in `work/results_70tab/prs_final/SBayesRC/`);
`$GENO`, `$STRATA`, `$ANNOT`, `$SBRC_EIGEN` from `genetic_analysis/setup/paths.sh`;
`/rds/user/rajd2/hpc-work/magma/gene_locations/NCBI37.3.gene.loc`; the
single-LMM exports `results_70tab{,_hcp}/prs_final_1lmm/pheno/`. A missing
weights/profile pair skips that arm with a warning. A missing gene.loc or
genotype file is fatal. If the MDD/EA `prs_final` directory is under a
different `$RES`, set it in `genetic_analysis/config.local.sh`, not here.

**Gates to check before reading any result:** `table_d7_score_gate.tsv`
`r_profile` > 0.9999 for every arm; the SUMMARY step-3 gate |Δβ| < 0.005 for
SCZ25_META, SCZ25_EUR, MDD_pooled and MDD_eur on `global_slope_1lmm`; K1c
behaving as the positive control. If K1c does not deplete, stop and diagnose
before reading K2/K3.

**Recording.** Done 2026-10-01: D7 block in
`genetic_analysis/build_current_results.py`, README_HPC §2.6, result above.

## Validation done on the laptop

- `code/test_partition_core.py` (5 tests): additivity of s and f; recovery of
  a planted enrichment; the matched null holding f; the linear `_zanc`
  transform adding up; clustered SE and bootstrap.
- Synthetic PLINK end-to-end (`plink --dummy`, 1,600 samples, 30,000 SNPs, 2 %
  missing calls, weights on either allele, real set Entrez ids on synthetic
  genes). Step 2 reproduced `plink --score … sum` at r = 1.000000 (max |Δ| 4e-6,
  float32), with mean imputation matching PLINK. With the signal planted in
  SynGO genes only, step 3 gave SynGO ER_rel 2.4, p_ER 0.004 (Bonferroni 0.028)
  and SynGO − oligodendrocyte p 0.004. Every other set sat inside its null
  (p 0.07–0.96).

## Deviations from the DIRECTIONS.md D7 text

- **Brain-expressed genes** are defined by the developmental snRNA-seq
  expression filter. They cover about 90 % of the genome, so the positive
  control is tested on the complement, as a depletion.
- **Symptoms** are tested in the imaged children only (8,716, of whom 8,596 are
  genotyped), not the full ~11k. This reuses the Figure-1 outcome construction
  verbatim. Extending to all children would need that construction rebuilt.
- **PQ-BC** (`mh_y_pps`) is not vendored yet. Add it to the symptom export when
  it is fetched.
