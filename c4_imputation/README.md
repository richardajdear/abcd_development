# c4_imputation — does imputed C4A expression predict the adolescent thinning rate?

**Status, 2026-09-24: run, null, closed.** Imputed C4A expression does not
predict the thinning rate (EUR arm, n 4,308; M1 β −0.011, 95 % CI −0.039 to
+0.017, p 0.45); the secondary models and the other two cells agree. Imputation
QC passed (4,837 array SNPs aligned; EUR structure frequencies match the
panel). Tables in `results/c4_assoc_*.tsv`; the record is in
[`genetic_analysis/README_HPC.md`](../genetic_analysis/README_HPC.md) §3.
The plan below is kept as it was run.

---

## 1. Question

Complement component 4A (C4A) is the best-characterised schizophrenia risk
mechanism. Sekar et al. (2016, *Nature* 530:177) showed that much of the MHC
schizophrenia association comes from structural variation at *C4*: copy number
of *C4A* and *C4B*, and a HERV insertion (the "long" form) that raises
expression. Risk rises with predicted *C4A* expression. In mice, C4 marks
synapses for microglial removal. The hypothesis is that excess complement-tagged
pruning in adolescence contributes to schizophrenia.

Adolescent cortical thinning is the MRI-visible process most often attributed
to pruning. The test here: **does higher genetically predicted C4A expression
predict a faster thinning rate** (`global_slope`, the project's primary trait,
genetic_analysis/README_HPC.md §1)? A positive result would move the project's
headline from "schizophrenia polygenic risk predicts thinning" to a named gene
and pathway.

## 2. What is already known

- **Hernandez et al. 2023** (*Genome Biology* 24:42) imputed C4 in ABCD (release
  3.0, n 7,789, multi-ancestry panel of Kamitaki et al. 2020). C4A expression
  was **not associated with global cortical thickness, surface area or volume at
  baseline** (ages 9–10, cross-sectional). It was associated with smaller
  entorhinal surface area, independently of the schizophrenia polygenic score.
- A UK Biobank study (*Psychological Medicine*, 2021; Sekar panel) reported
  associations of C4A expression with cognition and with thickness and area in
  selected regions in adults.

**What is new here** is the *rate*: the per-child slope from 2–4 scans across
ages 9–16 on release 7.0. The Hernandez null is for the level at one age. The
baseline-thickness control (M5) should replicate it.

## 3. Expected effect and power

The minimum detectable effect at 80 % power (α 0.05) in the EUR arm (n 4,308)
is **|β| ≈ 0.048 SD of thinning rate per SD of true C4A expression**. This
allows for imputation error (C4A GREx r = 0.89 with the truth,
results/power.tsv). With perfect genotypes it would be 0.043.

That is about the size of the entire schizophrenia polygenic-score effect
(β −0.035 to −0.044 in this arm). If C4 acted only through its share of
schizophrenia liability, the expected β would be far smaller: C4 carries a
small fraction of the SNP-based liability that the polygenic score does.
**The test is therefore informative only for a C4 effect on the thinning rate
that is disproportionate to its effect on schizophrenia risk.** That is the
mechanistic hypothesis, but it means a null result is weak evidence against C4
involvement. Report the confidence interval, not just p.

## 4. Analysis plan (pre-specified before any ABCD genotype is touched)

Model (`tools/prs_assoc.R`'s model, predictor swapped):

    scale(trait) ~ scale(C4A_GREx) + sex + age_c + PC1..PC10 + (1 | family_id)

| model | what | role |
|:--|:--|:--|
| **M1_primary** | EUR arm, all children, C4A GREx from expected allele dosages; trait = HCP-MMP single-LMM `global_slope_1lmm` | **the test**: one test, two-sided α 0.05, hypothesis β < 0 |
| M2_joint_C4B | + C4B GREx | specificity (C4B is not a schizophrenia risk gene) |
| M3_sex_* | C4A × sex; female and male fits | Kamitaki et al. 2020 report stronger C4 effects in men |
| M4_conditional_PRS | + SCZ 2025 EUR SBayesRC score; M4_PRS_alone alongside | is C4A independent of genome-wide risk? Read with step 5 (MHC content of the score) |
| M5_control_baseline | same predictor, baseline thickness | replication of the Hernandez null for the level |
| M6_common5 | children whose structures are all among AL, AL-AL, AL-BL, AL-BS, BS | Hernandez et al. filter |
| M7_copy_numbers | C4A, C4B and HERV copy numbers jointly | free of the expression weights |
| M8_posterior_filter | M1 restricted to structure posterior ≥ 0.7 | Hernandez et al. filter |
| cells | `hcp_1lmm` (primary), `dk_1lmm`, `hcp_perregion` | parcellation and slope-construction sensitivity |

Secondary models are reported whatever M1 shows and are not used to rescue it.

**Why no posterior filter in M1.** In the leave-out test (§5), Beagle's posteriors
were poorly calibrated: 85 % of people had a structure posterior ≥ 0.7, but
only 74 % of those were called correctly, against 70 % for everyone. The
posterior-weighted dosage needs no filter and keeps n.

**GREx weights** as printed by Hernandez et al. (citing Sekar 2016):
`C4A = 0.47·AL + 0.47·AS + 0.20·BL`, `C4B = 1.03·BL + 0.88·BS`. **Open item:**
check these against Sekar 2016 Extended Data before the cluster run. The BL
term in the C4A equation is unexpected. The weights are one constant in
`03_c4_grex.py`, and M7 is weight-free in any case.

## 5. Method decisions and the local results behind them

**Panel.** The public Sekar 2016 panel: 111 HapMap3 CEU individuals (222
haplotypes), 7,752 MHC markers plus the multi-allelic C4 marker (25 structural
alleles), GRCh37. It is **European only**, so the analysis is the EUR arm only.
The multi-ancestry panel (Kamitaki et al. 2020; WGS-derived, 2,530 haplotypes;
dbGaP **phs001992**, controlled access) is what Hernandez used, and it would add
the pooled arm. **Upgrade path:** a dbGaP data-access request (or ask the
McCarroll lab). The pipeline takes the panel as `PANEL=`, so only the panel
and its build need changing.

**Input genotypes.** The 7.0 Smokescreen array hard calls (hg19), not the
imputed dosages. The panel is dense only at HapMap3 sites, and imputation error
in the MHC would propagate into the C4 call. conform-gt aligns strand and
allele orientation to the panel by genotype correlation.

**Imputation.** Beagle 5.5 (27Feb25) with ap=true. Sekar and Hernandez used
Beagle 4.1; the imputec4 protocol uses Beagle 5.1. The genetic map is linear,
cM = bp × 1e-6. This was chosen by cross-validation
(results/smoketest_map_scale_grid.tsv): 1e-6 gave the best C4A GREx accuracy
at both SNP densities (r 0.89–0.90). The imputec4 default of 1e-7 gave
0.62–0.81, and 1e-5 gave 0.50–0.70. The scale was chosen on the same 111 people
it is scored on (5 candidate values), so its accuracy is mildly optimistic.

**Leave-fold-out accuracy** (results/smoketest_accuracy.tsv). 5 folds; each
fold's ~22 people were imputed from the other ~89 (a lower bound on the full
panel). Targets were array-like: C4 removed, unphased, 5 % strand-flipped, 50
decoy SNPs.

| SNP density | structure-pair concordance | haplotype concordance | r C4A GREx | r C4B GREx | r HERV copies |
|:--|--:|--:|--:|--:|--:|
| all panel SNPs | 0.71 | 0.84 | 0.90 | 0.63 | 0.90 |
| every third SNP (~2,500, array-like) | 0.70 | 0.83 | 0.89 | 0.60 | 0.90 |

C4A expression is imputed well, and C4B less well. The ~70 % structure-pair
concordance matches the accuracy usually quoted for the Sekar panel.

**Panel composition** (results/reference_panel_c4.tsv). Structures AL-BL 41 %,
AL-BS 31 %, AL-AL 11 %, BS 7 %, AL 3 % of haplotypes. 86 % of people carry only
the five common structures. C4A and C4B GREx correlate −0.43, so M2 is not
badly collinear.

**Association code.** A planted β = −0.100 on 3,000 synthetic children, with
C4 alleles drawn from the panel, was recovered as −0.109 (95 % CI −0.147,
−0.071). The null control phenotype stayed null. Both CSD3 sbatch scripts were
run locally as plain bash with every path overridden, and all outputs were
produced.

## 6. Running it

**Laptop** (tests only; no ABCD data involved):

```bash
bash c4_imputation/00_fetch_resources.sh     # panel + Beagle + conform-gt, sha256-checked
PLINK=<plink 1.9> JAVA=<java 17> BCFTOOLS=bcftools RSCRIPT=Rscript PY=python \
  bash c4_imputation/smoketest/run_smoketest.sh    # ~2 min; ends "smoketest PASSED"
python c4_imputation/smoketest/panel_summary.py    # panel table + power table
```

On the laptop, the Claude Science env `c4-imputation` has bcftools and java
(`$CONDA_PREFIX/lib/jvm/bin/java`). PLINK 1.9 is in `abcd-genetics`; lme4 and
lmerTest are in `r`.

**CSD3**:

```bash
ssh -o ControlPersist=4h login-q-1.hpc.cam.ac.uk
bash c4_imputation/run_on_csd3.sh setup    # git pull; fetch resources; build ~/rds/hpc-work/envs/c4 (bcftools, java)
bash c4_imputation/run_on_csd3.sh submit   # c4_impute.sbatch (~30 min) -> c4_assoc.sbatch (~10 min)
bash c4_imputation/run_on_csd3.sh pull     # results/*.tsv
```

The cluster paths are defaults inside the two sbatch scripts, and every one can
be overridden by an environment variable:

- **Genotypes:** `GENO_ARRAY` from `genetic_analysis/config.local.sh`.
- **Phenotypes:** `results_70tab{,_hcp}/prs_final_1lmm/pheno/`, plus the HCP
  per-region export.
- **EUR keep list:** `legacy/hpc/work/results/ancestry/eur_anchor.keep`.
- **Score:** `genetic_analysis/work/scores_scz2025/SBayesRC/SCZ25_EUR`.
- **R:** the GENESIS R.

`setup_csd3_env.sh` looks for micromamba, mamba or conda. If none is on the
PATH, `module load` one first.

**Check before reading M1** (and record it here):

1. `results/conform_marker_tally.txt`: how many array SNPs were aligned to the
   panel (expect ~2–3 k), and how many were strand-flipped or removed.
2. `results/c4_imputation_summary.tsv`, `keep` group (EUR): structure
   frequencies should resemble the panel rows in the same file. A large
   mismatch means the array does not tag the panel haplotypes, and M1 is not
   interpretable.
3. `results/mhc_in_scz_scores.tsv`: the number of MHC SNPs in each SCZ score.
   If they are many, M4 partly conditions C4 on itself; say so.

## 7. Files

| file | what |
|:--|:--|
| `00_fetch_resources.sh` | downloads + checksums the panel, Beagle 5.5, conform-gt into `resources/` (gitignored) |
| `01_extract_mhc.sh` | PLINK 1.9 → MHC VCF (GRCh37 6:24.89–33.89 Mb; SNPs, geno < 5 %, MAF ≥ 0.5 %) |
| `02_impute_c4.sh` | conform-gt + Beagle; `MAP_CM_PER_BP` (default 1e-6) |
| `03_c4_grex.py` | C4 marker → expected AL/AS/BL/BS, copy numbers, GREx, QC; per-subject table (never committed) + aggregate summary |
| `04_c4_assoc.R` | models M1–M8 for one phenotype cell |
| `05_mhc_in_scores.py` | MHC SNP counts in the SCZ score weight files |
| `c4_impute.sbatch`, `c4_assoc.sbatch` | CSD3 jobs (steps 1–3; steps 4–5) |
| `run_on_csd3.sh` | setup / submit / status / pull from the laptop |
| `env/` | CSD3 env spec + builder |
| `smoketest/` | leave-fold-out imputation test, planted-effect fixture, map-scale grid, panel/power summary |
| `results/` | committed summary tables only |

Per-subject C4 calls (`work/c4_calls.tsv`), the MHC extracts and imputed VCFs
are individual-level genotype-derived data. They stay in `work/` (gitignored)
and are never committed (genetic_analysis/README_HPC.md rule 16).

## 8. Caveats to carry into the write-up

- **MHC confounding.** C4 alleles are in LD with HLA alleles and other MHC
  variants, and MHC haplotype frequencies differ steeply by ancestry. The EUR
  restriction and 10 PCs address stratification; LD with HLA is not addressed.
  A positive M1 is "the C4 haplotype", not necessarily C4 expression.
- **Weights.** The GREx weights are from post-mortem adult brain. M7
  (weight-free copy numbers) is the check.
- **Imputation.** The panel is small (222 haplotypes), and C4B is imputed
  poorly (r 0.60), so M2's C4B coefficient is attenuated in power.
- **Non-EUR children** are imputed (all 11,670 genotyped) but not analysed; the
  CEU panel is not valid for them.
