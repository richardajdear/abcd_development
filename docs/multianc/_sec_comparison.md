Searched August 2026. The relevant comparison set is ABCD genetics
methodology papers and the release's own documentation.

### What others did

| study | ancestry handling | relatedness | GRM | MAF |
|---|---|---|---|---|
| **ABCD genotype resource** (Loughnan et al., *Behav Genet* 2023) | **PC-AiR** PCs; **1000G Phase 3 (2,504) projected** for AFR/AMR/EAS/EUR/SAS reference | **PC-AiR** → 8,005 unrelated, 3,384 projected | **PC-Relate** | not stated |
| **Cognitive heritability in ABCD** (Loughnan et al. 2023) | pooled multi-ancestry, 10 PCs as covariates; **not stratified** | relatives **retained** and modelled (FEMA sparse block) | **PC-Relate**, 158,103 pruned SNPs | not stated |
| **Subcortical volumes, multi-ethnic** (2024) | pooled multi-ethnic (n=7,234) **plus** a EUR sensitivity arm (n=3,669) | family structure removed (threshold unstated) | PLINK GRM, AdjHE-RE | **0.05** |
| **Most ABCD imaging-genetics** | European-ancestry subset only, PCA outlier exclusion | varies | GCTA | typically 0.01 |
| **This analysis** | pooled multi-ancestry **plus** stratified arms; in-sample PCs; k-means strata validated against self-report | **PC-AiR** (the release's own set) | GCTA GREML | 0.01, **plus a 0.001 arm** |

### Where we agree with the field

- **PC-AiR for relatedness in ABCD is what the release itself uses**, and what
  the PAGE consortium recommends for diverse-ancestry studies. Our 8,181-child
  unrelated set is essentially the release's own (they report 8,005 at an earlier
  data freeze). **We arrived at it by finding `--grm-cutoff` broken rather than
  by following the documentation — but we landed in the right place.**
- **Pooled-plus-EUR-sensitivity** is exactly the design of the subcortical-volume
  paper. Our §8.9(d) EUR arm serves the same role.
- **The problem is documented.** `FastSparseGRM` (2024) exists because "in
  multi-ancestry heterogeneous populations, GRMs are confounded with population
  structure," and GCTA's own documentation warns that sub-setting a GRM after
  relatedness pruning changes the sample allele frequencies. **Our finding is not
  novel; it is a known trap that the conventional recipe walks straight into.**

### Where we differ, and whether it matters

| difference | ours | theirs | assessment |
|---|---|---|---|
| **ancestry labels** | k-means on in-sample PCs + `abcd_eur` anchor + self-report validation | 1000 Genomes projection | **theirs is better.** We could not do it — no multi-ancestry 1000G panel on this account. Ours is validated but unanchored to an external reference. **Worth fixing if a panel becomes available.** |
| **GRM estimator** | GCTA additive GRM | PC-Relate | **theirs is more correct for multi-ancestry**, computing individual-specific allele frequencies from PCs. We use GCTA for comparability with this project's published EUR numbers; switching would make every h² incomparable with §2. A deliberate trade, and the compromise flagged in §7. |
| **relatives** | pruned (PC-AiR set) | retained and modelled (FEMA) | Retaining relatives keeps N and is more efficient, but requires a mixed model that separates shared environment from additive genetics. Pruning is the conservative choice and is what GCTA-GREML expects. Ours costs sample size; theirs costs an assumption. |
| **MAF** | 0.01 (+0.001 arm) | 0.05 (subcortical paper) | 0.05 is stricter than convention. Our 0.001 arm exists because this project's own published estimate used 0.001 — a comparison that would otherwise confound frequency floor with variant set. |
| **stratified heritability** | run, and reported as uninformative outside EUR | generally not run | We report SEs of 0.23–0.45 in the non-European strata rather than omitting them. **Reporting an uninformative estimate honestly is better than reporting only the arms that worked.** |

### The tension the field has not resolved, and where we sit

The ABCD behaviour-genetics review argues **against** European-only restriction —
that polygenic scores derived in European samples transfer poorly and that
excluding non-European participants "may exacerbate health disparities." Most
ABCD imaging-genetics papers nonetheless restrict to European ancestry, because
population stratification biases the estimators.

Both are right. Restriction is *statistically* safer and *scientifically*
narrower. Our position: **include everyone, and be explicit about which
conclusions the inclusion actually supports.** Concretely — the pooled
heritability estimate gains real precision from the extra 1,523 unrelated
children, while the non-European stratified estimates remain uninformative
(SE 0.23–0.45) because ABCD does not contain enough children per group. Doubling
the sample does not by itself make the analysis equitable; it makes the pooled
estimate better and leaves the group-specific questions where they were.

**Sources:**
- [Genotype data and derived genetic instruments of ABCD](https://pmc.ncbi.nlm.nih.gov/articles/PMC10635818/)
- [Heritability Estimation of Cognitive Phenotypes in the ABCD Study Using Mixed Models](https://pmc.ncbi.nlm.nih.gov/articles/PMC10154273/)
- [Heritability estimation of subcortical volumes in a multi-ethnic multi-site cohort study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10802572/)
- [ABCD Behavior Genetics: Twin, Family, and Genomic Studies](https://pmc.ncbi.nlm.nih.gov/articles/PMC10833231/)
- [Scalable analysis of large multi-ancestry biobanks (FastSparseGRM)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11601839/)
- [GENESIS PC-AiR / PC-Relate vignette](https://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html)
- [ABCD Study genetics documentation](https://docs.abcdstudy.org/latest/documentation/non_imaging/gn.html)
- [GCTA documentation](https://yanglab.westlake.edu.cn/software/gcta/)
