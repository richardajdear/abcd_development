**Figure 1 | Polygenic risk for schizophrenia predicts the rate of adolescent
cortical thinning.** All panels use the HCP-MMP1.0 parcellation (358 parcels)
and ABCD release 7.0.

**a**, Age at scan by study visit for the 8,716 children with at least two
usable scans (26,946 scans; 2,205 / 3,508 / 3,003
children with 2 / 3 / 4 scans). **b**, Mean cortical thickness across parcels
per scan against age. Thin lines join the scans of 250 randomly drawn children;
three children with four scans are highlighted; the black line is the ordinary
least-squares population trend (-19 µm per year). Each child's
age slope from a linear mixed model is the trait analysed below. **c**,
Reliability of the slope estimate by scans per child. Boxes: model-based
reliability (variance of the true slope over variance of the estimate) of a
single parcel's slope across the 358 parcels (median 0.15 / 0.21 / 0.25
for 2 / 3 / 4 scans). Squares: split-half consistency of the cortex-wide mean
slope (left- vs right-hemisphere means, Spearman–Brown corrected; 0.91 / 0.89 / 0.89).
**d**, Association of polygenic scores with the thinning rate (single linear
mixed model on the per-scan cortical mean; β per SD of score, 95 % CI). Four
scoring methods per trait (C+T, PRS-CS, SBayesR, SBayesRC); circles, European
ancestry arm scored with a European discovery GWAS (n = 4,308); diamonds, all
ancestries scored with a multi-ancestry discovery GWAS and standardised within
ancestry cluster (n = 8,596); filled, p < 0.05 after correction across C+T
thresholds. Schizophrenia (2025 multi-ancestry GWAS): 4/4 methods
significant in the European arm (β -0.044 to -0.035) and 4/4 in the pooled
arm; depression 3/4 and 3/4; Alzheimer's disease with the APOE region
excluded 0/4; educational attainment 3/4 (opposite sign); autism 0/4.
**e**, Genome-wide association of the thinning rate, European arm (n = 4,205;
7.5 M imputed variants; λGC 0.99; 0 loci at p < 5 × 10⁻⁸;
SNP heritability from GREML 0.16 ± 0.04). *Points are placeholders until
the summary statistics are exported.* **f**, MAGMA competitive gene-set
enrichment of the thinning-rate GWAS (red) and of baseline thickness (grey)
in schizophrenia gene sets: the curated locus pool of the 2022 PGC GWAS
(Trubetskoy et al. Supplementary Table 12; 455 genes) and its
non-prioritised (354) and prioritised (101) subsets, genes
under genome-wide-significant peaks of the 2022 and 2025 GWAS, and genes
significant in the 2025 gene-level analysis. The curated pool is enriched for
the thinning rate (p = 0.0005) and for baseline thickness
(p = 0.0001) alike, carried by its non-prioritised genes
(p = 0.0003) rather than the prioritised ones (p = 0.33).

Methods notes
- Thickness parsed from release FreeSurfer surfaces into HCP-MMP parcels; scans passing the release QC code; children with ≥ 2 scans.
- Mixed model per parcel and for the cortex-wide mean: thickness ~ age + sex + (1 + age | child) + (1 | site); the trait is the child's age slope.
- Polygenic scores from four methods; discovery GWAS matched to the target arm; pooled-arm scores z-scored within genetic-ancestry cluster; association model score + age + sex + 10 PCs + (1 | family).
- Sensitivity analyses (DK parcellation, mean of per-parcel slopes as the trait, PGC3 2022 discovery GWAS, per-ancestry strata) in Supplementary Figures.
