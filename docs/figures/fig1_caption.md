**Figure 1 | Polygenic risk for schizophrenia predicts the rate, not the baseline level, of adolescent cortical thinning.**

- **a** Group maps per parcel (bilateral mean, left hemisphere shown): baseline thickness and mean per-child thinning rate (-33 to -0.1 µm / year).
- **b** Age at scan by visit; 8,716 children, 26,946 scans; 2 / 3 / 4 scans per child for 2,205 / 3,508 / 3,003 children.
- **c** Per-scan cortical mean vs age for 250 random children (grey), three children with four scans (colour), population OLS trend -19 µm / year (black). Side densities: model-implied distribution of child-level thickness at ages 9 (left) and 17 (right) from the single LMM (mean 2.82 and 2.67 mm, SD 0.070 and 0.072 mm; site variance excluded).
- **d** Slope reliability 1 − v/τ², where v is a child's conditional (posterior) variance of the slope random effect and τ² the between-child slope variance: boxes, per-parcel LMMs (358 parcels; each value the mean over children; medians 0.15 / 0.21 / 0.25 for 2 / 3 / 4 scans); squares, the single LMM on the cortical mean (mean over children 0.17 / 0.32 / 0.38).
- **e, f** PRS association with thinning rate (e) and baseline thickness (f); β per SD with 95% CI. Schizophrenia (2025 GWAS): thinning 3/3 EUR / 3/3 pooled, baseline 0/3 / 0/3. Alzheimer's 3/3 with APOE, 0/3 without; depression 2/3 / 2/3; education 2/3 (opposite sign); autism 0/3.
- **g** MAGMA enrichment in SCZ curated loci (Trubetskoy 2022 ST12, 455 genes) and MDD high-confidence genes (MDD2025 prioritised; 213 genes), p thinning / baseline; p printed on points with p < 0.05. SCZ: EUR 0.005 / 1e-04, pooled 0.13 / 0.24. MDD: EUR 0.59 / 0.08, pooled 0.39 / 0.37. EUR arm on 1000 Genomes EUR LD; pooled arm on the ABCD analysis sample as its own LD reference.

Methods
- ABCD release 7.0, HCP-MMP1.0 (358 parcels; 1 lh polygon without a value drawn grey); 8,716 children with >= 2 QC-passing scans.
- Slope model: thickness ~ age + sex + (1 + age | child) + (1 | site), per parcel (a, d) and on the per-scan cortical mean (c, e–g); the trait is the child's age slope.
- PRS: C+T, PRS-CS, SBayesRC; EUR arm scored with European discovery GWAS, pooled arm with multi-ancestry GWAS, z-scored within ancestry cluster; model score + age + sex + 10 PCs + (1 | family); filled = p < 0.05 (C+T corrected over its thresholds).
- MAGMA competitive gene-set test; EUR arm on 1000 Genomes EUR LD; pooled arm on the ABCD analysis sample as LD reference.
- Sensitivity analyses (DK parcellation, mean of per-parcel slopes, PGC3 2022 GWAS, SBayesR, per-ancestry strata) in Supplementary Information.
