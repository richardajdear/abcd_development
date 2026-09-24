**Figure 1 | Polygenic risk for schizophrenia predicts the rate, not the baseline level, of adolescent cortical thinning.**

- **a** Group maps per parcel (bilateral mean, left hemisphere shown): baseline thickness and mean per-child thinning rate (-33 to -0.1 µm / year).
- **b** Age at scan by visit; 8,716 children, 26,946 scans; 2 / 3 / 4 scans per child for 2,205 / 3,508 / 3,003 children.
- **c** Per-scan cortical mean vs age for 250 random children (grey), three with four scans (blue), population OLS trend -19 µm / year (black).
- **d** Slope reliability: single-parcel model-based reliability across 358 parcels (boxes, medians 0.15 / 0.21 / 0.25 for 2 / 3 / 4 scans) vs split-half consistency of the cortex-wide slope (squares; 0.91 / 0.89 / 0.89).
- **e, f** PRS association with thinning rate (e) and baseline thickness (f); β per SD with 95% CI. Schizophrenia (2025 GWAS): thinning 3/3 EUR / 3/3 pooled, baseline 0/3 / 0/3. Alzheimer's 3/3 with APOE, 0/3 without; depression 2/3 / 2/3; education 2/3 (opposite sign); autism 0/3.
- **g** MAGMA enrichment in SCZ curated loci (Trubetskoy 2022 ST12, 455 genes) and MDD loci (1,881 genes), p thinning / baseline. SCZ: EUR 0.005 / 1e-04, pooled 0.13 / 0.24. MDD: EUR 0.34 / 0.06, pooled 0.39 / 0.49. EUR arm on 1000 Genomes EUR LD; pooled arm on the ABCD analysis sample as its own LD reference.

Methods
- ABCD release 7.0, HCP-MMP1.0 (358 parcels; 1 lh polygon without a value drawn grey); 8,716 children with >= 2 QC-passing scans.
- Slope model: thickness ~ age + sex + (1 + age | child) + (1 | site), per parcel (a, d) and on the per-scan cortical mean (c, e–g); the trait is the child's age slope.
- PRS: C+T, PRS-CS, SBayesRC; EUR arm scored with European discovery GWAS, pooled arm with multi-ancestry GWAS, z-scored within ancestry cluster; model score + age + sex + 10 PCs + (1 | family); filled = p < 0.05 (C+T corrected over its thresholds).
- MAGMA competitive gene-set test; EUR arm on 1000 Genomes EUR LD; pooled arm on the ABCD analysis sample as LD reference.
- Sensitivity analyses (DK parcellation, mean of per-parcel slopes, PGC3 2022 GWAS, SBayesR, per-ancestry strata) in Supplementary Information.
