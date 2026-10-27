# D1 step 9 -- per-child phenotypes for the geometry GWAS on CSD3 (CSD3_GWAS_INSTRUCTIONS.md).
# Same construction as the thinning GWAS phenotype (global_slope_1lmm): one LMM on the Figure-1 scan set,
#   y ~ age_c + sex + (1 + age_c | subject) + (1 | site), age centre 12.797, REML,
# per-child random slope and intercept (BLUPs). Five traits:
#   thk          = whole-cortex mean thickness (reproduction control; must give r ~ 1 with global_slope_1lmm)
#   area, vol    = HCP-MMP total surface area / GM volume, z-scored over scans (alone)
#   thk_area, thk_vol_area = per-scan z-sum composites (step 8), as used in the Figure-1 variants
# thk / thk_area / thk_vol_area BLUPs are read from the Figure-1-variant fits (genetic_analysis/fig1_variants);
# area and vol are fitted here with the identical model.
# Output (individual-level, gitignored): directions/d1_multimodal/work/gwas_pheno_geom.tsv
# Run from the repo root: Rscript directions/d1_multimodal/13_export_gwas_pheno.R   (env r)
suppressMessages({library(data.table); library(lme4)})
W <- "directions/d1_multimodal/work"; V <- "genetic_analysis/work/fig1_variants"
keep <- fread(file.path(V, "thk", "scan_means_cov.csv"))[, .(subject, visit)]
S <- fread(file.path(W, "scan_means_smri.csv.gz"))[metric %in% c("smri__area_hcp", "smri__vol_hcp")]
S <- merge(S, keep, by = c("subject", "visit"))
S[, y := (mean_val - mean(mean_val)) / sd(mean_val), by = metric]
fit <- function(d) {
  m <- lmer(y ~ age_c + sex + (1 + age_c | subject) + (1 | site), data = d, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  re <- ranef(m)$subject
  data.table(subject = rownames(re), re_intercept = re[["(Intercept)"]], re_slope = re[["age_c"]])
}
B <- list(area = fit(S[metric == "smri__area_hcp"]), vol = fit(S[metric == "smri__vol_hcp"]))
for (t in c("thk", "thk_area", "thk_vol_area")) B[[t]] <- fread(file.path(V, t, "blups.csv"))
cat(sprintf("scans %d; children per trait: %s\n", nrow(S) / 2,
            paste(sprintf("%s %d", names(B), sapply(B, nrow)), collapse = ", ")))
out <- NULL
for (t in c("thk", "area", "vol", "thk_area", "thk_vol_area")) {
  b <- copy(B[[t]]); setnames(b, c("re_intercept", "re_slope"), paste0(c("int_", "slope_"), t))
  out <- if (is.null(out)) b else merge(out, b, by = "subject")
}
sl <- grep("^slope_", names(out), value = TRUE)
cat("child-level r between slope BLUPs:\n"); print(round(cor(out[, ..sl]), 3))
fwrite(out, file.path(W, "gwas_pheno_geom.tsv"), sep = "\t")
cat(sprintf("wrote %s: %d children x %d traits (slope + intercept)\n", file.path(W, "gwas_pheno_geom.tsv"),
            nrow(out), length(sl)))
