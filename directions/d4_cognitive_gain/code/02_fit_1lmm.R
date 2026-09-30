#!/usr/bin/env Rscript
# Single LMM on the per-scan cortical mean -> per-child random intercept and slope.
# Same model as genetic_analysis/fig1_prep_1lmm.R (the primary trait, global_slope_1lmm):
#   mean_ct ~ age_c + sex + (1 + age_c | subject) + (1 | site)
# Writes work/hcp70_1lmm_blups.csv (individual-level; never committed).
# Run from repo root:  Rscript directions/d4_cognitive_gain/code/02_fit_1lmm.R   (env r)
suppressMessages({library(lme4); library(data.table)})
W <- "directions/d4_cognitive_gain/work"
d <- fread(file.path(W, "hcp70_scan_means_cov.csv"))
m <- lmer(mean_ct ~ age_c + sex + (1 + age_c | subject) + (1 | site), data = d, REML = TRUE,
          control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
print(as.data.frame(VarCorr(m))); print(fixef(m))
re <- ranef(m, condVar = TRUE)$subject
cv <- attr(re, "postVar")
tau1 <- as.data.frame(VarCorr(m))$vcov[2]
fwrite(data.table(subject = rownames(re), re_intercept = re[["(Intercept)"]], re_slope = re[["age_c"]],
                  slope_reliability = 1 - cv[2, 2, ] / tau1),
       file.path(W, "hcp70_1lmm_blups.csv"))
cat("wrote", nrow(re), "children\n")
