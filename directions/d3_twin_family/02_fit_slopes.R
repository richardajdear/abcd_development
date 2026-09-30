#!/usr/bin/env Rscript
# 02_fit_slopes.R -- single LMM per whole-cortex measure (the Figure-1 / genetics trait).
#   mean ~ age_c + sex + (1 + age_c | subject) + (1 | site)   (REML; no family effect)
# Same model as genetic_analysis/fig1_prep_1lmm.R. The family effect is deliberately
# absent (README "Design rules"): it would shrink co-twin BLUPs towards each other.
# Input : work/scan_means_{ct,t1t2}.csv   (01_scan_means.py)
# Output: work/blups_{ct,t1t2}.csv        per-child random intercept and slope (gitignored)
#         results/lmm_variance_components.tsv, results/slope_reliability.tsv  (aggregate)
# Run from the repo root: Rscript directions/d3_twin_family/02_fit_slopes.R   (env r)
suppressMessages({library(lme4); library(data.table)})
args <- commandArgs(FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
W <- file.path(HERE, "work"); RES <- file.path(HERE, "results")
vcs <- list(); rels <- list()
for (tag in c("ct", "t1t2")) {
  d <- fread(file.path(W, sprintf("scan_means_%s.csv", tag)))
  m <- lmer(mean_val ~ age_c + sex + (1 + age_c | subject) + (1 | site), data = d, REML = TRUE,
            control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  vc <- as.data.table(VarCorr(m)); vc[, measure := tag]; vcs[[tag]] <- vc
  re <- ranef(m, condVar = TRUE)$subject
  pv <- attr(re, "postVar")                       # 2 x 2 x n conditional covariances
  tau1 <- vc[grp == "subject" & var1 == "age_c" & is.na(var2), vcov]
  rel <- 1 - pv[2, 2, ] / tau1                    # per-child slope reliability
  nv <- d[, .(n_visits = .N), by = subject][match(rownames(re), subject), n_visits]
  rels[[tag]] <- data.table(measure = tag, n_visits = nv, rel = rel)[
    , .(n_children = .N, mean_rel = mean(rel), median_rel = median(rel)), by = .(measure, n_visits)]
  fwrite(data.table(subject = rownames(re), re_intercept = re[["(Intercept)"]], re_slope = re[["age_c"]],
                    slope_rel = rel), file.path(W, sprintf("blups_%s.csv", tag)))
  cat(tag, "fixed age_c =", fixef(m)[["age_c"]], "; n =", nrow(re), "\n")
}
fwrite(rbindlist(vcs), file.path(RES, "lmm_variance_components.tsv"), sep = "\t")
fwrite(rbindlist(rels)[order(measure, n_visits)], file.path(RES, "slope_reliability.tsv"), sep = "\t")
