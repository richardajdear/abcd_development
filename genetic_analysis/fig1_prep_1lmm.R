#!/usr/bin/env Rscript
# fig1_prep_1lmm.R -- single-LMM quantities for Figure 1 panels c and d.
#
# Input (individual-level, NEVER committed): per-scan HCP-MMP cortical mean with
# covariates, written by
#   python -c "import pandas as pd; d='out/thickness_hcp_70_aa6e91efba82';
#     mt=pd.read_parquet(f'{d}/model_table.parquet',columns=['subject','visit','value','age','age_c','sex','site','n_visits']);
#     mt.groupby(['subject','visit'],observed=True).agg(mean_ct=('value','mean'),age=('age','first'),age_c=('age_c','first'),
#       sex=('sex','first'),site=('site','first'),n_visits=('n_visits','first')).reset_index()
#       .to_csv('/tmp/hcp70_scan_means_cov.csv',index=False)"
# Model: mean_ct ~ age_c + sex + (1 + age_c | subject) + (1 | site)  (the trait's LMM; age centre 12.797)
# Also writes /tmp/hcp70_1lmm_blups.csv (per-child REs; never committed) for fig1_prep_cbcl.py.
# Outputs (aggregates only, committed):
#   fig1_inputs/hcp70_child_density.csv, hcp70_child_summary.csv   per-child phenotype
#       densities (thickness at the age centre; slope) for panel b
#   fig1_inputs/hcp70_ct_at_age.csv                     model-implied child-level CT
#       distribution at ages 9 and 17: mean = fixed part at the sample sex mix,
#       sd = sqrt(tau0 + 2 t c01 + t^2 tau1); site variance excluded (scanner, not child)
#   fig1_inputs/hcp70_global_slope_reliability_1lmm.csv  per-child slope reliability
#       1 - condVar(slope RE) / tau1 (same formula as phenotype.slope_reliability),
#       summarised by scans per child
# Run: Rscript genetic_analysis/fig1_prep_1lmm.R   (env r; needs lme4, data.table)
suppressMessages({library(lme4); library(data.table)})
d <- fread("/tmp/hcp70_scan_means_cov.csv")
m <- lmer(mean_ct ~ age_c + sex + (1 + age_c | subject) + (1 | site), data = d, REML = TRUE,
          control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
vc <- as.data.frame(VarCorr(m)); print(vc)
fx <- fixef(m); print(fx)
CEN <- 12.797
# per-child slope reliability, same formula as the parcel boxes: 1 - condvar / tau2
re <- ranef(m, condVar = TRUE)$subject
pv <- attr(re, "postVar")                     # 2 x 2 x n
tau1 <- vc$vcov[vc$grp == "subject" & vc$var1 == "age_c" & is.na(vc$var2)]
tau0 <- vc$vcov[vc$grp == "subject" & vc$var1 == "(Intercept)" & is.na(vc$var2)]
c01  <- vc$vcov[vc$grp == "subject" & !is.na(vc$var2)]
rel  <- 1 - pv[2, 2, ] / tau1
nv <- d[, .(n_visits = first(n_visits)), by = subject]
r <- data.table(subject = rownames(re), reliability = rel)[nv, on = "subject"]
rs <- r[, .(n = .N, mean = mean(reliability), median = median(reliability), q25 = quantile(reliability, .25),
            q75 = quantile(reliability, .75)), by = n_visits][order(n_visits)]
print(rs)
# model-implied distribution of child-level CT at ages 9 and 17 (fixed part at the
# sample sex mix; child random effects; site variance excluded = scanner, not biology)
psex <- mean(d$sex == sort(unique(d$sex))[2])
at <- rbindlist(lapply(c(9, 17), function(a) {
  t <- a - CEN
  data.table(age = a, mean_mm = fx[["(Intercept)"]] + fx[["age_c"]] * t + fx[[3]] * psex,
             sd_mm = sqrt(tau0 + 2 * t * c01 + t^2 * tau1))
}))
print(at)
fwrite(at, "genetic_analysis/fig1_inputs/hcp70_ct_at_age.csv")
# per-child random effects for fig1_prep_cbcl.py -- individual-level, /tmp ONLY
fwrite(data.table(subject = rownames(re), re_intercept = re[["(Intercept)"]],
                  re_slope = re[["age_c"]]), "/tmp/hcp70_1lmm_blups.csv")
# per-child phenotype distributions for Figure 1b (aggregate: density grid only).
# thickness = fixed intercept + sex + child intercept, i.e. thickness at the
# age centre (the 'baseline_thickness' phenotype); slope = fixed + child slope.
sx <- d[, .(sex = first(sex)), by = subject][match(rownames(re), subject)]
ct_c  <- fx[["(Intercept)"]] + fx[[3]] * (sx$sex == names(fx)[3] |> sub(pattern = "^sex", replacement = "")) + re[["(Intercept)"]]
sl_c  <- (fx[["age_c"]] + re[["age_c"]]) * 1000                  # um / yr
dens1 <- function(v, name) { k <- density(v, n = 256); data.table(var = name, x = k$x, density = k$y) }
# per-child CT / ΔCT for the Figure 1e scatter: individual-level -> gitignored
# sid = 0-based rank in the codepoint-sorted subject list, the key hcp70_scans.csv uses
sids <- match(rownames(re), sort(unique(rownames(re)), method = "radix")) - 1L
fwrite(data.table(sid = sids, ct_mm = round(ct_c, 4), dct_um_per_yr = round(sl_c, 3)),
       "genetic_analysis/fig1_inputs/hcp70_child_traits.csv")
# committed fallback: 2-D counts, cells with n < 10 suppressed
bx <- cut(ct_c, seq(2.4, 3.1, 0.01), labels = FALSE); by <- cut(sl_c, seq(-32, -8, 0.4), labels = FALSE)
h2 <- data.table(ct_bin = 2.4 + (bx - 0.5) * 0.01, dct_bin = -32 + (by - 0.5) * 0.4)[, .N, by = .(ct_bin, dct_bin)][N >= 10]
fwrite(h2[order(ct_bin, dct_bin)], "genetic_analysis/fig1_inputs/hcp70_child_hist2d.csv")
cat("r(CT, dCT) =", round(cor(ct_c, sl_c), 3), "  hist2d cells kept:", nrow(h2), "children in kept cells:", sum(h2$N), "\n")
fwrite(rbind(dens1(ct_c, "thickness_mm"), dens1(sl_c, "slope_um_per_yr")),
       "genetic_analysis/fig1_inputs/hcp70_child_density.csv")
fwrite(data.table(var = c("thickness_mm", "slope_um_per_yr"), n = length(ct_c),
                  mean = c(mean(ct_c), mean(sl_c)), sd = c(sd(ct_c), sd(sl_c)),
                  model_sd = c(sqrt(tau0), sqrt(tau1) * 1000),
                  frac_thinning = c(NA, mean(sl_c < 0)), age_centre = CEN),
       "genetic_analysis/fig1_inputs/hcp70_child_summary.csv")
fwrite(rs, "genetic_analysis/fig1_inputs/hcp70_global_slope_reliability_1lmm.csv")
cat("slope sd (mm/yr):", sqrt(tau1), " r(int,slope):", c01 / sqrt(tau0 * tau1), "\n")
