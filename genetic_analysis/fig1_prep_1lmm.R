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
# Outputs (aggregates only, committed):
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
fwrite(rs, "genetic_analysis/fig1_inputs/hcp70_global_slope_reliability_1lmm.csv")
cat("slope sd (mm/yr):", sqrt(tau1), " r(int,slope):", c01 / sqrt(tau0 * tau1), "\n")
