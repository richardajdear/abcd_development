#!/usr/bin/env Rscript
# D1 step 7b -- do time-varying scanner terms improve per-child slopes? Applied to the HCP-MMP structural measures
# (thickness, area, volume, GWC, T1 gm50, T1 wm1), T1w/T2w (DK) and 8 diffusion measures. Input: 09_smri_scan_means.py.
#
#   M0  y ~ age_c + sex + (1 + age_c | subject) + (1 | site)                 (Figure 1)
#   M1  M0 + (1 | serial)                                                     scanner unit
#   M2  M0 + (1 | epoch)                                                      scanner unit x software version
# y is z-scored over scans within measure, so slopes are in SD units per year.
# Epochs are partly confounded with age (calendar time ~ age in a single-age cohort; within-epoch age SD 0.88 y vs
# 2.45 y overall). An epoch term can therefore absorb part of the true mean trend; the fixed age slope is reported
# under every model as the diagnostic.
# Per measure x model: fixed age slope, between-child slope SD, residual SD, scanner variance share, mean per-child
# slope reliability (all, 2/3/4 visits), R2 of site on child slopes, within-site r of the child slope with the M0
# thinning slope, r of the slope with the same measure's M0 slope. PRS: Figure-1 model + site (as 03), PRS-CS and
# SBayesRC, EUR and pooled.
# Outputs: tables/scanner_terms_summary.tsv, tables/scanner_terms_prs.tsv; work/blups_scanner_terms.csv.gz.
# Run (env r, repo root).
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
W <- "directions/d1_multimodal/work"; OUT <- "directions/d1_multimodal/tables"
ctrl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
S <- fread(file.path(W, "scan_means_smri.csv.gz"))
S[, y := (mean_val - mean(mean_val)) / sd(mean_val), by = metric]
RHS <- c(M0 = "", M1 = " + (1 | serial)", M2 = " + (1 | epoch)")
metrics <- unique(S$metric)
SUM <- list(); BL <- list()
for (mt in metrics) for (mn in names(RHS)) {
  d <- S[metric == mt]
  m <- lmer(as.formula(paste0("y ~ age_c + sex + (1 + age_c | subject) + (1 | site)", RHS[[mn]])), data = d,
            REML = TRUE, control = ctrl)
  vc <- as.data.table(VarCorr(m))
  tau0 <- vc[grp == "subject" & var1 == "(Intercept)" & is.na(var2), vcov]
  tau1 <- vc[grp == "subject" & var1 == "age_c" & is.na(var2), vcov]
  s_sc <- if (mn == "M0") 0 else vc[grp %in% c("serial", "epoch"), vcov]
  s_site <- vc[grp == "site", vcov]; s_res <- vc[grp == "Residual", vcov]
  re <- ranef(m, condVar = TRUE)$subject; rel <- 1 - attr(re, "postVar")[2, 2, ] / tau1
  b <- data.table(subject = rownames(re), metric = mt, model = mn, slope = re[["age_c"]], rel = rel)
  b <- d[, .(nv = .N, site = site[1]), by = subject][b, on = "subject"]
  BL[[length(BL) + 1]] <- b
  SUM[[length(SUM) + 1]] <- data.table(metric = mt, model = mn, n_children = nrow(re), singular = isSingular(m),
    age_slope_sd_per_yr = fixef(m)[["age_c"]], child_slope_sd = sqrt(tau1), resid_sd = sqrt(s_res),
    scanner_var_share = s_sc / (s_sc + s_site + tau0 + s_res), site_var_share = s_site / (s_sc + s_site + tau0 + s_res),
    rel_mean = mean(rel), rel2 = mean(rel[b$nv == 2]), rel3 = mean(rel[b$nv == 3]), rel4 = mean(rel[b$nv == 4]),
    site_r2_slope = summary(lm(slope ~ factor(site), b))$r.squared)
  cat(sprintf("%-18s %s  age %+.3f SD/yr  slopeSD %.3f  resid %.3f  scanner %.3f  rel %.3f (%.2f/%.2f/%.2f)  siteR2 %.3f%s\n",
              mt, mn, fixef(m)[["age_c"]], sqrt(tau1), sqrt(s_res), SUM[[length(SUM)]]$scanner_var_share,
              mean(rel), SUM[[length(SUM)]]$rel2, SUM[[length(SUM)]]$rel3, SUM[[length(SUM)]]$rel4,
              SUM[[length(SUM)]]$site_r2_slope, ifelse(isSingular(m), " SINGULAR", "")), sep = "")
}
B <- rbindlist(BL); fwrite(B, file.path(W, "blups_scanner_terms.csv.gz"))
SM <- rbindlist(SUM)
thin0 <- B[metric == "smri__thk_hcp" & model == "M0", .(subject, thin = -slope, site)]
thin0[, thin_rs := thin - mean(thin), by = site]
for (i in seq_len(nrow(SM))) {
  b <- B[metric == SM$metric[i] & model == SM$model[i]]; b[, s_rs := slope - mean(slope), by = site]
  x <- thin0[b, on = "subject", nomatch = 0]
  b0 <- B[metric == SM$metric[i] & model == "M0", .(subject, s0 = slope)]
  SM[i, `:=`(r_thinfast_within_site = cor(x$s_rs, x$thin_rs), r_with_own_M0 = cor(b0[b, on = "subject"]$s0, b$slope))]
}
fwrite(SM, file.path(OUT, "scanner_terms_summary.tsv"), sep = "\t")

ch <- fread("genetic_analysis/work/results_70tab_hcp/d9_sex_moderation/children.csv")
ch[, `:=`(sex = factor(sex, levels = c("F", "M")), FID = factor(FID), site = factor(site), age_c = baseline_age)]
cov0 <- paste(c("sex", "age_c", paste0("PC", 1:10), "site"), collapse = " + ")
grid <- CJ(disorder = c("scz", "mdd"), method = c("prscs", "sbrc"), arm = c("eur", "pooled"))
P <- list()
for (mt in metrics) for (mn in names(RHS)) {
  y <- B[metric == mt & model == mn, .(subject, yv = slope)]
  if (mt == "smri__thk_hcp") y[, yv := -yv]                 # thickness reported as faster thinning
  dd <- y[ch, on = "subject", nomatch = 0]
  for (i in seq_len(nrow(grid))) {
    g <- grid[i]; col <- sprintf("%s_%s_%s", g$disorder, g$method, g$arm)
    d <- dd[!is.na(get(col)) & !is.na(yv) & !is.na(PC1)]
    if (g$arm == "eur") d <- d[eur == 1]
    d[, `:=`(prs = as.numeric(scale(get(col))), yz = as.numeric(scale(yv)))]
    s <- coef(summary(lmer(as.formula(paste("yz ~ prs +", cov0, "+ (1 | FID)")), d, control = lmerControl(calc.derivs = FALSE))))
    P[[length(P) + 1]] <- data.table(metric = mt, model = mn, disorder = toupper(g$disorder), method = g$method, arm = g$arm,
                                     est = s["prs", 1], se = s["prs", 2], p = s["prs", 5], n = nrow(d))
  }
  cat(".")
}
P <- rbindlist(P); fwrite(P, file.path(OUT, "scanner_terms_prs.tsv"), sep = "\t")
cat("\ndone\n")
