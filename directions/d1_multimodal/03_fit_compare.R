#!/usr/bin/env Rscript
# D1 step 3 -- single-LMM slopes for every DTI/RSI whole-cortex measure, compared with thickness and T1w/T2w.
#
#   mean ~ age_c + sex + (1 + age_c | subject) + (1 | site)      (REML; exactly the Figure-1 / D3 model)
#
# Inputs : work/scan_means_dmri.csv.gz (02_dmri_scan_means.py); D3 work/scan_means_{ct,t1t2}.csv; d9 children.csv
# Outputs (aggregates only):
#   tables/slope_summary.tsv     per metric: n, fixed age slope (raw, %/yr of the value at 12.8 y, and in units of
#                                the between-child SD at 12.8 y), between-child slope SD, residual SD, scan-level
#                                site ICC, R2 of site on child slopes, mean slope reliability at 2/3/4 visits,
#                                r of the child slope with faster thinning and with the T1w/T2w slope (all and
#                                within site), r of the child intercept with baseline thickness
#   tables/slope_cor.tsv         Spearman matrix of child slopes (all metrics, children with every metric)
#   tables/dmri_prs.tsv          Figure-1 PRS model on each child slope, base and + site (as 01_t1t2_prs.R)
# Individual-level BLUPs -> work/blups_dmri.csv.gz (gitignored).
# Run (env r, repo root):  Rscript directions/d1_multimodal/03_fit_compare.R
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
W <- "directions/d1_multimodal/work"; OUT <- "directions/d1_multimodal/tables"; D3W <- "directions/d3_twin_family/work"
ctrl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)

S <- fread(file.path(W, "scan_means_dmri.csv.gz"))
for (tag in c("ct", "t1t2")) {
  x <- fread(file.path(D3W, sprintf("scan_means_%s.csv", tag)))
  x[, metric := c(ct = "smri__thk_hcp", t1t2 = "smri__t1t2__gm")[[tag]]]
  S <- rbind(S, x[, .(subject, visit, metric, mean_val, age, age_c, sex, site, family_id, n_visits)])
}
metrics <- unique(S$metric)
cat(length(metrics), "metrics\n")

SUM <- list(); BL <- list()
for (mt in metrics) {
  d <- S[metric == mt]
  m <- lmer(mean_val ~ age_c + sex + (1 + age_c | subject) + (1 | site), data = d, REML = TRUE, control = ctrl)
  vc <- as.data.table(VarCorr(m))
  tau0 <- vc[grp == "subject" & var1 == "(Intercept)" & is.na(var2), vcov]
  tau1 <- vc[grp == "subject" & var1 == "age_c" & is.na(var2), vcov]
  s_site <- vc[grp == "site", vcov]; s_res <- vc[grp == "Residual", vcov]
  fe <- coef(summary(m))
  re <- ranef(m, condVar = TRUE)$subject; pv <- attr(re, "postVar")
  rel <- 1 - pv[2, 2, ] / tau1
  b <- data.table(subject = rownames(re), metric = mt, int = re[["(Intercept)"]], slope = re[["age_c"]], rel = rel)
  nv <- d[, .(nv = .N, site = site[1]), by = subject]; b <- nv[b, on = "subject"]
  BL[[mt]] <- b
  mu <- fixef(m)[["(Intercept)"]] + 0.5 * fixef(m)[["sexM"]]
  SUM[[mt]] <- data.table(metric = mt, n_children = nrow(re), n_scans = nrow(d), singular = isSingular(m),
    mean_at_12.8 = mu, age_slope = fe["age_c", "Estimate"], age_slope_t = fe["age_c", "t value"],
    pct_per_yr = 100 * fe["age_c", "Estimate"] / mu, sd_units_per_yr = fe["age_c", "Estimate"] / sqrt(tau0),
    child_slope_sd = sqrt(tau1), child_slope_cv = sqrt(tau1) / abs(fe["age_c", "Estimate"]),
    resid_sd = sqrt(s_res), site_icc = s_site / (s_site + tau0 + s_res),
    site_r2_slope = summary(lm(slope ~ factor(site), b))$r.squared,
    rel2 = mean(rel[b$nv == 2]), rel3 = mean(rel[b$nv == 3]), rel4 = mean(rel[b$nv == 4]))
  cat(sprintf("%-20s n=%5d  slope %+.3g (t %+.1f, %+.2f%%/yr, %+.2f SD/yr)  rel %.2f/%.2f/%.2f  siteICC %.2f  siteR2slope %.2f%s\n",
              mt, nrow(re), fe["age_c", "Estimate"], fe["age_c", "t value"], SUM[[mt]]$pct_per_yr, SUM[[mt]]$sd_units_per_yr,
              SUM[[mt]]$rel2, SUM[[mt]]$rel3, SUM[[mt]]$rel4, SUM[[mt]]$site_icc, SUM[[mt]]$site_r2_slope,
              ifelse(isSingular(m), "  SINGULAR", "")))
}
B <- rbindlist(BL)
fwrite(B, file.path(W, "blups_dmri.csv.gz"))

# ---- cross-modal relations of child slopes ---------------------------------------------------------
ws <- dcast(B, subject ~ metric, value.var = "slope")
wi <- dcast(B, subject ~ metric, value.var = "int")
site1 <- unique(B[, .(subject, site)], by = "subject")
ws <- site1[ws, on = "subject"]; wi <- site1[wi, on = "subject"]
ws[, thin_fast := -smri__thk_hcp]
resid_site <- function(v, s) { ok <- !is.na(v); r <- rep(NA_real_, length(v)); r[ok] <- resid(lm(v[ok] ~ factor(s[ok]))); r }
for (mt in metrics) {
  a <- ws[[mt]]
  ok1 <- !is.na(a) & !is.na(ws$thin_fast); ok2 <- !is.na(a) & !is.na(ws$smri__t1t2__gm)
  ra <- resid_site(a, ws$site); rt <- resid_site(ws$thin_fast, ws$site); r12 <- resid_site(ws$smri__t1t2__gm, ws$site)
  oki <- !is.na(wi[[mt]]) & !is.na(wi$smri__thk_hcp)
  SUM[[mt]][, `:=`(r_slope_thinfast = cor(a[ok1], ws$thin_fast[ok1]),
                   r_slope_thinfast_within_site = cor(ra[ok1], rt[ok1], use = "complete.obs"),
                   r_slope_t1t2 = cor(a[ok2], ws$smri__t1t2__gm[ok2]),
                   r_slope_t1t2_within_site = cor(ra[ok2], r12[ok2], use = "complete.obs"),
                   r_int_thickness = cor(wi[[mt]][oki], wi$smri__thk_hcp[oki]))]
}
SM <- rbindlist(SUM)
fwrite(SM, file.path(OUT, "slope_summary.tsv"), sep = "\t")
cc <- ws[complete.cases(ws[, ..metrics])]
C <- cor(as.matrix(cc[, ..metrics]), method = "spearman")
fwrite(data.table(metric = rownames(C), C), file.path(OUT, "slope_cor.tsv"), sep = "\t")
cat(sprintf("slope correlation matrix on %d children with every metric\n", nrow(cc)))

# ---- PRS -> each slope (Figure-1 model; as 01_t1t2_prs.R) ------------------------------------------
ch <- fread("genetic_analysis/work/results_70tab_hcp/d9_sex_moderation/children.csv")
ch[, `:=`(sex = factor(sex, levels = c("F", "M")), FID = factor(FID), site = factor(site), age_c = baseline_age)]
pcs <- paste0("PC", 1:10); cov0 <- paste(c("sex", "age_c", pcs), collapse = " + ")
ctrlp <- lmerControl(calc.derivs = FALSE)
grid <- CJ(disorder = c("scz", "mdd"), method = c("prscs", "sbrc"), arm = c("eur", "pooled"))
P <- list()
for (mt in metrics) {
  y <- B[metric == mt, .(subject, yv = slope)]
  if (mt == "smri__thk_hcp") y[, yv := -yv]             # report thickness as faster thinning
  dd <- y[ch, on = "subject", nomatch = 0]
  for (i in seq_len(nrow(grid))) {
    g <- grid[i]; col <- sprintf("%s_%s_%s", g$disorder, g$method, g$arm)
    d <- dd[!is.na(get(col)) & !is.na(yv) & !is.na(PC1)]
    if (g$arm == "eur") d <- d[eur == 1]
    d[, `:=`(prs = as.numeric(scale(get(col))), y = as.numeric(scale(yv)))]
    for (mn in c("base", "site")) {
      rhs <- if (mn == "base") cov0 else paste(cov0, "+ site")
      m <- lmer(as.formula(paste("y ~ prs +", rhs, "+ (1 | FID)")), data = d, REML = TRUE, control = ctrlp)
      s <- coef(summary(m))
      P[[length(P) + 1]] <- data.table(metric = mt, disorder = toupper(g$disorder), method = g$method, arm = g$arm,
                                       model = mn, est = s["prs", "Estimate"], se = s["prs", "Std. Error"],
                                       p = s["prs", "Pr(>|t|)"], n = nrow(d))
    }
  }
  cat(".")
}
P <- rbindlist(P); P[, mde80 := 2.8 * se]
fwrite(P, file.path(OUT, "dmri_prs.tsv"), sep = "\t")
cat("\ndone\n")
