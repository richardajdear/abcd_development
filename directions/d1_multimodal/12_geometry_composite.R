#!/usr/bin/env Rscript
# D1 step 8 -- (a) composites of thickness, GM volume and surface area (HCP-MMP, Figure-1 scans);
# (b) how much would less scan-to-scan noise (e.g. the FreeSurfer longitudinal stream) raise thinning-slope reliability?
#
# (a) Each measure is z-scored over scans; a per-scan composite is the sum of the z-scores with weights w. The composite
#     is fitted with the Figure-1 LMM, so its reliability is estimated on the composite itself (thickness, area and
#     volume share surfaces, hence errors; a child-level average of slopes would hide that). Slopes are reported as
#     faster LOSS (sign flipped), so a positive PRS beta = faster loss. PRS: Figure-1 model + site, 8 cells.
#     The weight grid is fixed in advance; picking the best cell afterwards is selection, so read the whole grid.
# (b) For the thickness LMM, the per-child posterior slope variance is recomputed with the residual variance scaled by
#     f = 1, 0.75, 0.5, 0.33, 0.25, holding the between-child (co)variances fixed and using each child's own scan ages:
#     V_i = (G^-1 + Z_i'Z_i / (f sigma2))^-1, reliability_i = 1 - V_i[2,2] / tau1. f = 1 is checked against lme4.
#     The expected PRS z scales with sqrt(reliability) (corr(PRS, BLUP) = corr(PRS, true slope) x sqrt(reliability)).
# Outputs: tables/geom_composite.tsv, tables/geom_composite_prs.tsv, tables/noise_scenarios.tsv. Run (env r, repo root).
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
W <- "directions/d1_multimodal/work"; OUT <- "directions/d1_multimodal/tables"
ctrl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
S <- fread(file.path(W, "scan_means_smri.csv.gz"))[metric %in% c("smri__thk_hcp", "smri__vol_hcp", "smri__area_hcp")]
S[, z := (mean_val - mean(mean_val)) / sd(mean_val), by = metric]
X <- dcast(S, subject + visit + age_c + sex + site ~ metric, value.var = "z")
setnames(X, c("smri__thk_hcp", "smri__vol_hcp", "smri__area_hcp"), c("thk", "vol", "area"))
cat(sprintf("scans with all three measures: %d (children %d); scan-level r: thk-vol %.2f, thk-area %.2f, vol-area %.2f\n",
            nrow(X), uniqueN(X$subject), cor(X$thk, X$vol), cor(X$thk, X$area), cor(X$vol, X$area)))
DEF <- list(thk = c(1, 0, 0), vol = c(0, 1, 0), area = c(0, 0, 1),
            thk_vol_0.5 = c(1, 0.5, 0), thk_vol_1 = c(1, 1, 0), thk_vol_2 = c(1, 2, 0),
            thk_area_0.5 = c(1, 0, 0.5), thk_area_1 = c(1, 0, 1), thk_vol_area = c(1, 1, 1))
fitc <- function(nm, w) {
  X[, y := w[1] * thk + w[2] * vol + w[3] * area]; X[, y := (y - mean(y)) / sd(y)]
  m <- lmer(y ~ age_c + sex + (1 + age_c | subject) + (1 | site), data = X, REML = TRUE, control = ctrl)
  vc <- as.data.table(VarCorr(m)); tau1 <- vc[grp == "subject" & var1 == "age_c" & is.na(var2), vcov]
  re <- ranef(m, condVar = TRUE)$subject; rel <- 1 - attr(re, "postVar")[2, 2, ] / tau1
  list(m = m, b = data.table(subject = rownames(re), loss = -re[["age_c"]], rel = rel), tau1 = tau1, vc = vc)
}
nv <- X[, .(nv = .N), by = subject]
FITS <- lapply(names(DEF), function(k) fitc(k, DEF[[k]])); names(FITS) <- names(DEF)
REL <- rbindlist(lapply(names(FITS), function(k) {
  b <- nv[FITS[[k]]$b, on = "subject"]
  data.table(series = k, w_thk = DEF[[k]][1], w_vol = DEF[[k]][2], w_area = DEF[[k]][3], rel_mean = mean(b$rel),
             rel2 = mean(b$rel[b$nv == 2]), rel3 = mean(b$rel[b$nv == 3]), rel4 = mean(b$rel[b$nv == 4]),
             r_with_thk_slope = cor(b$loss, FITS$thk$b$loss[match(b$subject, FITS$thk$b$subject)]))
}))
fwrite(REL, file.path(OUT, "geom_composite.tsv"), sep = "\t"); print(REL, digits = 3)

ch <- fread("genetic_analysis/work/results_70tab_hcp/d9_sex_moderation/children.csv")
ch[, `:=`(sex = factor(sex, levels = c("F", "M")), FID = factor(FID), site = factor(site), age_c = baseline_age)]
cov0 <- paste(c("sex", "age_c", paste0("PC", 1:10), "site"), collapse = " + ")
grid <- CJ(disorder = c("scz", "mdd"), method = c("prscs", "sbrc"), arm = c("eur", "pooled"))
P <- list()
for (k in names(FITS)) {
  dd <- FITS[[k]]$b[ch, on = "subject", nomatch = 0]
  for (i in seq_len(nrow(grid))) {
    g <- grid[i]; col <- sprintf("%s_%s_%s", g$disorder, g$method, g$arm)
    d <- dd[!is.na(get(col)) & !is.na(loss) & !is.na(PC1)]
    if (g$arm == "eur") d <- d[eur == 1]
    d[, `:=`(prs = as.numeric(scale(get(col))), yz = as.numeric(scale(loss)))]
    s <- coef(summary(lmer(as.formula(paste("yz ~ prs +", cov0, "+ (1 | FID)")), d, control = lmerControl(calc.derivs = FALSE))))
    P[[length(P) + 1]] <- data.table(series = k, disorder = toupper(g$disorder), method = g$method, arm = g$arm,
                                     est = s["prs", 1], se = s["prs", 2], z = s["prs", 1] / s["prs", 2], p = s["prs", 5], n = nrow(d))
  }
}
P <- rbindlist(P); fwrite(P, file.path(OUT, "geom_composite_prs.tsv"), sep = "\t")
P[, cell := sprintf("%+.3f (%+.1f)", est, z)]
print(dcast(P, series ~ disorder + method + arm, value.var = "cell"))
print(P[, .(mean_z_SCZ = mean(z[disorder == "SCZ"]), mean_z_MDD = mean(z[disorder == "MDD"]), mean_z_all = mean(z)), by = series])

# ---- (b) noise scenarios for the thickness slope ------------------------------------------------------
f0 <- FITS$thk; vc <- f0$vc
G <- matrix(c(vc[grp == "subject" & var1 == "(Intercept)" & is.na(var2), vcov],
              vc[grp == "subject" & var1 == "(Intercept)" & var2 == "age_c", vcov],
              vc[grp == "subject" & var1 == "(Intercept)" & var2 == "age_c", vcov], f0$tau1), 2)
s2 <- vc[grp == "Residual", vcov]; Gi <- solve(G)
Zs <- split(X[, .(subject, age_c)], by = "subject", keep.by = FALSE)
relf <- function(f) sapply(Zs, function(z) { Z <- cbind(1, z$age_c); V <- solve(Gi + crossprod(Z) / (f * s2)); 1 - V[2, 2] / f0$tau1 })
r1 <- relf(1)
cat(sprintf("\ncheck f = 1: analytic mean reliability %.3f vs lme4 %.3f (r over children %.3f)\n", mean(r1), mean(f0$b$rel),
            cor(r1, f0$b$rel[match(names(Zs), f0$b$subject)])))
NS <- rbindlist(lapply(c(1, 0.75, 0.5, 1/3, 0.25), function(f) {
  r <- relf(f); nvv <- sapply(Zs, nrow)
  data.table(resid_var_factor = f, resid_sd_factor = sqrt(f), rel_mean = mean(r), rel2 = mean(r[nvv == 2]),
             rel3 = mean(r[nvv == 3]), rel4 = mean(r[nvv == 4]), prs_z_multiplier = sqrt(mean(r) / mean(r1)))
}))
fwrite(NS, file.path(OUT, "noise_scenarios.tsv"), sep = "\t"); print(NS, digits = 3)
