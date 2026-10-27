#!/usr/bin/env Rscript
# D1 step 5b -- slope reliability of the cross-fitted diffusion composite (05_pls_composite.py), alone and mixed
# with thickness, and PRS on every derived score (composites, out-of-fold PLS prediction, PC1-PC5).
#
# Per-scan series (all on scans present in both the dMRI and the thickness sample, >= 2 per child):
#   thin    = -HCP whole-cortex thickness (slope positive = faster thinning)
#   comp    = per-scan PLS composite with the child's own out-of-fold weights
#   mix_w   = (1 - w) * thin / sd_slope(thin) + w * comp / sd_slope(comp),  w = 0.25, 0.5, 0.75
# Each is fitted with the Figure-1 LMM  y ~ age_c + sex + (1 + age_c | subject) + (1 | site).
# Reported: mean per-child slope reliability (overall and at 2/3/4 visits), R2 of site on child slopes, and the
# correlation of each child slope with the thickness-only slope.
# PRS: Figure-1 model (+ site), PRS-CS and SBayesRC, EUR and pooled arms, on the slopes above, on the child-level
# out-of-fold PLS prediction and on PC1-PC5 (signs arbitrary).
# Outputs: tables/pls_reliability.tsv, tables/pls_prs.tsv. Run (env r, repo root).
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
W <- "directions/d1_multimodal/work"; OUT <- "directions/d1_multimodal/tables"
ctrl <- lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
S <- fread(file.path(W, "pls_scan.csv.gz")); S[, thin := -thk]

fit <- function(v) {
  m <- lmer(as.formula(paste(v, "~ age_c + sex + (1 + age_c | subject) + (1 | site)")), data = S, REML = TRUE, control = ctrl)
  vc <- as.data.table(VarCorr(m)); tau1 <- vc[grp == "subject" & var1 == "age_c" & is.na(var2), vcov]
  re <- ranef(m, condVar = TRUE)$subject; rel <- 1 - attr(re, "postVar")[2, 2, ] / tau1
  list(tau1 = tau1, b = data.table(subject = rownames(re), slope = re[["age_c"]], rel = rel), singular = isSingular(m))
}
F <- list(thin = fit("thin"), comp = fit("comp"))
for (w in c(0.25, 0.5, 0.75)) {
  nm <- sprintf("mix_%02d", round(100 * w))
  S[, (nm) := (1 - w) * thin / sqrt(F$thin$tau1) + w * comp / sqrt(F$comp$tau1)]
  F[[nm]] <- fit(nm)
}
nv <- S[, .(nv = .N, site = site[1]), by = subject]
REL <- rbindlist(lapply(names(F), function(k) {
  b <- nv[F[[k]]$b, on = "subject"]; bt <- F$thin$b
  data.table(series = k, n_children = nrow(b), singular = F[[k]]$singular, rel_mean = mean(b$rel),
             rel2 = mean(b$rel[b$nv == 2]), rel3 = mean(b$rel[b$nv == 3]), rel4 = mean(b$rel[b$nv == 4]),
             site_r2_slope = summary(lm(slope ~ factor(site), b))$r.squared,
             r_with_thin_slope = cor(b$slope, bt$slope[match(b$subject, bt$subject)]))
}))
fwrite(REL, file.path(OUT, "pls_reliability.tsv"), sep = "\t"); print(REL, digits = 3)

# ---- PRS ----
ch <- fread("genetic_analysis/work/results_70tab_hcp/d9_sex_moderation/children.csv")
ch[, `:=`(sex = factor(sex, levels = c("F", "M")), FID = factor(FID), site = factor(site), age_c = baseline_age)]
pc <- fread(file.path(W, "pls_child.csv.gz"))
Y <- Reduce(function(a, b) merge(a, b, by = "subject", all = TRUE),
            c(lapply(names(F), function(k) setnames(copy(F[[k]]$b[, .(subject, slope)]), "slope", paste0("slope_", k))),
              list(pc[, .(subject, pls_oof, dPC1 = PC1, dPC2 = PC2, dPC3 = PC3, dPC4 = PC4, dPC5 = PC5)])))
ch <- Y[ch, on = "subject"]
OUTC <- c(paste0("slope_", names(F)), "pls_oof", paste0("dPC", 1:5))
cov0 <- paste(c("sex", "age_c", paste0("PC", 1:10), "site"), collapse = " + ")
grid <- CJ(disorder = c("scz", "mdd"), method = c("prscs", "sbrc"), arm = c("eur", "pooled"))
P <- list()
for (yv in OUTC) for (i in seq_len(nrow(grid))) {
  g <- grid[i]; col <- sprintf("%s_%s_%s", g$disorder, g$method, g$arm)
  d <- ch[!is.na(get(col)) & !is.na(get(yv)) & !is.na(PC1)]
  if (g$arm == "eur") d <- d[eur == 1]
  d[, `:=`(prs = as.numeric(scale(get(col))), y = as.numeric(scale(get(yv))))]
  # diffusion PCs are dPC1..dPC5; PC1..PC10 are the genetic covariates
  s <- coef(summary(lmer(as.formula(paste("y ~ prs +", cov0, "+ (1 | FID)")), d, control = lmerControl(calc.derivs = FALSE))))
  P[[length(P) + 1]] <- data.table(outcome = yv, disorder = toupper(g$disorder), method = g$method, arm = g$arm,
                                   est = s["prs", 1], se = s["prs", 2], p = s["prs", 5], n = nrow(d))
}
P <- rbindlist(P); fwrite(P, file.path(OUT, "pls_prs.tsv"), sep = "\t")
P[, cell := sprintf("%+.3f%s", est, ifelse(p < 0.05, "*", " "))]
print(dcast(P[method == "prscs"], outcome ~ disorder + arm, value.var = "cell"))
print(dcast(P[method == "sbrc"], outcome ~ disorder + arm, value.var = "cell"))
