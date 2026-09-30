#!/usr/bin/env Rscript
# 02_models.R -- D5: does the environment act on the same pace variable as polygenic
# risk, and which way does the thinning <-> symptom association run?
#
#   Rscript directions/d5_adversity_direction/code/02_models.R        # repo root, env r
#
# Input : out/d5_adversity_direction/{child,waves}.parquet  (01_build_tables.py;
#         individual-level, gitignored)
# Output: directions/d5_adversity_direction/results/table_d5_*.tsv (coefficients only)
#
# Conventions (all models)
#   * Outcome and every continuous predictor standardised within the analysis
#     sample, so beta is SD per SD.  global_slope_1lmm is the standardised random
#     slope: NEGATIVE = FASTER thinning.
#   * Genetics-arm covariates, as in tools/prs_assoc.R: sex + age_c (baseline age)
#     + PC1..PC10 + (1 | family_id), lmer / REML.  "site" is NOT in the PRS model
#     because the slope already carries a (1 | site) term from the trait LMM; the
#     environment models are repeated with site fixed effects (variant "+site")
#     because SES differs sharply between the 21 sites.
#   * Arms: pooled (all 8,596 genotyped children, within-ancestry z scores) and
#     EUR (4,308, raw scores), each PRS in its matched arm only.
#
# Sections
#   A0  reproduction: PRS -> slope must match genetic_analysis/current_results.tsv
#   A1  environment -> slope (and -> baseline thickness, for contrast)
#   A2  gene-environment correlation: environment ~ PRS
#   A3  PRS -> slope with SES and adversity in the model (attenuation; family-
#       cluster bootstrap of the change, 1,000 reps, OLS with the same fixed part)
#   A4  PRS x environment interaction on the slope
#   B1  reverse direction: slope ~ symptoms at baseline (ses-00A, before any change)
#   B2  forward direction: symptoms at 15-17 ~ slope + symptoms at baseline
#       (Figure 1g model: OLS, site FE, family-clustered SE), then + SES + adversity
#   B3  random-intercept cross-lagged panel model (Hamaker et al. 2015) on the four
#       imaging waves: within-child deviations of whole-cortex thickness and CBCL,
#       2-year lags, lags constrained equal over time, MLR with family clustering.
#       Each variable is first residualised within wave on age, sex and site.
suppressMessages({
  library(arrow); library(data.table); library(lme4); library(lmerTest)
  library(sandwich); library(lmtest); library(lavaan)
})
set.seed(5)
options(warn = 1)
REPO <- normalizePath(".")
stopifnot(file.exists(file.path(REPO, "abcd-data-release-7.0")))
IN  <- file.path(REPO, "out/d5_adversity_direction")
RES <- file.path(REPO, "directions/d5_adversity_direction/results")
C <- as.data.table(read_parquet(file.path(IN, "child.parquet")))
W <- as.data.table(read_parquet(file.path(IN, "waves.parquet")))
C[, age_c := baseline_age - mean(baseline_age)]
C[, sex := factor(sex)]; C[, site := factor(site)]
PCS <- paste0("PC", 1:10)
C[, (PCS) := lapply(.SD, function(x) as.numeric(scale(x))), .SDcols = PCS]   # scale only; betas unchanged
COV <- c("sex", "age_c", PCS)
ENV <- c("income", "parent_edu", "adi", "conflict_y", "conflict_p", "bad_events_y", "bad_events_p",
         "ses", "adversity")
ARMS <- list(SCZ_pooled = "pooled", MDD_pooled = "pooled", SCZ_eur = "EUR", MDD_eur = "EUR")
# EA scores are added by 01_build_tables.py only when copied from CSD3 (README, "EA step")
if ("EA_eur" %in% names(C) && sum(!is.na(C$EA_eur)) > 1000) ARMS$EA_eur <- "EUR"
cat("scores in this run:", paste(names(ARMS), collapse = ", "), "\n")
OUTC <- c("pfactor", "internal", "external", "depress")
zs <- function(x) as.numeric(scale(x))
arm_rows <- function(arm) if (arm == "EUR") C[eur == TRUE] else C

fit_lmer <- function(d, y, x, extra = character()) {
  vars <- unique(c(y, x, extra, COV, "family_id"))
  d <- na.omit(d[, ..vars])
  # standardise within the analysis sample (covariates too: in the EUR arm the
  # ancestry PCs have a tiny spread, which trips lme4's scaling check)
  for (v in c(y, x, extra, "age_c", PCS)) if (is.numeric(d[[v]])) d[[v]] <- zs(d[[v]])
  f <- as.formula(sprintf("%s ~ %s + (1 | family_id)", y, paste(c(x, extra, COV), collapse = " + ")))
  m <- lmer(f, data = d, REML = TRUE, control = lmerControl(calc.derivs = FALSE))
  s <- coef(summary(m))
  list(m = m, n = nrow(d), s = s, d = d)
}
row_of <- function(fit, term, ...) {
  data.table(..., term = term, beta = fit$s[term, "Estimate"], se = fit$s[term, "Std. Error"],
             p = fit$s[term, "Pr(>|t|)"], n = fit$n)
}

# ---------------------------------------------------------------- A0 -------------
# cluster betas, genetic_analysis/current_results.tsv (hcp, 1lmm, global_slope, PRSCS)
REF <- c(SCZ_pooled = -0.029276, SCZ_eur = -0.036081, MDD_pooled = -0.024396, MDD_eur = -0.02939,
         EA_eur = 0.039505)
A0 <- rbindlist(lapply(names(ARMS), function(k) {
  f <- fit_lmer(arm_rows(ARMS[[k]]), "global_slope_1lmm", k)
  row_of(f, k, score = k, arm = ARMS[[k]])[, cluster_beta := REF[[k]]]
}))
A0[, abs_diff := abs(beta - cluster_beta)]
print(A0)
stopifnot(max(A0$abs_diff) < 0.003)
fwrite(A0, file.path(RES, "table_d5_A0_reproduction.tsv"), sep = "\t")

# ---------------------------------------------------------------- A1 -------------
A1 <- rbindlist(lapply(c("global_slope_1lmm", "baseline_thickness_1lmm"), function(y)
  rbindlist(lapply(ENV, function(e) rbindlist(list(
    row_of(fit_lmer(C, y, e), e, outcome = y, variant = "base"),
    row_of(fit_lmer(C, y, e, "site"), e, outcome = y, variant = "+site"),
    row_of(fit_lmer(C, y, e, c("site", "qc_defects")), e, outcome = y, variant = "+site+qc")))))))
jt <- fit_lmer(C, "global_slope_1lmm", c("ses", "adversity"))
A1 <- rbind(A1, row_of(jt, "ses", outcome = "global_slope_1lmm", variant = "joint ses+adversity"),
            row_of(jt, "adversity", outcome = "global_slope_1lmm", variant = "joint ses+adversity"))
fwrite(A1, file.path(RES, "table_d5_A1_env_slope.tsv"), sep = "\t")
print(A1[outcome == "global_slope_1lmm" & variant != "+site"])

# ---------------------------------------------------------------- A2 -------------
A2 <- rbindlist(lapply(names(ARMS), function(k) rbindlist(lapply(c("ses", "adversity", ENV[1:7]), function(e)
  row_of(fit_lmer(arm_rows(ARMS[[k]]), e, k), k, score = k, arm = ARMS[[k]], env = e)))))
fwrite(A2, file.path(RES, "table_d5_A2_rGE.tsv"), sep = "\t")
print(A2[env %in% c("ses", "adversity")])

# ---------------------------------------------------------------- A3 -------------
ols_beta <- function(d, x, extra) {
  f <- as.formula(sprintf("global_slope_1lmm ~ %s", paste(c(x, extra, COV), collapse = " + ")))
  coef(lm(f, data = d))[[x]]
}
A3 <- rbindlist(lapply(names(ARMS), function(k) {
  vars <- unique(c("global_slope_1lmm", k, "ses", "adversity", "site", COV, "family_id"))
  d <- na.omit(arm_rows(ARMS[[k]])[, ..vars])
  for (v in c("global_slope_1lmm", k, "ses", "adversity")) d[[v]] <- zs(d[[v]])
  sets <- list(`PRS alone` = character(), `+ses` = "ses", `+adversity` = "adversity",
               `+ses+adversity` = c("ses", "adversity"), `+ses+adversity+site` = c("ses", "adversity", "site"))
  out <- rbindlist(lapply(names(sets), function(s) {
    f <- fit_lmer(d, "global_slope_1lmm", k, sets[[s]])
    row_of(f, k, score = k, arm = ARMS[[k]], model = s)
  }))
  # family-cluster bootstrap of beta(alone) - beta(+ses+adversity), OLS fixed part
  fam <- split(seq_len(nrow(d)), d$family_id)
  diffs <- replicate(1000, {
    b <- d[unlist(fam[sample(length(fam), replace = TRUE)], use.names = FALSE)]
    ols_beta(b, k, character()) - ols_beta(b, k, c("ses", "adversity"))
  })
  b0 <- out[model == "PRS alone", beta]; b1 <- out[model == "+ses+adversity", beta]
  out[, `:=`(attenuation_pct = 100 * (b0 - beta) / b0)]
  out[model == "+ses+adversity", `:=`(diff = b0 - b1, diff_lo = quantile(diffs, .025), diff_hi = quantile(diffs, .975))]
  out
}), fill = TRUE)
fwrite(A3, file.path(RES, "table_d5_A3_attenuation.tsv"), sep = "\t")
print(A3[, .(score, model, beta, p, n, attenuation_pct, diff, diff_lo, diff_hi)])

# ---------------------------------------------------------------- A4 -------------
A4 <- rbindlist(lapply(names(ARMS), function(k) rbindlist(lapply(c("ses", "adversity"), function(e) {
  vars <- unique(c("global_slope_1lmm", k, e, COV, "family_id"))
  d <- copy(na.omit(arm_rows(ARMS[[k]])[, ..vars]))
  for (v in c("global_slope_1lmm", k, e)) set(d, j = v, value = zs(d[[v]]))
  d[, gxe := get(k) * get(e)]
  f <- fit_lmer(d, "global_slope_1lmm", c(k, e, "gxe"))
  rbind(row_of(f, k, score = k, arm = ARMS[[k]], env = e, role = "PRS"),
        row_of(f, e, score = k, arm = ARMS[[k]], env = e, role = "environment"),
        row_of(f, "gxe", score = k, arm = ARMS[[k]], env = e, role = "PRS x environment"))
}))))
fwrite(A4, file.path(RES, "table_d5_A4_gxe.tsv"), sep = "\t")
print(A4[role == "PRS x environment"])

# ---------------------------------------------------------------- B1, B2 ---------
B1 <- rbindlist(lapply(OUTC, function(o) rbindlist(list(
  row_of(fit_lmer(C, "global_slope_1lmm", paste0(o, "_base")), paste0(o, "_base"), outcome = o, model = "base"),
  row_of(fit_lmer(C, "global_slope_1lmm", paste0(o, "_base"), c("ses", "adversity")), paste0(o, "_base"),
         outcome = o, model = "+ses+adversity")))))
B1[, direction := "baseline symptoms -> subsequent slope"]

fig1g <- function(o, extra = character()) {
  y <- paste0(o, "_late"); yb <- paste0(o, "_base")
  vars <- unique(c(y, yb, "global_slope_1lmm", "age_late", "sex", "site", "family_id", extra))
  d <- na.omit(C[, ..vars])
  sdy <- sd(d[[y]])
  f <- as.formula(sprintf("%s ~ global_slope_1lmm + %s + age_late + sex + site%s", y, yb,
                          paste0(if (length(extra)) " + " else "", paste(extra, collapse = " + "))))
  m <- lm(f, data = d)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = ~family_id))
  data.table(outcome = o, term = "global_slope_1lmm", beta = ct["global_slope_1lmm", 1] / sdy,
             se = ct["global_slope_1lmm", 2] / sdy, p = ct["global_slope_1lmm", 4], n = nrow(d))
}
B2 <- rbindlist(lapply(OUTC, function(o) rbind(fig1g(o)[, model := "base"],
                                               fig1g(o, c("ses", "adversity"))[, model := "+ses+adversity"])))
B2[, direction := "slope -> symptoms at 15-17 | baseline"]
B12 <- rbind(B1, B2, fill = TRUE)
fwrite(B12, file.path(RES, "table_d5_B12_direction.tsv"), sep = "\t")
print(B12)

# ---------------------------------------------------------------- B3 -------------
W <- merge(W, C[, .(guid8, sex, site, family_id)], by = "guid8")
WAV <- c("v0", "v2", "v4", "v6")
resid_wave <- function(v, extra = NULL) W[, {
  r <- rep(NA_real_, .N); dd <- data.frame(y = get(v), age, sex, site)
  if (!is.null(extra)) dd$q <- get(extra)
  ok <- complete.cases(dd)
  f <- if (is.null(extra)) y ~ age + sex + site else y ~ age + sex + site + q
  r[ok] <- residuals(lm(f, data = dd[ok, ])); r
}, by = visit]$V1
setorder(W, visit)
for (v in c("mean_ct", "depress", "internal", "external")) W[, (paste0(v, "_r")) := resid_wave(v)]
W[, mean_ct_rq := resid_wave("mean_ct", "qc_defects")]      # thickness net of scan quality
riclpm <- function(o, xvar = "mean_ct_r") {
  wide <- dcast(W, guid8 + family_id ~ visit, value.var = c(xvar, paste0(o, "_r")))
  setnames(wide, gsub(paste0("^", o, "_r_"), "y_", gsub(paste0("^", xvar, "_"), "x_", names(wide))))
  sdx <- sd(unlist(wide[, paste0("x_", WAV), with = FALSE]), na.rm = TRUE)
  sdy <- sd(unlist(wide[, paste0("y_", WAV), with = FALSE]), na.rm = TRUE)
  for (w in WAV) { wide[[paste0("x_", w)]] <- wide[[paste0("x_", w)]] / sdx
                   wide[[paste0("y_", w)]] <- wide[[paste0("y_", w)]] / sdy }
  X <- paste0("x_", WAV); Y <- paste0("y_", WAV); wx <- paste0("wx", 1:4); wy <- paste0("wy", 1:4)
  mod <- c(sprintf("RIx =~ %s", paste0("1*", X, collapse = " + ")),
           sprintf("RIy =~ %s", paste0("1*", Y, collapse = " + ")),
           sprintf("%s =~ 1*%s", wx, X), sprintf("%s =~ 1*%s", wy, Y),
           sprintf("%s ~~ 0*%s", c(X, Y), c(X, Y)),
           sprintf("%s ~ ax*%s + cyx*%s", wx[-1], wx[-4], wy[-4]),     # CBCL_t -> CT_t+1
           sprintf("%s ~ ay*%s + cxy*%s", wy[-1], wy[-4], wx[-4]),     # CT_t -> CBCL_t+1
           sprintf("%s ~~ %s", wx, wy), sprintf("%s ~~ %s", wx, wx), sprintf("%s ~~ %s", wy, wy),
           "RIx ~~ RIx", "RIy ~~ RIy", "RIx ~~ RIy",
           sprintf("RIx ~~ 0*%s", c(wx[1], wy[1])), sprintf("RIy ~~ 0*%s", c(wx[1], wy[1])))
  fit <- lavaan(paste(mod, collapse = "\n"), data = as.data.frame(wide), missing = "fiml",
                estimator = "MLR", cluster = "family_id", meanstructure = TRUE, int.ov.free = TRUE)
  pe <- as.data.table(parameterEstimates(fit))
  fm <- fitMeasures(fit, c("cfi.robust", "rmsea.robust", "srmr"))
  lab <- c(cxy = "CT_t -> symptoms_t+2y", cyx = "symptoms_t -> CT_t+2y", ax = "CT autoregression",
           ay = "symptoms autoregression")
  r <- unique(pe[label %in% names(lab), .(label, est, se, pvalue)], by = "label")
  ri <- pe[lhs == "RIx" & op == "~~" & rhs == "RIy", .(label = "RI covariance", est, se, pvalue)]
  # standardised (std.all): lags are constrained equal on the raw scale, so the
  # standardised value varies slightly by wave; report the mean over the 3 lags
  ss <- as.data.table(standardizedSolution(fit))
  if ("label" %in% names(ss)) ss[, label := NULL]
  ss <- merge(ss, pe[, .(lhs, op, rhs, label)], by = c("lhs", "op", "rhs"))
  std <- ss[label %in% names(lab), .(std_all = mean(est.std)), by = label]
  std <- rbind(std, data.table(label = "RI covariance",
                               std_all = ss[lhs == "RIx" & op == "~~" & rhs == "RIy", est.std]))
  merge(rbind(r, ri), std, by = "label")[, `:=`(
    path = c(lab, `RI covariance` = "between-child (RI) correlation")[label], outcome = o,
    thickness = if (xvar == "mean_ct_r") "age+sex+site" else "age+sex+site+scan quality",
    n = lavInspect(fit, "ntotal"), cfi = fm[["cfi.robust"]], rmsea = fm[["rmsea.robust"]], srmr = fm[["srmr"]])]
}
B3 <- rbindlist(lapply(c("depress", "internal", "external"), function(o)
  rbind(riclpm(o), riclpm(o, "mean_ct_rq"))))
fwrite(B3, file.path(RES, "table_d5_B3_riclpm.tsv"), sep = "\t")
print(B3)
