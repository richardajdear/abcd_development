#!/usr/bin/env Rscript
# D1 step 1 -- T1w/T2w check: do the SCZ / MDD polygenic scores predict the T1w/T2w slope,
# and is the ABCD T1w/T2w measure good enough to carry the myelination question?
#
# Phenotypes (one single LMM per whole-cortex measure, as Figure 1):
#   thin_fast   = -global_slope_1lmm   (HCP-MMP thickness, Figure 1 trait; HIGHER = FASTER thinning)
#   t1t2_slope  = re_slope of directions/d3_twin_family/02_fit_slopes.R on the DK whole-cortex mean
#                 T1w/T2w ratio (out/t1t2_ratio_dsk_70_9a62dde44370); HIGHER = FASTER rise
#   t1t2_int    = re_intercept of the same model (level at the centring age)
# HCP-MMP T1w/T2w does not exist yet (needs the CSD3 job in CSD3_INSTRUCTIONS.md), so DK is used.
#
# Prediction under the myelination account of thinning (DIRECTIONS.md 2.2): a score that
# speeds thinning (positive beta on thin_fast) also speeds the T1w/T2w rise (positive beta on
# t1t2_slope). Under the neuropil account, t1t2_slope beta ~ 0.
#
# Model template (tools/prs_assoc.R, as directions/d9_sex_moderation/02_fit_prs_dct.R):
#   base      scale(y) ~ scale(PRS) + sex + age_c + PC1..10 + (1 | family)    (REML, lmerTest)
#   site      base + site (fixed)            -- T1w/T2w is not bias-field corrected, so site matters
#   inv       base on the rank-inverse-normal y   -- the T1w/T2w slope has |z| up to ~28
#   site_inv  site + inv
# Arms: EUR (EUR weights, eur_anchor.keep, z within arm) and pooled (_zanc multi-ancestry
# weights, z over the analysis sample), PRS-CS and SBayesRC.
# Gate: base model on thin_fast must reproduce -beta_1lmm of fig1_inputs/hcp70_prs_key_arms.tsv.
#
# Inputs : gitignored per-child tables (d9 children.csv; d3 blups_t1t2.csv, scan_means_*.csv)
# Outputs: tables/t1t2_prs.tsv, tables/t1t2_gate.tsv, tables/t1t2_measure_qc.tsv (aggregates only)
# Run (env r, repo root):  Rscript directions/d1_multimodal/01_t1t2_prs.R
suppressMessages({library(data.table); library(lme4); library(lmerTest)})
WORK <- "genetic_analysis/work/results_70tab_hcp/d9_sex_moderation"
D3W  <- "directions/d3_twin_family/work"
D3R  <- "directions/d3_twin_family/results"
OUT  <- "directions/d1_multimodal/tables"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
pcs  <- paste0("PC", 1:10)
ctrl <- lmerControl(calc.derivs = FALSE)

ch <- fread(file.path(WORK, "children.csv"))
bt <- fread(file.path(D3W, "blups_t1t2.csv"))[, .(subject, t1t2_slope = re_slope, t1t2_int = re_intercept,
                                                   t1t2_rel = slope_rel)]
nv <- fread(file.path(D3W, "scan_means_t1t2.csv"))[, .(n_visits_t1t2 = .N), by = subject]
ch <- nv[bt[ch, on = "subject"], on = "subject"]
ch[, `:=`(sex = factor(sex, levels = c("F", "M")), FID = factor(FID), site = factor(site),
          age_c = baseline_age, thin_fast = -global_slope_1lmm)]
inv_norm <- function(x) { r <- rank(x, na.last = "keep"); qnorm((r - 0.5) / sum(!is.na(x))) }
cat(sprintf("children %d; with thickness slope %d; with T1w/T2w slope %d\n",
            nrow(ch), sum(!is.na(ch$thin_fast)), sum(!is.na(ch$t1t2_slope))))

# ---- 1. measure quality -----------------------------------------------------------------------
vc <- fread(file.path(D3R, "lmm_variance_components.tsv"))
rel <- fread(file.path(D3R, "slope_reliability.tsv"))
QC <- list()
for (tag in c("ct", "t1t2")) {
  v <- vc[measure == tag]
  s_int <- v[grp == "subject" & var1 == "(Intercept)" & (is.na(var2) | var2 == ""), vcov]
  s_sl  <- v[grp == "subject" & var1 == "age_c" & (is.na(var2) | var2 == ""), vcov]
  s_site <- v[grp == "site", vcov]; s_res <- v[grp == "Residual", vcov]
  QC[[length(QC) + 1]] <- data.table(measure = tag, quantity = c(
    "scan-level site ICC (site / (site + child intercept + residual))",
    "residual SD / between-child slope SD (years of follow-up to reach SNR 1)"),
    value = c(s_site / (s_site + s_int + s_res), sqrt(s_res) / sqrt(s_sl)))
  r <- rel[measure == tag]
  QC[[length(QC) + 1]] <- data.table(measure = tag, quantity = sprintf("mean per-child slope reliability, %d visits", r$n_visits),
                                     value = r$mean_rel)
}
for (y in c("thin_fast", "t1t2_slope")) {
  d <- ch[!is.na(get(y))]; z <- as.numeric(scale(d[[y]]))
  f <- lm(d[[y]] ~ d$site); a <- anova(f)
  QC[[length(QC) + 1]] <- data.table(measure = ifelse(y == "thin_fast", "ct", "t1t2"), quantity = c(
    "slope: R2 of site (child BLUP ~ site)", "slope: site F-test p", "slope: n children |z| > 4", "slope: n children |z| > 8",
    "slope: n children"), value = c(summary(f)$r.squared, a$`Pr(>F)`[1], sum(abs(z) > 4), sum(abs(z) > 8), nrow(d)))
}
d <- ch[!is.na(thin_fast) & !is.na(t1t2_slope)]
rw <- cor(resid(lm(thin_fast ~ site, d)), resid(lm(t1t2_slope ~ site, d)))
QC[[length(QC) + 1]] <- data.table(measure = "ct x t1t2", quantity = c(
  "r(faster thinning, faster T1w/T2w rise), Pearson", "same, Spearman", "same, within site (site-residualised)", "n children"),
  value = c(cor(d$thin_fast, d$t1t2_slope), cor(d$thin_fast, d$t1t2_slope, method = "spearman"), rw, nrow(d)))
QC <- rbindlist(QC)
fwrite(QC, file.path(OUT, "t1t2_measure_qc.tsv"), sep = "\t")
print(QC, digits = 3)

# ---- 2. PRS -> slope ---------------------------------------------------------------------------
key <- fread("genetic_analysis/fig1_inputs/hcp70_prs_key_arms.tsv")[phenotype == "global_slope"]
KEYARM  <- c(scz_eur = "SCZ25_EUR", scz_pooled = "SCZ25_META", mdd_eur = "MDD_eur", mdd_pooled = "MDD_pooled")
KEYMETH <- c(prscs = "PRSCS", sbrc = "SBayesRC")
cov0 <- paste(c("sex", "age_c", pcs), collapse = " + ")
fitm <- function(rhs, d) lmer(as.formula(paste("y ~ prs +", rhs, "+ (1 | FID)")), data = d, REML = TRUE, control = ctrl)
co <- function(m) { s <- coef(summary(m)); data.table(est = s["prs", "Estimate"], se = s["prs", "Std. Error"],
                                                       p = s["prs", "Pr(>|t|)"]) }
MODELS <- list(base = list(rhs = cov0, inv = FALSE), site = list(rhs = paste(cov0, "+ site"), inv = FALSE),
               inv = list(rhs = cov0, inv = TRUE), site_inv = list(rhs = paste(cov0, "+ site"), inv = TRUE))
ROWS <- list(); GATE <- list()
grid <- CJ(disorder = c("scz", "mdd"), method = c("prscs", "sbrc"), arm = c("eur", "pooled"))
for (i in seq_len(nrow(grid))) {
  g <- grid[i]; col <- sprintf("%s_%s_%s", g$disorder, g$method, g$arm)
  k <- data.table(disorder = toupper(g$disorder), method = g$method, arm = g$arm)
  for (yv in c("thin_fast", "t1t2_slope", "t1t2_int")) {
    d <- ch[!is.na(get(col)) & !is.na(get(yv)) & !is.na(PC1)]
    if (g$arm == "eur") d <- d[eur == 1]
    d[, prs := as.numeric(scale(get(col)))]
    for (mn in names(MODELS)) {
      M <- MODELS[[mn]]
      d[, y := as.numeric(scale(if (M$inv) inv_norm(get(yv)) else get(yv)))]
      r <- cbind(k, outcome = yv, model = mn, co(fitm(M$rhs, d)), n = nrow(d))
      ROWS[[length(ROWS) + 1]] <- r
      if (yv == "thin_fast" && mn == "base") {
        kr <- key[trait_arm == KEYARM[[paste(g$disorder, g$arm, sep = "_")]] & method == KEYMETH[[g$method]]]
        GATE[[length(GATE) + 1]] <- cbind(k, est = r$est, se = r$se, n = r$n,
                                          committed_faster = -kr$beta_1lmm, committed_se = kr$se_1lmm, committed_n = kr$n)
      }
    }
  }
  cat(".", if (i %% 4 == 0) "\n" else "")
}
R <- rbindlist(ROWS); R[, mde80 := 2.80 * se]
G <- rbindlist(GATE); G[, diff := est - committed_faster]
fwrite(R, file.path(OUT, "t1t2_prs.tsv"), sep = "\t")
fwrite(G, file.path(OUT, "t1t2_gate.tsv"), sep = "\t")
cat(sprintf("\ngate: max |beta - committed| = %.2g; n mismatches = %d\n", max(abs(G$diff)), sum(G$n != G$committed_n)))
w <- dcast(R[method == "prscs"], disorder + arm + outcome ~ model, value.var = c("est", "p"))
print(w, digits = 2)
