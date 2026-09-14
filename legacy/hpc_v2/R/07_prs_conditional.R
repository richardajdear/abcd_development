#!/usr/bin/env Rscript
# Is the disorder association independent of polygenic EDUCATIONAL ATTAINMENT?
#
# Motivation.  EA is not inert here: under SBayesR it associates with
# global_slope at +0.030 (p_adj = .050), POSITIVE -- higher polygenic education,
# slower cortical thinning.  Because the sign is opposite to the disorder
# scores, an EA-mediated confound cannot manufacture a negative disorder
# coefficient directly.  But "cannot manufacture it" is weaker than "the
# disorder coefficient survives adjustment for it", and only the joint model
# answers the latter.
#
# MODEL, per disorder x method x stratum:
#
#   marginal:    phenotype ~ scale(PRS_disorder) + covars + (1 | family_id)
#   conditional: phenotype ~ scale(PRS_disorder) + scale(PRS_EA) + covars + (1 | family_id)
#
# and the quantity of interest is the CHANGE in the disorder coefficient.
#
# HOW TO READ IT, because the naive reading is wrong in a specific way.  The two
# scores are correlated, so conditioning will shrink the disorder coefficient
# somewhat whatever the truth is.  What distinguishes confounding from
# collinearity:
#   * attenuation toward zero with the EA term significant, on a score pair
#     whose correlation is substantial  -> EA is plausibly carrying the signal
#   * coefficient essentially unchanged                       -> independent
#   * BOTH coefficients inflating with inflated SEs           -> collinearity
#     artefact, not a result; the r between the scores is reported so this is
#     checkable rather than guessed at
# The score-score correlation within the stratum is printed for every row, since
# without it the attenuation is uninterpretable.
suppressPackageStartupMessages({
  library(data.table); library(lme4); library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--final",       type = "character"),
  make_option("--pheno",       type = "character"),
  make_option("--covar-quant", type = "character"),
  make_option("--covar-cat",   type = "character"),
  make_option("--manifest",    type = "character"),
  make_option("--eur-ids",     type = "character"),
  make_option("--out",         type = "character")
)))

ph   <- fread(opt$pheno); qcov <- fread(opt$`covar-quant`)
ccov <- fread(opt$`covar-cat`); man <- fread(opt$manifest)
pheno_names <- man$name[man$name %in% names(ph)]
pc_cols <- grep("^PC[0-9]+$", names(qcov), value = TRUE)
d <- merge(ph, qcov, by = c("FID", "IID"))
d <- merge(d, ccov[, c("FID", "IID", "sex", "site"), with = FALSE], by = c("FID", "IID"))
d[, family_id := FID]
# Same age fix as prs_assoc.R: the covariate export names the column
# baseline_age, so asking for age_c silently gave an age-unadjusted model.
if (!"age_c" %in% names(d)) {
  if ("baseline_age" %in% names(d)) {
    d[, age_c := baseline_age - mean(baseline_age, na.rm = TRUE)]
    message("age_c derived by centring baseline_age")
  } else stop("no age column (age_c / baseline_age)")
}
eur <- fread(opt$`eur-ids`, header = FALSE)$V1
d[, is_eur := IID %in% eur]

# Pick, for each (method, trait), the single score file the main table reports:
# the best threshold for C+T, the single score for everything else.
pick <- function(meth, arm) {
  dir <- file.path(opt$final, meth, arm)
  f <- list.files(dir, pattern = "^score_.*\\.profile$", full.names = TRUE)
  if (!length(f)) return(NULL)
  if (meth != "CT") return(f[1])
  # C+T: the threshold the main table selected, read back from it so the
  # conditional test is run on the same score and not a differently chosen one
  tm <- file.path(opt$final, "table_main.tsv")
  if (!file.exists(tm)) return(NULL)
  t <- fread(tm)
  t <- t[method == meth & trait_arm == arm & phenotype == "global_slope" &
         matched == "yes" & score == "raw"]
  if (!nrow(t)) return(NULL)
  hit <- grep(paste0("_", t$threshold[1], "\\.profile$"), f, value = TRUE)
  if (length(hit)) hit[1] else NULL
}
readscore <- function(f, nm) {
  s <- fread(f)
  col <- if ("SCORESUM" %in% names(s)) "SCORESUM" else "SCORE"
  s <- s[, c("IID", col), with = FALSE]; setnames(s, col, nm); s
}

DISORDERS <- list(c("SCZ_eur", "EUR"), c("SCZ_pooled", "full"),
                  c("MDD_eur", "EUR"), c("MDD_pooled", "full"),
                  c("ALZ_noAPOE", "EUR"), c("ALZ", "EUR"),
                  c("ALZ_IGAP", "EUR"), c("ALZ_IGAP_noAPOE", "EUR"), c("ASD", "EUR"))
METHODS <- c("CT", "PRSCS", "SBayesR", "SBayesRC")

rows <- list()
for (meth in METHODS) {
  fe <- pick(meth, "EA"); if (is.null(fe)) next
  ea <- readscore(fe, "PRS_EA")
  for (spec in DISORDERS) {
    arm <- spec[1]; stratum <- spec[2]
    fd <- pick(meth, arm); if (is.null(fd)) next
    dd <- merge(merge(d, readscore(fd, "PRS_D"), by = "IID"), ea, by = "IID")
    sub <- if (stratum == "EUR") dd[is_eur == TRUE] else dd
    for (p in pheno_names) {
      s <- sub[is.finite(get(p)) & is.finite(PRS_D) & is.finite(PRS_EA)]
      if (nrow(s) < 100 || uniqueN(s$family_id) < 50) next
      cv <- paste(c("sex", "age_c", pc_cols), collapse = " + ")
      f1 <- as.formula(sprintf("scale(%s) ~ scale(PRS_D) + %s + (1|family_id)", p, cv))
      f2 <- as.formula(sprintf("scale(%s) ~ scale(PRS_D) + scale(PRS_EA) + %s + (1|family_id)", p, cv))
      m1 <- try(lmer(f1, s, REML = TRUE, control = lmerControl(calc.derivs = FALSE)), silent = TRUE)
      m2 <- try(lmer(f2, s, REML = TRUE, control = lmerControl(calc.derivs = FALSE)), silent = TRUE)
      if (inherits(m1, "try-error") || inherits(m2, "try-error")) next
      c1 <- summary(m1)$coefficients; c2 <- summary(m2)$coefficients
      if (!("scale(PRS_D)" %in% rownames(c1)) || !("scale(PRS_D)" %in% rownames(c2))) next
      b1 <- c1["scale(PRS_D)", "Estimate"]; b2 <- c2["scale(PRS_D)", "Estimate"]
      rows[[length(rows) + 1]] <- data.table(
        method = meth, trait_arm = arm, stratum = stratum, phenotype = p,
        n = nrow(s), score_file = basename(fd),
        r_scores = cor(s$PRS_D, s$PRS_EA),
        beta_marginal = b1, se_marginal = c1["scale(PRS_D)", "Std. Error"],
        p_marginal = c1["scale(PRS_D)", "Pr(>|t|)"],
        beta_cond = b2, se_cond = c2["scale(PRS_D)", "Std. Error"],
        p_cond = c2["scale(PRS_D)", "Pr(>|t|)"],
        beta_EA = c2["scale(PRS_EA)", "Estimate"],
        p_EA = c2["scale(PRS_EA)", "Pr(>|t|)"],
        pct_attenuation = 100 * (1 - b2 / b1))
    }
  }
}
res <- rbindlist(rows)
if (!nrow(res)) stop("no models fitted")
fwrite(res, opt$out, sep = "\t")
cat(sprintf("wrote %s (%d rows)\n", opt$out, nrow(res)))

cat("\nglobal_slope: disorder coefficient before and after adjusting for EA\n")
cat("r = score-score correlation in this stratum; att% = attenuation\n\n")
g <- res[phenotype == "global_slope"][order(trait_arm, method)]
cat(sprintf("%-17s %-9s %-5s %6s %18s %18s %8s %8s\n", "trait_arm", "method",
            "strat", "r", "marginal", "conditional", "att%", "p_EA"))
for (i in seq_len(nrow(g))) with(g[i], cat(sprintf(
  "%-17s %-9s %-5s %+6.3f  %+7.4f (%.3f)  %+7.4f (%.3f) %7.1f %8.3f\n",
  trait_arm, method, stratum, r_scores, beta_marginal, p_marginal,
  beta_cond, p_cond, pct_attenuation, p_EA)))
