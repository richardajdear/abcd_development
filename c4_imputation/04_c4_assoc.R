#!/usr/bin/env Rscript
# Step 4 (CSD3): imputed C4A expression -> adolescent cortical thinning rate.
#
#   Rscript 04_c4_assoc.R --c4 work/c4_calls.tsv --pheno-dir <dir> --eur-ids <keep> \
#     [--prs-dir <SCZ score dir> | --prs-file <FID IID score table>] \
#     [--slope-col global_slope_1lmm --baseline-col baseline_thickness_1lmm] \
#     [--post-min 0.7] [--label hcp_1lmm] --out results/c4_assoc_<label>.tsv
#
# The model is tools/prs_assoc.R's, with the C4 predictor in place of the score,
# so a C4 beta and a PRS beta are on the same scale (SD of trait per SD of
# predictor) and directly comparable:
#
#   scale(trait) ~ scale(C4A_GREx) + sex + age_c + PC1..PC10 + (1 | family_id)
#
# PRE-SPECIFIED (README.md, "Analysis plan"):
#   M1_primary            EUR arm, every child, C4A GREx from expected (posterior-
#                         weighted) allele dosages; trait = the primary thinning
#                         rate (single-LMM slope).  One test, alpha 0.05, two-sided.
#                         Hypothesis: beta < 0.  No posterior filter: on the
#                         leave-out test the Beagle posteriors were poorly
#                         calibrated and filtering at 0.7 barely raised accuracy,
#                         while expected dosages need no filter.
# Secondary, reported whatever M1 shows:
#   M2_joint_C4B          + scale(C4B_GREx)  (C4B is not a SCZ risk gene; its
#                         coefficient is a specificity check)
#   M3_sex_interaction    C4A x sex (Kamitaki 2020: stronger C4 effects in men),
#                         plus M3_female / M3_male fits
#   M4_conditional_PRS    + scale(SCZ PRS): is C4A independent of genome-wide risk?
#                         The PRS term's beta is reported too (with vs without C4A)
#   M5_control_baseline   same predictor, baseline thickness (the intercept)
#   M6_common5            only children whose two structures are among the five
#                         common ones (Hernandez et al. 2023 filter)
#   M7_copy_numbers       C4A copies + C4B copies + HERV copies, weight-free
#   M8_posterior_filter   M1 restricted to structure posterior >= post-min (0.7;
#                         the Hernandez et al. 2023 filter)
# mde80 = 2.80 * se: the |beta| this fit had 80 % power to detect at alpha 0.05.
suppressPackageStartupMessages({ library(data.table); library(lme4); library(lmerTest) })

args <- commandArgs(trailingOnly = TRUE)
opt <- list(`post-min` = "0.7", label = "c4", `slope-col` = NA, `baseline-col` = NA)
i <- 1
while (i <= length(args)) {
  stopifnot(startsWith(args[i], "--"), i < length(args))
  opt[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2
}
stopifnot(!is.null(opt$c4), !is.null(opt$`pheno-dir`), !is.null(opt$out), !is.null(opt$`eur-ids`))
post_min <- as.numeric(opt$`post-min`)
token <- function(x) sub("^.*?([A-Z0-9]{8})$", "\\1", x)

# ---- phenotype + covariates (aligned export; FID = family id) ----------------
pd   <- opt$`pheno-dir`
ph   <- fread(file.path(pd, "phenotypes_gcta.txt"))
qcov <- fread(file.path(pd, "covar_quant.txt"))
ccov <- fread(file.path(pd, "covar_categorical.txt"))
pick <- function(given, choices) {
  if (!is.na(given)) return(given)
  hit <- intersect(choices, names(ph)); if (!length(hit)) stop("none of ", paste(choices, collapse = "/"), " in phenotypes"); hit[1]
}
slope_col <- pick(opt$`slope-col`, c("global_slope_1lmm", "global_slope"))
base_col  <- pick(opt$`baseline-col`, c("baseline_thickness_1lmm", "baseline_thickness"))
d <- merge(ph[, c("FID", "IID", slope_col, base_col), with = FALSE], qcov, by = c("FID", "IID"))
d <- merge(d, ccov[, .(FID, IID, sex)], by = c("FID", "IID"))
setnames(d, c(slope_col, base_col), c("slope", "baseline"))
d[, family_id := as.character(FID)]
d[, sex := factor(sex)]
d[, age_c := baseline_age - mean(baseline_age, na.rm = TRUE)]
d[, tok := token(IID)]
pc_cols <- grep("^PC[0-9]+$", names(qcov), value = TRUE)
stopifnot(length(pc_cols) == 10)

# ---- C4 calls ------------------------------------------------------------------
c4 <- fread(opt$c4)
c4[, tok := token(IID)]
stopifnot(!anyDuplicated(c4$tok))
n_ph <- nrow(d)
d <- merge(d, c4[, .(tok, C4A_GREx, C4B_GREx, C4A_copies, C4B_copies, HERV_copies, post_mean, common5)], by = "tok")
cat(sprintf("joined %d of %d phenotyped children to C4 calls (%d calls)\n", nrow(d), n_ph, nrow(c4)))
if (nrow(d) < 0.8 * n_ph) stop("fewer than 80% of phenotyped children have C4 calls -- ID mismatch (rule 1)?")

eur <- token(fread(opt$`eur-ids`, header = FALSE)[[ncol(fread(opt$`eur-ids`, header = FALSE, nrows = 1))]])
d <- d[tok %in% eur]
cat(sprintf("EUR arm: %d children\n", nrow(d)))

# ---- optional SCZ score ------------------------------------------------------
has_prs <- FALSE
if (!is.null(opt$`prs-dir`) || !is.null(opt$`prs-file`)) {
  if (!is.null(opt$`prs-dir`)) {
    prof <- list.files(opt$`prs-dir`, pattern = "^score_.*\\.profile$", full.names = TRUE)
    if (length(prof) != 1) stop("expected one score_*.profile in ", opt$`prs-dir`, ", found ", length(prof))
    s <- fread(prof); col <- if ("SCORESUM" %in% names(s)) "SCORESUM" else "SCORE"
  } else {
    s <- fread(opt$`prs-file`); col <- setdiff(names(s), c("FID", "IID"))[1]
  }
  s <- s[, .(tok = token(IID), PRS = get(col))]
  d <- merge(d, s, by = "tok", all.x = TRUE)
  has_prs <- TRUE
}

# ---- models --------------------------------------------------------------------
covs <- c("sex", "age_c", pc_cols)
fit <- function(data, y, preds, model, extra_covs = covs) {
  data <- copy(data)
  data[, y__ := as.numeric(scale(get(y)))]
  for (p in preds) data[, (paste0(p, "_z")) := as.numeric(scale(get(p)))]
  rhs <- c(paste0(preds, "_z"), extra_covs)
  f <- as.formula(paste("y__ ~", paste(rhs, collapse = " + "), "+ (1 | family_id)"))
  m <- lmer(f, data = data, REML = TRUE, control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  s <- summary(m)$coefficients
  keep <- grep("C4A|C4B|HERV|PRS", rownames(s), value = TRUE)
  data.table(label = opt$label, model = model, trait = y, term = keep,
             beta = s[keep, "Estimate"], se = s[keep, "Std. Error"], p = s[keep, "Pr(>|t|)"],
             n = nobs(m), n_families = length(unique(data$family_id)))
}
fit_int <- function(data, model) {   # C4A x sex
  data <- copy(data)
  data[, y__ := as.numeric(scale(slope))]; data[, C4A_GREx_z := as.numeric(scale(C4A_GREx))]
  f <- as.formula(paste("y__ ~ C4A_GREx_z * sex +", paste(setdiff(covs, "sex"), collapse = " + "), "+ (1 | family_id)"))
  m <- lmer(f, data = data, REML = TRUE, control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE))
  s <- summary(m)$coefficients; keep <- grep("C4A", rownames(s), value = TRUE)
  data.table(label = opt$label, model = model, trait = "slope", term = keep,
             beta = s[keep, "Estimate"], se = s[keep, "Std. Error"], p = s[keep, "Pr(>|t|)"],
             n = nobs(m), n_families = length(unique(data$family_id)))
}

q <- d
res <- list(
  fit(q, "slope", "C4A_GREx", "M1_primary"),
  fit(q, "slope", c("C4A_GREx", "C4B_GREx"), "M2_joint_C4B"),
  fit_int(q, "M3_sex_interaction"))
for (lv in levels(q$sex))
  res[[length(res) + 1]] <- fit(q[sex == lv], "slope", "C4A_GREx", paste0("M3_sex_", lv), setdiff(covs, "sex"))
if (has_prs) {
  qp <- q[is.finite(PRS)]
  res[[length(res) + 1]] <- fit(qp, "slope", "PRS", "M4_PRS_alone")
  res[[length(res) + 1]] <- fit(qp, "slope", c("C4A_GREx", "PRS"), "M4_conditional_PRS")
}
res[[length(res) + 1]] <- fit(q, "baseline", "C4A_GREx", "M5_control_baseline")
res[[length(res) + 1]] <- fit(q[common5 == 1], "slope", "C4A_GREx", "M6_common5")
res[[length(res) + 1]] <- fit(q, "slope", c("C4A_copies", "C4B_copies", "HERV_copies"), "M7_copy_numbers")
res[[length(res) + 1]] <- fit(d[post_mean >= post_min], "slope", "C4A_GREx", "M8_posterior_filter")
out <- rbindlist(res)
out[, trait := fifelse(trait == "slope", slope_col, fifelse(trait == "baseline", base_col, trait))]
out[, mde80 := 2.80 * se]
out[, post_min := post_min]
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(out, opt$out, sep = "\t")
p1 <- out[model == "M1_primary" & term == "C4A_GREx_z"]
cat(sprintf("M1 primary: beta %.4f (se %.4f), p %.3g, n %d, mde80 %.3f -> %s\n",
            p1$beta, p1$se, p1$p, p1$n, p1$mde80, opt$out))
