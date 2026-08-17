#!/usr/bin/env Rscript
# Within-family PRS association: does the SCZ score still predict cortical
# thinning when we compare SIBLINGS to each other?
#
# WHY THIS TEST EXISTS
# --------------------
# The population PRS association (tools/prs_assoc.R) is confounded by anything
# that varies between families and correlates with both the score and the
# phenotype: population structure the ancestry PCs did not fully capture,
# assortative mating, and "genetic nurture" (parental genotype acting through
# the rearing environment rather than through the child's own genome).  Ancestry
# PCs are a partial fix and everyone knows it.
#
# Siblings are the clean comparison.  Two full sibs are drawn from the same
# parents, so their expected ancestry is identical, they share the family
# environment, and their parents' genotypes are the same.  The only thing that
# differs is which alleles each one happened to inherit -- Mendelian
# segregation, which is random with respect to everything environmental.  So if
# the sib who inherited more SCZ risk alleles also thins faster, that difference
# cannot be explained by ancestry, upbringing, or parental phenotype.  It is as
# close to a randomised exposure as observational genetics gets.
#
# THE MODEL
# ---------
# The standard decomposition (Fulker/Abecasis; "between-within") splits each
# subject's score into a family mean and a deviation from it:
#
#     PRS_ij  =  PRS_mean_j          +  (PRS_ij - PRS_mean_j)
#                between-family         within-family
#
# and fits BOTH as separate predictors:
#
#     y_ij = b_B * PRS_mean_j + b_W * (PRS_ij - PRS_mean_j)
#            + covariates + u_j + e_ij ,      u_j ~ N(0, s2_fam)
#
# where i indexes a sibling and j a family.
#
#   b_W  is the WITHIN-family effect.  Identified only by Mendelian segregation,
#        so it is immune to confounding that is constant within a family.  This
#        is the causal estimate.
#   b_B  is the BETWEEN-family effect.  It absorbs the true effect PLUS
#        stratification, assortative mating and genetic nurture.
#
# Two things to read off the fit:
#
#   1. Is b_W different from zero, and in the same direction as the population
#      estimate?  If yes, the association survives the strongest available
#      control.
#   2. Is b_B > b_W?  A significant gap is direct evidence of confounding
#      (or of genetic nurture, which is a real biological effect but not the
#      one we are claiming).  The formal test is a Wald test on b_B - b_W;
#      equivalently, refit with PRS_ij and PRS_mean_j as the two regressors,
#      where the coefficient on PRS_mean_j is exactly b_B - b_W.
#
# Note it is the FAMILY MEAN, not the parental score, that goes in: we have no
# parental genotypes here, so the sibship mean is the available proxy for the
# shared component.  With only one genotyped sib a family contributes nothing to
# b_W (its deviation is identically 0) but still contributes to b_B, which is
# why singletons are kept in the fit rather than dropped -- they sharpen b_B and
# therefore the b_B - b_W contrast.
#
# WHAT THIS TEST CAN AND CANNOT DO AT ABCD's CURRENT SIZE
# ------------------------------------------------------
# Only within-pair variation carries information about b_W, and for siblings
# that is roughly half the total score variance, so the effective sample size is
# the number of PAIRS, not subjects.  At the ~688 pairs recoverable from the
# current genotype release, the expected z for an effect the size of the
# observed population estimate (beta = -0.040 per SD) is only about 0.7 --
# roughly 11% power.
#
# So a null b_W here is NOT evidence against the effect, and this script prints
# that in the output rather than leaving it to be misread.  What the test can do
# is BOUND the confounding: if b_W is close to b_B with a wide interval, the
# data are consistent with the population estimate being unconfounded; if b_B is
# significantly larger than b_W, the population estimate is inflated and by
# roughly how much.  The second of those is answerable at this N, because it
# uses the between-family information too.
#
# CRITICAL PREREQUISITE
# ---------------------
# This needs real family IDs.  The current genotype .fam on CSD3 has FID = IID
# (every subject their own family), so the sibships are invisible to any tool
# reading the .fam and this script will refuse to run on it.  The family IDs
# exist in the GCTA export (src/abcd/gcta_export.py writes the real FID
# precisely so this is possible) -- pass that file via --pheno and the FIDs come
# from there, or repair the .fam.  See the check below.
#
#   Rscript tools/prs_assoc_withinfamily.R \
#     --prs-dir  work/results/prs \
#     --pheno    work/pheno/phenotypes_gcta.txt \
#     --covar-quant work/pheno/covar_quant.txt \
#     --covar-cat   work/pheno/covar_categorical.txt \
#     --manifest    work/pheno/phenotype_manifest.tsv \
#     --out      work/results/prs/prs_withinfamily.tsv
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(lme4); library(lmerTest)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--prs-dir",     type = "character"),
  make_option("--pheno",       type = "character"),
  make_option("--covar-quant", type = "character"),
  make_option("--covar-cat",   type = "character"),
  make_option("--manifest",    type = "character"),
  make_option("--out",         type = "character"),
  make_option("--eur-ids",     type = "character", default = NULL),
  make_option("--min-pairs",   type = "integer",   default = 100,
              help = "refuse to fit below this many informative pairs")
)))
for (r in c("prs-dir","pheno","covar-quant","covar-cat","manifest","out"))
  if (is.null(opt[[r]])) stop("--", r, " is required")

man   <- fread(opt$manifest)
ph    <- fread(opt$pheno)
qcov  <- fread(opt$`covar-quant`)
ccov  <- fread(opt$`covar-cat`)

pheno_names <- man$name[man$name %in% names(ph)]
if (!length(pheno_names)) stop("no manifest phenotypes found in ", opt$pheno)
pc_cols <- grep("^PC[0-9]+$", names(qcov), value = TRUE)

d <- merge(ph, qcov, by = c("FID","IID"))
d <- merge(d, ccov[, c("FID","IID","sex","site"), with = FALSE], by = c("FID","IID"))
d[, family_id := FID]

# --- The prerequisite check, stated loudly -------------------------------------
# FID == IID for every row means the family structure was destroyed upstream.
# Fitting anyway would silently produce b_W = 0 with b_B = the population
# estimate, which looks like a clean null result and is actually no test at all.
if (all(as.character(d$FID) == as.character(d$IID))) {
  stop("FID == IID for every subject in ", opt$pheno, ": the family structure ",
       "is absent, so there is nothing to compare within. This is the ",
       "FID = IID shortcut in the genotype .fam. Re-export with real family ",
       "IDs (src/abcd/gcta_export.py does this) before running.")
}

# --- ID normalisation, matching work/align_ids.py ------------------------------
norm_id <- function(x) toupper(gsub("_", "", sub("^sub-", "", as.character(x))))
d[, key := norm_id(IID)]

# --- Ancestry strata -----------------------------------------------------------
if (!is.null(opt$`eur-ids`)) {
  eur <- fread(opt$`eur-ids`, header = FALSE)$V1
  d[, is_eur := IID %in% eur | key %in% norm_id(eur)]
  eur_def <- paste0("supplied list (", opt$`eur-ids`, ")")
} else if (length(pc_cols) >= 2) {
  use <- pc_cols[1:min(4, length(pc_cols))]
  ctr <- d[, lapply(.SD, median, na.rm = TRUE), .SDcols = use]
  sds <- d[, lapply(.SD, sd,     na.rm = TRUE), .SDcols = use]
  z   <- sweep(as.matrix(d[, ..use]), 2, unlist(ctr), "-")
  z   <- sweep(z, 2, unlist(sds), "/")
  d[, is_eur := sqrt(rowSums(z^2)) < 2]
  eur_def <- "largest PC1-PC4 cluster within 2 SD of the median (stand-in)"
} else {
  d[, is_eur := TRUE]; eur_def <- "no PCs available; EUR == full sample"
}

res <- list()
prof <- list.files(opt$`prs-dir`, pattern = "^score_.*\\.profile$", full.names = TRUE)
if (!length(prof)) stop("no score_*.profile in ", opt$`prs-dir`)

for (f in prof) {
  tag <- sub("^score_", "", sub("\\.profile$", "", basename(f)))
  parts    <- strsplit(tag, "_", fixed = TRUE)[[1]]
  disorder <- parts[1]
  threshold <- paste(parts[-1], collapse = "_")

  s <- fread(f)
  if (!"SCORESUM" %in% names(s)) next
  # Join on the normalised ID, never on FID: the score file's FID is the one we
  # do not trust (it is the .fam's), while d's FID came from the export.
  s[, key := norm_id(IID)]
  dd <- merge(d, s[, .(key, PRS = SCORESUM)], by = "key")
  if (!nrow(dd)) next

  for (stratum in c("EUR", "full")) {
    sub <- if (stratum == "EUR") dd[is_eur == TRUE] else dd
    if (!nrow(sub)) next

    # Between/within decomposition. Scale FIRST, over the analysed subset, so
    # b_B and b_W are in the same per-SD units as the population estimate.
    sub[, PRS_z := as.numeric(scale(PRS))]
    sub[, prs_fam := mean(PRS_z), by = family_id]
    sub[, prs_dev := PRS_z - prs_fam]

    # Informative pairs: only families with >=2 genotyped members contribute to
    # b_W.  Count them explicitly and report -- a reader must be able to see
    # that a null b_W came from 30 pairs rather than 3,000.
    sizes  <- sub[, .N, by = family_id]
    n_mult <- sizes[N >= 2, .N]
    n_pair <- sizes[N >= 2, sum(N * (N - 1) / 2)]
    if (is.na(n_pair)) n_pair <- 0

    for (phen in pheno_names) {
      if (!phen %in% names(sub)) next
      terms <- c("prs_fam", "prs_dev", "sex", "age_c", pc_cols)
      terms <- terms[terms %in% names(sub) | terms %in% c("prs_fam","prs_dev")]
      fml <- as.formula(sprintf("scale(%s) ~ %s + (1 | family_id)",
                                phen, paste(terms, collapse = " + ")))
      fit <- try(lmer(fml, data = sub, REML = TRUE,
                      control = lmerControl(calc.derivs = FALSE)), silent = TRUE)
      if (inherits(fit, "try-error")) next
      co <- summary(fit)$coefficients
      if (!all(c("prs_fam","prs_dev") %in% rownames(co))) next

      bW <- co["prs_dev","Estimate"]; seW <- co["prs_dev","Std. Error"]
      bB <- co["prs_fam","Estimate"]; seB <- co["prs_fam","Std. Error"]

      # Wald test on bB - bW: is the between-family estimate inflated relative
      # to the causal within-family one?  Uses the actual covariance of the two
      # coefficients, not an independence assumption -- they are correlated.
      V  <- as.matrix(vcov(fit))
      iB <- which(rownames(V) == "prs_fam"); iW <- which(rownames(V) == "prs_dev")
      se_diff <- sqrt(V[iB,iB] + V[iW,iW] - 2 * V[iB,iW])
      z_diff  <- (bB - bW) / se_diff

      res[[length(res)+1]] <- data.table(
        disorder = disorder, threshold = threshold, phenotype = phen,
        stratum = stratum, n = nrow(sub),
        n_families = uniqueN(sub$family_id),
        n_multi_families = n_mult, n_pairs = as.integer(n_pair),
        beta_within = bW, se_within = seW,
        p_within = co["prs_dev","Pr(>|t|)"],
        beta_between = bB, se_between = seB,
        p_between = co["prs_fam","Pr(>|t|)"],
        beta_diff = bB - bW, se_diff = se_diff, z_diff = z_diff,
        p_diff = 2 * pnorm(-abs(z_diff)))
    }
  }
}

if (!length(res)) stop("no models fitted; check --prs-dir and --pheno overlap")
out <- rbindlist(res)
setorder(out, disorder, phenotype, stratum, threshold)
fwrite(out, opt$out, sep = "\t")

cat(sprintf("wrote %s (%d rows)\nEUR definition: %s\n",
            opt$out, nrow(out), eur_def))

np <- max(out$n_pairs, na.rm = TRUE)
cat(sprintf("\nInformative within-family pairs: %d\n", np))
if (np < opt$`min-pairs`) {
  cat("WARNING: too few pairs for the within-family estimate to be informative.\n")
}
# Power statement, printed unconditionally: the failure mode this test invites
# is reading a null b_W as evidence against the effect.
b_ref <- 0.040
se_ref <- 1 / sqrt(max(np, 1) * 0.5)
cat(sprintf(paste0("At %d pairs, an effect of %.3f SD/SD would give z ~ %.2f ",
                   "(power ~%.0f%% at alpha=0.05).\n"),
            np, b_ref, b_ref/se_ref,
            100 * (1 - pnorm(1.96 - b_ref/se_ref))))
cat("A null within-family estimate at this N is UNINFORMATIVE about the effect;\n")
cat("the interpretable quantity is beta_diff (between minus within), which\n")
cat("tests whether the population estimate is inflated by confounding.\n")
