#!/usr/bin/env Rscript
# Ancestry-stratified polygenic-score association, and a test of whether the
# pooled arm is entitled to pool.
#
# WHY THIS EXISTS (step 9, 2026-09-17).  The pooled arm reports one beta for
# 8,596 children of four ancestry clusters, with the score standardised within
# cluster so that between-cluster differences in score cannot masquerade as an
# effect.  That still ASSUMES the within-cluster effect is the same in every
# cluster.  With ancestry-specific discovery GWAS now available (EUR, AFR, EAS
# and their meta), the assumption can be checked directly and the transfer of
# each weight set into each cluster can be read off:
#
#   1. WITHIN-STRATUM fits: for each k-means cluster (EURlike, cluster1 =
#      African-American-like, cluster2 = Hispanic-like, cluster3 = mixed/Asian)
#      and for the EUR anchor set,
#          scale(y) ~ scale(PRS) + sex + age_c + PC1..PC10 + (1 | family_id)
#      the same model as tools/prs_assoc.R, fitted inside the stratum.
#
#   2. POOLED WITH STRATUM FIXED EFFECTS: PRS z-scored within stratum, stratum
#      as a factor, so the beta is the within-stratum effect averaged over
#      strata -- the pooled arm's estimand, stated explicitly.
#
#   3. HETEROGENEITY: likelihood-ratio test of PRS x stratum against the model
#      in 2 (ML fits).  A significant p_het says the clusters do not share one
#      slope and the pooled beta is an average of different things; a null
#      p_het at these n is NOT proof of homogeneity (cluster3 has ~380
#      children), so read the per-stratum SEs too.
#
# Strata come from strata_k4.tsv (legacy/hpc/README_HPC.md validates them
# blind to self-report as 91 % White / 80 % Black / 95 % Hispanic / 45 % Asian).
# The EUR anchor set (eur_anchor.keep) is the stricter EUR definition the EUR
# arm uses; it is a subset of EURlike and is reported as its own row.
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(lme4); library(lmerTest)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--prs-dir",     type = "character"),
  make_option("--pheno",       type = "character"),
  make_option("--covar-quant", type = "character"),
  make_option("--covar-cat",   type = "character"),
  make_option("--manifest",    type = "character"),
  make_option("--strata",      type = "character", help = "strata_k4.tsv (IID, stratum)"),
  make_option("--eur-ids",     type = "character", default = NULL),
  make_option("--out",         type = "character"),
  make_option("--min-n",       type = "integer", default = 150),
  make_option("--phenotypes",  type = "character", default = NULL,
              help = "comma-separated subset of manifest phenotypes [default all]")
)))
for (r in c("prs-dir", "pheno", "covar-quant", "covar-cat", "manifest", "strata", "out"))
  if (is.null(opt[[r]])) stop("--", r, " is required")

ph   <- fread(opt$pheno); qcov <- fread(opt$`covar-quant`)
ccov <- fread(opt$`covar-cat`); man <- fread(opt$manifest)
pheno_names <- man$name[man$name %in% names(ph)]
if (!is.null(opt$phenotypes)) pheno_names <- intersect(pheno_names, strsplit(opt$phenotypes, ",")[[1]])
if (!length(pheno_names)) stop("no phenotypes to fit")
pc_cols <- grep("^PC[0-9]+$", names(qcov), value = TRUE)

d <- merge(ph, qcov, by = c("FID", "IID"))
d <- merge(d, ccov[, .(FID, IID, sex, site)], by = c("FID", "IID"))
d[, family_id := as.character(FID)]
if (!"age_c" %in% names(d)) {
  if (!"baseline_age" %in% names(d)) stop("no baseline_age / age_c covariate")
  d[, age_c := baseline_age - mean(baseline_age, na.rm = TRUE)]
}
norm_id <- function(x) toupper(gsub("_", "", sub("^sub-", "", x)))
st <- fread(opt$strata)
st[, token := norm_id(IID)]
d[, token := norm_id(IID)]
n0 <- nrow(d)
d <- merge(d, unique(st[, .(token, stratum)]), by = "token")
message(sprintf("strata: %d of %d subjects assigned; %s", nrow(d), n0,
                paste(sprintf("%s=%d", names(table(d$stratum)), table(d$stratum)), collapse = ", ")))
d[, stratum := factor(stratum)]
has_eur <- !is.null(opt$`eur-ids`)
if (has_eur) {
  eur <- fread(opt$`eur-ids`, header = FALSE)$V1
  d[, eur_anchor := IID %in% eur]
  message(sprintf("EUR anchor: %d subjects", sum(d$eur_anchor)))
}
ctrl <- lmerControl(calc.derivs = FALSE)
covs <- paste(c("sex", "age_c", pc_cols), collapse = " + ")

prof <- list.files(opt$`prs-dir`, pattern = "^score_.*\\.profile$", full.names = TRUE)
if (!length(prof)) stop("no score_*.profile files in ", opt$`prs-dir`)
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.table(...)

for (f in prof) {
  tag <- sub("^score_", "", sub("\\.profile$", "", basename(f)))
  parts <- strsplit(tag, "_", fixed = TRUE)[[1]]
  disorder <- parts[1]; threshold <- paste(parts[-1], collapse = "_")
  s <- fread(f)
  col <- if ("SCORESUM" %in% names(s)) "SCORESUM" else "SCORE"
  s <- s[, .(IID, PRS = get(col))]
  dd <- merge(d, s, by = "IID")
  if (!nrow(dd)) { warning("no ID overlap for ", tag); next }
  dd <- dd[is.finite(PRS)]
  dd[, PRSz := as.numeric(scale(PRS)), by = stratum]

  for (p in pheno_names) {
    groups <- split(dd, by = "stratum", keep.by = TRUE)
    if (has_eur) groups[["EURanchor"]] <- dd[eur_anchor == TRUE]
    for (g in names(groups)) {
      sub <- groups[[g]][is.finite(get(p))]
      if (nrow(sub) < opt$`min-n` || uniqueN(sub$family_id) < 50 || sd(sub$PRS) == 0) next
      fml <- as.formula(sprintf("scale(%s) ~ scale(PRS) + %s + (1 | family_id)", p, covs))
      fit <- try(lmer(fml, data = sub, REML = TRUE, control = ctrl), silent = TRUE)
      if (inherits(fit, "try-error")) next
      co <- summary(fit)$coefficients
      if (!"scale(PRS)" %in% rownames(co)) next
      add(disorder = disorder, threshold = threshold, phenotype = p,
          analysis = "within_stratum", stratum = g,
          n = nrow(sub), n_families = uniqueN(sub$family_id),
          beta = co["scale(PRS)", "Estimate"], se = co["scale(PRS)", "Std. Error"],
          p = co["scale(PRS)", "Pr(>|t|)"], chisq = NA_real_, df = NA_integer_)
    }
    # pooled, stratum fixed effects, PRS standardised within stratum; then the
    # PRS x stratum interaction by LRT (both ML).
    sub <- dd[is.finite(get(p))]
    if (nrow(sub) < opt$`min-n` || uniqueN(sub$stratum) < 2) next
    f0 <- as.formula(sprintf("scale(%s) ~ PRSz + stratum + %s + (1 | family_id)", p, covs))
    f1 <- as.formula(sprintf("scale(%s) ~ PRSz * stratum + %s + (1 | family_id)", p, covs))
    m0 <- try(lmer(f0, data = sub, REML = FALSE, control = ctrl), silent = TRUE)
    m1 <- try(lmer(f1, data = sub, REML = FALSE, control = ctrl), silent = TRUE)
    if (inherits(m0, "try-error") || inherits(m1, "try-error")) next
    co0 <- summary(m0)$coefficients
    add(disorder = disorder, threshold = threshold, phenotype = p,
        analysis = "pooled_stratumFE", stratum = "all",
        n = nrow(sub), n_families = uniqueN(sub$family_id),
        beta = co0["PRSz", "Estimate"], se = co0["PRSz", "Std. Error"],
        p = co0["PRSz", "Pr(>|t|)"], chisq = NA_real_, df = NA_integer_)
    a <- anova(m0, m1)
    add(disorder = disorder, threshold = threshold, phenotype = p,
        analysis = "heterogeneity_PRSxStratum", stratum = "all",
        n = nrow(sub), n_families = uniqueN(sub$family_id),
        beta = NA_real_, se = NA_real_, p = a[2, "Pr(>Chisq)"],
        chisq = a[2, "Chisq"], df = as.integer(a[2, "Df"]))
  }
  message(sprintf("  %s: %d rows so far", tag, length(rows)))
}
res <- rbindlist(rows)
if (!nrow(res)) stop("no models fitted")
# Bonferroni over thresholds within disorder x phenotype x analysis x stratum
# (only C+T has more than one); mirrors tools/prs_assoc.R.
res[, n_thresholds := uniqueN(threshold), by = .(disorder, phenotype, analysis, stratum)]
res[, p_adj := pmin(1, p * n_thresholds)]
setorder(res, phenotype, analysis, stratum, disorder, threshold)
fwrite(res, opt$out, sep = "\t")
cat(sprintf("wrote %s (%d rows)\n", opt$out, nrow(res)))
