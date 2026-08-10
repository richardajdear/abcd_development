#!/usr/bin/env Rscript
# Associate polygenic scores with the ABCD developmental phenotypes.
#
# Called by hpc/06_prs.sbatch; also runnable standalone:
#
#   Rscript tools/prs_assoc.R --prs-dir <dir> --pheno phenotypes_gcta.txt \
#     --covar-quant covar_quant.txt --covar-cat covar_categorical.txt \
#     --manifest phenotype_manifest.tsv --out prs_association.tsv
#
# MODEL
#
#   phenotype ~ scale(PRS) + sex + age_c + PC1..PC10 + (1 | family_id)
#
# Three parts of that are load-bearing:
#
# * (1 | family_id) -- ABCD contains twin and sibling pairs.  Fitting OLS over
#   related subjects treats each sibling as independent evidence and shrinks the
#   standard error toward a false positive.  The random intercept is the cheap
#   correct thing; it is also what the phenotype model itself uses.
#
# * PC1..PC10 -- adjust for population structure.  Necessary but NOT sufficient
#   for cross-ancestry validity, which is why the EUR subset is primary.
#
# * scale(PRS) -- the score's raw units are arbitrary (a sum of allele dosages
#   weighted by log-odds), so beta is only interpretable per SD of score.
#
# ANCESTRY STRATA
#
# "EUR" is defined by proximity to the EUR centroid in ancestry-PC space rather
# than by self-reported race, because the score's transferability depends on
# genetic ancestry.  With no labelled reference panel here, the operational
# definition is a Mahalanobis-style cut on the first few PCs around the largest
# cluster; the cut and the resulting N are recorded in the output so the
# definition is auditable rather than implicit.  If a labelled reference is
# available on the cluster, prefer projecting onto it -- see --eur-ids.
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
  make_option("--eur-ids",     type = "character", default = NULL,
              help = "optional file of IIDs defining the EUR subset; overrides the PC-based cut"),
  make_option("--eur-sd",      type = "double", default = 3,
              help = "PC-space cut for the EUR subset, in SDs [default 3]")
)))

stopifnot(!is.null(opt$`prs-dir`), !is.null(opt$out))

ph   <- fread(opt$pheno)
qcov <- fread(opt$`covar-quant`)
ccov <- fread(opt$`covar-cat`)
man  <- fread(opt$manifest)

# Phenotype file is FID IID <phenotypes...>; the manifest names them.
pheno_names <- man$name[man$name %in% names(ph)]
if (!length(pheno_names)) stop("no manifest phenotypes found in ", opt$pheno)

pc_cols <- grep("^PC[0-9]+$", names(qcov), value = TRUE)
if (!length(pc_cols)) {
  warning("no ancestry PCs in ", opt$`covar-quant`,
          " -- association tests will be unadjusted for population structure")
}

d <- merge(ph, qcov, by = c("FID", "IID"))
d <- merge(d, ccov[, c("FID", "IID", "sex", "site"), with = FALSE],
           by = c("FID", "IID"))
# FID is the family identifier in the GCTA export, so it IS family_id.
d[, family_id := FID]

# ---------------------------------------------------------------------------
# Ancestry strata.
# ---------------------------------------------------------------------------
if (!is.null(opt$`eur-ids`)) {
  eur <- fread(opt$`eur-ids`, header = FALSE)$V1
  d[, is_eur := IID %in% eur]
  eur_def <- paste0("supplied list (", opt$`eur-ids`, ")")
} else if (length(pc_cols) >= 2) {
  # Largest cluster on PC1/PC2 by a robust centre, then an SD cut.  This is a
  # pragmatic stand-in for reference projection, and it is reported as such.
  use <- pc_cols[1:min(4, length(pc_cols))]
  ctr <- d[, lapply(.SD, median, na.rm = TRUE), .SDcols = use]
  sds <- d[, lapply(.SD, sd, na.rm = TRUE), .SDcols = use]
  z <- as.matrix(d[, ..use])
  z <- sweep(z, 2, unlist(ctr), "-")
  z <- sweep(z, 2, unlist(sds), "/")
  dist <- sqrt(rowSums(z^2, na.rm = TRUE))
  d[, is_eur := dist <= opt$`eur-sd`]
  eur_def <- sprintf("PC-space cut: within %.1f SD of median on %s",
                     opt$`eur-sd`, paste(use, collapse = "/"))
} else {
  d[, is_eur := NA]
  eur_def <- "undetermined (no ancestry PCs)"
}
message(sprintf("EUR subset: %d of %d subjects [%s]",
                sum(d$is_eur, na.rm = TRUE), nrow(d), eur_def))

# ---------------------------------------------------------------------------
# Load scores.  PLINK 1.9 --score writes SCORESUM when 'sum' is requested.
# ---------------------------------------------------------------------------
prof <- list.files(opt$`prs-dir`, pattern = "^score_.*\\.profile$", full.names = TRUE)
if (!length(prof)) stop("no score_*.profile files in ", opt$`prs-dir`)

rows <- list()
for (f in prof) {
  tag <- sub("^score_", "", sub("\\.profile$", "", basename(f)))
  parts <- strsplit(tag, "_", fixed = TRUE)[[1]]
  disorder <- parts[1]
  threshold <- paste(parts[-1], collapse = "_")

  s <- fread(f)
  col <- if ("SCORESUM" %in% names(s)) "SCORESUM" else "SCORE"
  s <- s[, c("IID", col), with = FALSE]
  setnames(s, col, "PRS")
  n_snp_file <- file.path(opt$`prs-dir`,
                          paste0(disorder, "_", threshold, ".weights"))
  n_snp <- if (file.exists(n_snp_file)) nrow(fread(n_snp_file, header = FALSE)) else NA_integer_

  dd <- merge(d, s, by = "IID")
  if (!nrow(dd)) { warning("no ID overlap for ", tag); next }

  for (p in pheno_names) {
    for (stratum in c("EUR", "full")) {
      sub <- if (stratum == "EUR") dd[is_eur == TRUE] else dd
      sub <- sub[is.finite(get(p)) & is.finite(PRS)]
      # A random intercept needs replicated families to be identifiable.
      if (nrow(sub) < 100 || uniqueN(sub$family_id) < 50) next

      terms <- c("scale(PRS)", "sex", "age_c", pc_cols)
      terms <- terms[terms %in% names(sub) | terms == "scale(PRS)"]
      fml <- as.formula(sprintf("scale(%s) ~ %s + (1 | family_id)",
                                p, paste(terms, collapse = " + ")))
      fit <- try(lmer(fml, data = sub, REML = TRUE,
                      control = lmerControl(calc.derivs = FALSE)), silent = TRUE)
      if (inherits(fit, "try-error")) next
      co <- summary(fit)$coefficients
      if (!"scale(PRS)" %in% rownames(co)) next

      rows[[length(rows) + 1]] <- data.table(
        disorder = disorder, threshold = threshold, phenotype = p,
        stratum = stratum, n = nrow(sub), n_families = uniqueN(sub$family_id),
        n_snps = n_snp,
        # beta is in SD of phenotype per SD of score, so it is comparable
        # across phenotypes measured in different units.
        beta = co["scale(PRS)", "Estimate"],
        se   = co["scale(PRS)", "Std. Error"],
        p    = co["scale(PRS)", "Pr(>|t|)"])
    }
  }
}

res <- rbindlist(rows)
if (!nrow(res)) stop("no models fitted -- check ID overlap between scores and phenotypes")

# Bonferroni across thresholds within each disorder x phenotype x stratum:
# the threshold is a nuisance choice, not a hypothesis, so a hit that survives
# only at one threshold must be corrected for having looked at all of them.
res[, n_thresholds := uniqueN(threshold), by = .(disorder, phenotype, stratum)]
res[, p_adj := pmin(1, p * n_thresholds)]
res[, r2_partial := beta^2]   # phenotype is scaled, so beta^2 ~ variance explained
setorder(res, stratum, disorder, phenotype, threshold)

attr(res, "eur_definition") <- eur_def
fwrite(res, opt$out, sep = "\t")
cat(sprintf("wrote %s (%d rows)\nEUR definition: %s\n",
            opt$out, nrow(res), eur_def))
