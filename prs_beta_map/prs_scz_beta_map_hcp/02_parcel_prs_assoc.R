#!/usr/bin/env Rscript
# Step 2 (CSD3): regress every HCP-MMP parcel's thinning slope on one polygenic
# score.  Output is a 358-row beta map, plus a whole-cortex row that must
# reproduce the tabled global_slope beta for the same cell.
#
# MODEL (tools/prs_assoc.R, unchanged -- so a parcel beta and the tabled
# whole-cortex beta are the same estimator on different outcomes):
#
#   scale(slope_parcel) ~ scale(PRS) + sex + age_c + PC1..PC10 + (1 | family_id)
#
# A second fit per parcel adds scale(cortex_mean) -- the whole-cortex mean of the
# 358 slopes -- as a covariate (beta_cond): the parcel-specific association once
# the global one is held fixed.  r_global is each parcel's correlation with the
# cortex mean (the 'loading' map a purely global effect would reproduce).
#
# The outcome is standardised within parcel, so beta is SD/SD, the unit of
# every PRS table in genetic_analysis/work/results_70tab_hcp/.  The pooled arm
# reads the within-ancestry-standardised score (_zanc); the EUR arm reads the
# raw score on the EUR anchor set -- the score directory decides, as in
# step9_scz2025_assoc.sbatch.
#
#   Rscript 02_parcel_prs_assoc.R --slopes work/parcel_slopes_hcp70.tsv.gz \
#     --pheno-dir $PHENO_DIR --prs-dir <score dir> [--eur-ids <keep>] \
#     --stratum full|EUR --cell SCZ25_META_SBayesRC_zanc --out results/beta_map_<cell>.tsv
#
# --prs-file <file> [--prs-col <col>] accepts one score table instead of a
# PLINK profile directory (the laptop pipeline test uses it).
suppressPackageStartupMessages({ library(data.table); library(lme4); library(lmerTest) })

# --key value argument parser (base R; the local R has no optparse)
args <- commandArgs(trailingOnly = TRUE)
opt <- list(stratum = "full", cell = "cell", `max-parcels` = NA_integer_)
i <- 1
while (i <= length(args)) {
  stopifnot(startsWith(args[i], "--"), i < length(args))
  opt[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2
}
opt$`max-parcels` <- as.integer(opt$`max-parcels`)
stopifnot(!is.null(opt$slopes), !is.null(opt$`pheno-dir`), !is.null(opt$out),
          opt$stratum %in% c("full", "EUR"))

token <- function(x) sub("^.*?([A-Z0-9]{8})$", "\\1", x)   # align_export.py::token

# ---- covariates: the aligned 7.0 export (FID = family id) -------------------
pd   <- opt$`pheno-dir`
ph   <- fread(file.path(pd, "phenotypes_gcta.txt"))
qcov <- fread(file.path(pd, "covar_quant.txt"))
ccov <- fread(file.path(pd, "covar_categorical.txt"))
d <- merge(ph[, .(FID, IID, global_slope)], qcov, by = c("FID", "IID"))
d <- merge(d, ccov[, .(FID, IID, sex)], by = c("FID", "IID"))
d[, family_id := as.character(FID)]
d[, age_c := baseline_age - mean(baseline_age, na.rm = TRUE)]
d[, tok := token(IID)]
pc_cols <- grep("^PC[0-9]+$", names(qcov), value = TRUE)
stopifnot(length(pc_cols) == 10)

# ---- score ------------------------------------------------------------------
if (!is.null(opt$`prs-dir`)) {
  prof <- list.files(opt$`prs-dir`, pattern = "^score_.*\\.profile$", full.names = TRUE)
  if (length(prof) != 1)
    stop("expected exactly one score_*.profile in ", opt$`prs-dir`, ", found ", length(prof),
         ": ", paste(basename(prof), collapse = " "), " -- pass the C+T threshold dir, not its parent")
  s <- fread(prof)
  col <- if ("SCORESUM" %in% names(s)) "SCORESUM" else "SCORE"
  s <- s[, .(tok = token(IID), PRS = get(col))]
  score_src <- prof
} else if (!is.null(opt$`prs-file`)) {
  s <- fread(opt$`prs-file`)
  col <- if (!is.null(opt$`prs-col`)) opt$`prs-col` else setdiff(names(s), c("FID", "IID"))[1]
  s <- s[, .(tok = token(IID), PRS = get(col))]
  score_src <- opt$`prs-file`
} else stop("give --prs-dir or --prs-file")
d <- merge(d, s, by = "tok")

if (opt$stratum == "EUR") {
  stopifnot(!is.null(opt$`eur-ids`))
  eur <- token(fread(opt$`eur-ids`, header = FALSE)[[1]])
  d <- d[tok %in% eur]
}
d <- d[is.finite(PRS)]

# ---- per-parcel slopes ------------------------------------------------------
# gz via a pipe: the GENESIS R on CSD3 has data.table without R.utils, so
# fread() cannot open a .gz path directly
sl <- if (grepl("\\.gz$", opt$slopes)) fread(cmd = paste("gzip -dc", shQuote(opt$slopes))) else fread(opt$slopes)
sl[, tok := token(IID)]; sl[, IID := NULL]
parcels <- setdiff(names(sl), "tok")
if (!is.na(opt$`max-parcels`)) parcels <- parcels[seq_len(opt$`max-parcels`)]
d <- merge(d, sl, by = "tok")
message(sprintf("[%s] stratum %s: n = %d children, %d families, %d parcels, score %s",
                opt$cell, opt$stratum, nrow(d), uniqueN(d$family_id), length(parcels),
                basename(score_src)))
stopifnot(nrow(d) > 1000)

rhs <- paste("scale(PRS) + sex + age_c +", paste(pc_cols, collapse = " + "), "+ (1 | family_id)")
lmer_q <- function(fml) suppressMessages(suppressWarnings(lmer(as.formula(fml), data = d, REML = TRUE)))
fit_one <- function(y) {
  d[, .y := scale(get(y))[, 1]]
  m <- lmer_q(paste(".y ~", rhs))
  co <- summary(m)$coefficients["scale(PRS)", ]
  out <- data.table(label = y, beta = co[["Estimate"]], se = co[["Std. Error"]],
                    t = co[["t value"]], p = co[["Pr(>|t|)"]], n = nobs(m))
  if (y == "cortex_mean") {
    out[, `:=`(beta_cond = NA_real_, se_cond = NA_real_, t_cond = NA_real_, p_cond = NA_real_,
               r_global = 1)]
  } else {
    # the parcel-specific part: the same model with the whole-cortex mean slope
    # as a covariate.  A score that only moves the whole cortex still produces a
    # structured beta map (each parcel's loading on the mean), so this is the
    # column that says whether the map has spatial content of its own.
    m2 <- lmer_q(paste(".y ~ scale(cortex_mean) +", rhs))
    co2 <- summary(m2)$coefficients["scale(PRS)", ]
    out[, `:=`(beta_cond = co2[["Estimate"]], se_cond = co2[["Std. Error"]],
               t_cond = co2[["t value"]], p_cond = co2[["Pr(>|t|)"]],
               r_global = cor(d$.y, d$cortex_mean))]
  }
  out
}

# the whole-cortex row first: mean of the 358 parcel slopes = the per-region
# global_slope, so this must match the tabled beta for the same cell
d[, cortex_mean := rowMeans(.SD), .SDcols = setdiff(names(sl), "tok")]
r_check <- cor(d$cortex_mean, d$global_slope)
message(sprintf("cortex mean of the parcel slopes vs export global_slope: r = %.4f", r_check))
res <- rbindlist(c(list(fit_one("cortex_mean")), lapply(parcels, fit_one)))
res[, `:=`(cell = opt$cell, stratum = opt$stratum, score = basename(score_src),
           r_cortexmean_globalslope = r_check)]
setcolorder(res, c("cell", "stratum", "label"))
fwrite(res, opt$out, sep = "\t")
message(sprintf("global (cortex_mean) beta = %.4f (se %.4f, p %.3g); %d parcel rows -> %s",
                res$beta[1], res$se[1], res$p[1], nrow(res) - 1, opt$out))
