#!/usr/bin/env Rscript
#' Fit the regional mixed-effects models for one assembled run.
#'
#' Reads   out/<run_id>/model_table.parquet
#' Writes  out/<run_id>/fits/{blups,fixed,varcomp,diagnostics}.parquet
#'
#' Usage (the run comes from $ABCD_CONFIG; no run_id to type):
#'   export ABCD_CONFIG=ct_70_genetic
#'   Rscript R/fit_lmm.R --cores 8
#'   Rscript R/fit_lmm.R --config ct_70_baseline        # override the export
#'   Rscript R/fit_lmm.R --regions lh_insula,lh_cuneus  # subset, for debugging
#'   Rscript R/fit_lmm.R --run-dir out/<run_id>         # explicit, still works
#'
#' What this does differently from the 5.1 pipeline
#' ------------------------------------------------
#' 1. **Singular fits are kept and flagged, never dropped.** The old code
#'    returned NULL for singular regions, so the region set silently varied
#'    between runs and downstream maps had holes. A singular fit is a variance
#'    component estimated at zero -- the fixed effects and BLUPs are still
#'    usable, and the flag travels with them.
#' 2. **BLUPs carry their conditional SD.** The random slope for a subject with
#'    two visits is shrunk far harder than for one with three. Without the
#'    conditional SD you cannot tell a genuinely flat trajectory from an
#'    uninformative one, and that difference correlates with visit count --
#'    i.e. with attrition, which is not random.
#' 3. **Nothing is predicted out of sample.** The old code differenced model
#'    predictions at ages 9 and 16 when the oldest scan is 15.75 and most
#'    subjects span ~2 years. The slope itself is the phenotype.

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lme4)
  library(yaml)
})

here <- function(...) file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ...)
source(here("model_spec.R"))

# --------------------------------------------------------------------------

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run-dir", type = "character", default = NULL,
              help = "run directory; omit to resolve from $ABCD_CONFIG"),
  make_option("--config", type = "character", default = NULL,
              help = "config name, overriding $ABCD_CONFIG"),
  make_option("--regions", type = "character", default = NULL, help = "comma-separated subset of labels (debugging)"),
  make_option("--cores", type = "integer", default = 1L, help = "parallel workers"),
  make_option("--out-name", type = "character", default = "fits", help = "subdirectory for outputs")
)))

#' Resolve the run directory, preferring an explicit --run-dir.
#'
#' Delegates to `python -m abcd.run_dir` rather than recomputing the config
#' hash here: the hash is the run's identity, and a second implementation in R
#' would drift silently the moment a config field is added, pointing these fits
#' at a stale run directory.
resolve_run_dir <- function(opt) {
  if (!is.null(opt$`run-dir`)) return(opt$`run-dir`)
  args <- c("-m", "abcd.run_dir")
  if (!is.null(opt$config)) args <- c(args, opt$config)
  # The interpreter that has abcd's dependencies; `make` and tools/rerun_local.sh
  # export PY.  A bare "python" resolved to a shim without pandas and failed here.
  py <- Sys.getenv("PY", unset = "python")
  res <- suppressWarnings(system2(py, args, stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
    stop("could not resolve a run directory. Either pass --run-dir, or set\n",
         "  export ABCD_CONFIG=ct_70_genetic\n",
         "python said: ", paste(res, collapse = " "), call. = FALSE)
  }
  tail(res[nzchar(res)], 1)
}

run_dir <- normalizePath(resolve_run_dir(opt), mustWork = TRUE)
cfg <- yaml::read_yaml(file.path(run_dir, "config.yaml"))
out_dir <- file.path(run_dir, opt$`out-name`)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- arrow::read_parquet(file.path(run_dir, "model_table.parquet")) |>
  mutate(
    subject   = factor(subject),
    site      = factor(site),
    family_id = factor(family_id),
    sex       = factor(sex, levels = c("F", "M"))
  )

form <- build_formula(cfg)
message("run:     ", basename(run_dir))
message("formula: ", deparse1(form))
message("data:    ", nrow(dat), " rows, ", nlevels(dat$subject), " subjects, ",
        length(unique(dat$label)), " regions")

labels <- sort(unique(dat$label))
if (!is.null(opt$regions)) labels <- intersect(labels, strsplit(opt$regions, ",")[[1]])

# --------------------------------------------------------------------------

#' Fit one region and return a list of tidy frames.
#'
#' Returns even when the fit is singular or fails to converge; the caller
#' decides what to do about it. Only a hard error yields status "error".
fit_region <- function(lab) {
  d <- dat[dat$label == lab, , drop = FALSE]
  d <- droplevels(d)

  t0 <- Sys.time()
  m <- tryCatch(
    withCallingHandlers(
      lme4::lmer(form, data = d, control = lmer_control(), REML = TRUE),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(m, "error")) {
    return(list(diagnostics = tibble::tibble(
      label = lab, status = "error", message = conditionMessage(m),
      singular = NA, converged = NA, n_obs = nrow(d), n_subjects = nlevels(d$subject),
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )))
  }

  # --- random effects with conditional variances --------------------------
  re <- lme4::ranef(m, condVar = TRUE)
  rs <- re$subject
  pv <- attr(rs, "postVar")               # k x k x n_subjects
  cond_sd <- t(apply(pv, 3, function(x) sqrt(diag(as.matrix(x)))))
  colnames(cond_sd) <- paste0("se_", colnames(rs))

  blups <- tibble::as_tibble(rs, rownames = "subject") |>
    bind_cols(tibble::as_tibble(cond_sd)) |>
    mutate(label = lab, .before = 1)
  names(blups) <- sub("^\\(Intercept\\)$", "re_intercept", names(blups))
  names(blups) <- sub("^age_c$",           "re_slope",     names(blups))
  names(blups) <- sub("^se_\\(Intercept\\)$", "se_re_intercept", names(blups))
  names(blups) <- sub("^se_age_c$",           "se_re_slope",     names(blups))

  # --- fixed effects ------------------------------------------------------
  cf <- summary(m)$coefficients
  fixed <- tibble::tibble(
    label = lab, term = rownames(cf),
    estimate = cf[, "Estimate"], std_error = cf[, "Std. Error"],
    statistic = cf[, "t value"]
  )

  # --- variance components ------------------------------------------------
  vc <- as.data.frame(lme4::VarCorr(m))
  varcomp <- tibble::tibble(
    label = lab, grp = vc$grp, var1 = vc$var1, var2 = vc$var2,
    vcov = vc$vcov, sdcor = vc$sdcor
  )

  # --- diagnostics --------------------------------------------------------
  conv <- m@optinfo$conv$lme4
  diagnostics <- tibble::tibble(
    label = lab,
    status = "ok",
    message = if (length(warns)) paste(unique(warns), collapse = " | ") else NA_character_,
    singular = lme4::isSingular(m),
    converged = length(conv$messages) == 0,
    n_obs = stats::nobs(m),
    n_subjects = nlevels(d$subject),
    logLik = as.numeric(stats::logLik(m)),
    AIC = stats::AIC(m),
    BIC = stats::BIC(m),
    sigma = stats::sigma(m),
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )

  list(blups = blups, fixed = fixed, varcomp = varcomp, diagnostics = diagnostics)
}

# a warning collector visible to fit_region's handler
warns <- character(0)
fit_one <- function(lab) { warns <<- character(0); fit_region(lab) }

# --------------------------------------------------------------------------

message("fitting ", length(labels), " regions on ", opt$cores, " core(s) ...")
t_all <- Sys.time()

if (opt$cores > 1L) {
  suppressPackageStartupMessages(library(parallel))
  res <- parallel::mclapply(labels, fit_one, mc.cores = opt$cores)
} else {
  res <- lapply(seq_along(labels), function(i) {
    r <- fit_one(labels[i])
    message(sprintf("  [%2d/%2d] %-22s %s", i, length(labels), labels[i],
                    if (isTRUE(r$diagnostics$singular)) "singular" else "ok"))
    r
  })
}

gather <- function(key) bind_rows(lapply(res, `[[`, key))
diagnostics <- gather("diagnostics")

arrow::write_parquet(gather("blups"),   file.path(out_dir, "blups.parquet"))
arrow::write_parquet(gather("fixed"),   file.path(out_dir, "fixed.parquet"))
arrow::write_parquet(gather("varcomp"), file.path(out_dir, "varcomp.parquet"))
arrow::write_parquet(diagnostics,       file.path(out_dir, "diagnostics.parquet"))

writeLines(yaml::as.yaml(list(
  formula = deparse1(form),
  n_regions = length(labels),
  n_singular = sum(diagnostics$singular, na.rm = TRUE),
  n_nonconverged = sum(!diagnostics$converged, na.rm = TRUE),
  n_error = sum(diagnostics$status == "error"),
  elapsed_seconds = round(as.numeric(difftime(Sys.time(), t_all, units = "secs")), 1),
  fitted_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  R_version = R.version.string,
  lme4_version = as.character(utils::packageVersion("lme4"))
)), file.path(out_dir, "fit_manifest.yaml"))

message(sprintf(
  "done in %.1f min | singular %d/%d | non-converged %d | errors %d -> %s",
  as.numeric(difftime(Sys.time(), t_all, units = "mins")),
  sum(diagnostics$singular, na.rm = TRUE), nrow(diagnostics),
  sum(!diagnostics$converged, na.rm = TRUE),
  sum(diagnostics$status == "error"), out_dir
))
