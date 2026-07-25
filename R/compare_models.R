#!/usr/bin/env Rscript
#' Compare model specifications on a common subset of regions.
#'
#' Every specification is fitted to the SAME rows, so AIC/BIC are comparable.
#' ML (not REML) is used whenever fixed effects differ between specifications,
#' because REML likelihoods of models with different fixed effects are not
#' comparable -- a mistake that is easy to make and hard to notice.
#'
#' Usage:
#'   Rscript R/compare_models.R --run-dir out/<run_id> --n-regions 12
#'
#' Writes out/<run_id>/comparison/{model_comparison,slope_agreement}.parquet

suppressPackageStartupMessages({
  library(optparse); library(arrow); library(dplyr); library(tidyr)
  library(lme4); library(purrr)
})

here <- function(...) file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ...)
source(here("model_spec.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run-dir", type = "character"),
  make_option("--n-regions", type = "integer", default = 12L),
  make_option("--cores", type = "integer", default = 6L),
  make_option("--seed", type = "integer", default = 1L)
)))
run_dir <- normalizePath(opt$`run-dir`, mustWork = TRUE)
out_dir <- file.path(run_dir, "comparison"); dir.create(out_dir, showWarnings = FALSE)

dat <- read_parquet(file.path(run_dir, "model_table.parquet")) |>
  mutate(across(c(subject, site, family_id), factor),
         sex = factor(sex, levels = c("F", "M")))

set.seed(opt$seed)
labels <- sample(sort(unique(dat$label)), min(opt$`n-regions`, n_distinct(dat$label)))
message("comparing on ", length(labels), " regions: ", paste(labels, collapse = ", "))

#' Specifications to compare. Each is a full formula string; naming them here
#' keeps the comparison honest -- no specification can be quietly dropped
#' because it looked bad.
SPECS <- list(
  `1 intercept only`          = "value ~ sex + age_c + (1 | subject)",
  `2 slope, uncorrelated`     = "value ~ sex + age_c + (1 | subject) + (0 + age_c | subject)",
  `3 slope, correlated`       = "value ~ sex + age_c + (1 + age_c | subject)",
  `4 + site`                  = "value ~ sex + age_c + (1 + age_c | subject) + (1 | site)",
  `5 + site + family`         = "value ~ sex + age_c + (1 + age_c | subject) + (1 | site) + (1 | family_id)",
  `6 + global (split)`        = "value ~ sex + age_c + global_between_c + global_within + (1 + age_c | subject) + (1 | site) + (1 | family_id)",
  `7 quadratic age`           = "value ~ sex + age_c + I(age_c^2) + global_between_c + global_within + (1 + age_c | subject) + (1 | site) + (1 | family_id)",
  `8 spline age (df=3)`       = "value ~ sex + splines::ns(age_c, df = 3) + global_between_c + global_within + (1 + age_c | subject) + (1 | site) + (1 | family_id)",
  `9 sex x age`               = "value ~ sex * age_c + global_between_c + global_within + (1 + age_c | subject) + (1 | site) + (1 | family_id)"
)

fit_spec <- function(lab, nm) {
  d <- droplevels(dat[dat$label == lab, , drop = FALSE])
  t0 <- Sys.time()
  m <- tryCatch(suppressMessages(suppressWarnings(
    lmer(as.formula(SPECS[[nm]]), data = d, control = lmer_control(), REML = FALSE)
  )), error = function(e) e)
  if (inherits(m, "error")) {
    return(tibble::tibble(label = lab, spec = nm, status = "error",
                          AIC = NA_real_, BIC = NA_real_, logLik = NA_real_,
                          df = NA_integer_, singular = NA, sigma = NA_real_,
                          seconds = NA_real_))
  }
  tibble::tibble(
    label = lab, spec = nm, status = "ok",
    AIC = AIC(m), BIC = BIC(m), logLik = as.numeric(logLik(m)),
    df = attr(logLik(m), "df"), singular = isSingular(m), sigma = sigma(m),
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

grid <- expand.grid(label = labels, spec = names(SPECS), stringsAsFactors = FALSE)
res <- parallel::mcmapply(fit_spec, grid$label, grid$spec,
                          SIMPLIFY = FALSE, mc.cores = opt$cores)
comp <- bind_rows(res)
write_parquet(comp, file.path(out_dir, "model_comparison.parquet"))

#' Does the specification change the phenotype, or only the fit statistics?
#' If the subject slopes agree at r > 0.99 the choice is inconsequential
#' downstream; if they do not, the choice IS the phenotype definition.
slope_of <- function(lab, nm) {
  d <- droplevels(dat[dat$label == lab, , drop = FALSE])
  m <- tryCatch(suppressMessages(suppressWarnings(
    lmer(as.formula(SPECS[[nm]]), data = d, control = lmer_control(), REML = TRUE)
  )), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  re <- ranef(m)$subject
  if (!"age_c" %in% names(re)) return(NULL)
  tibble::tibble(label = lab, spec = nm, subject = rownames(re), slope = re$age_c)
}

with_slope <- setdiff(names(SPECS), "1 intercept only")
ref <- "6 + global (split)"
agree <- map_dfr(labels, function(lab) {
  sl <- map_dfr(with_slope, ~ slope_of(lab, .x))
  if (!nrow(sl)) return(NULL)
  w <- pivot_wider(sl, names_from = spec, values_from = slope)
  map_dfr(setdiff(with_slope, ref), function(nm) {
    tibble::tibble(label = lab, spec = nm,
                   r_vs_ref = suppressWarnings(cor(w[[nm]], w[[ref]],
                                                   use = "complete.obs")))
  })
})
write_parquet(agree, file.path(out_dir, "slope_agreement.parquet"))

message("\n== mean AIC across regions (lower is better) ==")
print(comp |> group_by(spec) |>
        summarise(mean_AIC = mean(AIC, na.rm = TRUE),
                  n_singular = sum(singular, na.rm = TRUE),
                  median_sec = median(seconds, na.rm = TRUE)) |>
        arrange(mean_AIC) |> as.data.frame(), digits = 6)

message("\n== subject-slope agreement with '", ref, "' ==")
print(agree |> group_by(spec) |>
        summarise(mean_r = mean(r_vs_ref, na.rm = TRUE),
                  min_r = min(r_vs_ref, na.rm = TRUE)) |>
        arrange(desc(mean_r)) |> as.data.frame(), digits = 4)
message("\nwritten -> ", out_dir)
