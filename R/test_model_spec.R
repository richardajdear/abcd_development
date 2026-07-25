#!/usr/bin/env Rscript
#' Tests for the formula builder.
#'
#' Deliberately dependency-free (no testthat) so it runs anywhere R runs,
#' including a bare HPC node. Usage:  Rscript R/test_model_spec.R
#'
#' The critical test here is `test_config_schema_in_step`: the formula builder
#' reads field names owned by src/abcd/config.py. If a field is renamed on the
#' Python side, R would silently fall back to a default and fit the wrong
#' model -- e.g. dropping the family random effect from a genetics phenotype.
#' That test reads the real config files and fails loudly instead.

here <- function(...) file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ...)
source(here("model_spec.R"))
suppressPackageStartupMessages(library(yaml))

failures <- 0L
check <- function(desc, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    message("  ERROR ", desc, ": ", conditionMessage(e)); FALSE
  })
  cat(sprintf("  [%s] %s\n", if (ok) "ok  " else "FAIL", desc))
  if (!ok) failures <<- failures + 1L
  invisible(ok)
}

cat("fixed_terms\n")
check("linear age", fixed_terms("linear", "none") == "sex + age_c")
check("quadratic age adds a squared term",
      grepl("I(age_c^2)", fixed_terms("quadratic", "none"), fixed = TRUE))
check("spline age uses ns() with the requested df",
      fixed_terms("spline", "none", spline_df = 4) == "sex + splines::ns(age_c, df = 4)")
check("global covariate enters split within/between",
      fixed_terms("linear", "observed_mean") ==
        "sex + age_c + global_between_c + global_within")
check("global covariate 'none' omits both global terms",
      !grepl("global", fixed_terms("linear", "none")))
check("sex interaction multiplies through the age basis",
      grepl("sex \\* \\(age_c\\)", fixed_terms("linear", "none", sex_interaction = TRUE)))
check("unknown age basis errors",
      inherits(try(fixed_terms("cubic-ish", "none"), silent = TRUE), "try-error"))

cat("random_terms\n")
check("correlated slope is the default structure",
      random_terms("slope_correlated", FALSE, FALSE) == "(1 + age_c | subject)")
check("independent slope splits the subject terms",
      random_terms("slope_independent", FALSE, FALSE) ==
        "(1 | subject) + (0 + age_c | subject)")
check("intercept-only omits the slope",
      random_terms("intercept_only", FALSE, FALSE) == "(1 | subject)")
check("site and family enter as crossed intercepts",
      random_terms("slope_correlated", TRUE, TRUE) ==
        "(1 + age_c | subject) + (1 | site) + (1 | family_id)")
check("unknown re_structure errors",
      inherits(try(random_terms("wishful", TRUE, TRUE), silent = TRUE), "try-error"))

cat("build_formula\n")
cfg <- list(age_basis = "linear", global_covariate = "observed_mean",
            re_structure = "slope_correlated", site_effect = TRUE,
            family_effect = TRUE, covariates = list("sex"))
f <- build_formula(cfg)
check("returns a formula", inherits(f, "formula"))
check("response is value", all.vars(f)[1] == "value")
check("full baseline formula is exactly as intended",
      deparse1(f) == paste("value ~ sex + age_c + global_between_c + global_within +",
                           "(1 + age_c | subject) + (1 | site) + (1 | family_id)"))
cfg_int <- cfg
cfg_int$covariates <- list("sex", "sex:age")   # not modifyList(): it recurses into unnamed lists
check("sex:age in covariates turns on the interaction",
      grepl("sex \\* ", deparse1(build_formula(cfg_int))))
check("plain sex covariates leave the interaction off",
      !grepl("sex \\* ", deparse1(build_formula(cfg))))

cat("config schema in step with Python\n")
cfg_files <- list.files(here("..", "configs"), pattern = "\\.yaml$", full.names = TRUE)
check("config files found", length(cfg_files) > 0)
required <- c("age_basis", "spline_df", "re_structure", "global_covariate",
              "site_effect", "family_effect", "covariates")
for (p in cfg_files) {
  y <- yaml::read_yaml(p)
  missing <- setdiff(required, names(y))
  check(sprintf("%s has every field build_formula() reads", basename(p)),
        length(missing) == 0)
  if (length(missing)) message("    missing: ", paste(missing, collapse = ", "))
  check(sprintf("%s builds a valid formula", basename(p)),
        inherits(build_formula(y), "formula"))
}

cat(sprintf("\n%s: %d failure(s)\n", if (failures == 0L) "PASS" else "FAIL", failures))
quit(status = if (failures == 0L) 0L else 1L)
