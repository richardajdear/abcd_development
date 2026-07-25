#' Model specification for the ABCD longitudinal pipeline.
#'
#' Formulas are *built*, not written out per analysis, so that a config field
#' maps to exactly one formula and there is no risk of two scripts fitting
#' subtly different models under the same name.
#'
#' Design decisions this file encodes, and why:
#'
#'  * **Random slope AND intercept by subject.** The scientific target is the
#'    subject-specific rate of change, so it must be a modelled random effect,
#'    not a difference between two fitted values.
#'  * **Site as a random intercept.** 22 sites, near-constant within subject.
#'    Modelling it as a fixed effect burns 21 df and treats the sites as the
#'    population of interest, which they are not.
#'  * **Family as a crossed random intercept.** ~5,989 families for ~6,937
#'    subjects, including twins. A genetics phenotype that ignores family
#'    structure inherits shared-environment and relatedness variance directly
#'    into the trait. The 5.1 pipeline omitted this entirely.
#'  * **Global covariate split within/between subject.** A subject whose whole
#'    cortex is thin is a different phenomenon from a subject whose cortex
#'    thinned at this visit; a single global term conflates them.
#'  * **Sex-by-age interaction** is available but off by default -- puberty
#'    timing differs by sex, so it matters, but it changes what the random
#'    slope means (deviation from a sex-specific mean), which is a deliberate
#'    choice rather than a default.

suppressPackageStartupMessages({
  library(lme4)
})

#' Build the fixed-effect part of the formula.
#'
#' @param age_basis "linear", "quadratic", or "spline"
#' @param global_covariate "observed_mean" or "none"
#' @param sex_interaction logical
#' @param spline_df df for the natural spline basis
fixed_terms <- function(age_basis = "linear",
                        global_covariate = "observed_mean",
                        sex_interaction = FALSE,
                        spline_df = 3) {
  age <- switch(
    age_basis,
    linear    = "age_c",
    quadratic = "age_c + I(age_c^2)",
    spline    = sprintf("splines::ns(age_c, df = %d)", spline_df),
    stop("unknown age_basis: ", age_basis)
  )
  terms <- c("sex", age)
  if (sex_interaction) {
    terms <- c("sex", sprintf("sex * (%s)", age))
    terms <- unique(terms)
  }
  if (global_covariate != "none") {
    terms <- c(terms, "global_between_c", "global_within")
  }
  paste(terms, collapse = " + ")
}

#' Build the random-effect part of the formula.
#'
#' @param re_structure one of "intercept_only", "slope_correlated",
#'   "slope_independent"
#' @param include_site include a site random intercept
#' @param include_family include a family random intercept
random_terms <- function(re_structure = "slope_correlated",
                         include_site = TRUE,
                         include_family = TRUE) {
  subj <- switch(
    re_structure,
    intercept_only    = "(1 | subject)",
    slope_correlated  = "(1 + age_c | subject)",
    slope_independent = "(1 | subject) + (0 + age_c | subject)",
    stop("unknown re_structure: ", re_structure)
  )
  parts <- subj
  if (include_site)   parts <- c(parts, "(1 | site)")
  if (include_family) parts <- c(parts, "(1 | family_id)")
  paste(parts, collapse = " + ")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Assemble the full model formula for one region from a RunConfig.
#'
#' Field names here must match src/abcd/config.py exactly. The Python side
#' owns the schema; this function is the only place R interprets it, and
#' test_model_spec.R asserts the two stay in step.
build_formula <- function(cfg) {
  sex_int <- "sex:age" %in% (cfg$covariates %||% character(0))
  f <- sprintf(
    "value ~ %s + %s",
    fixed_terms(cfg$age_basis, cfg$global_covariate, sex_int, cfg$spline_df %||% 3),
    random_terms(cfg$re_structure,
                 isTRUE(cfg$site_effect), isTRUE(cfg$family_effect))
  )
  stats::as.formula(f)
}

#' Standard lmer control: bobyqa converges more reliably on these data, and
#' we do NOT want lmer's default convergence warnings suppressed.
lmer_control <- function(maxfun = 2e5) {
  lme4::lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = maxfun),
    calc.derivs = TRUE
  )
}
