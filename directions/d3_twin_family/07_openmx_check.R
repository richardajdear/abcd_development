#!/usr/bin/env Rscript
# 07_openmx_check.R -- cross-check ace.py's univariate ACE fits against OpenMx.
# OPTIONAL, NOT RUN on the laptop: OpenMx has no conda build for osx-arm64, and building
# it from CRAN source fails in the conda R toolchain (no arm64-apple-darwin20 clang
# wrapper for mvtnorm, no cmake for RcppParallel; tried 2026-09-30). Run it wherever
# OpenMx is installed (e.g. an R on CSD3) after copying work/ there; it stops otherwise.
# ace.py is validated instead by 00_ace_selftest.py (parameter recovery on simulated twins).
# Input : work/traits.parquet, work/pairs.csv
# Output: results/openmx_check.tsv   (OpenMx vs ace.py h2 / c2 / -2LL per trait)
# Run from the repo root: Rscript directions/d3_twin_family/07_openmx_check.R   (env r)
args <- commandArgs(FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
if (!requireNamespace("OpenMx", quietly = TRUE)) stop("OpenMx is not installed in this R; see header")
suppressMessages({library(OpenMx); library(arrow); library(data.table)})
mxOption(NULL, "Default optimizer", "SLSQP")
T <- as.data.table(read_parquet(file.path(HERE, "work", "traits.parquet")), keep.rownames = FALSE)
if (!"subject" %in% names(T)) T[, subject := read_parquet(file.path(HERE, "work", "traits.parquet"),
                                                           as_data_frame = TRUE) |> rownames()]
P <- fread(file.path(HERE, "work", "pairs.csv"))
ours <- fread(file.path(HERE, "results", "ace_univariate.tsv"))[model == "ACE"]

fit_ace <- function(tr) {
  get <- function(ty) {
    p <- P[type == ty]
    d <- data.frame(t1 = T[[tr]][match(p$id1, T$subject)], t2 = T[[tr]][match(p$id2, T$subject)])
    d[!(is.na(d$t1) & is.na(d$t2)), ]
  }
  mzd <- get("MZ"); dzd <- get("DZ")
  sv <- sqrt(var(c(mzd$t1, mzd$t2, dzd$t1, dzd$t2), na.rm = TRUE) / 3)
  base <- list(
    mxMatrix("Full", 1, 1, TRUE, sv, "a", name = "a"), mxMatrix("Full", 1, 1, TRUE, sv, "c", name = "c"),
    mxMatrix("Full", 1, 1, TRUE, sv, "e", name = "e"), mxMatrix("Full", 1, 2, TRUE, 0, "mu", name = "M"),
    mxAlgebra(a %*% t(a), name = "A"), mxAlgebra(c %*% t(c), name = "C"), mxAlgebra(e %*% t(e), name = "E"))
  grp <- function(nm, k, dat) mxModel(nm, base,
    mxAlgebra(rbind(cbind(A + C + E, k * A + C), cbind(k * A + C, A + C + E)), name = "S"),
    mxData(dat, type = "raw"), mxExpectationNormal("S", "M", dimnames = c("t1", "t2")), mxFitFunctionML())
  m <- mxModel("ACE", grp("MZ", 1, mzd), grp("DZ", 0.5, dzd), mxFitFunctionMultigroup(c("MZ", "DZ")))
  r <- mxRun(m, silent = TRUE)
  va <- mxEval(MZ.A, r)[1]; vc <- mxEval(MZ.C, r)[1]; ve <- mxEval(MZ.E, r)[1]; vt <- va + vc + ve
  data.table(trait = tr, openmx_h2 = va / vt, openmx_c2 = vc / vt, openmx_minus2ll = r$output$Minus2LogLikelihood,
             openmx_status = r$output$status$code)
}
traits <- c("dCT_r", "dCT_rq", "CT0_r", "dT1T2_r", "pub_timing_r", "pub_tempo_r", "depress_chg_r", "pfactor_chg_r")
out <- rbindlist(lapply(traits, fit_ace))
out <- merge(out, ours[, .(trait, acepy_h2 = h2, acepy_c2 = c2, acepy_minus2ll = minus2ll)], by = "trait")
out[, `:=`(d_h2 = openmx_h2 - acepy_h2, d_minus2ll = openmx_minus2ll - acepy_minus2ll)]
fwrite(out, file.path(HERE, "results", "openmx_check.tsv"), sep = "\t")
print(out[, .(trait, openmx_h2, acepy_h2, openmx_c2, acepy_c2, d_minus2ll)], digits = 4)
