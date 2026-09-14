#!/usr/bin/env Rscript
# Convert a PLINK1 binary fileset to SNPRelate GDS.
#
# GENESIS reads GDS, not .bed.  snpgdsBED2GDS is lossless for hard calls and
# keeps sample IDs verbatim, so the ID conventions checked in 00_check_inputs
# survive the conversion unchanged.
#
#   Rscript 01_make_gds.R --bfile PREFIX --out FILE.gds [--force]
#
# Idempotent: refuses to overwrite an existing .gds unless --force, because a
# half-written GDS from a killed job looks like a complete file to `ls`.
suppressPackageStartupMessages({
  library(optparse)
  library(SNPRelate)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bfile", type = "character", help = "PLINK fileset prefix"),
  make_option("--out",   type = "character", help = "output .gds path"),
  make_option("--force", action = "store_true", default = FALSE)
)))
for (r in c("bfile", "out")) if (is.null(opt[[r]])) stop("--", r, " is required")

for (ext in c(".bed", ".bim", ".fam")) {
  f <- paste0(opt$bfile, ext)
  if (!file.exists(f)) stop("missing ", f)
}

if (file.exists(opt$out) && !opt$force) {
  # Validate rather than skip blindly: a truncated GDS must not pass as done.
  g <- try(snpgdsOpen(opt$out), silent = TRUE)
  if (inherits(g, "try-error")) {
    stop(opt$out, " exists but does not open as GDS (truncated write?). ",
         "Remove it or pass --force.")
  }
  n <- length(read.gdsn(index.gdsn(g, "sample.id")))
  m <- length(read.gdsn(index.gdsn(g, "snp.id")))
  snpgdsClose(g)
  n_fam <- length(readLines(paste0(opt$bfile, ".fam")))
  m_bim <- length(readLines(paste0(opt$bfile, ".bim")))
  if (n == n_fam && m == m_bim) {
    cat(sprintf("EXISTS  %s (n=%d, m=%d) matches %s -- nothing to do\n",
                opt$out, n, m, opt$bfile))
    quit(status = 0)
  }
  stop(opt$out, " exists with n=", n, ", m=", m, " but ", opt$bfile,
       " has n=", n_fam, ", m=", m_bim, ". Stale conversion; pass --force.")
}

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
snpgdsBED2GDS(
  bed.fn = paste0(opt$bfile, ".bed"),
  bim.fn = paste0(opt$bfile, ".bim"),
  fam.fn = paste0(opt$bfile, ".fam"),
  out.gdsfn = opt$out,
  cvt.chr = "int",
  verbose = TRUE
)

# Verify the conversion round-trips the counts (the check 00 does for .bed).
g <- snpgdsOpen(opt$out)
n <- length(read.gdsn(index.gdsn(g, "sample.id")))
m <- length(read.gdsn(index.gdsn(g, "snp.id")))
snpgdsClose(g)
n_fam <- length(readLines(paste0(opt$bfile, ".fam")))
m_bim <- length(readLines(paste0(opt$bfile, ".bim")))
stopifnot(n == n_fam, m == m_bim)
cat(sprintf("WROTE  %s  (n=%d samples, m=%d variants)\n", opt$out, n, m))
