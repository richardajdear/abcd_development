#!/usr/bin/env Rscript
# Concatenate per-chromosome GENESIS association output into one sumstats file
# per phenotype, apply MAF/MAC filters, and emit the v1-style summary row
# (lambda_GC, hit counts, N) so v1 and v2 GWAS are directly comparable.
#
#   Rscript 05_collect_assoc.R --assoc-dir DIR --pheno NAME --out-dir DIR \
#       [--maf 0.01] [--mac 20]
#
# Output:
#   <out-dir>/<pheno>.sumstats.tsv.gz   SNP CHR POS A1 A2 AF1 N BETA SE P
#   <out-dir>/gwas_summary.tsv          appended/updated row for this phenotype
#
# BETA here is the score-test effect estimate (Est), in phenotype-SD units per
# allele; P is the score p-value.  lambda_GC from the median score chi-square.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--assoc-dir", type = "character"),
  make_option("--pheno",     type = "character"),
  make_option("--out-dir",   type = "character"),
  make_option("--maf",       type = "double",  default = 0.01),
  make_option("--mac",       type = "double",  default = 20)
)))
for (r in c("assoc-dir", "pheno", "out-dir")) if (is.null(opt[[r]])) stop("--", r, " is required")
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$`assoc-dir`,
                    pattern = sprintf("^%s_chr[0-9]+\\.tsv\\.gz$", opt$pheno),
                    full.names = TRUE)
if (!length(files)) stop("no ", opt$pheno, "_chr*.tsv.gz under ", opt$`assoc-dir`)
chrs <- as.integer(sub(".*_chr([0-9]+)\\.tsv\\.gz$", "\\1", files))
cat(sprintf("collecting %d chromosome files (chr %s)\n",
            length(files), paste(sort(chrs), collapse = ",")))

# fread() on a .gz path needs the R.utils package; piping through gzip -dc does
# not, and gzip is present everywhere this runs.  Avoiding the dependency keeps
# the cluster environment one package smaller.
read_gz <- function(f) fread(cmd = sprintf("gzip -dc %s", shQuote(f)))
dt <- rbindlist(lapply(files, read_gz), fill = TRUE)

# GENESIS 'freq' is the frequency of the FIRST/counted allele (A1 here).
dt[, maf := pmin(freq, 1 - freq)]
n_total <- nrow(dt)
dt <- dt[maf >= opt$maf & MAC >= opt$mac]
cat(sprintf("variants: %d total, %d after MAF>=%g & MAC>=%g\n",
            n_total, nrow(dt), opt$maf, opt$mac))

out <- dt[, .(SNP, CHR = chr, POS = pos, A1, A2, AF1 = freq, N = n.obs,
              BETA = Est, SE = Est.SE, P = Score.pval)]
setorder(out, CHR, POS)
ss_path <- file.path(opt$`out-dir`, sprintf("%s.sumstats.tsv.gz", opt$pheno))
# data.table's own gzip writer (no R.utils needed for WRITING).
fwrite(out, ss_path, sep = "\t", compress = "gzip")

# lambda_GC from the score statistic (chi2, 1 df).
chi2 <- (out$BETA / out$SE)^2
lambda <- median(chi2, na.rm = TRUE) / qchisq(0.5, df = 1)

row <- data.table(
  phenotype = opt$pheno,
  n_snps = nrow(out),
  n_mean = as.integer(mean(out$N)),
  lambda_gc = round(lambda, 4),
  n_p5e8 = sum(out$P < 5e-8, na.rm = TRUE),
  n_p1e5 = sum(out$P < 1e-5, na.rm = TRUE),
  min_p = suppressWarnings(min(out$P, na.rm = TRUE)),
  n_chr = length(files)
)

# Update-in-place: re-running one phenotype must not duplicate its row.
sum_path <- file.path(opt$`out-dir`, "gwas_summary.tsv")
if (file.exists(sum_path)) {
  old <- fread(sum_path)
  old <- old[phenotype != opt$pheno]
  row <- rbind(old, row, fill = TRUE)
}
setorder(row, phenotype)
fwrite(row, sum_path, sep = "\t")

cat(sprintf("WROTE  %s\nWROTE  %s\n", ss_path, sum_path))
print(row[phenotype == opt$pheno], row.names = FALSE)
if (length(files) < 22) {
  cat(sprintf("WARN: only %d chromosomes present -- a partial-genome GWAS.\n", length(files)))
  cat("      (Expected on the fixture; a FAILURE on real data.)\n")
}
