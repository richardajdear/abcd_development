#!/usr/bin/env Rscript
# Local genetic correlation between our cortical phenotypes and SCZ / MDD.
#
# The targeted test is the point: LD blocks containing SCZ locus-pool genes are
# selected from MAGMA's result, which was computed WITHOUT reference to any
# local rg, so the selection is not circular.
suppressPackageStartupMessages({
  library(optparse); library(data.table)
  .libPaths(c("/rds/user/rajd2/hpc-work/lava_local", .libPaths()))
  library(R.utils)   # LAVA reads .sumstats.gz via fread, which needs this
  library(LAVA)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--ref",       type = "character"),
  make_option("--loci",      type = "character"),
  make_option("--gene-loc",  type = "character"),
  make_option("--geneset",   type = "character"),
  make_option("--set-name",  type = "character", default = "SCZ_locus_pool"),
  make_option("--abcd-dir",  type = "character"),
  make_option("--scz",       type = "character"),
  make_option("--mdd",       type = "character"),
  make_option("--out-dir",   type = "character"),
  make_option("--n-abcd",    type = "integer", default = 4017)
)))
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

loci <- fread(opt$loci)
setnames(loci, toupper(names(loci)))
cat(sprintf("loci: %d\n", nrow(loci)))

# --- which blocks carry SCZ locus-pool genes ---------------------------------
gl  <- fread(opt$`gene-loc`, header = FALSE,
             col.names = c("GENE","CHR","START","STOP","STRAND","SYMBOL"))
gs  <- fread(opt$geneset, header = FALSE, col.names = c("GENE","SET"))
pool <- gl[GENE %in% gs[SET == opt$`set-name`, GENE]]
cat(sprintf("%s: %d genes, %d with locations\n", opt$`set-name`,
            gs[SET == opt$`set-name`, .N], nrow(pool)))
pool[, CHR := as.character(CHR)]
loci[, CHR := as.character(CHR)]
hit <- loci[pool, on = .(CHR), allow.cartesian = TRUE][
             START <= i.STOP & STOP >= i.START, unique(LOC)]
cat(sprintf("blocks containing at least one pool gene: %d of %d\n",
            length(hit), nrow(loci)))
fwrite(data.table(LOC = hit), file.path(opt$`out-dir`, "targeted_blocks.txt"))

# --- input files, in LAVA's expected shape -----------------------------------
phenos <- c("baseline_thickness","global_slope","slope_PC1","slope_PC2","slope_PC3")
files  <- c(setNames(file.path(opt$`abcd-dir`, paste0(phenos, ".sumstats.gz")), phenos),
            SCZ = opt$scz, MDD = opt$mdd)
files  <- files[file.exists(files)]
cat("inputs:\n"); for (n in names(files)) cat(sprintf("  %-20s %s\n", n, files[[n]]))

# LAVA's info file is exactly four whitespace-separated columns and will not
# accept NA in cases/controls -- read.table counts fields and errors out.
# Quantitative traits take cases/controls = 0; the disorders take their real
# counts, which LAVA uses only to convert h2 to the liability scale.
ncas <- c(SCZ = 55085, MDD = 412305); ncon <- c(SCZ = 78957, MDD = 1588397)
info <- data.table(
  phenotype = names(files),
  cases     = ifelse(names(files) %in% names(ncas), ncas[names(files)], 0),
  controls  = ifelse(names(files) %in% names(ncon), ncon[names(files)], 0),
  filename  = unname(files))
info[is.na(cases), cases := 0][is.na(controls), controls := 0]
info_f <- file.path(opt$`out-dir`, "input.info.txt")
write.table(info, info_f, sep = "\t", quote = FALSE, row.names = FALSE)

# Sample overlap between our five phenotypes is complete (same subjects), and
# between ours and PGC is nil.  LAVA takes that as a genetic-covariance
# intercept matrix; with no overlap against the disorders the off-diagonal is 0.
inp <- process.input(input.info.file = info_f, sample.overlap.file = NULL,
                     ref.prefix = opt$ref, phenos = names(files))

run_pairs <- function(block_ids, tag) {
  res <- list(); uni <- list()
  for (i in seq_along(block_ids)) {
    lid <- block_ids[i]
    lrow <- loci[LOC == lid]
    if (!nrow(lrow)) next
    lo <- try(process.locus(as.data.frame(lrow)[1, ], inp), silent = TRUE)
    if (inherits(lo, "try-error") || is.null(lo)) next
    u <- try(run.univ(lo), silent = TRUE)
    if (inherits(u, "try-error") || is.null(u) || !nrow(u)) next
    u$LOC <- lid; uni[[length(uni)+1]] <- u
    # Bivariate only where BOTH traits clear the univariate filter, LAVA's rule.
    for (d in intersect(c("SCZ","MDD"), u$phen[u$p < 0.05])) {
      for (p in intersect(phenos, u$phen[u$p < 0.05])) {
        b <- try(run.bivar(lo, phenos = c(p, d)), silent = TRUE)
        if (inherits(b, "try-error") || is.null(b) || !nrow(b)) next
        b$LOC <- lid; res[[length(res)+1]] <- b
      }
    }
    if (i %% 100 == 0) cat(sprintf("  [%s] %d/%d blocks\n", tag, i, length(block_ids)))
  }
  if (length(uni)) fwrite(rbindlist(uni, fill = TRUE),
                          file.path(opt$`out-dir`, paste0("univ_", tag, ".tsv")), sep = "\t")
  if (length(res)) {
    out <- rbindlist(res, fill = TRUE)
    fwrite(out, file.path(opt$`out-dir`, paste0("bivar_", tag, ".tsv")), sep = "\t")
    cat(sprintf("[%s] %d bivariate tests written\n", tag, nrow(out)))
    print(head(out[order(p)], 12))
  } else {
    cat(sprintf("[%s] no bivariate tests -- no block cleared the univariate filter in both traits\n", tag))
  }
  invisible(NULL)
}

cat("\n=== TARGETED: blocks with SCZ locus-pool genes ===\n")
run_pairs(hit, "targeted")
cat("\n=== GENOME-WIDE ===\n")
run_pairs(loci$LOC, "genomewide")
cat("done\n")
