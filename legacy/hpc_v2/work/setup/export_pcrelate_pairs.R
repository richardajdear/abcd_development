#!/usr/bin/env Rscript
# Export the PC-Relate related pairs that define the Zaitlen second component.
#
# WHY THIS EXISTS: step 05 as shipped builds its bK GRM with GCTA's
# --make-bK on the POOLED multi-ancestry GRM, i.e. it decides who is a close
# relative from that GRM's off-diagonals.  Those off-diagonals confound kinship
# with ancestry -- the single most important finding in v1
# (hpc/README_HPC.md 8.6) and the reason hpc_v2 exists.  Measured on this run:
# --make-bK 0.05 kept 3,837,371 pairs, 5.6% of all pairs, against ~1,339 true
# sibling pairs in the phenotyped sample.  So G2 was an ancestry component, not
# a pedigree one, and it swallowed the additive variance.
#
# PC-Relate estimates kinship CONDITIONAL on the ancestry PCs, which is exactly
# the separation that fixes this.  Output: ID1 ID2 kin for pairs above the
# threshold, in GCTA .grm.id spelling.
suppressPackageStartupMessages({library(GENESIS); library(data.table); library(optparse)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--pcrelate", type = "character"),
  make_option("--out",      type = "character"),
  make_option("--kin-thresh", type = "double", default = 0.025,
              help = paste("raw PC-Relate kinship floor.  GCTA GRM values are",
                           "~2*kinship, so the pipeline's BK_THRESH on GRM",
                           "values corresponds to BK_THRESH/2 here."))
)))
pcrel <- readRDS(opt$pcrelate)
kb <- as.data.table(pcrel$kinBtwn)
cat(sprintf("kinBtwn pairs: %d\n", nrow(kb)))
keep <- kb[kin >= opt$`kin-thresh`, .(ID1, ID2, kin)]
cat(sprintf("pairs with kinship >= %.4f: %d\n", opt$`kin-thresh`, nrow(keep)))
cat(sprintf("  1st-degree (kin>=0.177): %d\n", nrow(kb[kin >= 0.177])))
cat(sprintf("  2nd-degree (0.0884-0.177): %d\n", nrow(kb[kin >= 0.0884 & kin < 0.177])))
cat(sprintf("  3rd-degree (0.0442-0.0884): %d\n", nrow(kb[kin >= 0.0442 & kin < 0.0884])))
fwrite(keep, opt$out, sep = "\t")
cat(sprintf("wrote %s\n", opt$out))
