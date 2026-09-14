#!/usr/bin/env Rscript
# KING-robust -> PC-AiR -> PC-Relate on the array GDS.
#
# WHY THIS SEQUENCE (the fix for hpc/README_HPC.md §8.6)
# ------------------------------------------------------
# A pooled GCTA GRM scores same-ancestry pairs against pooled allele
# frequencies, so its off-diagonals confound ancestry with kinship;
# `--grm-cutoff 0.05` on it returned a 94%-European "unrelated" set while
# reporting nothing wrong.  The estimators here separate the two axes:
#
#   KING-robust   -- kinship robust to population structure (but biased by
#                    admixture); used only to SEED the partition.
#   PC-AiR        -- PCs computed on an unrelated subset chosen with KING,
#                    then related subjects projected in.  Ancestry axes
#                    uncontaminated by family structure.
#   PC-Relate     -- kinship CONDITIONAL on the PC-AiR PCs: recent-kinship
#                    estimates with the ancestry component removed.
#
# The final sparse kinship matrix is the mixed-model random effect for the
# GENESIS null model; the PC-AiR PCs are its fixed-effect covariates; the
# PC-AiR unrelated partition is reported for comparison with the ABCD-shipped
# genesis/unrelateds_individuals.txt.
#
#   Rscript 02_kinship.R --gds FILE.gds --out-dir DIR \
#       [--maf 0.05] [--ld-r2 0.1] [--kin-thresh 0.02209709] \
#       [--sparse-thresh 0.04419417] [--max-density 0.10] \
#       [--n-pcs 32] [--threads 8] [--seed 20260823]
#
# TWO thresholds, doing different jobs -- do not conflate them:
#   --kin-thresh    (default 2^(-11/2) ~ 0.02210) is RAW kinship, and sets
#                   PC-AiR's unrelated/related partition (3rd-degree boundary).
#   --sparse-thresh (default 2^(-9/2)  ~ 0.04419) is compared against
#                   2*kinship, because pcrelateToMatrix applies it AFTER
#                   scaleKin=2.  It therefore corresponds to raw kinship
#                   ~0.02210 -- numerically the same boundary, reached through
#                   a different scaling.  See config.sh's SPARSE_KIN_THRESH.
#
# Outputs (all under --out-dir):
#   pruned_snps.rds       LD-pruned SNP set used by every estimator
#   king_matrix.rds       KING-robust kinship (dense, snpgdsIBDKING)
#   pcair.rds             pcair object: PCs ($vectors), unrels, rels
#   pcair_pcs.tsv         FID IID PC1..PCk   (FID = IID here; join on IID)
#   pcair_unrelated.txt   IIDs of the PC-AiR unrelated partition
#   pcrelate.rds          pcrelate object (kinBtwn/kinSelf tables)
#   kinship_sparse.rds    block-diagonal symmetric Matrix from
#                         pcrelateToMatrix -- the null model's covariance
#                         input.  NOTE it clusters by transitive closure at
#                         --sparse-thresh rather than zeroing entries
#                         element-wise, so it can come out DENSE; the summary
#                         reports sparse_density / sparse_largest_block and a
#                         warning fires above --max-density.
#   kinship_summary.tsv   counts + degree tallies, the "did it work" file
suppressPackageStartupMessages({
  library(optparse)
  library(SNPRelate)
  library(GWASTools)
  library(GENESIS)
  library(Matrix)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gds",        type = "character"),
  make_option("--out-dir",    type = "character"),
  make_option("--maf",        type = "double",  default = 0.05),
  make_option("--ld-r2",      type = "double",  default = 0.1),
  make_option("--kin-thresh", type = "double",  default = 2^(-11/2),
              help = "PC-AiR unrelated/related partition threshold (raw kinship)"),
  make_option("--sparse-thresh", type = "double", default = 2^(-9/2),
              help = paste("pcrelateToMatrix clustering threshold, compared against",
                           "2*kinship (see config.sh); NOT element-wise zeroing")),
  make_option("--max-density", type = "double", default = 0.10,
              help = "warn if the 'sparse' kinship matrix exceeds this density"),
  make_option("--n-pcs",      type = "integer", default = 32),
  make_option("--threads",    type = "integer", default = 4),
  make_option("--seed",       type = "integer", default = 20260823)
)))
for (r in c("gds", "out-dir")) if (is.null(opt[[r]])) stop("--", r, " is required")
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)
outp <- function(...) file.path(opt$`out-dir`, ...)

gds <- snpgdsOpen(opt$gds)
sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))
cat(sprintf("samples: %d\n", length(sample.id)))

# --- 1. LD pruning ------------------------------------------------------------
# One pruned set, reused by KING, PC-AiR and PC-Relate, so their inputs are
# identical and differences between estimators are method, not SNP set.
snpset <- snpgdsLDpruning(
  gds, maf = opt$maf, missing.rate = 0.02,
  method = "corr", slide.max.bp = 500e3, ld.threshold = sqrt(opt$`ld-r2`),
  num.thread = opt$threads, verbose = FALSE
)
pruned <- unlist(snpset, use.names = FALSE)
saveRDS(pruned, outp("pruned_snps.rds"))
cat(sprintf("LD-pruned SNPs: %d (maf>=%g, r2<=%g)\n", length(pruned), opt$maf, opt$`ld-r2`))
if (length(pruned) < 5000) {
  cat("WARN: fewer than 5,000 pruned SNPs; kinship estimates will be noisy.\n",
      "     Fine on a synthetic fixture, alarming on real data.\n")
}

# --- 2. KING-robust -----------------------------------------------------------
king <- snpgdsIBDKING(gds, snp.id = pruned, num.thread = opt$threads, verbose = FALSE)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
saveRDS(kingMat, outp("king_matrix.rds"))
snpgdsClose(gds)

# --- 3. PC-AiR ----------------------------------------------------------------
# kinobj = divobj = KING: kinship seeds the unrelated partition, and KING's
# negative values between diverged pairs serve as the divergence measure.
geno <- GdsGenotypeReader(opt$gds)
genoData <- GenotypeData(geno)

pca <- pcair(
  genoData,
  kinobj = kingMat, divobj = kingMat,
  kin.thresh = opt$`kin-thresh`, div.thresh = -opt$`kin-thresh`,
  snp.include = pruned,
  num.cores = opt$threads
)
saveRDS(pca, outp("pcair.rds"))

npc <- min(opt$`n-pcs`, ncol(pca$vectors))
pcs <- data.frame(
  FID = rownames(pca$vectors), IID = rownames(pca$vectors),
  pca$vectors[, seq_len(npc)], check.names = FALSE
)
colnames(pcs)[-(1:2)] <- paste0("PC", seq_len(npc))
write.table(pcs, outp("pcair_pcs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(pca$unrels, outp("pcair_unrelated.txt"))
cat(sprintf("PC-AiR: %d unrelated / %d related, %d PCs written\n",
            length(pca$unrels), length(pca$rels), npc))

# --- 4. PC-Relate ---------------------------------------------------------------
# Kinship conditional on the ancestry PCs.  training.set = the PC-AiR unrelated
# partition, per the GENESIS vignette: ISAF betas estimated on unrelateds.
# Iterating pcair<->pcrelate once more changes ABCD estimates negligibly and
# costs another full pass; one round is the published default.
genoIt <- GenotypeBlockIterator(genoData, snpBlock = 10000, snpInclude = pruned)
pcrel <- pcrelate(
  genoIt,
  pcs = pca$vectors[, 1:2, drop = FALSE],
  training.set = pca$unrels,
  BPPARAM = BiocParallel::SerialParam()
)
saveRDS(pcrel, outp("pcrelate.rds"))

# Sparse kinship for the null model.
#
# CAUTION: `thresh` CLUSTERS samples by transitive closure -- it does not zero
# individual entries.  Any pair above the threshold joins a cluster, every
# within-cluster pair is then retained, and only BETWEEN-cluster entries become
# zero.  So a chain of noise-level pairs can merge the cohort into one block and
# produce a DENSE matrix, which defeats the sparse design silently.  Hence the
# density report and warning below -- see config.sh's SPARSE_KIN_THRESH note.
kinM <- pcrelateToMatrix(pcrel, thresh = opt$`sparse-thresh`, scaleKin = 2)
kinM <- Matrix(kinM, sparse = TRUE)
saveRDS(kinM, outp("kinship_sparse.rds"))

# Density and block structure of what we just wrote.
n_k <- nrow(kinM)
density <- Matrix::nnzero(kinM) / (as.numeric(n_k) * n_k)
# Largest block = largest connected component of the supra-threshold graph.
hi <- pcrel$kinBtwn[2 * pcrel$kinBtwn$kin > opt$`sparse-thresh`, c("ID1", "ID2")]
largest_block <- 1L
if (nrow(hi)) {
  adj <- split(c(as.character(hi$ID2), as.character(hi$ID1)),
               c(as.character(hi$ID1), as.character(hi$ID2)))
  seen <- new.env(hash = TRUE, parent = emptyenv())
  for (s in names(adj)) {
    if (!is.null(seen[[s]])) next
    comp <- s; queue <- s
    while (length(queue)) {
      v <- queue[1]; queue <- queue[-1]
      nb <- setdiff(adj[[v]], comp)
      comp <- c(comp, nb); queue <- c(queue, nb)
    }
    for (m in comp) assign(m, TRUE, envir = seen)
    largest_block <- max(largest_block, length(comp))
  }
}
cat(sprintf("\nsparse kinship: %d x %d, density %.4f, largest block %d (%.1f%% of sample)\n",
            n_k, n_k, density, largest_block, 100 * largest_block / n_k))
if (density > opt$`max-density`) {
  cat("\n*** WARNING: the 'sparse' kinship matrix is not sparse. ***\n")
  cat(sprintf("    density %.3f exceeds --max-density %.3f, and %d of %d subjects\n",
              density, opt$`max-density`, largest_block, n_k))
  cat("    sit in ONE block.  pcrelateToMatrix clusters by TRANSITIVE CLOSURE, so\n")
  cat("    weak noise-level pairs can chain the whole cohort together.  Consequences:\n")
  cat("    fitNullModel will hold a dense matrix (~1.1 GB at n=11,670) and may OOM,\n")
  cat("    and the association scan loses the speed the sparse design buys.\n")
  cat(sprintf("    FIX: raise SPARSE_KIN_THRESH (currently %.5f, i.e. raw kinship > %.5f);\n",
              opt$`sparse-thresh`, opt$`sparse-thresh` / 2))
  cat("    0.0884 keeps only 2nd-degree-and-closer pairs.  Check n_pruned_snps first:\n")
  cat("    too few pruned SNPs is the usual cause of noisy kinship estimates.\n\n")
}

# --- 5. Summary -----------------------------------------------------------------
kb <- pcrel$kinBtwn
deg <- function(lo, hi) sum(kb$kin >= lo & kb$kin < hi, na.rm = TRUE)
# Counts are written as plain INTEGERS and the two fractional rows with fixed
# decimals.  A single numeric `value` column would coerce everything to double
# and print counts as 1.500000e+03 -- which reads fine to a human but breaks
# every downstream `[[ $n -gt 0 ]]`, since bash cannot compare scientific
# notation.  Formatting per row is the fix.
int_row <- function(x) formatC(as.numeric(x), format = "d")
summary_df <- data.frame(
  quantity = c("n_samples", "n_pruned_snps", "n_pcair_unrelated", "n_pcair_related",
               "pairs_kin_ge_0.354 (dup/MZ)", "pairs_kin_0.177_0.354 (1st deg)",
               "pairs_kin_0.088_0.177 (2nd deg)", "pairs_kin_0.044_0.088 (3rd deg)",
               "sparse_nonzero_offdiag", "sparse_density", "sparse_largest_block",
               "sparse_thresh_used"),
  value = c(int_row(length(sample.id)), int_row(length(pruned)),
            int_row(length(pca$unrels)), int_row(length(pca$rels)),
            int_row(deg(0.3536, Inf)), int_row(deg(0.1768, 0.3536)),
            int_row(deg(0.0884, 0.1768)), int_row(deg(0.0442, 0.0884)),
            int_row((Matrix::nnzero(kinM) - nrow(kinM)) / 2),
            formatC(density, format = "f", digits = 5),
            int_row(largest_block),
            formatC(opt$`sparse-thresh`, format = "f", digits = 8)),
  stringsAsFactors = FALSE
)
write.table(summary_df, outp("kinship_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n"); print(summary_df, row.names = FALSE)
cat(sprintf("\nWROTE  %s\n", outp("kinship_summary.tsv")))
