#!/usr/bin/env Rscript
# Single-variant association for one chromosome against a fitted null model.
#
#   Rscript 04_assoc.R --gds chrN.gds --null-model F.rds --out chrN.tsv.gz \
#       [--block 5000]
#
# assocTestSingle scores every variant against the SAME null model, so the
# per-chromosome scans are embarrassingly parallel and their concatenation is
# one GWAS.  Output columns follow GENESIS (variant.id chr pos allele freq MAC
# n.obs Est Est.SE Score Score.SE Score.Stat Score.pval); the collector renames
# to the v1 fastGWA-like layout so downstream tooling (Manhattan, MAGMA prep)
# reads either pipeline's output.
suppressPackageStartupMessages({
  library(optparse)
  library(GENESIS)
  library(GWASTools)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gds",        type = "character"),
  make_option("--null-model", type = "character"),
  make_option("--out",        type = "character"),
  make_option("--block",      type = "integer", default = 5000)
)))
for (r in c("gds", "null-model", "out")) if (is.null(opt[[r]])) stop("--", r, " is required")

nullmod <- readRDS(opt$`null-model`)

geno <- GdsGenotypeReader(opt$gds)
genoData <- GenotypeData(geno)
iter <- GenotypeBlockIterator(genoData, snpBlock = opt$block)

assoc <- assocTestSingle(
  iter, null.model = nullmod,
  test = "Score",
  BPPARAM = BiocParallel::SerialParam()
)
close(genoData)

# Attach rsIDs: GWASTools carries snpID as integer index; map back to the
# .gds snp.rs.id if present so downstream tools can match on rsID.
g <- SNPRelate::snpgdsOpen(opt$gds)
rsid <- tryCatch(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.rs.id")),
                 error = function(e) NULL)
snpid <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.id"))
alleles <- tryCatch(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.allele")),
                    error = function(e) NULL)
SNPRelate::snpgdsClose(g)

dt <- as.data.table(assoc)
if (!is.null(rsid)) {
  map <- data.table(variant.id = snpid, SNP = rsid)
  dt <- merge(dt, map, by = "variant.id", all.x = TRUE, sort = FALSE)
} else {
  dt[, SNP := as.character(variant.id)]
}
# A1 = the allele GENESIS scored (effect allele = GDS first allele);
# snp.allele is "A/B" with A the counted allele.
if (!is.null(alleles)) {
  al <- tstrsplit(alleles, "/", fixed = TRUE)
  amap <- data.table(variant.id = snpid, A1 = al[[1]], A2 = al[[2]])
  dt <- merge(dt, amap, by = "variant.id", all.x = TRUE, sort = FALSE)
}

setorder(dt, chr, pos)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(dt, opt$out, sep = "\t", compress = if (grepl("\\.gz$", opt$out)) "gzip" else "none")
cat(sprintf("WROTE  %s  (%d variants, n.obs median %d)\n",
            opt$out, nrow(dt), as.integer(median(dt$n.obs, na.rm = TRUE))))
