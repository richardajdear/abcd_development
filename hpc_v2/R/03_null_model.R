#!/usr/bin/env Rscript
# GENESIS null model: phenotype ~ covariates + PC-AiR PCs + (kinship random
# effect).  One phenotype per invocation; the sbatch array loops the manifest.
#
# This replaces BOTH v1 exclusions at once: relatives stay in the sample
# (the sparse PC-Relate kinship models them) and the multi-ancestry sample is
# analysed pooled (PC-AiR PCs are valid ancestry axes in the presence of
# relatives, which GCTA's PCA was not).
#
#   Rscript 03_null_model.R --pheno-file F --pheno NAME --covar-quant F \
#       --covar-cat F --pcs F --kinship F --out F.rds \
#       [--n-pcs 10] [--strata F] [--inv-norm]
#
# ID JOIN: on the normalised 8-char NDAR token (align_ids.py convention),
# because the .fam/GDS IDs and the export IDs spell subjects differently.
# The output null model carries the GDS spelling, which is what
# assocTestSingle will see.
suppressPackageStartupMessages({
  library(optparse)
  library(GENESIS)
  library(Matrix)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--pheno-file",  type = "character"),
  make_option("--pheno",       type = "character"),
  make_option("--covar-quant", type = "character"),
  make_option("--covar-cat",   type = "character"),
  make_option("--pcs",         type = "character"),
  make_option("--kinship",     type = "character"),
  make_option("--out",         type = "character"),
  make_option("--n-pcs",       type = "integer", default = 10),
  make_option("--strata",      type = "character", default = NULL,
              help = "FID IID stratum file; enables heteroscedastic residuals (group.var)"),
  make_option("--inv-norm",    action = "store_true", default = FALSE,
              help = "rank-inverse-normalise the phenotype before fitting")
)))
for (r in c("pheno-file","pheno","covar-quant","covar-cat","pcs","kinship","out"))
  if (is.null(opt[[r]])) stop("--", r, " is required")

norm_id <- function(x) toupper(gsub("_", "", sub("^sub-", "", as.character(x))))

ph   <- fread(opt$`pheno-file`)
if (!opt$pheno %in% names(ph)) stop("phenotype '", opt$pheno, "' not in ", opt$`pheno-file`)
qcov <- fread(opt$`covar-quant`)
ccov <- fread(opt$`covar-cat`)
pcs  <- fread(opt$pcs)
kin  <- readRDS(opt$kinship)

# Assemble the analysis table on the token; keep the KINSHIP/GDS spelling as
# the sample.id GENESIS matches on.
kin_ids <- rownames(kin)
kd <- data.table(sample.id = kin_ids, key_tok = norm_id(kin_ids))

for (dt in list(ph, qcov, ccov, pcs)) dt[, key_tok := norm_id(IID)]
d <- merge(kd, ph, by = "key_tok")
# Quantitative covariates: drop the export's release PCs (PC1..PC10 columns);
# the PC-AiR PCs replace them.  Everything else (age, n_visits, ...) stays.
qkeep <- setdiff(names(qcov), c(grep("^PC[0-9]+$", names(qcov), value = TRUE)))
d <- merge(d, qcov[, ..qkeep], by = "key_tok", suffixes = c("", ".q"))
d <- merge(d, ccov[, .(key_tok, sex, site)], by = "key_tok", suffixes = c("", ".c"))
pc_cols <- paste0("PC", seq_len(opt$`n-pcs`))
missing_pcs <- setdiff(pc_cols, names(pcs))
if (length(missing_pcs)) stop("PCs file lacks ", paste(missing_pcs, collapse = ","))
d <- merge(d, pcs[, c("key_tok", pc_cols), with = FALSE], by = "key_tok")

# One row per subject; a token collision would silently duplicate rows.
if (anyDuplicated(d$sample.id)) stop("duplicate sample.id after joins -- token collision?")

group_var <- NULL
if (!is.null(opt$strata)) {
  st <- fread(opt$strata)
  stopifnot(ncol(st) >= 3)
  setnames(st, 1:3, c("FID", "IID", "stratum"))
  st[, key_tok := norm_id(IID)]
  d <- merge(d, st[, .(key_tok, stratum)], by = "key_tok")
  d[, stratum := as.factor(stratum)]
  group_var <- "stratum"
}

y <- d[[opt$pheno]]
if (opt$`inv-norm`) {
  # Rank-based inverse normal (Blom).  For BLUP phenotypes already ~normal this
  # is a no-op in practice; the flag exists for sensitivity runs.
  y <- qnorm((rank(y, na.last = "keep") - 3/8) / (sum(!is.na(y)) + 1/4))
}
d[, y__ := y]

covars <- c("sex", "site", intersect(c("baseline_age", "age_c", "n_visits"), names(d)), pc_cols)
d[, sex := as.factor(sex)]
d[, site := as.factor(site)]

df <- as.data.frame(d)
rownames(df) <- df$sample.id
cat(sprintf("analysis sample: %d subjects (%d in kinship, %d phenotyped)\n",
            nrow(df), nrow(kin), nrow(ph)))
if (nrow(df) == 0) stop("empty join: phenotype IDs and kinship IDs do not intersect")

# Subset the kinship to the analysed subjects, in df order.
kin_sub <- kin[df$sample.id, df$sample.id]

annot <- Biobase::AnnotatedDataFrame(df)
nullmod <- fitNullModel(
  annot,
  outcome = "y__",
  covars = covars,
  cov.mat = kin_sub,
  group.var = group_var,
  family = "gaussian"
)
saveRDS(nullmod, opt$out)

# Variance components + a fit summary a human can read without loading the rds.
vc <- nullmod$varComp
h2_kin <- if (!is.null(vc) && length(vc) >= 2) vc[["V_A"]] / sum(unlist(vc)) else NA
smry <- data.frame(
  phenotype = opt$pheno,
  n = nrow(df),
  n_covars = length(covars),
  group_var = ifelse(is.null(group_var), "none", group_var),
  inv_norm = opt$`inv-norm`,
  varcomp_kin = if (!is.null(vc)) vc[[grep("^V_A|A$", names(vc))[1]]] else NA,
  varcomp_resid = if (!is.null(vc)) vc[[grep("V_E|resid", names(vc))[1]]] else NA,
  prop_kin = h2_kin
)
out_tsv <- sub("\\.rds$", "_summary.tsv", opt$out)
write.table(smry, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("WROTE  %s\nWROTE  %s\n", opt$out, out_tsv))
print(smry, row.names = FALSE)
