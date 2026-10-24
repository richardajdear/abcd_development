# SCTWAS stage 2: signed gene-level association for one weights tag x all traits x both models.
# usage: Rscript 03_assoc.R <tag>
#   z_g = w' z / sqrt(w' R w); w = model weights (standardised scale, effect allele = panel A1),
#   z = GWAS BETA/SE aligned to panel A1 (stage 0), R = 1000G EUR LD correlation of the model SNPs; p = 2*pnorm(-|z_g|).
#   Raw z only; no trait sign conventions are applied here.
# Inputs : weights/<tag>.{susie,top1}.tsv.gz, gwas/<trait>.z.tsv.gz, ld/panel_keep
# Outputs: assoc/<trait>.<tag>.<model>.tsv.gz  (gene n_model_snps z p n_cs s_rss w_flag)
#   w_flag (from 02b_collect.R): implausible SuSiE weights, max|w| > 3 |z_top|/sqrt(n); flagged, not dropped. Always FALSE in the
#   top1 model (a single marginal weight cannot be implausible).
args <- commandArgs(trailingOnly = TRUE); tag <- args[1]
W <- "/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas"
suppressPackageStartupMessages(library(data.table))
source("/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/sctwas/plink_bed.R")
TRAITS <- c("ABCD_thin_EUR", "ABCD_base_EUR", "SCZ25_EUR", "MDD_EUR", "PGC3_EUR")
dir.create(file.path(W, "assoc"), showWarnings = FALSE)

gw <- lapply(setNames(TRAITS, TRAITS), function(t) { x <- fread(cmd = paste("zcat", shQuote(file.path(W, "gwas", paste0(t, ".z.tsv.gz"))))); setNames(x$z, x$rsid) })
gs <- fread(cmd = paste("zcat", shQuote(file.path(W, "weights", paste0(tag, ".genes.tsv.gz")))))
snp <- fread(file.path(W, "ld", "snps_keep.tsv"), colClasses = list(character = c("rsid", "chr")))

for (model in c("susie", "top1")) {
  wt <- fread(cmd = paste("zcat", shQuote(file.path(W, "weights", sprintf("%s.%s.tsv.gz", tag, model)))), colClasses = list(character = "rsid"))
  wt[, chr := snp$chr[match(rsid, snp$rsid)]]
  stopifnot(!anyNA(wt$chr))
  Z <- sapply(TRAITS, function(t) unname(gw[[t]][match(wt$rsid, names(gw[[t]]))]))   # one lookup, rows = wt rows
  rows <- split(seq_len(nrow(wt)), wt$gene)
  res <- list()
  for (ch in unique(wt$chr)) {
    B <- open_bed(file.path(W, "ld", "panel_keep"), ch)
    for (g in unique(wt[chr == ch, gene])) {
      ix <- rows[[g]]; d <- wt[ix]
      R <- if (nrow(d) > 1) ld_cor(B, d$rsid) else matrix(1)
      for (t in TRAITS) {
        z <- Z[ix, t]; ok <- !is.na(z)
        if (!any(ok)) next
        w <- d$weight[ok]; Rk <- R[ok, ok, drop = FALSE]
        v <- drop(t(w) %*% Rk %*% w)
        if (!is.finite(v) || v <= 0) next
        zg <- sum(w * z[ok]) / sqrt(v)
        res[[length(res) + 1]] <- list(trait = t, gene = g, n_model_snps = sum(ok), z = zg)
      }
    }
  }
  res <- rbindlist(res)
  res[, p := 2 * pnorm(-abs(z))]
  res <- merge(res, gs[, .(gene, n_cs, s_rss, w_flag)], by = "gene", all.x = TRUE)
  for (t in TRAITS) {
    o <- res[trait == t, .(gene, n_model_snps, z, p, n_cs, s_rss, w_flag)][order(p)]
    if (model == "top1") o[, w_flag := FALSE]
    fwrite(o, file.path(W, "assoc", sprintf("%s.%s.%s.tsv.gz", t, tag, model)), sep = "\t")
    k <- !o$w_flag
    cat(sprintf("%s %s %s: %d genes, lambda %.3f, p<2.5e-6: %d | unflagged %d genes, lambda %.3f, p<2.5e-6: %d\n", t, tag, model,
                nrow(o), median(o$z^2) / 0.455, sum(o$p < 2.5e-6), sum(k), median(o$z[k]^2) / 0.455, sum(o$p[k] < 2.5e-6)))
  }
}
