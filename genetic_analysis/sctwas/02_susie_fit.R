# SCTWAS stage 1b: SuSiE-RSS expression weights for one cell type x chromosome (SCTWAS_INSTRUCTIONS.md, Stage 1).
# usage: Rscript 02_susie_fit.R <CT> <chr> <n> [<out tag, default CT>]
# Inputs : stage1/<CT>/zin_chr<chr>.tsv.gz (gene_version rsid z; z = fixed_z aligned to panel A1)
#          ld/panel_keep (1000G EUR restricted to the stage-0 SNP set; bim col5 = A1)
# Per gene: R = panel LD of the gene's cis SNPs (same order as z);
#           fit = susie_rss(z, R, n, L = 5, estimate_residual_variance = FALSE), other arguments default;
#           w = coef(fit)[-1] (posterior means, standardised scale); s = estimate_s_rss(z, R, n).
# Outputs: stage1/<tag>/chr<chr>.{genes,susie,top1}.tsv.gz
#   genes : one row per eGene (diagnostics, including failures)
#   susie : SNPs in any 95% credible set, genes with >= 1 CS (primary model)
#   top1  : lead |z| SNP per eGene, w = z_top / sqrt(n) (secondary SMR-like model)
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]; chr <- args[2]; n <- as.numeric(args[3]); tag <- if (length(args) > 3) args[4] else ct
W <- "/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas"
.libPaths(c(file.path(W, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(susieR) })
source("/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/sctwas/plink_bed.R")

out <- file.path(W, "stage1", tag); dir.create(out, showWarnings = FALSE, recursive = TRUE)
fin <- file.path(out, sprintf("chr%s.genes.tsv.gz", chr))
if (file.exists(fin)) { cat("exists:", fin, "\n"); quit(save = "no") }
zin <- fread(cmd = paste("zcat", shQuote(file.path(W, "stage1", ct, sprintf("zin_chr%s.tsv.gz", chr)))))
B <- open_bed(file.path(W, "ld", "panel_keep"), chr)
bim <- B$bim
zin <- zin[rsid %in% bim$rsid]
zin[, ord := B$pos[rsid]]
setorder(zin, gene_version, ord)
genes <- unique(zin$gene_version)
cat(sprintf("%s chr%s n=%g: %d genes, %d rows  %s\n", tag, chr, n, length(genes), nrow(zin), format(Sys.time())))

G <- vector("list", length(genes)); S <- vector("list", length(genes)); T1 <- vector("list", length(genes))
t0 <- Sys.time()
for (i in seq_along(genes)) {
  g <- genes[i]; d <- zin[gene_version == g]
  z <- d$z; rs <- d$rsid; p <- length(z)
  gid <- sub("\\..*$", "", g)
  j <- which.max(abs(z))
  T1[[i]] <- data.table(gene = gid, gene_version = g, rsid = rs[j], weight = z[j] / sqrt(n), n_snps_locus = p, top_z = z[j])
  row <- list(gene = gid, gene_version = g, n_snps_locus = p, n_cs = NA_integer_, cs_purity_min = NA_real_,
              cs_purity_mean = NA_real_, cs_size = NA_character_, s_rss = NA_real_, converged = NA, niter = NA_integer_,
              max_pip = NA_real_, top_abs_z = abs(z[j]), status = "ok", warn = "")
  wmsg <- character(0)
  res <- tryCatch(withCallingHandlers({
    R <- ld_cor(B, rs)
    row$s_rss <- estimate_s_rss(z, R, n)
    fit <- susie_rss(z = z, R = R, n = n, L = 5, estimate_residual_variance = FALSE)
    fit
  }, warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") }),
  error = function(e) { row$status <<- paste("error:", conditionMessage(e)); NULL })
  if (length(wmsg)) row$warn <- paste(unique(substr(wmsg, 1, 120)), collapse = " | ")
  if (!is.null(res)) {
    fit <- res
    cs <- fit$sets$cs
    row$converged <- isTRUE(fit$converged); row$niter <- as.integer(fit$niter); row$max_pip <- max(fit$pip)
    row$n_cs <- length(cs)
    if (length(cs)) {
      pur <- fit$sets$purity
      row$cs_purity_min <- min(pur$min.abs.corr); row$cs_purity_mean <- mean(pur$mean.abs.corr)
      row$cs_size <- paste(lengths(cs), collapse = ",")
      w <- coef(fit)[-1]
      csid <- rep(NA_character_, p)
      for (k in seq_along(cs)) {
        v <- cs[[k]]; lab <- names(cs)[k]
        csid[v] <- ifelse(is.na(csid[v]), lab, paste(csid[v], lab, sep = ";"))
      }
      keep <- which(!is.na(csid))
      S[[i]] <- data.table(gene = gid, gene_version = g, rsid = rs[keep], weight = w[keep], cs_id = csid[keep],
                           pip = fit$pip[keep], n_snps_locus = p, n_cs = length(cs), s_rss = row$s_rss, converged = row$converged)
    }
  }
  G[[i]] <- as.data.table(row)
  if (i %% 25 == 0) cat(sprintf("  %d/%d genes  %.1f min\n", i, length(genes), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
al <- bim[, .(rsid, a1, a2)]
add_alleles <- function(x) { x <- merge(rbindlist(x), al, by = "rsid", sort = FALSE); x }
S <- if (length(Filter(Negate(is.null), S))) add_alleles(S) else data.table()
T1 <- add_alleles(T1)
fwrite(S, file.path(out, sprintf("chr%s.susie.tsv.gz", chr)), sep = "\t")
fwrite(T1, file.path(out, sprintf("chr%s.top1.tsv.gz", chr)), sep = "\t")
fwrite(rbindlist(G), fin, sep = "\t")
cat(sprintf("done %s chr%s: %d genes, %d with CS, %d errors, %.1f min\n", tag, chr, length(genes),
            sum(rbindlist(G)$n_cs > 0, na.rm = TRUE), sum(rbindlist(G)$status != "ok"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
