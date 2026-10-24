# SCTWAS stage 1c: gather per-chromosome fits into the weights tables (SCTWAS_INSTRUCTIONS.md, Stage 1 step 6).
# usage: Rscript 02b_collect.R <tag>
# Outputs: weights/<tag>.susie.tsv.gz, weights/<tag>.top1.tsv.gz  (gene gene_version rsid a1 a2 weight cs_id pip
#          n_snps_locus n_cs s_rss converged), weights/<tag>.genes.tsv.gz (per-eGene diagnostics).
# The SNP lists in these tables are what stage 2 uses to rebuild the LD matrices from ld/panel_keep.
args <- commandArgs(trailingOnly = TRUE); tag <- args[1]
W <- "/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/sctwas"
suppressPackageStartupMessages(library(data.table))
d <- file.path(W, "stage1", tag)
chrs <- sub("^chr(.*)\\.genes\\.tsv\\.gz$", "\\1", list.files(d, "^chr.*\\.genes\\.tsv\\.gz$"))
want <- sub("^zin_chr(.*)\\.tsv\\.gz$", "\\1", list.files(file.path(W, "stage1", sub("_n[0-9]+$", "", tag)), "^zin_chr"))
miss <- setdiff(want, chrs)
if (length(miss)) stop("missing chromosomes: ", paste(miss, collapse = ","))
rd <- function(k) rbindlist(lapply(chrs, function(c) { f <- file.path(d, sprintf("chr%s.%s.tsv.gz", c, k)); if (file.size(f) > 30) fread(cmd = paste("zcat", shQuote(f)), colClasses = list(character = "rsid")) }), fill = TRUE)
G <- rd("genes"); S <- rd("susie"); T1 <- rd("top1")
# Implausible-weight flag (added after the Ext1 run; genes are flagged, never dropped): on the standardised scale a
# posterior-mean effect should not exceed the marginal |z_top|/sqrt(n) several-fold. Out-of-sample LD occasionally
# makes susie_rss fit large opposite-signed effects on near-collinear SNPs (|w| up to ~30), which blows up w'z/sqrt(w'Rw).
# w_ratio = max|w| / (|z_top|/sqrt(n)); w_flag = w_ratio > 3.
n <- if (grepl("_n[0-9]+$", tag)) as.numeric(sub(".*_n", "", tag)) else 757
if (nrow(S)) {
  mw <- S[, .(max_abs_w = max(abs(weight))), by = gene_version]
  G <- merge(G, mw, by = "gene_version", all.x = TRUE)
} else G[, max_abs_w := NA_real_]
G[, w_ratio := max_abs_w / (top_abs_z / sqrt(n))]
G[, w_flag := !is.na(w_ratio) & w_ratio > 3]
S <- merge(S, G[, .(gene_version, w_flag)], by = "gene_version")
cols <- c("gene", "gene_version", "rsid", "a1", "a2", "weight", "cs_id", "pip", "n_snps_locus", "n_cs", "s_rss", "converged")
T1 <- merge(T1, G[, .(gene_version, n_cs, s_rss, converged)], by = "gene_version")
T1[, `:=`(cs_id = NA_character_, pip = NA_real_)]
dir.create(file.path(W, "weights"), showWarnings = FALSE)
fwrite(S[, c(cols, "w_flag"), with = FALSE][order(gene)], file.path(W, "weights", paste0(tag, ".susie.tsv.gz")), sep = "\t")
fwrite(T1[, ..cols][order(gene)], file.path(W, "weights", paste0(tag, ".top1.tsv.gz")), sep = "\t")
fwrite(G[order(gene)], file.path(W, "weights", paste0(tag, ".genes.tsv.gz")), sep = "\t")
cat(sprintf("%s: %d eGenes fitted, %d errors, %d with >=1 CS, %d susie SNP rows, s_rss>0.3: %.3f, non-converged: %d, w_flag: %d\n",
            tag, nrow(G), sum(G$status != "ok"), sum(G$n_cs > 0, na.rm = TRUE), nrow(S), mean(G$s_rss > 0.3, na.rm = TRUE),
            sum(!G$converged, na.rm = TRUE), sum(G$w_flag)))
