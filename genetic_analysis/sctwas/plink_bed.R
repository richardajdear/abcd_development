# Minimal PLINK .bed reader shared by SCTWAS stage 1 and stage 2 (no bed-reading package in the abcdR env).
# Dosage counts the bim column-5 allele (panel A1): 00 -> 2, 10 -> 1, 11 -> 0, 01 -> missing.
# open_bed() reads one chromosome's contiguous block of SNP-major bytes; ld_cor() returns the LD correlation matrix
# of the requested SNPs (mean-imputed, sample correlation over the 503 1000G EUR individuals).

open_bed <- function(prefix, chr) {
  bim <- data.table::fread(paste0(prefix, ".bim"), header = FALSE,
                           col.names = c("chr", "rsid", "cm", "pos", "a1", "a2"), colClasses = list(character = c(1, 2, 5, 6)))
  N <- nrow(data.table::fread(paste0(prefix, ".fam"), header = FALSE))
  nb <- ceiling(N / 4)
  idx <- which(bim$chr == as.character(chr))
  stopifnot(length(idx) > 0, all(diff(idx) == 1))
  con <- file(paste0(prefix, ".bed"), "rb"); on.exit(close(con))
  magic <- readBin(con, "raw", 3)
  stopifnot(identical(magic, as.raw(c(0x6c, 0x1b, 0x01))))
  seek(con, 3 + (idx[1] - 1) * nb)
  blk <- matrix(readBin(con, "raw", length(idx) * nb), nrow = nb)
  code <- c(2, NA, 1, 0)                       # 2-bit value 0,1,2,3
  lut <- t(sapply(0:255, function(b) code[bitwAnd(bitwShiftR(b, c(0, 2, 4, 6)), 3L) + 1]))   # 256 x 4
  list(bim = bim[idx], blk = blk, N = N, nb = nb, lut = lut, pos = setNames(seq_along(idx), bim$rsid[idx]))
}

geno <- function(B, rsids) {
  j <- B$pos[rsids]
  stopifnot(!anyNA(j))
  b <- B$blk[, j, drop = FALSE]
  g <- B$lut[as.integer(b) + 1L, , drop = FALSE]          # (nb*k) x 4, rows = byte-major within SNP
  X <- matrix(as.vector(t(g)), nrow = B$nb * 4)[seq_len(B$N), , drop = FALSE]
  colnames(X) <- rsids
  X
}

ld_cor <- function(B, rsids) {
  X <- geno(B, rsids)
  mu <- colMeans(X, na.rm = TRUE)
  miss <- which(is.na(X), arr.ind = TRUE)
  if (nrow(miss)) X[miss] <- mu[miss[, 2]]
  X <- sweep(X, 2, mu)
  s <- sqrt(colSums(X^2))
  X <- sweep(X, 2, s, "/")
  R <- crossprod(X)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R
}
