# Thinned GWAS summary statistics

Full summary statistics are ~9.08M SNPs and ~700 MB per phenotype per arm, which
neither GitHub nor a plotting script wants. These keep everything that shows in a
figure and sample the rest.

**Thinning rule:** every SNP with **p < 1e-3** is kept; the remainder is sampled
deterministically at **1 in 50**. The `kept_frac` column records which rule
applied to each row (`1` or `0.02`).

**For a Manhattan plot** this is lossless where it matters — no point above
−log10(p) = 3 has been dropped.

**For a QQ plot you must correct the expected quantiles**, or the null will look
inflated: the sampled rows each stand for 50 SNPs. Reconstruct the expected
distribution using `kept_frac` (i.e. weight each thinned row by `1/kept_frac`)
rather than ranking the rows as if they were the whole genome.

Columns: `CHR POS SNP A1 A2 BETA SE P kept_frac`. Positions are **GRCh38**
(the imputed genotypes' build); `SNP` is an rsID, which is build-independent.

Arms: `gwas_imp__*` = pooled multi-ancestry (n = 7,932);
`gwas_imp_eur__*` = EUR-stratified (n = 4,039).
