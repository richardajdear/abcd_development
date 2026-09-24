#!/bin/bash
# Step 1: pull the extended MHC out of a PLINK fileset as a bgzipped VCF.
#
#   PLINK=... BCFTOOLS=... bash 01_extract_mhc.sh <bfile prefix> <out prefix> [keep file]
#
# Input must be GRCh37 (the ABCD 7.0 Smokescreen array fileset is hg19).  The
# region is the span of the Sekar 2016 C4 reference panel (imputec4 README):
# 6:24894177-33890574.  Hard-called array genotypes are used, not imputed
# dosages: the panel is dense only at HapMap3 sites and imputation error in
# the MHC would propagate into the C4 call.
#
# Filters: biallelic A/C/G/T SNPs, per-variant missingness < 5 %, MAF >= 0.5 %.
# Strand and allele orientation are NOT fixed here -- 02_impute_c4.sh runs
# conform-gt against the panel, which does that by genotype correlation.
# Output: <out>.vcf.gz (+ .tbi), sample IDs = the .fam IID column.
set -euo pipefail
BFILE="$1"; OUT="$2"; KEEP="${3:-}"
PLINK="${PLINK:-plink}"; BCFTOOLS="${BCFTOOLS:-bcftools}"
REGION_FROM=24894177; REGION_TO=33890574
for x in bed bim fam; do [[ -s "$BFILE.$x" ]] || { echo "missing $BFILE.$x" >&2; exit 2; }; done
"$PLINK" --version | grep -q "v1.9" || { echo "PLINK must be 1.9 (got: $("$PLINK" --version))" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"
"$PLINK" --bfile "$BFILE" --chr 6 --from-bp "$REGION_FROM" --to-bp "$REGION_TO" \
  --snps-only just-acgt --geno 0.05 --maf 0.005 ${KEEP:+--keep "$KEEP"} \
  --recode vcf-iid bgz --out "$OUT" > "$OUT.plink.stdout"
"$BCFTOOLS" index -f -t "$OUT.vcf.gz"
n_var=$("$BCFTOOLS" index -n "$OUT.vcf.gz"); n_smp=$("$BCFTOOLS" query -l "$OUT.vcf.gz" | wc -l | tr -d ' ')
echo "extract_mhc: $n_var variants x $n_smp samples in chr6:$REGION_FROM-$REGION_TO -> $OUT.vcf.gz"
[[ "$n_var" -gt 500 ]] || { echo "only $n_var MHC variants -- wrong build or chromosome coding?" >&2; exit 3; }
