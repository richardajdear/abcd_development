#!/bin/bash
#SBATCH --job-name=v2_eurdown
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#
# EUR-stratified arm, downstream: MAGMA and LDSC on the v2 EUR sumstats.
#
# Deliberately runs v1's OWN 04_magma.sbatch and 05_ldsc_rg.sbatch, with
# GWAS_DIR / MAGMA_DIR / LDSC_DIR redirected, rather than reimplementing them.
# The point of this run is a v1-vs-v2 comparison; using the same code on both
# sides means any difference is the GWAS, not the downstream analysis.
#
# Staging: MAGMA cannot read gzip and wants SNP/P/N by name; LDSC's munge wants
# A1/A2 as well.  One uncompressed file per phenotype serves both.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/genetic_analysis/config.sh"
STAGE="$OUT_V2/gwas_eur_staged"
ensure_dirs "$STAGE"

for pheno in $(phenotype_names); do
  src="$OUT_V2/assoc_eur/$pheno.sumstats.tsv.gz"
  dst="$STAGE/$pheno.fastGWA"
  [[ -f "$src" ]] || { echo "FATAL: $src missing" >&2; exit 3; }
  if [[ ! -s "$dst" ]]; then
    echo "=== stage $pheno  $(date) ==="
    # v2 sumstats columns: SNP CHR POS A1 A2 AF1 N BETA SE P
    zcat "$src" | awk 'BEGIN{OFS="\t"} NR==1{print "SNP","A1","A2","AF1","N","BETA","SE","P"; next}
                       {print $1,$4,$5,$6,$7,$8,$9,$10}' > "$dst"
  fi
  echo "  $pheno: $(( $(wc -l < "$dst") - 1 )) SNPs"
done

export GWAS_DIR="$STAGE"
export MAGMA_DIR="$OUT_V2/magma_eur"
export LDSC_DIR="$OUT_V2/ldsc_eur"
export DRY_RUN=0
mkdir -p "$MAGMA_DIR" "$LDSC_DIR"

echo "=== MAGMA (v1 04_magma.sbatch, GWAS_DIR redirected)  $(date) ==="
bash "$REPO/legacy/hpc/04_magma.sbatch"

echo "=== LDSC  (v1 05_ldsc_rg.sbatch, GWAS_DIR redirected)  $(date) ==="
bash "$REPO/legacy/hpc/05_ldsc_rg.sbatch"

echo "=== done  $(date) ==="
