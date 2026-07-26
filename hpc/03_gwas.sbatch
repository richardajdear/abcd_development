#!/bin/bash
#SBATCH -A VERTES-SL2-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --job-name=abcd-gwas
#SBATCH --array=1-2
#SBATCH -o slurm/abcd-gwas-%A_%a.out
#SBATCH -e slurm/abcd-gwas-%A_%a.err
source "$(dirname "$0")/config.sh"

PHENO="${PHENOTYPES[$((SLURM_ARRAY_TASK_ID-1))]}"
PHENO_FILE="$PHENO_DIR/phenotypes_gcta.txt"
COL=$(head -1 "$PHENO_FILE" | tr ' ' '\n' | grep -n "^${PHENO}$" | cut -d: -f1)
COL=$((COL-2))

# fastGWA mixed model with the sparse GRM as random effect. This keeps all
# related individuals in the analysis (dropping them would discard roughly a
# third of ABCD) while controlling relatedness properly.
"$GCTA" --bfile "$GENO" --grm-sparse "${GRM}_sparse" \
        --fastGWA-mlm \
        --pheno "$PHENO_FILE" --mpheno "$COL" \
        --qcovar "$PHENO_DIR/covar_quant.txt" \
        --covar "$PHENO_DIR/covar_categorical.txt" \
        --maf 0.01 --geno 0.05 \
        --thread-num "$SLURM_CPUS_PER_TASK" \
        --out "$OUT/gwas/${PHENO}"

# Sanity: lambda_GC. A slope GWAS on 5.1 should show lambda ~1 and no hits;
# if lambda is well above 1 with this few effective observations, suspect
# residual population structure or a phenotype outlier, not real signal.
awk 'NR>1{print $10}' "$OUT/gwas/${PHENO}.fastGWA" | sort -g \
  | awk '{a[NR]=$1} END{printf "lambda_GC = %.4f\n", \
      (a[int(NR/2)]>0 ? 0.4549364/a[int(NR/2)] : -1)}'
