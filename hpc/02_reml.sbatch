#!/bin/bash
#SBATCH -A VERTES-SL2-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --job-name=abcd-reml
#SBATCH --array=1-2
#SBATCH -o slurm/abcd-reml-%A_%a.out
#SBATCH -e slurm/abcd-reml-%A_%a.err
source "$(dirname "$0")/config.sh"

PHENO="${PHENOTYPES[$((SLURM_ARRAY_TASK_ID-1))]}"
PHENO_FILE="$PHENO_DIR/phenotypes_gcta.txt"
# GCTA takes the phenotype by column position (1-based, after FID/IID).
COL=$(head -1 "$PHENO_FILE" | tr ' ' '\n' | grep -n "^${PHENO}$" | cut -d: -f1)
COL=$((COL-2))
echo "phenotype $PHENO -> --mpheno $COL"

# Unrelated GRM: SNP-h2 without shared-environment inflation.
"$GCTA" --grm "${GRM}_unrel" --pheno "$PHENO_FILE" --mpheno "$COL" \
        --qcovar "$PHENO_DIR/covar_quant.txt" \
        --covar "$PHENO_DIR/covar_categorical.txt" \
        --reml --reml-no-constrain \
        --thread-num "$SLURM_CPUS_PER_TASK" \
        --out "$OUT/reml/${PHENO}_unrel"

# Full GRM for comparison. If h2(full) >> h2(unrelated), shared environment
# or assortative mating is inflating the family-based estimate -- report both.
"$GCTA" --grm "$GRM" --pheno "$PHENO_FILE" --mpheno "$COL" \
        --qcovar "$PHENO_DIR/covar_quant.txt" \
        --covar "$PHENO_DIR/covar_categorical.txt" \
        --reml --thread-num "$SLURM_CPUS_PER_TASK" \
        --out "$OUT/reml/${PHENO}_full"

grep -H "^V(G)/Vp" "$OUT/reml/${PHENO}"_*.hsq
