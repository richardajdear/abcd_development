#!/bin/bash
#SBATCH --job-name=v2_cond
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
"$RSCRIPT" "$V2ROOT/R/07_prs_conditional.R" \
  --final "$FINAL" --pheno "$PHENO" --covar-quant "$COVQ" --covar-cat "$COVC" \
  --manifest "$MANIF" --eur-ids "$EURKEEP" \
  --out "$FINAL/table_conditional_EA.tsv"
