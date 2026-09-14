#!/bin/bash
#SBATCH --job-name=v2_renorm
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=1:00:00
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
"$PYTHON" "$V2ROOT/work/setup/normalise_gwas.py" "$GWAS_DIR" "$FINAL/gwas" \
  "$IN/ref/maf_union.tsv" ASD,ALZ
wc -l "$FINAL"/gwas/*.ma
