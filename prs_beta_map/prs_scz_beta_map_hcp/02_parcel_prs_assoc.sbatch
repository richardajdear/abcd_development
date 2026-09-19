#!/bin/bash
#SBATCH --job-name=prs_beta_map
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=1-2
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 2 on CSD3: the per-parcel SCZ-PRS beta map on HCP-MMP, two cells:
#   1  SCZ25_META  SBayesRC  _zanc  -> full (pooled) arm   [primary]
#   2  SCZ25_EUR   SBayesRC  raw    -> EUR anchor set
# (the primary pooled and primary EUR cells of README_HPC.md step 9; the HCP
# global_slope betas to reproduce are -0.0276 and -0.0418, per-region slope,
# results_70tab_hcp/prs_scz2025/table_scz2025_main.tsv)
#
#   cd ~/rds/hpc-work/abcd_development
#   PARC=hcp sbatch --export=ALL,PARC=hcp prs_beta_map/prs_scz_beta_map_hcp/02_parcel_prs_assoc.sbatch
#
# Needs work/parcel_slopes_hcp70.tsv.gz rsynced from the laptop first (run_on_csd3.sh).
set -euo pipefail
export PARC=hcp
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/genetic_analysis/config.sh"
HERE="$REPO/prs_beta_map/prs_scz_beta_map_hcp"
ROOT="${SCZ25_ROOT:-$ABCD_HPC_ROOT/scores_scz2025}"
SLOPES="$HERE/work/parcel_slopes_hcp70.tsv.gz"
require_paths PHENO_DIR EUR_IDS RSCRIPT SLOPES
[[ "$(basename "$PHENO_DIR")" == *hcp* ]] || { echo "PHENO_DIR=$PHENO_DIR is not the HCP export; PARC=hcp not honoured by config.local.sh" >&2; exit 2; }
mkdir -p "$HERE/results" slurm

case "$SLURM_ARRAY_TASK_ID" in
  1) DIR="$ROOT/SBayesRC/SCZ25_META/_zanc"; STRATUM=full; CELL=SCZ25_META_SBayesRC_zanc ;;
  2) DIR="$ROOT/SBayesRC/SCZ25_EUR";        STRATUM=EUR;  CELL=SCZ25_EUR_SBayesRC_raw ;;
esac
[[ -d "$DIR" ]] || { echo "score dir $DIR missing; ls $ROOT/SBayesRC:"; ls "$ROOT/SBayesRC"; exit 3; }
echo "$CELL: $DIR  stratum=$STRATUM  PHENO_DIR=$PHENO_DIR  $(date)"
run "$RSCRIPT" "$HERE/02_parcel_prs_assoc.R" --slopes "$SLOPES" --pheno-dir "$PHENO_DIR" \
  --prs-dir "$DIR" --eur-ids "$EUR_IDS" --stratum "$STRATUM" --cell "$CELL" \
  --out "$HERE/results/beta_map_${CELL}.tsv"
echo "done $(date)"
