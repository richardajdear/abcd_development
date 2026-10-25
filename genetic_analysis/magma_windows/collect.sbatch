#!/bin/bash
#SBATCH --job-name=magmawin_collect
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/%x_%A.log
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/magma_windows/common.sh
"$PY" "$REPO/genetic_analysis/magma_windows/collect.py"
rmdir "$W/tmp" 2>/dev/null || true
