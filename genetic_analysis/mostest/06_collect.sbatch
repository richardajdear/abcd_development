#!/bin/bash
#SBATCH --job-name=m16_collect
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=slurm/%x_%j.log
# Step 16.7: tracked tables + gate verdicts (06_collect.py).
set -eo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/genetic_analysis/mostest/paths.sh"
"$PY" "$REPO/genetic_analysis/mostest/06_collect.py" --root "$M16"
