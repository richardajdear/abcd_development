#!/bin/bash
#SBATCH --job-name=minp70
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=genetic_analysis/work/results_70tab/logs/%x_%j.log
# README_HPC.md 3 step 6, min-p permutation (C+T only).  See the .py docstring.
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
legacy/hpc/work/envs/abcd/bin/python genetic_analysis/step6_minp_permutation.py "$@"
