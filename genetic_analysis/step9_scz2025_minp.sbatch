#!/bin/bash
#SBATCH --job-name=scz25_minp
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-2
#SBATCH --output=slurm/%x_%A_%a.log
# Step 9 min-p permutation; array task 1 = DK, 2 = HCP-MMP.
set -eo pipefail
cd /home/rajd2/rds/hpc-work/abcd_development
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
PARCS=(dsk hcp); P="${PARCS[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"
echo "=== min-p permutation, PARC=$P  $(date)"
legacy/hpc/work/envs/abcd/bin/python genetic_analysis/step9_scz2025_minp.py --parc "$P" "$@"
echo "=== done $(date)"
