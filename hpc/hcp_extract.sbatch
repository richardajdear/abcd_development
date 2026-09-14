#!/bin/bash
#SBATCH --job-name=hcp_extract
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=hpc/slurm/%x_%j.log
#
# Parse the HCP-MMP1.0 mris_anatomical_stats tables for every ABCD session
# into release-style tables under abcd-7.0/processed/hcp/.  See
# src/abcd/hcp_stats.py for what is read and written.
#
# Submit from the repo root:   sbatch hpc/hcp_extract.sbatch
# Roughly 30k sessions x 3 small files; ~10 min on 16 cores.

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p hpc/slurm

# conda.sh references PS1, so -u only after activation
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hpc/work/envs/abcd
set -u

export PYTHONPATH=src:${PYTHONPATH:-}
python -m abcd.hcp_stats --workers "${SLURM_CPUS_PER_TASK:-8}" "$@"
