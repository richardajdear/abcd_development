#!/bin/bash
#SBATCH --job-name=dk_extract
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=slurm/%x_%j.log
#
# Parse FreeSurfer's own Desikan tables (stats/?h.aparc.stats) for every ABCD
# session into release-shaped tables under abcd-data-release-7.0/processed/dsk_local/.
# Comparison object for the published DK tables, not a replacement.
#
# Submit from the repo root:   sbatch tools/dk_extract.sbatch

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p slurm
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ABCD_CONDA_ENV:-legacy/hpc/work/envs/abcd}"
set -u
export PYTHONPATH=src:${PYTHONPATH:-}
python -m abcd.dk_stats --workers "${SLURM_CPUS_PER_TASK:-8}" "$@"
