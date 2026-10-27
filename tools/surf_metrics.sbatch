#!/bin/bash
#SBATCH --job-name=surf_metrics
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/%x_%A_%a.log
#
# D1: sample GWC and T1w intensity on the FreeSurfer surfaces, per HCP-MMP and DK parcel,
# for one chunk of sessions per array task (src/abcd/surf_metrics.py).  Chunks are 200
# sessions in abcd.hcp_stats discovery order; a task whose part file exists is skipped,
# so the array can simply be resubmitted after failures.
#
# Submit from the repo root:
#   N=$(PYTHONPATH=src python -m abcd.surf_metrics --list-chunks)
#   sbatch --array=0-$((N-1)) tools/surf_metrics.sbatch
# then:  PYTHONPATH=src python tools/surf_metrics_collect.py   (on a compute node)
# ~3 s per session per core: 200 sessions on 8 cores is ~80 s of sampling plus ~1 min discovery.

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p slurm

# conda.sh references PS1, so -u only after activation
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ABCD_CONDA_ENV:-legacy/hpc/work/envs/abcd}"
set -u

export PYTHONPATH=src:${PYTHONPATH:-}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m abcd.surf_metrics --chunk "${SLURM_ARRAY_TASK_ID}" \
    --workers "${SLURM_CPUS_PER_TASK:-8}" "$@"
