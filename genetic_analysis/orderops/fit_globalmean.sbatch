#!/bin/bash
#SBATCH --job-name=fit_gmean
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=genetic_analysis/work/results_70tab_hcp/logs/%x_%j.log
# One LMM on the per-scan whole-cortex mean, per atlas (order-of-operations check).
set -eo pipefail
export TMPDIR="${SLURM_SUBMIT_DIR}/genetic_analysis/work/tmp"   # node-local /tmp was full (job 35616285)
cd "${SLURM_SUBMIT_DIR}"
export PYTHONPATH=src PY=$PWD/legacy/hpc/work/envs/abcd/bin/python
R=$PWD/legacy/hpc/work/envs/abcdR/bin/Rscript
for d in out/thickness_dsk_70_139406217085_globalmean out/thickness_hcp_70_aa6e91efba82_globalmean; do
  echo "=== $d  $(date -u +%T)"; "$R" R/fit_lmm.R --run-dir "$d" --cores 4
  "$PY" genetic_analysis/work/pc_compare/check_blups.py "$d"
done
echo "=== done $(date -u +%T)"
