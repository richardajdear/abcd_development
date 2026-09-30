#!/bin/bash
#SBATCH --job-name=m16_most
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --array=1-4
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 16.4: MOSTest per family (tasks 1-3 = FAMILIES in paths.sh: slope, ct,
# slopeols) and, as task 4, the engine check (gate G1, 04_engine_check.py).
set -eo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/genetic_analysis/mostest/paths.sh"
i=${SLURM_ARRAY_TASK_ID:-1}
cd "$REPO/genetic_analysis/mostest"          # mostest_core import
if (( i <= ${#FAMILIES[@]} )); then
  FAM="${FAMILIES[$((i-1))]}"
  echo "=== 16.4 MOSTest  $FAM  $(date -u +%T)"
  "$PY" 03_mostest.py --family "$FAM" --root "$M16"
  [[ -s "$SS/mostest_${FAM}.sumstats.tsv.gz" ]] || { echo "FATAL: no sumstats" >&2; exit 4; }
else
  echo "=== 16.5 engine check  $(date -u +%T)"
  "$PY" 04_engine_check.py --root "$M16" --genesis-dir "$GENESIS_1LMM"
fi
echo "=== done  $(date -u +%T)"
