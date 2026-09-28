#!/bin/bash
#SBATCH --job-name=c3axis_build
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=genetic_analysis/work/results_70tab_hcp/logs/%x_%j.log
# README_HPC.md section 8, step C3-A: C3-axis phenotypes.
#   1. per-scan transcriptional-axis scores   (c3axis/build_c3axis_model_table.py)
#   2. one LMM per score, same formula as the primary trait   (R/fit_lmm.R)
#   3. GCTA/GENESIS/prs_assoc export on the 8,596-child HCP analysis set   (c3axis/build_c3axis_pheno.py)
# Submit from the repo root: sbatch genetic_analysis/c3axis/run_c3axis_build.sbatch
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
export TMPDIR="$PWD/genetic_analysis/work/tmp"; mkdir -p "$TMPDIR"     # rule 10
export PYTHONPATH=src
PY=$PWD/legacy/hpc/work/envs/abcd/bin/python
R=$PWD/legacy/hpc/work/envs/abcdR/bin/Rscript
RUN=out/thickness_hcp_70_aa6e91efba82
PH=genetic_analysis/work/pheno_70tab_hcp
OUT=genetic_analysis/work/results_70tab_hcp/c3axis
for f in "$RUN/fits/blups.parquet" "$RUN/fits/fixed.parquet" "$RUN/model_table.parquet" "$PH/phenotypes_gcta.txt"; do
  [[ -s "$f" ]] || { echo "FATAL: missing $f" >&2; exit 2; }
done
mkdir -p "$OUT"
echo "=== 1 build  $(date -u +%T)"
"$PY" genetic_analysis/c3axis/build_c3axis_model_table.py --run "$RUN"
git diff --exit-code genetic_analysis/c3axis/c3axis_weights_check.tsv \
  && echo "weights check identical to the committed laptop build" \
  || echo "WARNING: c3axis_weights_check.tsv differs from the committed version -- inspect before step C3-B"
echo "=== 2 fit  $(date -u +%T)"
"$R" R/fit_lmm.R --run-dir "${RUN}_c3axis" --cores 7
echo "=== 3 export  $(date -u +%T)"
"$PY" genetic_analysis/c3axis/build_c3axis_pheno.py "${RUN}_c3axis" "$PH" "$OUT/pheno"
# the global_mean label must reproduce the primary single-LMM trait
"$PY" -c "
import pandas as pd
a = pd.read_csv('$OUT/pheno/phenotypes_gcta.txt', sep=' ', dtype=str).set_index('IID')['global_slope_c3axis'].astype(float)
b = pd.read_csv('genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno/phenotypes_gcta.txt', sep=' ', dtype=str).set_index('IID')['global_slope_1lmm'].astype(float)
r = a.corr(b.reindex(a.index)); print(f'r(global_slope_c3axis, global_slope_1lmm) = {r:.4f}')
assert r > 0.99, 'global_mean label does not reproduce the primary trait -- stop'
"
echo "=== done $(date -u +%T)"
