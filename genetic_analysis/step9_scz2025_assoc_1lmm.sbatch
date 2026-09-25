#!/bin/bash
#SBATCH --job-name=scz25_assoc_1lmm
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-2
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 9, order-of-operations arm: the 2025-SCZ scores against the slope of
# ONE LMM fitted to the per-scan whole-cortex mean (global_slope_1lmm), the
# phenotype the 2026-09-15 order-of-operations check built
# (results_70tab*/prs_final_1lmm/pheno/, aligned to the same 8,596 children),
# instead of the pipeline's mean of per-region slope BLUPs.  Same
# tools/prs_assoc.R, same covariates, raw and within-cluster-standardised
# scores; task 1 = DK, task 2 = HCP-MMP.  Output:
# results_70tab*/prs_scz2025_1lmm/assoc_*.tsv, then
# step9_scz2025_orderops.py writes prs_scz2025/table_order_of_operations_scz2025.tsv.
set -eo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; cd "$REPO"
export TMPDIR="$REPO/genetic_analysis/work/tmp"
R=/home/rajd2/rds/hpc-work/envs/genesis/bin/Rscript
PY="$REPO/legacy/hpc/work/envs/abcd/bin/python"
SUF=""; [[ "${SLURM_ARRAY_TASK_ID:-1}" == 2 ]] && SUF="_hcp"
ROOT="$REPO/genetic_analysis/work/results_70tab$SUF"
PH="$ROOT/prs_final_1lmm/pheno"; OUT="$ROOT/${PRS_TAG:-prs_scz2025}_1lmm"; mkdir -p "$OUT"
ARM_GLOB="${ARM_GLOB:-SCZ25_*}"   # step 15 reuse; defaults are step 9
SC="${SCZ25_ROOT:-$REPO/genetic_analysis/work/scores_scz2025}"; EUR="$REPO/legacy/hpc/work/results/ancestry/eur_anchor.keep"
[[ -s "$PH/phenotypes_gcta.txt" ]] || { echo "FATAL: $PH/phenotypes_gcta.txt missing (run the order-of-operations fit first)" >&2; exit 2; }
n=0
for d in "$SC"/{CT,PRSCS,SBayesR,SBayesRC,PRSCSX}/$ARM_GLOB; do
  [[ -d "$d" ]] || continue; ls "$d"/score_*.profile >/dev/null 2>&1 || continue
  M="$(basename "$(dirname "$d")")"; T="$(basename "$d")"
  for z in "" _zanc; do
    dd="$d${z:+/_zanc}"; [[ -d "$dd" ]] || continue
    echo "=== $M/$T$z  $(date +%T)"
    "$R" tools/prs_assoc.R --prs-dir "$dd" --pheno "$PH/phenotypes_gcta.txt" --covar-quant "$PH/covar_quant.txt" \
      --covar-cat "$PH/covar_categorical.txt" --manifest "$PH/phenotype_manifest.tsv" --eur-ids "$EUR" \
      --out "$OUT/assoc_${M}_${T}${z}.tsv" >/dev/null 2>&1 && n=$((n+1)) || echo "  (failed: $M/$T$z)"
  done
done
echo "=== $n association tables in $OUT  $(date +%T)"
[[ -z "${PRS_TAG:-}" ]] && "$PY" genetic_analysis/step9_scz2025_orderops.py || true
