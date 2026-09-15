#!/bin/bash
#SBATCH --job-name=prs_1lmm
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-2
#SBATCH --output=genetic_analysis/work/results_70tab_hcp/logs/%x_%A_%a.log
# Order-of-operations check: PRS associations on the slope of ONE LMM fitted to
# the per-scan whole-cortex mean, vs the pipeline's mean of 68/358 per-region
# slope BLUPs.  Task 1 = DK, task 2 = HCP.  Matched SCZ/MDD cells only.
set -eo pipefail
export TMPDIR="${SLURM_SUBMIT_DIR}/genetic_analysis/work/tmp"
cd "${SLURM_SUBMIT_DIR}"
export PYTHONPATH=src
PY=$PWD/legacy/hpc/work/envs/abcd/bin/python
R=/home/rajd2/rds/hpc-work/envs/genesis/bin/Rscript
if [[ "$SLURM_ARRAY_TASK_ID" == 1 ]]; then RUN=out/thickness_dsk_70_139406217085_globalmean; PH=genetic_analysis/work/pheno_70tab; OUT=genetic_analysis/work/results_70tab/prs_final_1lmm
else RUN=out/thickness_hcp_70_aa6e91efba82_globalmean; PH=genetic_analysis/work/pheno_70tab_hcp; OUT=genetic_analysis/work/results_70tab_hcp/prs_final_1lmm; fi
mkdir -p "$OUT/pheno"
# --- phenotype file from the single-LMM BLUPs, aligned to the atlas's analysis set ---
"$PY" genetic_analysis/work/pc_compare/build_1lmm_pheno.py "$RUN" "$PH" "$OUT/pheno"
cp "$PH/covar_quant.txt" "$PH/covar_categorical.txt" "$OUT/pheno/"
EUR=legacy/hpc/work/results/ancestry/eur_anchor.keep; SC=legacy/hpc_v2/work/results_v2/prs_final
for m in CT PRSCS SBayesR SBayesRC; do for t in SCZ_pooled SCZ_eur MDD_pooled MDD_eur; do
  for z in "" _zanc; do
    d="$SC/$m/$t${z:+/_zanc}"
    echo "=== $m/$t$z  $(date -u +%T)"
    "$R" tools/prs_assoc.R --prs-dir "$d" --pheno "$OUT/pheno/phenotypes_gcta.txt" --covar-quant "$OUT/pheno/covar_quant.txt" \
      --covar-cat "$OUT/pheno/covar_categorical.txt" --manifest "$OUT/pheno/phenotype_manifest.tsv" --eur-ids "$EUR" \
      --out "$OUT/assoc_${m}_${t}${z}.tsv" >/dev/null 2>&1 || echo "  (failed: $m/$t$z)"
  done
done; done
echo "=== done $(date -u +%T): $(ls $OUT/assoc_*.tsv | wc -l)/32 outputs"
