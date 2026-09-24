#!/bin/bash
#SBATCH --job-name=c4_assoc
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --output=slurm/%x_%j.log
#
# Step 4 on CSD3: C4A GREx -> thinning rate, EUR arm, three phenotype cells:
#   hcp_1lmm       HCP-MMP single-LMM slope      [PRIMARY: M1 in this cell is the test]
#   dk_1lmm        DK single-LMM slope           (sensitivity: parcellation)
#   hcp_perregion  HCP mean of parcel BLUPs      (sensitivity: slope construction)
# plus step 5, the MHC content of the SCZ scores (does the PRS already carry C4?).
#
#   sbatch --dependency=afterok:<c4_impute job> c4_imputation/c4_assoc.sbatch
set -euo pipefail
REPO="${REPO:-/home/rajd2/rds/hpc-work/abcd_development}"; cd "$REPO"
C4="$REPO/c4_imputation"; WORK="${C4_WORK:-$C4/work}"; RES="${C4_RESULTS:-$C4/results}"; mkdir -p "$WORK" "$RES" slurm
export TMPDIR="${TMPDIR_OVERRIDE:-$WORK/tmp}"; mkdir -p "$TMPDIR"
RSCRIPT="${RSCRIPT:-/home/rajd2/rds/hpc-work/envs/genesis/bin/Rscript}"
PY="${PY:-$REPO/legacy/hpc/work/envs/abcd/bin/python}"
EUR_KEEP="${EUR_KEEP:-$REPO/legacy/hpc/work/results/ancestry/eur_anchor.keep}"
PRS_DIR="${PRS_DIR:-$REPO/genetic_analysis/work/scores_scz2025/SBayesRC/SCZ25_EUR}"
CALLS="${CALLS:-$WORK/c4_calls.tsv}"
GA="$REPO/genetic_analysis/work"
# (plain variables, not an associative array: macOS bash 3.2 runs this as a local test)
PH_hcp_1lmm="${PH_HCP_1LMM:-$GA/results_70tab_hcp/prs_final_1lmm/pheno}"
PH_dk_1lmm="${PH_DK_1LMM:-$GA/results_70tab/prs_final_1lmm/pheno}"
PH_hcp_perregion="${PH_HCP_PERREGION:-}"
if [[ -z "$PH_hcp_perregion" && -f genetic_analysis/config.sh ]]; then
  PH_hcp_perregion="$(PARC=hcp bash -c 'source genetic_analysis/config.sh >/dev/null 2>&1; echo "$PHENO_DIR"')"; fi
[[ -s "$CALLS" ]] || { echo "FATAL: $CALLS missing -- run c4_impute.sbatch first" >&2; exit 2; }
PRS_ARG=(); [[ -d "$PRS_DIR" ]] && PRS_ARG=(--prs-dir "$PRS_DIR") || echo "WARN: $PRS_DIR missing; M4 skipped"
for cell in hcp_1lmm dk_1lmm hcp_perregion; do
  v="PH_$cell"; d="${!v:-}"
  [[ -n "$d" && -s "$d/phenotypes_gcta.txt" ]] || { echo "skip $cell: no phenotypes at '${d}'"; continue; }
  echo "=== $cell  $d  $(date +%T)"
  "$RSCRIPT" "$C4/04_c4_assoc.R" --c4 "$CALLS" --pheno-dir "$d" --eur-ids "$EUR_KEEP" "${PRS_ARG[@]}" \
    --label "$cell" --out "$RES/c4_assoc_${cell}.tsv"
done
if [[ -z "${MHC_BIM:-}" && -f genetic_analysis/config.sh ]]; then
  MHC_BIM="$(bash -c 'source genetic_analysis/config.sh >/dev/null 2>&1; echo "$GENO_IMP_DIR"')/imp_union_chr6.bim"; fi
if [[ -z "${SKIP_MHC_CHECK:-}" && -s "${MHC_BIM:-}" ]]; then
  "$PY" "$C4/05_mhc_in_scores.py" --score-root "${SCORE_ROOT:-$GA/scores_scz2025}" \
    --bim "$MHC_BIM" \
    --out "$RES/mhc_in_scz_scores.tsv" || echo "WARN: MHC check failed (non-fatal)"
fi
echo "done $(date)"
