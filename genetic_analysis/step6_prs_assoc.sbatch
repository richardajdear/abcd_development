#!/bin/bash
#SBATCH --job-name=prs70
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-40
#
# README_HPC.md 3 step 6: the PRS ASSOCIATION LAYER ONLY.
#
# 3 is explicit that scoring is phenotype-independent and must be reused: "a
# score never touches our phenotype".  All 40 cells (4 methods x 10 trait-arms)
# already exist under legacy/hpc_v2/work/results_v2/prs_final/, with their
# _zanc within-ancestry-standardised copies -- verified 40/40 present, C+T with
# 8 thresholds each, the Bayesian methods with 1.  So this does NOT run
# setup/prs_final.sbatch, which would re-clump, re-run PRS-CS and re-run GCTB
# for days to reproduce files we already have.  6 Step 6 allows exactly this:
# "otherwise call prs_assoc.R directly over the existing score directories".
#
# --prs-dir therefore points at the LEGACY score directory (read-only) while
# every output goes to results_70tab/prs_final/ (rule 9, fresh dirs).
#
# $PHENO_DIR (FID = family_id) is the right export here, NOT the _fidiid
# sibling: prs_assoc.R and 06_prs_family.R both fit `(1 | family_id)` and read
# the family from the FID column (rule 2).
#
# Four calls per cell: population association and the Fulker between/within
# decomposition, each on the raw score and on the within-ancestry-standardised
# score.  `|| true` mirrors the original -- a cell with too few pairs should
# not kill the array.
#
#   sbatch genetic_analysis/step6_prs_assoc.sbatch
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate genetic_analysis/config.sh" >&2; exit 2; }
source "$_cfg"

SCORE_ROOT="${SCORE_ROOT:-$V2_LEGACY/prs_final}"
FINAL="${FINAL:-$OUT_V2/prs_final}"
ASSOC_R="${ASSOC_R:-$REPO_ROOT/tools/prs_assoc.R}"
FAM_R="${FAM_R:-$RSCRIPT_DIR/06_prs_family.R}"

PHENO="$PHENO_DIR/phenotypes_gcta.txt"
COVQ="$PHENO_DIR/covar_quant.txt"
COVC="$PHENO_DIR/covar_categorical.txt"
MANIF="$PHENO_DIR/phenotype_manifest.tsv"
require_paths PHENO COVQ COVC MANIF EUR_IDS RSCRIPT ASSOC_R FAM_R
ensure_dirs "$FINAL"

METHODS=(CT PRSCS SBayesR SBayesRC)
TRAITS=(SCZ_pooled SCZ_eur MDD_pooled MDD_eur ASD ALZ ALZ_noAPOE
        ALZ_IGAP ALZ_IGAP_noAPOE EA)
TODO=()
for m in "${METHODS[@]}"; do for t in "${TRAITS[@]}"; do TODO+=("$m:$t"); done; done

# A cell list may be passed POSITIONALLY (never via --export: sbatch splits
# that on commas, which once silently re-ran the last cell three times).
CELL_LIST="${1:-}"
if [[ -n "$CELL_LIST" ]]; then
  IFS=, read -r -a CELLS <<< "$CELL_LIST"
  IDX=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
  (( IDX >= 0 && IDX < ${#CELLS[@]} )) \
    || { echo "FATAL: task ${SLURM_ARRAY_TASK_ID:-1} has no entry in '$CELL_LIST'" >&2; exit 2; }
  SEL="${CELLS[$IDX]}"; T=-1
  for i in "${!TODO[@]}"; do [[ "${TODO[$i]}" == "$SEL" ]] && { T=$i; break; }; done
  (( T >= 0 )) || { echo "FATAL: '$SEL' is not a cell" >&2; printf '  %s\n' "${TODO[@]}" >&2; exit 2; }
else
  T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
fi
(( T >= 0 && T < ${#TODO[@]} )) || { echo "FATAL: cell index $((T+1)) outside 1-${#TODO[@]}" >&2; exit 2; }

SPEC="${TODO[$T]}"; METH="${SPEC%%:*}"; TRAIT="${SPEC##*:}"
OUTD="$SCORE_ROOT/$METH/$TRAIT"
ZD="$OUTD/_zanc"
[[ -d "$OUTD" ]] || { echo "FATAL: no score dir $OUTD" >&2; exit 3; }
n_prof=$(ls "$OUTD"/*.profile 2>/dev/null | wc -l)
(( n_prof > 0 )) || { echo "FATAL: no .profile in $OUTD" >&2; exit 3; }
echo "=== step6 $METH / $TRAIT : $n_prof score profiles  $(date) ==="

run "$RSCRIPT" "$ASSOC_R" --prs-dir "$OUTD" --pheno "$PHENO" --covar-quant "$COVQ" \
  --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
  --out "$FINAL/assoc_${METH}_${TRAIT}.tsv" || true
run "$RSCRIPT" "$FAM_R" --prs-dir "$OUTD" --pheno "$PHENO" --covar-quant "$COVQ" \
  --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
  --out "$FINAL/fam_${METH}_${TRAIT}.tsv" || true

if [[ -d "$ZD" ]]; then
  run "$RSCRIPT" "$ASSOC_R" --prs-dir "$ZD" --pheno "$PHENO" --covar-quant "$COVQ" \
    --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
    --out "$FINAL/assoc_${METH}_${TRAIT}_zanc.tsv" || true
  run "$RSCRIPT" "$FAM_R" --prs-dir "$ZD" --pheno "$PHENO" --covar-quant "$COVQ" \
    --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EUR_IDS" \
    --out "$FINAL/fam_${METH}_${TRAIT}_zanc.tsv" || true
else
  echo "WARN: no _zanc dir for $METH/$TRAIT"
fi
echo "=== step6 $METH / $TRAIT done  $(date) ==="
