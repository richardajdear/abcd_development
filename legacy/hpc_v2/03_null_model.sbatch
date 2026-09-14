#!/bin/bash
#SBATCH --job-name=v2_null
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-5
#
# Step 3: GENESIS null model, one phenotype per array task.
#
# The array index selects a manifest row (v1 convention: the phenotype set is
# defined by the export, not a hand-kept --array range; run_all.sh sizes the
# array from the manifest).  Outside SLURM, loops over all phenotypes.
#
# The null model is fitted ONCE per phenotype on the FULL sample -- relatives
# included, their covariance carried by the sparse PC-Relate kinship -- and
# reused by every per-chromosome association task in step 04.
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc_v2/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths PHENO COVAR_QUANT COVAR_CAT MANIFEST RSCRIPT
ensure_dirs "$NULL_DIR"

kin="$KIN_DIR/kinship_sparse.rds"
pcs="$KIN_DIR/pcair_pcs.tsv"
[[ -f "$kin" && -f "$pcs" ]] || { echo "FATAL: kinship outputs missing -- run 02_kinship first" >&2; exit 3; }

names=()   # not mapfile: bash 3.2 compatibility
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=$((SLURM_ARRAY_TASK_ID - 1))
  if (( idx >= ${#names[@]} )); then
    echo "array index $SLURM_ARRAY_TASK_ID exceeds ${#names[@]} phenotypes; nothing to do"
    exit 0
  fi
  todo=("${names[$idx]}")
else
  todo=("${names[@]}")
fi

# Empty-array expansion under `set -u`: bash 3.2 treats "${arr[@]}" on an EMPTY
# array as an unbound variable and aborts.  The ${arr[@]+...} guard expands to
# nothing when the array is empty and to the elements when it is not.  (bash 4.4+
# fixed this; CSD3 and macOS both predate that in places, so guard everywhere.)
strata_args=()
[[ -n "$STRATA_FILE" && -f "$STRATA_FILE" ]] && strata_args=(--strata "$STRATA_FILE")

for pheno in "${todo[@]}"; do
  echo "=== 03_null_model  $pheno  $(date) ==="
  run "$RSCRIPT" "$RSCRIPT_DIR/03_null_model.R" \
    --pheno-file "$PHENO" \
    --pheno "$pheno" \
    --covar-quant "$COVAR_QUANT" \
    --covar-cat "$COVAR_CAT" \
    --pcs "$pcs" \
    --kinship "$kin" \
    --n-pcs "$N_PCS" \
    ${strata_args[@]+"${strata_args[@]}"} \
    --out "$NULL_DIR/${pheno}_null.rds"
  if [[ "$DRY_RUN" != "1" ]]; then
    [[ -s "$NULL_DIR/${pheno}_null.rds" ]] \
      || { echo "FATAL: $NULL_DIR/${pheno}_null.rds missing/empty" >&2; exit 4; }
  fi
done

echo "=== 03_null_model done  $(date) ==="
