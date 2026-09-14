#!/bin/bash
#SBATCH --job-name=v2_assoc
#SBATCH --partition=icelake
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-110
#
# Step 4: GENESIS single-variant association, one (phenotype, chromosome) pair
# per array task: task = (pheno_idx - 1) * 22 + chr.  run_all.sh sizes the
# array as n_phenotypes * 22.  Outside SLURM, loops over everything that has
# both a null model and a GDS.
#
# After the last chromosome of each phenotype, the collector runs (cheap,
# idempotent, rebuilt from whatever chromosome files exist) producing
# <pheno>.sumstats.tsv.gz + the shared gwas_summary.tsv.
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc_v2/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths MANIFEST RSCRIPT
ensure_dirs "$ASSOC_DIR"

names=()   # not mapfile: bash 3.2 compatibility
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)

run_pair() {
  local pheno="$1" chr="$2"
  local nm="$NULL_DIR/${pheno}_null.rds"
  local gds="$GDS_DIR/imp_chr${chr}.gds"
  [[ -f "$nm" ]]  || { echo "FATAL: $nm missing -- run 03_null_model" >&2; return 4; }
  if [[ ! -f "$gds" ]]; then
    echo "SKIP: $gds not found (chromosome not converted)"
    return 0
  fi
  echo "=== 04_assoc  $pheno chr$chr  $(date) ==="
  run "$RSCRIPT" "$RSCRIPT_DIR/04_assoc.R" \
    --gds "$gds" \
    --null-model "$nm" \
    --out "$ASSOC_DIR/${pheno}_chr${chr}.tsv.gz"
}

collect() {
  local pheno="$1"
  run "$RSCRIPT" "$RSCRIPT_DIR/05_collect_assoc.R" \
    --assoc-dir "$ASSOC_DIR" \
    --pheno "$pheno" \
    --out-dir "$ASSOC_DIR" \
    --maf "$ASSOC_MAF" \
    --mac "$ASSOC_MAC"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  t=$((SLURM_ARRAY_TASK_ID - 1))
  pheno_idx=$((t / 22))
  chr=$((t % 22 + 1))
  if (( pheno_idx >= ${#names[@]} )); then
    echo "array index $SLURM_ARRAY_TASK_ID exceeds ${#names[@]} phenotypes x 22; nothing to do"
    exit 0
  fi
  pheno="${names[$pheno_idx]}"
  run_pair "$pheno" "$chr"
  # The last chromosome of each phenotype collects.  If chromosomes finish out
  # of order the collector may run on a partial set here; run_all.sh submits a
  # dedicated final collection pass with afterok on the whole array, which is
  # the one to trust.  This inline call is a convenience preview.
  if (( chr == 22 )); then collect "$pheno" || true; fi
else
  for pheno in "${names[@]}"; do
    for chr in $(seq 1 22); do run_pair "$pheno" "$chr"; done
    collect "$pheno"
  done
fi

echo "=== 04_assoc done  $(date) ==="
