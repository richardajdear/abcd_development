#!/bin/bash
#SBATCH --job-name=v2_kinship
#SBATCH --partition=icelake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=%x_%j.log
#
# Step 2: KING-robust -> PC-AiR -> PC-Relate on the array GDS.
#
# This is the replacement for BOTH v1 relatedness artefacts:
#   - GCTA --grm-cutoff on a pooled GRM kept a 94%-EUR "unrelated" set
#     (hpc/README_HPC.md §8.6); PC-AiR partitions evenly across ancestry.
#   - The GCTA PCA mixed family structure into the ancestry axes; PC-AiR
#     computes PCs on unrelateds and projects relatives in.
#
# ~12k samples x ~100k pruned SNPs: KING and PC-Relate are each a few
# CPU-hours; 64G covers the dense KING matrix (11670^2 doubles ~ 1.1 GB, plus
# GENESIS working copies).
#
# Sanity gate on completion: the PC-AiR unrelated count should be in the same
# range as ABCD's shipped genesis/unrelateds_individuals.txt (8,181 of 11,670).
# A materially smaller number means the KING seed or thresholds are off.
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc_v2/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths RSCRIPT
ensure_dirs "$KIN_DIR"

gds="$GDS_DIR/array.gds"
[[ -f "$gds" ]] || { echo "FATAL: $gds not found -- run 01_gds first" >&2; exit 3; }

run "$RSCRIPT" "$RSCRIPT_DIR/02_kinship.R" \
  --gds "$gds" \
  --out-dir "$KIN_DIR" \
  --maf "$KIN_MAF" \
  --ld-r2 "$KIN_LD_R2" \
  --kin-thresh "$KIN_THRESH" \
  --sparse-thresh "$SPARSE_KIN_THRESH" \
  --max-density "$SPARSE_KIN_MAX_DENSITY" \
  --threads "$THREADS"

# Verify outputs, not exit status (v1 lesson: 13 jobs COMPLETED with nothing).
if [[ "$DRY_RUN" != "1" ]]; then
  for f in kinship_sparse.rds pcair_pcs.tsv pcair_unrelated.txt kinship_summary.tsv; do
    [[ -s "$KIN_DIR/$f" ]] || { echo "FATAL: expected output $KIN_DIR/$f missing/empty" >&2; exit 4; }
  done
  echo
  column -t -s $'\t' "$KIN_DIR/kinship_summary.tsv" || cat "$KIN_DIR/kinship_summary.tsv"
fi

echo "=== 02_kinship done  $(date) ==="
