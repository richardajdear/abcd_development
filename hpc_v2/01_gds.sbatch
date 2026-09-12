#!/bin/bash
#SBATCH --job-name=v2_gds
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=0-22
#
# Step 1: PLINK -> GDS.
#
# Task 0 converts the merged ARRAY fileset (kinship input); tasks 1-22 convert
# the per-chromosome IMPUTED filesets (association input).  Missing imputed
# chromosomes are skipped with a note -- 00_check_inputs decides whether that
# is fatal.  Outside SLURM, converts the array fileset and every imputed
# fileset that exists, sequentially.
#
# Sized for backfill: conversion is I/O-bound, 4 CPUs is plenty (v1 lesson:
# ask for what the job needs, queue time dominates compute time).
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc_v2/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths GENO_ARRAY RSCRIPT
ensure_dirs "$GDS_DIR"

convert_one() {
  local bfile="$1" out="$2"
  if [[ ! -f "$bfile.bed" ]]; then
    echo "SKIP: $bfile.bed not found"
    return 0
  fi
  run "$RSCRIPT" "$RSCRIPT_DIR/01_make_gds.R" --bfile "$bfile" --out "$out"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  t="$SLURM_ARRAY_TASK_ID"
  if [[ "$t" == "0" ]]; then
    convert_one "$GENO_ARRAY" "$GDS_DIR/array.gds"
  else
    convert_one "$(geno_imp_prefix "$t")" "$GDS_DIR/imp_chr${t}.gds"
  fi
else
  convert_one "$GENO_ARRAY" "$GDS_DIR/array.gds"
  for chr in $(seq 1 22); do
    convert_one "$(geno_imp_prefix "$chr")" "$GDS_DIR/imp_chr${chr}.gds"
  done
fi

echo "=== 01_gds done  $(date) ==="
