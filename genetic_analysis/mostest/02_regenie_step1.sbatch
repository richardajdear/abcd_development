#!/bin/bash
#SBATCH --job-name=m16_rg1
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-7
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 16.2: REGENIE step 1 (whole-genome ridge regression -> LOCO polygenic
# predictions), one array task per run in paths.sh RUNS:
#   1 slope  2 slope_perm  3 ct  4 ct_perm  5 slopeols  6 slopeols_perm  7 global
# Relatives are KEPT: REGENIE's LOCO predictions absorb relatedness, the same
# role GENESIS's sparse kinship plays in steps 3-4.  Phenotypes are
# rank-inverse-normalised (--apply-rint, both steps), as in MOSTest.
set -eo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/genetic_analysis/mostest/paths.sh"
[[ -x "$REGENIE" ]] || { echo "FATAL: REGENIE=$REGENIE not executable (build envs/regenie_env.yml)" >&2; exit 2; }
RUN="${RUNS[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"; FAM="$(run_family "$RUN")"
SUF=""; run_is_perm "$RUN" && SUF="_perm"
PHENO="$PH/pheno_${FAM}${SUF}.txt"; COVAR="$PH/covar${SUF}.txt"
[[ -s "$PHENO" && -s "$COVAR" && -s "$GENO_DIR/step1_array.bed" ]] || { echo "FATAL: run 01_build first" >&2; exit 2; }
# Skip only a COMPLETE step 1: REGENIE writes pred.list before the last .loco,
# so a crashed run leaves a list that points at missing files.
k_exp=$(wc -l < "$PH/traits_${FAM}.txt")
if [[ -s "$RG/step1_${RUN}_pred.list" ]] && (( $(wc -l < "$RG/step1_${RUN}_pred.list") == k_exp )) \
   && awk '{print $2}' "$RG/step1_${RUN}_pred.list" | while read -r f; do [[ -s "$f" ]] || exit 1; done; then
  echo "exists: $RG/step1_${RUN}_pred.list"; exit 0
fi
rm -f "$RG/step1_${RUN}_pred.list" "$RG/step1_${RUN}"_*.loco
echo "=== 16.2 step 1  run=$RUN  $(date -u +%T)"
# Per-run level-0 directory: a shared prefix let run "slope" delete
# "slope_perm"'s files (l0_slope* matches l0_slope_perm_*) while it ran.
L0="$TMPDIR/l0_${RUN}"; rm -rf "$L0"; mkdir -p "$L0"
"$REGENIE" --step 1 \
  --bed "$GENO_DIR/step1_array" --keep "$PH/keep_imp.txt" \
  --phenoFile "$PHENO" --covarFile "$COVAR" \
  --catCovarList sex,site --maxCatLevels 30 \
  --qt --apply-rint --bsize 1000 \
  --lowmem --lowmem-prefix "$L0/l0" \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  --out "$RG/step1_${RUN}"
k_got=$(wc -l < "$RG/step1_${RUN}_pred.list")
(( k_got == k_exp )) || { echo "FATAL: pred.list has $k_got of $k_exp traits" >&2; exit 4; }
# REGENIE exits 0 even when a .loco write hits the quota (0-byte file seen 2026-09-30)
awk '{print $2}' "$RG/step1_${RUN}_pred.list" | while read -r f; do
  [[ -s "$f" ]] || { echo "FATAL: empty $f" >&2; exit 4; }; done
rm -rf "$L0"
echo "=== done  $(date -u +%T)"
