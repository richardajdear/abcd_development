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
if [[ -s "$RG/step1_${RUN}_pred.list" ]]; then echo "exists: $RG/step1_${RUN}_pred.list"; exit 0; fi
echo "=== 16.2 step 1  run=$RUN  $(date -u +%T)"
"$REGENIE" --step 1 \
  --bed "$GENO_DIR/step1_array" --keep "$PH/keep_imp.txt" \
  --phenoFile "$PHENO" --covarFile "$COVAR" \
  --catCovarList sex,site --maxCatLevels 30 \
  --qt --apply-rint --bsize 1000 \
  --lowmem --lowmem-prefix "$TMPDIR/l0_${RUN}" \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  --out "$RG/step1_${RUN}"
k_exp=$(wc -l < "$PH/traits_${FAM}.txt"); k_got=$(wc -l < "$RG/step1_${RUN}_pred.list")
(( k_got == k_exp )) || { echo "FATAL: pred.list has $k_got of $k_exp traits" >&2; exit 4; }
rm -f "$TMPDIR/l0_${RUN}"*
echo "=== done  $(date -u +%T)"
