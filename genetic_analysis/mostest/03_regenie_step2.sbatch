#!/bin/bash
#SBATCH --job-name=m16_rg2
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-154
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 16.3: REGENIE step 2 on the imputed PRS fileset (the same 7.07 M-SNP
# fileset that is the step-14 MAGMA LD reference, so SNP IDs join), then
# conversion of the k per-trait outputs to one z-matrix per (run, chr):
#   task = (run_idx - 1) * 22 + chr,  run_idx as in 02_regenie_step1.
# Permuted runs delete their per-trait REGENIE text after conversion (only z is
# needed); real runs keep it (per-region univariate sumstats, ~10 GB per run).
set -eo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/genetic_analysis/mostest/paths.sh"
t=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 )); RI=$(( t / 22 )); CHR=$(( t % 22 + 1 ))
RUN="${RUNS[$RI]}"; FAM="$(run_family "$RUN")"
SUF=""; run_is_perm "$RUN" && SUF="_perm"
PRED="$RG/step1_${RUN}_pred.list"
[[ -s "$PRED" ]] || { echo "FATAL: $PRED missing -- run 02_regenie_step1" >&2; exit 2; }
OUTP="$RG/step2_${RUN}_chr${CHR}"
if [[ -s "$ZM/${RUN}_chr${CHR}.npz" ]]; then echo "exists: $ZM/${RUN}_chr${CHR}.npz"; exit 0; fi
echo "=== 16.3 step 2  run=$RUN chr=$CHR  $(date -u +%T)"
"$REGENIE" --step 2 \
  --bed "$GENO_IMP_PRS" --chr "$CHR" --keep "$PH/keep_imp.txt" \
  --phenoFile "$PH/pheno_${FAM}${SUF}.txt" --covarFile "$PH/covar${SUF}.txt" \
  --catCovarList sex,site --maxCatLevels 30 \
  --qt --apply-rint --pred "$PRED" \
  --bsize 400 --minMAC 20 --gz \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --out "$OUTP"
"$PY" "$REPO/genetic_analysis/mostest/02_to_zmat.py" \
  --prefix "$OUTP" --traits "$PH/traits_${FAM}.txt" --chr "$CHR" --out "$ZM/${RUN}_chr${CHR}.npz"
if run_is_perm "$RUN"; then rm -f "${OUTP}"_*.regenie.gz; fi
echo "=== done  $(date -u +%T)"
