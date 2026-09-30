#!/bin/bash
#SBATCH --job-name=d7_annot
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=directions/d7_partitioned_prs/work/logs/%x_%A_%a.log
#
# D7 step 4b (secondary): refit SBayesRC for the four disorder arms with the D7
# primary sets appended to the baseline 2.2 annotation.  Same GCTB, LD eigen
# reference, imputed summary statistics, chain (3000/1000) and seed as the
# Figure-1 fits, so the ONLY change is the annotation.
#   step 4a first:  $D7PY code/04_build_annot.py ... --out $WORK/annot_d7.txt
#   sbatch --array=1-4 directions/d7_partitioned_prs/d7_annot.sbatch
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/directions/d7_partitioned_prs/config.sh
cd "$REPO"
AN="$WORK/annot_d7.txt"; [[ -s "$AN" ]] || { echo "FATAL: $AN missing (step 4a)" >&2; exit 2; }
T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
spec="${ARMS[$T]}"; name=${spec%%|*}; r=${spec#*|}; wf=${r%%|*}; cell=${spec##*|}
[[ "$name" == EA ]] && { echo "EA is not refitted (control trait)"; exit 0; }
SRC="$(dirname "$wf")"; IMP="$SRC/${name}_imp.imputed.ma"
[[ -s "$IMP" ]] || { echo "FATAL: $IMP missing (the Figure-1 fit's imputed sumstats)" >&2; exit 3; }
OUTD="$WORK/sbrc_annot/$name"; mkdir -p "$OUTD"; TAG="${name//_/}"
echo "=== $name  $(date +%T)"
"$GCTB" --sbayes RC --ldm-eigen "$SBRC_EIGEN" --gwas-summary "$IMP" --annot "$AN" \
  --chain-length 3000 --burn-in 1000 --seed 2026 --thread "${SLURM_CPUS_PER_TASK:-8}" \
  --out "$OUTD/$name" 2>&1 | tail -30
[[ -s "$OUTD/$name.snpRes" ]] || { echo "FATAL: no snpRes" >&2; exit 4; }
awk 'NR>1{print $2"\t"$5"\t"$8}' "$OUTD/$name.snpRes" > "$OUTD/${TAG}_sbrcannot.weights"
"$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_sbrcannot.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_sbrcannot" > /dev/null
ls -la "$OUTD"
echo "arm string for steps 2-3:  ${name}_annot|$OUTD/${TAG}_sbrcannot.weights|$OUTD/score_${TAG}_sbrcannot.profile|$cell"
echo "=== done $(date +%T)"
