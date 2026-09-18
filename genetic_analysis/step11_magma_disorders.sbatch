#!/bin/bash
#SBATCH --job-name=magma_disorders
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-6
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 11: the rest of the disorder panel through MAGMA on the step-10 union
# annotation, so every trait in the PRS panel has a gene-level result built the
# same way (35/10 kb, 1000G EUR panel, per-SNP N from the .ma):
#   1 MDD_EUR      PGC MDD2025 eur     x g1000_eur   the LD-matched MDD result
#   2 MDD_div      PGC MDD2025 div     x g1000_eur   step 7's choice (multi-ancestry
#                                                     file on a EUR panel): comparator
#   3 ASD          SPARK+iPSYCH+PGC    x g1000_eur
#   4 ALZ          Wightman (incl. APOE) x g1000_eur
#   5 ALZ_noAPOE   Wightman minus APOE
#   6 EA           Okbay 2016 EduYears
set -euo pipefail
MG=/rds/user/rajd2/hpc-work/magma; REF="$MG/reference_data"
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"
ROOT="${SCZ25_MAGMA_ROOT:-$REPO/genetic_analysis/work/magma_scz2025}"
ANNOT="$ROOT/ref/ncbi37.window35-10.union3.genes.annot"
GW3="$REPO/legacy/hpc_v2/work/results_v2/prs_final/gwas"
export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR" "$ROOT/genes"
CELLS=("MDD_EUR:$GW3/MDD_eur.ma" "MDD_div:$GW3/MDD_pooled.ma" "ASD:$GW3/ASD.ma"
       "ALZ:$GW3/ALZ.ma" "ALZ_noAPOE:$GW3/ALZ_noAPOE.ma" "EA:$GW3/EA.ma")
T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 )); SPEC="${CELLS[$T]}"; NAME="${SPEC%%:*}"; MA="${SPEC#*:}"
[[ -s "$MA" && -s "$ANNOT" ]] || { echo "FATAL: $MA or $ANNOT missing" >&2; exit 2; }
if [[ -s "$ROOT/genes/$NAME.genes.raw" ]]; then echo "exists: $NAME"; exit 0; fi
echo "=== $NAME : $MA x g1000_eur  $(date)"
"$MAGMA" --bfile "$REF/g1000_eur" --gene-annot "$ANNOT" --pval "$MA" use=SNP,p ncol=N --out "$ROOT/genes/$NAME" | tail -6
[[ -s "$ROOT/genes/$NAME.genes.raw" ]] || { echo "FATAL: no genes.raw for $NAME" >&2; exit 4; }
echo "  genes: $(( $(wc -l < "$ROOT/genes/$NAME.genes.out") - 1 ))"
echo "=== $NAME done  $(date)"
