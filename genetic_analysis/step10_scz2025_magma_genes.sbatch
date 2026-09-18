#!/bin/bash
#SBATCH --job-name=scz25_magma_genes
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-5
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 10, disorder-side gene analyses (one per array task), all on the union
# annotation from step10_scz2025_magma_prep.sbatch, per-SNP N from the .ma:
#   1 SCZ25_EUR        2025 EUR GWAS          x g1000_eur   } the three inputs
#   2 SCZ25_AFR        2025 AFR GWAS          x g1000_afr   } to the --meta
#   3 SCZ25_EAS        2025 EAS GWAS          x g1000_eas   } gene meta-analysis
#   4 SCZ25_METAnaive  2025 AFR+EUR+EAS meta  x g1000_eur   "normal MAGMA" on the
#                      meta file, LD-mismatched; kept to measure what the
#                      mismatch does against the proper meta
#   5 PGC3_EUR         PGC3 european          x g1000_eur   the like-for-like
#                      PGC3 comparator for the EUR arm (step 7 only ran the
#                      multi-ancestry PGC3 primary, against g1000_eur)
# Gene-gene correlations are kept (no --genes-only): the gene-set and
# gene-property tests downstream need them.
set -euo pipefail
MG=/rds/user/rajd2/hpc-work/magma; REF="$MG/reference_data"
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"
ROOT="${SCZ25_MAGMA_ROOT:-$REPO/genetic_analysis/work/magma_scz2025}"
ANNOT="$ROOT/ref/ncbi37.window35-10.union3.genes.annot"
GW25="$REPO/genetic_analysis/work/scores_scz2025/gwas"
GW3="$REPO/legacy/hpc_v2/work/results_v2/prs_final/gwas"
export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR" "$ROOT/genes"
CELLS=("SCZ25_EUR:$GW25/SCZ25_EUR.ma:eur" "SCZ25_AFR:$GW25/SCZ25_AFR.ma:afr" "SCZ25_EAS:$GW25/SCZ25_EAS.ma:eas"
       "SCZ25_METAnaive:$GW25/SCZ25_META.ma:eur" "PGC3_EUR:$GW3/SCZ_eur.ma:eur")
T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 )); SPEC="${CELLS[$T]}"
NAME="${SPEC%%:*}"; REST="${SPEC#*:}"; MA="${REST%%:*}"; POP="${REST##*:}"
[[ -s "$MA" && -s "$ANNOT" && -s "$REF/g1000_$POP.bed" ]] || { echo "FATAL: input missing ($MA / $ANNOT / g1000_$POP)" >&2; exit 2; }
echo "=== $NAME : $MA x g1000_$POP  $(date)"
"$MAGMA" --bfile "$REF/g1000_$POP" --gene-annot "$ANNOT" --pval "$MA" use=SNP,p ncol=N --out "$ROOT/genes/$NAME" | tail -8
[[ -s "$ROOT/genes/$NAME.genes.raw" ]] || { echo "FATAL: no genes.raw for $NAME" >&2; exit 4; }
echo "  genes: $(( $(wc -l < "$ROOT/genes/$NAME.genes.out") - 1 ))"
echo "=== $NAME done  $(date)"
