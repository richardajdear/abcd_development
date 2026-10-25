#!/bin/bash
#SBATCH --job-name=magmawin_genes
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-19
#SBATCH --output=slurm/%x_%A_%a.log
#
# MAGMA window robustness, step 2: gene analysis, one task per window x trait (default snp-wise=mean).
# Optional: MULTI=1 with --array=0-3 runs w35_10 with --gene-model multi=snp-wise (tag w35_10_multi).
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/magma_windows/common.sh
i=$SLURM_ARRAY_TASK_ID
if [[ "${MULTI:-0}" == 1 ]]; then t=w35_10_multi; atag=w35_10; tr=${TRAITS[$i]}; EXTRA=(--gene-model multi=snp-wise)
else t=${TAGS[$((i / 4))]}; atag=$t; tr=${TRAITS[$((i % 4))]}; EXTRA=(); fi
OUT="$W/genes/$t"; mkdir -p "$OUT"
[[ -s "$OUT/$tr.genes.raw" ]] && { echo "skip $t $tr"; exit 0; }
ANNOT="$W/ref/ncbi37.$atag.union3.genes.annot"; [[ -s "$ANNOT" ]] || { echo "FATAL: $ANNOT missing" >&2; exit 2; }
src="${SRC[$tr]}"; [[ -s "$src" ]] || { echo "FATAL: $src missing" >&2; exit 3; }
echo "=== $t $tr  $(date)"
if [[ "$tr" == ABCD_* ]]; then
  P="$W/tmp/${t}_$tr.snp_p.tsv"
  zcat "$src" | awk 'BEGIN{OFS="\t"} NR==1{if($1!="SNP"||$7!="N"||$10!="P"){print "bad header" > "/dev/stderr"; exit 5} print "SNP","P","N"; next}{print $1,$10,$7}' > "$P"
  "$MAGMA" --bfile "$REF/g1000_eur" --gene-annot "$ANNOT" --pval "$P" use=SNP,P ncol=N "${EXTRA[@]}" --out "$OUT/$tr" | tail -4
  rm -f "$P"
else
  "$MAGMA" --bfile "$REF/g1000_eur" --gene-annot "$ANNOT" --pval "$src" use=SNP,p ncol=N "${EXTRA[@]}" --out "$OUT/$tr" | tail -4
fi
echo "=== done $(date)"
