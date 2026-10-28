#!/bin/bash
#SBATCH --job-name=geom_magma
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#
# D1 step 9: MAGMA gene analysis of the scan_geom GWAS, one phenotype per array
# task (phenotypes passed positionally; rule 10).  Settings copied from
#   eur    -> step13_magma_panel.sbatch:      1000G EUR, ncbi37.window35-10.union3 annotation
#   pooled -> step14_magma_pooled_genes.sbatch: in-sample ABCD reference + its 35/10 annotation
# Default snp-wise=mean; SNP/P/N columns built as there.  Writes only under scan_geom/.
#   sbatch --array=1-5 --output=.../%x_%A_%a.log magma_genes.sbatch <eur|pooled> <pheno> ...
set -euo pipefail
ARM="$1"; shift; PHENOS=("$@")
P="${PHENOS[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"
W="$REPO/genetic_analysis/work/results_70tab_hcp/scan_geom"
case "$ARM" in
  eur)    BFILE=/rds/user/rajd2/hpc-work/magma/reference_data/g1000_eur
          ANNOT="$REPO/genetic_analysis/work/magma_scz2025/ref/ncbi37.window35-10.union3.genes.annot"
          SS="$W/assoc_eur/$P.sumstats.tsv.gz"; OUT="$W/magma_eur/genes/$P" ;;
  pooled) REF="$REPO/genetic_analysis/work/magma_scz2025/ref_abcd"; BFILE="$REF/abcd_analysis"
          ANNOT="$REF/ncbi37.window35-10.abcd.genes.annot"
          SS="$W/assoc/$P.sumstats.tsv.gz"; OUT="$W/magma_pooled/genes/$P" ;;
  *) echo "ARM must be eur|pooled" >&2; exit 2 ;;
esac
export TMPDIR="$W/tmp"; mkdir -p "$TMPDIR" "$(dirname "$OUT")"
[[ -s "$SS" && -s "$BFILE.bed" && -s "$ANNOT" ]] || { echo "FATAL: $SS / $BFILE / $ANNOT missing" >&2; exit 2; }
if [[ -s "$OUT.genes.raw" ]]; then echo "exists $OUT"; exit 0; fi
echo "=== $ARM $P  $(date)"
zcat "$SS" | awk 'BEGIN{OFS="\t"} NR==1{print "SNP","P","N"; next}{print $1,$10,$7}' > "$OUT.snp_p.tsv"
"$MAGMA" --bfile "$BFILE" --gene-annot "$ANNOT" --pval "$OUT.snp_p.tsv" use=SNP,P ncol=N --out "$OUT" | tail -4
rm -f "$OUT.snp_p.tsv"; [[ -s "$OUT.genes.raw" ]] || { echo "FATAL: no genes.raw" >&2; exit 4; }
echo "  genes: $(( $(wc -l < "$OUT.genes.out") - 1 ))  $(date)"
