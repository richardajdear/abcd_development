#!/bin/bash
#SBATCH --job-name=magma_pooled_genes
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-8
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 14: MAGMA gene analysis of the POOLED-arm phenotype GWAS (8,596 children)
# with the ABCD analysis sample as LD reference (step14_magma_pooled_prep).
# Cells: {dk,hcp} x {global_slope, baseline_thickness} x {perregion (assoc/), 1lmm (scan_1lmm/assoc/)}.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"
REF="$REPO/genetic_analysis/work/magma_scz2025/ref_abcd"
export TMPDIR="$REPO/genetic_analysis/work/magma_scz2025/tmp"
CELLS=()
for parc in dk hcp; do s=""; [[ $parc == hcp ]] && s="_hcp"; R="$REPO/genetic_analysis/work/results_70tab$s"
  for p in global_slope baseline_thickness; do
    CELLS+=("$parc:$p:$R/assoc/$p.sumstats.tsv.gz:$R/magma_pooled/genes/$p")
    CELLS+=("$parc:${p}_1lmm:$R/scan_1lmm/assoc/${p}_1lmm.sumstats.tsv.gz:$R/magma_pooled/genes/${p}_1lmm")
  done; done
SPEC="${CELLS[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"; PARC="${SPEC%%:*}"; R1="${SPEC#*:}"; P="${R1%%:*}"; R2="${R1#*:}"; SS="${R2%%:*}"; OUT="${R2#*:}"
[[ -s "$SS" && -s "$REF/abcd_analysis.bed" ]] || { echo "FATAL: $SS or reference missing" >&2; exit 2; }
mkdir -p "$(dirname "$OUT")"
if [[ -s "$OUT.genes.raw" ]]; then echo "exists $OUT"; exit 0; fi
echo "=== $PARC $P  $(date)"
zcat "$SS" | awk 'BEGIN{OFS="\t"} NR==1{print "SNP","P","N"; next}{print $1,$10,$7}' > "$OUT.snp_p.tsv"
"$MAGMA" --bfile "$REF/abcd_analysis" --gene-annot "$REF/ncbi37.window35-10.abcd.genes.annot" --pval "$OUT.snp_p.tsv" use=SNP,P ncol=N --out "$OUT" | tail -4
rm -f "$OUT.snp_p.tsv"; [[ -s "$OUT.genes.raw" ]] || { echo "FATAL: no genes.raw" >&2; exit 4; }
echo "  genes: $(( $(wc -l < "$OUT.genes.out") - 1 ))  $(date)"
