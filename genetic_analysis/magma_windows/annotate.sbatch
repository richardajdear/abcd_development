#!/bin/bash
#SBATCH --job-name=magmawin_annot
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm/%x_%A.log
#
# MAGMA window robustness, step 1: annotate the union-3 SNP locations under each window.
# w35_10 must reproduce work/magma_scz2025/ref/ncbi37.window35-10.union3.genes.annot.
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/magma_windows/common.sh
"$MAGMA" --version || true
for t in "${TAGS[@]}"; do
  [[ -s "$W/ref/ncbi37.$t.union3.genes.annot" ]] || "$MAGMA" --annotate window=$(win $t) --snp-loc "$UNION" --gene-loc "$GENELOC" --out "$W/ref/ncbi37.$t.union3" | tail -3
done
OLD="$REPO/genetic_analysis/work/magma_scz2025/ref/ncbi37.window35-10.union3.genes.annot"
if cmp "$OLD" "$W/ref/ncbi37.w35_10.union3.genes.annot"; then echo "REPRO_ANNOT: identical"
else echo "REPRO_ANNOT: DIFFERS" >&2; exit 4; fi
