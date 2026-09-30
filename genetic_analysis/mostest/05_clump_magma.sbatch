#!/bin/bash
#SBATCH --job-name=m16_magma
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-6
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 16.6: loci and genes for each MOSTest scan.
#   tasks 1-3: slope, ct, slopeols (real);  tasks 4-6: the permuted scans
#   (null calibration of loci and gene-level p, gates G2b and G2c).
# LD reference = the ABCD analysis sample itself (step-14 ref_abcd/abcd_analysis,
# 8,596 children; rule 19: the MOSTest scan IS that sample).
#   clump:  PLINK 1.9, p1 5e-8 / p2 1e-6 / r2 0.1 / 500 kb, and a suggestive
#           pass at p1 1e-6.  Loci are counted from the index SNPs.
#   MAGMA:  gene analysis on P (the MOSTest p; MAGMA needs only unsigned SNP p),
#           35/10 kb NCBI37.3 annotation, per-SNP N; then the step-15b gene
#           sets (direction=greater) on the genes.raw.
set -eo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/genetic_analysis/mostest/paths.sh"
SCANS=(slope ct slopeols slope_perm ct_perm slopeols_perm)
SC="${SCANS[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"
SSF="$SS/mostest_${SC}.sumstats.tsv.gz"
[[ -s "$SSF" && -s "$REF_ABCD/abcd_analysis.bed" ]] || { echo "FATAL: $SSF or $REF_ABCD missing" >&2; exit 2; }
mkdir -p "$M16/clump" "$MG/genes" "$MG/sets"
echo "=== 16.6 $SC  $(date -u +%T)"
zcat "$SSF" | awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; print "SNP","P","N"; next}
                   {print $c["SNP"],$c["P"],$c["N"]}' > "$TMPDIR/${SC}.snp_p.tsv"
for thr in 5e-8 1e-6; do
  "$PLINK1" --bfile "$REF_ABCD/abcd_analysis" --clump "$TMPDIR/${SC}.snp_p.tsv" \
    --clump-snp-field SNP --clump-field P --clump-p1 "$thr" --clump-p2 1e-4 \
    --clump-r2 0.1 --clump-kb 500 --threads 8 --memory 50000 \
    --out "$M16/clump/${SC}_p${thr}" > /dev/null || true
  n=0; [[ -s "$M16/clump/${SC}_p${thr}.clumped" ]] && n=$(awk 'NF>0 && NR>1' "$M16/clump/${SC}_p${thr}.clumped" | wc -l)
  echo "  loci at p < $thr: $n"
done
OUT="$MG/genes/mostest_${SC}"
if [[ ! -s "$OUT.genes.raw" ]]; then
  "$MAGMA" --bfile "$REF_ABCD/abcd_analysis" --gene-annot "$REF_ABCD/ncbi37.window35-10.abcd.genes.annot" \
    --pval "$TMPDIR/${SC}.snp_p.tsv" use=SNP,P ncol=N --out "$OUT" | tail -4
fi
[[ -s "$OUT.genes.raw" ]] || { echo "FATAL: no genes.raw for $SC" >&2; exit 4; }
echo "  genes: $(( $(wc -l < "$OUT.genes.out") - 1 ))"
for ss in "scz:$SETS_SCZ" "mdd:$SETS_MDD" "top:$SETS_TOP" "wes:$SETS_WES"; do
  [[ -s "${ss#*:}" ]] || { echo "  MISSING set file ${ss#*:} (skipped)"; continue; }
  "$MAGMA" --gene-results "$OUT.genes.raw" --set-annot "${ss#*:}" col=1,2 \
    --model direction=greater --out "$MG/sets/${SC}__${ss%%:*}" > /dev/null
done
rm -f "$TMPDIR/${SC}.snp_p.tsv"
echo "=== done  $(date -u +%T)"
