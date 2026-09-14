#!/bin/bash
#SBATCH --job-name=v2_eurkin
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#
# EUR-stratified arm, step 1: PC-AiR / PC-Relate computed WITHIN the European
# anchor set.
#
# WHY within-EUR rather than reusing the full-sample PC-AiR axes: inside a
# single continental group the leading full-sample PCs are near-constant, so
# they carry no information about the within-EUR structure that actually needs
# adjusting.  v1's EUR arm used EUR-specific PCs for the same reason
# (covar_quant_eurPC10.txt).  Kinship is re-estimated too, so PC-Relate's ISAF
# betas are fitted on the analysed population.
#
# The anchor set (5,656 genotyped) is used for structure/kinship; the null
# model then restricts to the 4,116 who are also phenotyped -- the same 4,116
# v1's EUR GWAS analysed.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/config.sh
GEUR="$ABCD_HPC_ROOT/genotype_eur/array_eur"
GDSEUR="$OUT_V2/gds_eur"
KINEUR="$OUT_V2/kinship_eur"
ensure_dirs "$(dirname "$GEUR")" "$GDSEUR" "$KINEUR"

if [[ ! -f "$GEUR.bed" ]]; then
  echo "=== subset array fileset to the EUR anchor  $(date) ==="
  run "$PLINK" --bfile "$GENO_ARRAY" --keep "$EUR_IDS" --make-bed --out "$GEUR"
fi
echo "EUR array fileset: $(wc -l < "$GEUR.fam") individuals, $(wc -l < "$GEUR.bim") SNPs"

if [[ ! -f "$GDSEUR/array.gds" ]]; then
  echo "=== PLINK -> GDS  $(date) ==="
  run "$RSCRIPT" "$RSCRIPT_DIR/01_make_gds.R" --bfile "$GEUR" --out "$GDSEUR/array.gds"
fi

echo "=== KING / PC-AiR / PC-Relate within EUR  $(date) ==="
run "$RSCRIPT" "$RSCRIPT_DIR/02_kinship.R" \
  --gds "$GDSEUR/array.gds" --out-dir "$KINEUR" \
  --maf "$KIN_MAF" --ld-r2 "$KIN_LD_R2" \
  --kin-thresh "$KIN_THRESH" --sparse-thresh "$SPARSE_KIN_THRESH" \
  --max-density "$SPARSE_KIN_MAX_DENSITY" --threads "$THREADS"

for f in kinship_sparse.rds pcair_pcs.tsv pcair_unrelated.txt kinship_summary.tsv; do
  [[ -s "$KINEUR/$f" ]] || { echo "FATAL: $KINEUR/$f missing/empty" >&2; exit 4; }
done
column -t -s $'\t' "$KINEUR/kinship_summary.tsv"
echo "=== done  $(date) ==="
