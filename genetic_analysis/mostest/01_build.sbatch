#!/bin/bash
#SBATCH --job-name=m16_build
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=slurm/%x_%j.log
#
# Step 16.1: regional phenotypes (01_build_pheno.py) + the REGENIE step-1
# genotype set (array SNPs, QC'd, restricted to the analysis set, IDs rewritten
# to the imputed .fam spelling so REGENIE steps 1 and 2 see the same sample).
#   sbatch genetic_analysis/mostest/01_build.sbatch      (from the repo root)
set -eo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/genetic_analysis/mostest/paths.sh"
export PYTHONPATH="$REPO/src:${PYTHONPATH:-}"
for f in "$DK_RUN/fits/blups.parquet" "$DK_RUN/model_table.parquet" "$REF_1LMM/phenotypes_gcta.txt" \
         "$GENO_IMP_PRS.fam" "$GENO_ARRAY.bed"; do
  [[ -s "$f" ]] || { echo "FATAL: missing $f" >&2; exit 2; }
done

echo "=== 16.1a phenotypes  $(date -u +%T)"
"$PY" "$REPO/genetic_analysis/mostest/01_build_pheno.py" \
  --run "$DK_RUN" --ref "$REF_1LMM" --geno-imp "$GENO_IMP_PRS" --geno-array "$GENO_ARRAY" --out "$PH"

echo "=== 16.1b step-1 genotypes  $(date -u +%T)"
# Autosomal array SNPs, MAF >= 1 %, call rate >= 98 %, HWE p >= 1e-15 (REGENIE
# guidance: ~500k common directly-genotyped SNPs; no LD pruning needed).
"$PLINK1" --bfile "$GENO_ARRAY" --keep "$PH/keep_array.txt" --chr 1-22 \
  --maf 0.01 --geno 0.02 --hwe 1e-15 --make-bed --threads 8 --memory 40000 \
  --out "$GENO_DIR/step1_tmp" > "$GENO_DIR/step1_tmp.plinklog"
"$PLINK1" --bfile "$GENO_DIR/step1_tmp" --update-ids "$PH/array_idmap.txt" --make-bed \
  --threads 8 --memory 40000 --out "$GENO_DIR/step1_array" > "$GENO_DIR/step1_array.plinklog"
rm -f "$GENO_DIR"/step1_tmp.{bed,bim,fam}
n_fam=$(wc -l < "$GENO_DIR/step1_array.fam"); n_snp=$(wc -l < "$GENO_DIR/step1_array.bim")
n_keep=$(wc -l < "$PH/keep_imp.txt")
echo "step-1 genotypes: $n_fam children, $n_snp SNPs (analysis set $n_keep)"
(( n_fam == n_keep )) || { echo "FATAL: step-1 fam has $n_fam of $n_keep children" >&2; exit 3; }
(( n_snp > 250000 )) || { echo "FATAL: only $n_snp array SNPs after QC" >&2; exit 3; }
# every step-1 ID must appear verbatim in the imputed .fam (rule 1)
awk 'NR==FNR{k[$1" "$2]=1; next} !(($1" "$2) in k){bad++} END{exit bad>0}' \
  "$GENO_IMP_PRS.fam" "$GENO_DIR/step1_array.fam" \
  || { echo "FATAL: step-1 IDs do not match the imputed .fam spelling" >&2; exit 3; }
printf 'n_children\t%s\nn_step1_snps\t%s\n' "$n_fam" "$n_snp" > "$GENO_DIR/step1_counts.tsv"
echo "=== done  $(date -u +%T)"
