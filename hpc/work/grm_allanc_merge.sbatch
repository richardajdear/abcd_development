#!/bin/bash
#SBATCH --job-name=abcd_grm_mrg
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --output=/home/rajd2/rds/hpc-work/abcd_development/hpc/work/logs/%x_%j.log
#
# Step 2 of 2: combine the 22 per-chromosome GRMs and derive everything the
# pipeline consumes, so that pointing config.local.sh at the all-ancestry set is
# a one-line change rather than another round of jobs.
#
#   abcd_all              genome-wide GRM         (02_reml, and --pca below)
#   abcd_all_sp           sparse GRM at 0.05      (03_gwas / fastGWA)
#   abcd_all.unrel        relatedness-pruned GRM  (02_reml's unrelated subset)
#   abcd_all_pca.eigenvec 20 in-sample ancestry PCs
#
# The PCA is not an extra: it is the fix for the defect this job was scoped
# alongside.  The covariates currently in use are the release's
# `ab_g_stc__gen_pc__01..10`, computed on the full multi-ancestry cohort, and
# within the EUR analysis sample PC1 has mean 0.00631 against SD 0.000523 --
# a mean 12x its own spread, which is what a PC that separates ancestry groups
# looks like when you restrict to one of them.  Such a PC absorbs almost no
# within-sample structure.  For an all-ancestry analysis, in-sample PCs are not
# merely preferable but required.
#
# himem: the merge holds the 10,072 x 10,072 matrix plus per-chromosome inputs;
# 16 CPUs x 6760 MB = ~108 GB leaves ample headroom.
set -euo pipefail

W=/home/rajd2/rds/hpc-work/abcd_development/hpc/work
OUT=/rds/user/rajd2/hpc-work/ABCD/grm_allanc
GCTA="$W/bin/gcta64"
T="${SLURM_CPUS_PER_TASK:-16}"

# Fail loudly and early if step 1 did not finish: --mgrm on a short list would
# silently build a GRM from fewer chromosomes, which is exactly the class of
# quiet-wrong result this pipeline has produced twice already (README_HPC §6d,
# §12b).
LIST="$OUT/grm_chr_list.txt"
: > "$LIST"
missing=0
for c in $(seq 1 22); do
  if [[ -f "$OUT/chr${c}.grm.bin" && -f "$OUT/chr${c}.grm.id" ]]; then
    echo "$OUT/chr${c}" >> "$LIST"
  else
    echo "MISSING: chr${c}" >&2; missing=1
  fi
done
n_listed=$(wc -l < "$LIST")
echo "GRMs to merge: $n_listed of 22"
[[ "$missing" == "0" ]] || { echo "FATAL: refusing to merge an incomplete set" >&2; exit 2; }

echo "=== merge  $(date) ==="
"$GCTA" --mgrm "$LIST" --make-grm --thread-num "$T" --out "$OUT/abcd_all"
echo "subjects in merged GRM: $(wc -l < "$OUT/abcd_all.grm.id")"

echo "=== sparse GRM (fastGWA)  $(date) ==="
"$GCTA" --grm "$OUT/abcd_all" --make-bK-sparse 0.05 \
        --thread-num "$T" --out "$OUT/abcd_all_sp"

echo "=== relatedness-pruned GRM (REML)  $(date) ==="
"$GCTA" --grm "$OUT/abcd_all" --grm-cutoff 0.05 --make-grm \
        --thread-num "$T" --out "$OUT/abcd_all.unrel"
echo "unrelated at 0.05: $(wc -l < "$OUT/abcd_all.unrel.grm.id") of $(wc -l < "$OUT/abcd_all.grm.id")"

echo "=== in-sample ancestry PCs  $(date) ==="
"$GCTA" --grm "$OUT/abcd_all" --pca 20 --thread-num "$T" --out "$OUT/abcd_all_pca"

echo
echo "=== variance explained by the first 10 in-sample PCs ==="
awk '{s+=$1; e[NR]=$1} END {
        for (i = 1; i <= 10 && i <= NR; i++)
          printf "  PC%-2d  %8.4f  (%5.2f%% of total)\n", i, e[i], 100*e[i]/s
      }' "$OUT/abcd_all_pca.eigenval"

echo
echo "=== done  $(date) ==="
ls -la "$OUT"/abcd_all*.grm.* "$OUT"/abcd_all_pca.* 2>/dev/null | awk '{print "  ", $5, $9}'
echo
echo "To switch the pipeline over, set in hpc/config.local.sh:"
echo "  GENO=<all-ancestry bfile>   # needs a merged fileset; per-chr for GWAS"
echo "  GRM=$OUT/abcd_all"
echo "  GRM_SPARSE=$OUT/abcd_all_sp"
echo "  GRM_UNREL=$OUT/abcd_all.unrel"
echo "and replace the release PCs in covar_quant.txt with abcd_all_pca.eigenvec."
