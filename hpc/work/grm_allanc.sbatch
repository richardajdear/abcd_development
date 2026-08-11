#!/bin/bash
#SBATCH --job-name=abcd_grm_anc
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --array=1-22%12
#SBATCH --output=/home/rajd2/rds/hpc-work/abcd_development/hpc/work/logs/%x_%A_%a.log
#
# Step 1 of 2: per-chromosome GRMs over the ALL-ANCESTRY genotypes.
#
# Why this exists: the analysis so far runs on `abcd_eur`, 5,678 European-only
# subjects, of whom 4,126 are also phenotyped.  The all-ancestry per-chromosome
# files hold 10,072 subjects, 7,111 of them phenotyped -- **+72 % N** (verified
# by joining on the 8-character NDAR token, not assumed).  N is the binding
# constraint on every result in README_HPC.md §12-§14, so this is the one change
# that alters what the pipeline can detect.
#
# Per chromosome rather than one merged fileset: the 22 autosome files total
# ~60 GB and merging them would cost disk and a long serial pass, whereas GCTA
# can build GRMs independently and combine them with --mgrm, which is both
# parallel and exact (the combined GRM is the SNP-count-weighted mean).
#
# %12 caps concurrency at 12 x 16 = 192 CPUs, inside the 448-CPU per-user limit.
# No --mem: icelake caps at 3370 MB/CPU, so 16 CPUs already implies ~53 GB, well
# clear of the ~0.8 GB the 10,072 x 10,072 matrix needs plus genotype blocks.
#
# X is excluded deliberately: a GRM is conventionally autosomal, and the X
# dosage convention differs between males and females.
set -euo pipefail

W=/home/rajd2/rds/hpc-work/abcd_development/hpc/work
SRC=/home/rajd2/rds/rds-genetics_hpc-Nl99R8pHODQ/ABCD
OUT=/rds/user/rajd2/hpc-work/ABCD/grm_allanc
mkdir -p "$OUT"

CHR=${SLURM_ARRAY_TASK_ID:?this script must run as an array job}
BFILE="$SRC/ABCD_chr${CHR}_hg19"

[[ -f "$BFILE.bed" ]] || { echo "FATAL: $BFILE.bed not found" >&2; exit 2; }

echo "=== chr$CHR  $(date) ==="
echo "input : $BFILE  ($(wc -l < "$BFILE.bim") SNPs, $(wc -l < "$BFILE.fam") subjects)"

# --maf 0.01 matches what 03_gwas applied, so the GRM and the association scan
# describe the same variant set.  Note the MAF is computed across all ancestries
# pooled, which is the standard choice but not a neutral one -- see the caveat
# in README_HPC.md: a single GRM across ancestries assumes a common allele
# frequency and LD structure that does not actually hold.
"$W/bin/gcta64" \
  --bfile "$BFILE" \
  --maf 0.01 \
  --make-grm \
  --thread-num "${SLURM_CPUS_PER_TASK:-16}" \
  --out "$OUT/chr${CHR}"

echo "=== chr$CHR done  $(date) ==="
ls -la "$OUT/chr${CHR}".grm.* | awk '{print "  ", $5, $9}'
