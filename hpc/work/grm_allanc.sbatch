#!/bin/bash
#SBATCH --job-name=abcd_grm_anc
#SBATCH --account=VERTES-SL3-CPU   # override at submit: sbatch -A <acct>
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --array=1-22%12
#SBATCH --output=%x_%A_%a.log
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
#
# PATHS COME FROM hpc/config.sh.  The first version of this script hardcoded
# /home/rajd2/... and rds-genetics_hpc-Nl99R8pHODQ, which are 4.0-era locations
# on one account; every other script in this pipeline is path-modular and this
# one now matches.  Override GENO_ALLANC_DIR / GENO_ALLANC_TPL in
# config.local.sh when the release tree moves or is renamed.
set -euo pipefail

# Locate the repository.  Under sbatch, $0 and ${BASH_SOURCE[0]} name a COPY of
# this script in /var/spool/slurm/slurmd/jobNNN/, so the obvious
# `dirname "${BASH_SOURCE[0]}"/../..` resolves to /var/spool/slurm and every
# array task dies in two seconds on
#   /var/spool/slurm/hpc/config.sh: No such file or directory
# 02_reml.sbatch and 06_prs.sbatch already carry this scar; these work/ scripts
# did not, only because they had always been run with REPO_ROOT already
# exported into the submitting environment, which hid it.  Observed here on
# 2026-08-17, job 33789921, all 22 tasks.  Walk up from the submit directory
# (and $PWD) looking for hpc/config.sh, which is what actually identifies the
# repo, rather than trusting any one path to be the right one.
if [[ -z "${REPO_ROOT:-}" || ! -f "${REPO_ROOT:-}/hpc/config.sh" ]]; then
  REPO_ROOT=""
  for _d in "$(dirname "${BASH_SOURCE[0]}")/../.." "${SLURM_SUBMIT_DIR:-}" "$PWD"; do
    [[ -n "$_d" ]] || continue
    _p="$(cd "$_d" 2>/dev/null && pwd)" || continue
    while [[ -n "$_p" && "$_p" != "/" ]]; do
      if [[ -f "$_p/hpc/config.sh" ]]; then REPO_ROOT="$_p"; break 2; fi
      _p="$(dirname "$_p")"
    done
  done
fi
[[ -n "$REPO_ROOT" ]] || {
  echo "FATAL: cannot find hpc/config.sh above $PWD or ${SLURM_SUBMIT_DIR:-unset}." >&2
  echo "       Submit from inside the repo, or export REPO_ROOT." >&2
  exit 2
}
# shellcheck source=/dev/null
source "$REPO_ROOT/hpc/config.sh"

ensure_dirs "$GRM_ALLANC_DIR" "$LOG_DIR"

CHR=${SLURM_ARRAY_TASK_ID:?this script must run as an array job}
BFILE="$(geno_allanc_prefix "$CHR")"
OUT="$GRM_ALLANC_DIR"

# Integrity check BEFORE GCTA touches the fileset.  On 4.0 all 22 tasks died in
# seconds with "Unexpected PLINK 1 .bed file size" because every .bim had been
# filtered to ~78 % without regenerating the .bed.  This costs about a second
# and reports which of the two files is wrong, instead of leaving that to be
# reconstructed from GCTA's error afterwards.
bash "$REPO_ROOT/hpc/work/check_bfile_integrity.sh" "$BFILE" || {
  echo "FATAL: chr$CHR failed the .bed/.bim/.fam integrity check (see above)." >&2
  echo "       Do NOT work around this by filtering: if the .bim was filtered" >&2
  echo "       without rebuilding the .bed, the variant identities are unknown." >&2
  exit 2
}

echo "=== chr$CHR  $(date) ==="
echo "input : $BFILE  ($(wc -l < "$BFILE.bim") SNPs, $(wc -l < "$BFILE.fam") subjects)"

# --maf 0.01 matches what 03_gwas applied, so the GRM and the association scan
# describe the same variant set.  Note the MAF is computed across all ancestries
# pooled, which is the standard choice but not a neutral one -- see the caveat
# in README_HPC.md: a single GRM across ancestries assumes a common allele
# frequency and LD structure that does not actually hold.
#
# GRM_ALLANC_KEEP is how you stop that being a caveat and make it a measurement.
# GCTA estimates allele frequencies from the samples actually included, so
# passing a single-ancestry keep-list here builds a GRM whose frequencies come
# from THAT group -- which is a different object from the pooled GRM restricted
# to the same people afterwards, and the difference is exactly the assumption
# under scrutiny.  Restricting after the fact keeps the pooled frequencies and
# so keeps the problem.
#
# This matters more than it sounds.  Measured on the pooled GRM here: within
# stratum the off-diagonal averages +0.0402, between strata -0.0348.  A 0.05
# relatedness cutoff therefore lands inside the within-stratum bulk and prunes
# the non-European strata to 1.7-22 % while keeping 75.6 % of the European one.
keep_args=()
if [[ -n "${GRM_ALLANC_KEEP:-}" ]]; then
  [[ -f "$GRM_ALLANC_KEEP" ]] || { echo "FATAL: GRM_ALLANC_KEEP=$GRM_ALLANC_KEEP not found" >&2; exit 3; }
  echo "keep  : $GRM_ALLANC_KEEP ($(wc -l < "$GRM_ALLANC_KEEP") subjects; allele frequencies come from THESE)"
  keep_args+=(--keep "$GRM_ALLANC_KEEP")
fi

"$GCTA" \
  --bfile "$BFILE" \
  "${keep_args[@]+"${keep_args[@]}"}" \
  --maf 0.01 \
  --make-grm \
  --thread-num "${SLURM_CPUS_PER_TASK:-16}" \
  --out "$OUT/chr${CHR}"

echo "=== chr$CHR done  $(date) ==="
ls -la "$OUT/chr${CHR}".grm.* | awk '{print "  ", $5, $9}'
