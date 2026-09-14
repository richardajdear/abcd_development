#!/bin/bash
#SBATCH --job-name=abcd_split70
#SBATCH --account=VERTES-SL3-CPU   # override at submit: sbatch -A <acct>
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --array=1-22
#SBATCH --output=%x_%A_%a.log
#
# Step 0 of the all-ancestry GRM: split the 7.0 curated genotype fileset into
# 22 per-chromosome PLINK filesets, because grm_allanc.sbatch is an array over
# chromosomes and the 7.0 release ships ONE merged fileset, not 22.
#
# WHY THE SMOKESCREEN ARRAY DATA AND NOT THE IMPUTED VCFs
# -------------------------------------------------------
# The 7.0 tree under Data_Genetics/genotype_microarray has three products:
#
#   smokescreen/merged_chroms.{bed,bim,fam}   1.5 GB   11,670 x 515,228   <- this
#   imputed/chr{1..22}_dose.vcf.gz            ~890 GB  TOPMed r3 dosages
#   genesis/pcrelate_relatedness.*            GENESIS PC-Relate kinship
#
# The imputed set is ~890 GB of VCF; converting it to PLINK would cost days of
# wall-clock and several hundred GB before a single GRM was built.  A GRM does
# not need imputed density -- GREML's SNP-heritability is estimated from the
# LD-tagging of ~500k common array SNPs, which is the standard input for it.
# What this DOES mean is recorded in README_HPC.md: the multi-ancestry GRM is
# built on a different variant set from the EUR GRM (13.7M imputed, MAF 0.001),
# so a pooled-vs-EUR h2 difference confounds ancestry with SNP set.  The fix is
# the EUR-stratified REML on THIS GRM, which holds the variant set fixed.
#
# The genesis/ PC-Relate matrix was not used either: it is a kinship estimator
# with a different definition (PC-adjusted, meant for association-model
# relatedness control), not a GCTA additive GRM, and substituting it would make
# every h2 here incomparable with the EUR numbers.  It is a useful cross-check
# on the relatedness pruning, not an input.
#
# By default no MAF or missingness filter is applied here: GCTA applies
# --maf 0.01 in grm_allanc.sbatch, so the filter lives in exactly one place and
# the per-chromosome files stay a faithful split of the source.
#
# GENO_ALLANC_KEEP and GENO_ALLANC_GENO exist because that default turned out to
# be inadequate for the 7.0 array data, and the reason is worth stating rather
# than burying in a flag.  Measured on chr1 (work/qc_missing_allanc.sbatch):
# per-sample call rate is not noise, it is BATCH -- 87.6 % of the variance in
# per-sample F_MISS is explained by the ten Smokescreen batches, and four of
# them (BATCH_5_Saliva, BATCH_3_WB, BATCH_2_WB, BATCH_6; 511 subjects) sit at
# 7-10 % missing while the other six sit at 0.4-1.7 %.
#
# What that does to a GRM -- checked against chr1.grm.N.bin rather than assumed,
# because the obvious guess is wrong.  GCTA does NOT mean-impute: it estimates
# each pair from the SNPs non-missing in BOTH members and records the count.
# On chr1 that count ranges 29,324-37,800 and averages 35,256 for the
# high-missing subjects against 37,683 for the clean ones.
#
# So the damage is not a shared imputation bias, it is ~6 % fewer SNPs behind
# every relatedness estimate involving a bad-batch subject.  Noisier relatedness
# attenuates GREML h2, and here the noise is structured BY BATCH rather than
# spread at random.  That is a milder problem than imputation would have been --
# milder enough that it justifies a sensitivity track, not distrust of the
# primary one.  Which is why both are built and compared.
#
# GENO_ALLANC_KEEP must be a GENOME-WIDE sample list, not a per-chromosome
# --mind: --mind applied inside this array would drop a DIFFERENT sample set on
# each chromosome, leaving 22 filesets with 22 different .fam files, and GCTA
# --mgrm would then combine GRMs over inconsistent samples.  work/
# qc_keeplist_allanc.sbatch computes the list once, genome-wide, for this reason.
#
# Both tracks are built and compared rather than one being chosen up front; see
# README_HPC.md.  Neither filter is applied unless the variable is set, so the
# unfiltered track is still what you get by default.
#
# Chromosomes 23/25/26 (X, PAR, MT) are present in the source .bim and are
# dropped here: a GRM is conventionally autosomal.
#
# PATHS COME FROM hpc/config.sh.  GENO_ALLANC_SRC is the merged 7.0 prefix.
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

: "${GENO_ALLANC_SRC:?set GENO_ALLANC_SRC in hpc/config.local.sh to the merged 7.0 bfile prefix}"

ensure_dirs "$GENO_ALLANC_DIR" "$LOG_DIR"

CHR=${SLURM_ARRAY_TASK_ID:?this script must run as an array job}
OUT="$(geno_allanc_prefix "$CHR")"

# Gate on the SOURCE before splitting.  A .bed/.bim mismatch in the source
# would propagate silently into all 22 outputs -- plink would happily emit
# filesets whose variant identities are wrong -- so check it 22 times rather
# than trust that task 1 checked it.
bash "$REPO_ROOT/hpc/work/check_bfile_integrity.sh" "$GENO_ALLANC_SRC" || {
  echo "FATAL: the 7.0 source fileset failed the integrity check (see above)." >&2
  exit 2
}

echo "=== split chr$CHR  $(date) ==="
echo "source: $GENO_ALLANC_SRC"
echo "out   : $OUT"

filters=()
if [[ -n "${GENO_ALLANC_KEEP:-}" ]]; then
  [[ -f "$GENO_ALLANC_KEEP" ]] || { echo "FATAL: GENO_ALLANC_KEEP=$GENO_ALLANC_KEEP not found" >&2; exit 5; }
  echo "keep  : $GENO_ALLANC_KEEP ($(wc -l < "$GENO_ALLANC_KEEP") subjects)"
  filters+=(--keep "$GENO_ALLANC_KEEP")
fi
if [[ -n "${GENO_ALLANC_GENO:-}" ]]; then
  echo "geno  : --geno $GENO_ALLANC_GENO"
  filters+=(--geno "$GENO_ALLANC_GENO")
fi

"$PLINK" \
  --bfile "$GENO_ALLANC_SRC" \
  --chr "$CHR" \
  "${filters[@]+"${filters[@]}"}" \
  --make-bed \
  --threads "${SLURM_CPUS_PER_TASK:-4}" \
  --out "$OUT"

# Verify the output too.  This is the cheap half of the argument for the gate:
# it costs a second and it is the only thing standing between a truncated
# write and a GRM built from garbage.
bash "$REPO_ROOT/hpc/work/check_bfile_integrity.sh" "$OUT" || {
  echo "FATAL: chr$CHR output failed its own integrity check." >&2
  exit 3
}

# The subject count must be exactly what was asked for -- the source count when
# no --keep was given, the keep-list count when one was.  plink keeps all
# samples unless told otherwise, so any other number means a filter arrived
# that nobody asked for, and a GRM built from 22 filesets with inconsistent
# .fam files is the quiet-wrong failure this whole gate exists to prevent.
if [[ -n "${GENO_ALLANC_KEEP:-}" ]]; then
  n_want=$(wc -l < "$GENO_ALLANC_KEEP")
else
  n_want=$(wc -l < "$GENO_ALLANC_SRC.fam")
fi
n_out=$(wc -l < "$OUT.fam")
[[ "$n_want" == "$n_out" ]] || {
  echo "FATAL: chr$CHR has $n_out subjects, expected $n_want" >&2
  exit 4
}

echo "=== chr$CHR done  $(date):  $(wc -l < "$OUT.bim") variants, $n_out subjects ==="
