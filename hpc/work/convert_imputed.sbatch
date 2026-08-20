#!/bin/bash
#SBATCH --job-name=abcd_impconv
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.log
#
# Phase 2: convert one imputed TOPMed VCF to a filtered PLINK fileset.
#
# The filters are applied in the SAME pass as the conversion.  The unfiltered
# form -- ~300M variants x 11,670 subjects, roughly 875 GB as PLINK 1 -- must
# never touch disk; there is 184 GB free.  plink2 reads the gzipped VCF once,
# drops variants on the way through, and writes only survivors.
#
# TWO FILESETS, AND THE REASON IS README section 8.6 FOR THE FOURTH TIME
# ---------------------------------------------------------------------
# The MAF in the INFO field is the POOLED cohort frequency.  Measured on chr22
# (job 34039199, 1-in-200 sample): of 1,117 variants common within at least one
# ancestry stratum, a pooled MAF>=0.01 filter keeps 700 and DISCARDS 417 -- and
# 95% of those discards are non-European variation (78.7% commonest in cluster1
# alone).  A pooled threshold is not ancestry-neutral; it is a European filter
# wearing a neutral name.
#
#   union   --mac $MAC   keeps anything common in ANY stratum.  A variant at
#                        MAF 0.01 in the smallest stratum (cluster3, n=523)
#                        carries >=10 minor alleles, so MAC 10 is the union
#                        rule expressed in a form plink2 can apply in one pass.
#   pooled  --maf 0.01   the conventional fileset, built for comparison so the
#                        cost of the convention is a measured number.
#
# The union fileset forecloses nothing, which is why it is preferred: GCTA
# recomputes allele frequencies on whichever samples --keep leaves (verified --
# the same chr1 fileset yields 37,807 SNPs over 11,670 subjects and 32,099 over
# the 5,656 EUR subset).  So a pooled analysis on the union fileset still gets a
# pooled MAF filter, while a stratified one automatically gets that stratum's.
# Choosing the union here costs disk and decides nothing downstream.
#
# R2>=0.8 is Minimac's imputation r-squared, conventional.  It costs ~10% of
# stratum-common variants (89.9% pass).  Quality is NOT the binding constraint
# here and is not ancestry-neutral either -- cluster3 runs 17.0% substantially
# uncertain against EURlike's 4.5% -- but cluster3 is 523 subjects and already
# below the stratified-REML floor, so this is recorded, not corrected.
#
# --snps-only --max-alleles 2 : a GRM is defined over biallelic SNPs; indels
# and multiallelics have no single dosage coding that GCTA will read correctly.
#
# --double-id sets FID = IID = the VCF sample name.  WITHOUT it plink2 writes
# FID = 0 for every sample, and that is the README section 6 ID gotcha arriving
# through a new door: GCTA identifies an individual by the FID+IID PAIR, so a
# keep-list written FID=IID (which is what every other fileset in this project
# and every keep list we hold uses) matches nothing.  Observed 2026-08-20, jobs
# 34058416/34058417: "Get 5656 samples from list ... After keeping individuals,
# 0 subjects remain."
#
# Those two failed LOUDLY only because they used --keep.  The pooled GRM has no
# --keep, so it would have completed normally and produced a GRM whose .grm.id
# carried FID=0 -- and every later --keep step, including the PC-AiR unrelated
# subset, would then have silently matched nothing.  The loud failure was luck.
set -euo pipefail

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
[[ -n "$REPO_ROOT" ]] || { echo "FATAL: cannot find hpc/config.sh" >&2; exit 2; }
# shellcheck source=/dev/null
source "$REPO_ROOT/hpc/config.sh"

IMP_DIR="${IMP_DIR:-/rds/project/rds-CeXlNYOYMxw/Data_Genetics/genotype_microarray/imputed}"
IMP_OUT="${IMP_OUT:-$ABCD_HPC_ROOT/genotype_imputed}"
IMP_R2="${IMP_R2:-0.8}"
IMP_MAC="${IMP_MAC:-10}"
IMP_MODE="${IMP_MODE:-union}"        # union | pooled | both
T="${SLURM_CPUS_PER_TASK:-8}"

CHR=${SLURM_ARRAY_TASK_ID:?run as an array job, or set SLURM_ARRAY_TASK_ID}
VCF="$IMP_DIR/chr${CHR}_dose.vcf.gz"
[[ -f "$VCF" ]] || { echo "FATAL: no $VCF" >&2; exit 2; }
ensure_dirs "$IMP_OUT"

echo "=== chr$CHR  $(date) ==="
echo "vcf   : $VCF ($(du -h "$VCF" | cut -f1))"
echo "mode  : $IMP_MODE   R2>=$IMP_R2   MAC>=$IMP_MAC"

convert() {   # $1 = tag, $2.. = frequency filter args
  local tag="$1"; shift
  local out="$IMP_OUT/${tag}_chr${CHR}"
  echo "--- $tag ---"
  "$PLINK2" \
    --vcf "$VCF" dosage=DS \
    --double-id \
    --extract-if-info "R2 >= $IMP_R2" \
    "$@" \
    --snps-only --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 100 missing \
    --make-bed \
    --threads "$T" \
    --out "$out"
  # PLINK 2 writes a PLINK 1 .bed here; gate it exactly as the array data was.
  bash "$REPO_ROOT/hpc/work/check_bfile_integrity.sh" "$out" || {
    echo "FATAL: chr$CHR/$tag failed the integrity check" >&2; exit 3; }
  echo "    $tag chr$CHR: $(wc -l < "$out.bim") variants, $(wc -l < "$out.fam") subjects, $(du -h "$out.bed" | cut -f1)"
}

case "$IMP_MODE" in
  union)  convert imp_union  --mac "$IMP_MAC" ;;
  pooled) convert imp_pooled --maf 0.01 ;;
  both)   convert imp_union  --mac "$IMP_MAC"; convert imp_pooled --maf 0.01 ;;
  *) echo "FATAL: IMP_MODE must be union|pooled|both" >&2; exit 2 ;;
esac

echo "=== chr$CHR done  $(date) ==="
