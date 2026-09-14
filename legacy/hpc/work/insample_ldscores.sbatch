#!/bin/bash
#SBATCH --job-name=abcd_ldscore
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --array=1-22%11
#SBATCH --output=%x_%A_%a.log
#
# Compute LD scores IN THE ANALYSIS SAMPLE, replacing the European reference.
#
# THE PROBLEM: LDSC's shipped LD scores (eur_w_ld_chr) are estimated in
# Europeans.  Our GWAS sample is 32 % non-European.  LDSC regresses chi-square on
# LD score, so if the LD scores describe a different population than the sample,
# the model is misspecified -- and LDSC's INTERCEPT is precisely the term that
# absorbs misspecification.  Measured: intercept 1.0995 (0.0085) for
# baseline_thickness, which reads as confounding but is equally consistent with
# the wrong LD scores.  It is not possible to tell those apart from the fit.
#
# THE FIX: we hold the genotypes.  LD scores computed from the analysis sample
# itself describe that sample's LD by construction -- no ancestry-matching
# question arises, and 11,670 individuals is >20x the 503 EUR individuals in
# 1000 Genomes, so they are also less noisy.
#
# WHAT THIS DOES NOT FIX, stated plainly: in an admixed or structured sample,
# admixture induces long-range LD that a 1 cM window does not capture, so
# in-sample LD scores are themselves biased.  cov-LDSC (Luo et al. 2021) is the
# method that addresses this, using a much larger window (20 cM) plus PC
# covariates; it is not implemented in this LDSC build.  So this is a strict
# improvement on a EUR panel, not a solution to admixture.  The EUR-stratified
# LDSC run remains the arm with no ancestry mismatch at all.
#
# --ld-wind-kb, not --ld-wind-cm: plink2 wrote cM = 0 for every variant when it
# converted the VCFs (there was no genetic map to carry over), so a cM window
# would silently include nothing.  1000 kb is the conventional stand-in for 1 cM.
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

IMP_OUT="${IMP_OUT:-$ABCD_HPC_ROOT/genotype_imputed}"
D="${INSAMPLE_LD_DIR:-$OUT/ld_insample}"
KEEP="${LDSCORE_KEEP:-}"          # optional: restrict to one ancestry stratum
ensure_dirs "$D"
CHR=${SLURM_ARRAY_TASK_ID:?run as an array job}

echo "=== chr$CHR  $(date) ==="
# HapMap3 SNPs only: that is the set LDSC's regression uses, so scoring the full
# 26M would be ~20x the work for variants munge_sumstats discards anyway.
#
# PER-CHROMOSOME path, not a shared one.  The first version wrote $D/hm3.snps
# from every array task at once; with 11 running concurrently, a task could read
# the file while another was truncating it.  Five chromosomes died on
# "--extract: 0 variants remaining" (job 34110235, chr4/6/7/8/9).  They failed
# loudly, but the same race could equally have delivered a PARTIAL list and
# produced LD scores over a silently truncated SNP set -- which nothing
# downstream would have flagged.  Verified afterwards that the surviving
# chromosomes scale correctly with chromosome size, so none were truncated.
SNPF="$D/hm3_chr${CHR}.snps"
awk 'NR>1{print $1}' "$HM3_SNPLIST" > "$SNPF"
keep_args=()
[[ -n "$KEEP" ]] && keep_args=(--keep "$KEEP")

"$PLINK2" --bfile "$IMP_OUT/imp_union_chr${CHR}" \
          --extract "$SNPF" "${keep_args[@]+"${keep_args[@]}"}" \
          --make-bed --threads "${SLURM_CPUS_PER_TASK:-4}" \
          --out "$D/hm3_chr${CHR}"

"$LDSC" --bfile "$D/hm3_chr${CHR}" --l2 --ld-wind-kb 1000 --out "$D/${CHR}"

echo "  chr$CHR: $(wc -l < "$D/hm3_chr${CHR}.bim") SNPs scored"
ls -la "$D/${CHR}".l2.* 2>/dev/null | awk '{print "   ", $5, $9}'
echo "=== chr$CHR done  $(date) ==="
