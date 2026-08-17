#!/bin/bash
#SBATCH --job-name=abcd_qckeep
#SBATCH --account=VERTES-SL3-CPU   # override at submit: sbatch -A <acct>
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:40:00
#SBATCH --output=%x_%j.log
#
# Compute GENOME-WIDE per-sample call rate on the 7.0 array data and write the
# sample keep-list that split_allanc_chroms.sbatch consumes as
# GENO_ALLANC_KEEP.
#
# WHY A SEPARATE JOB RATHER THAN `plink --mind` IN THE SPLIT
# ----------------------------------------------------------
# --mind inside the per-chromosome split would drop a different sample set on
# every chromosome -- a subject at 4 % genome-wide can be at 6 % on chr9 -- and
# the 22 filesets would then carry 22 different .fam files.  GCTA --mgrm would
# combine them anyway, over an inconsistent sample.  That is precisely the
# quiet-wrong class of failure this pipeline has already produced twice, so the
# threshold is applied ONCE, to one genome-wide number, and the resulting list
# is a file on disk that can be read, counted and committed to.
#
# WHY IT IS NEEDED AT ALL
# -----------------------
# Per-sample missingness in this fileset is not noise, it is batch: 87.6 % of
# its variance is explained by the ten Smokescreen batches, and four of them
# (511 subjects) sit at 7-10 % missing against 0.4-1.7 % for the other six.
#
# GCTA does NOT mean-impute -- verified against chr1.grm.N.bin, whose per-pair
# SNP count varies 29,324-37,800 rather than being constant.  It uses the SNPs
# non-missing in both members of a pair.  So the consequence is ~6 % fewer SNPs
# behind every estimate touching a bad-batch subject (35,256 vs 37,683 on
# chr1): batch-correlated noise in the GRM, which attenuates h2.  That is a
# sensitivity analysis, not an emergency, and this list exists so the
# comparison can be made.
#
# The threshold is reported alongside the counts at several cutoffs, so the
# choice is visible and someone disagreeing with 0.05 can see what 0.02 or 0.10
# would have cost before rerunning.
set -euo pipefail

# Locate the repository -- see the note in split_allanc_chroms.sbatch.
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
  exit 2
}
# shellcheck source=/dev/null
source "$REPO_ROOT/hpc/config.sh"

: "${GENO_ALLANC_SRC:?set GENO_ALLANC_SRC in hpc/config.local.sh}"
MIND="${GENO_ALLANC_MIND:-0.05}"
D="$OUT/qc_allanc"
ensure_dirs "$D"

echo "=== genome-wide --missing  $(date) ==="
"$PLINK" --bfile "$GENO_ALLANC_SRC" --autosome --missing \
         --out "$D/genomewide" --threads "${SLURM_CPUS_PER_TASK:-2}"

echo
echo "=== per-sample call rate, autosomes ==="
awk 'NR>1{n++; s+=$6; for(t=1;t<=4;t++){} \
          if($6>0.01)a++; if($6>0.02)b++; if($6>0.05)c++; if($6>0.10)d++} END{
  printf "  %d samples, mean F_MISS %.5f\n", n, s/n;
  printf "  would be dropped at --mind 0.01: %d   0.02: %d   0.05: %d   0.10: %d\n",
         a+0, b+0, c+0, d+0}' "$D/genomewide.imiss"

KEEP="$D/keep_mind${MIND}.txt"
awk -v m="$MIND" 'NR>1 && $6<=m {print $1"\t"$2}' "$D/genomewide.imiss" > "$KEEP"
echo
echo "threshold applied: --mind $MIND"
echo "wrote $KEEP  ($(wc -l < "$KEEP") of $(( $(wc -l < "$D/genomewide.imiss") - 1 )) subjects kept)"

# Per-variant missingness is reported but NOT turned into a list here: --geno is
# a plink flag the split can apply directly and, unlike --mind, applying it per
# chromosome is safe -- it changes the variant set, which is per-chromosome by
# construction, not the sample set, which must be identical across all 22.
echo
echo "=== per-variant call rate, autosomes ==="
awk 'NR>1{n++; s+=$5; if($5>0.01)a++; if($5>0.02)b++; if($5>0.05)c++} END{
  printf "  %d variants, mean F_MISS %.5f\n", n, s/n;
  printf "  would be dropped at --geno 0.01: %d (%.2f%%)   0.02: %d (%.2f%%)   0.05: %d (%.2f%%)\n",
         a+0, 100*(a+0)/n, b+0, 100*(b+0)/n, c+0, 100*(c+0)/n}' "$D/genomewide.lmiss"

echo
echo "=== done  $(date) ==="
