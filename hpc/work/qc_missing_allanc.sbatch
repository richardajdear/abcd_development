#!/bin/bash
#SBATCH --job-name=abcd_qcmiss
#SBATCH --account=VERTES-SL3-CPU   # override at submit: sbatch -A <acct>
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:25:00
#SBATCH --output=%x_%j.log
#
# A QC statistic, deliberately NOT a filter.
#
# config.sh's MAF note records that --maf is the only SNP filter GCTA
# implements: --geno and --hwe are PLINK options that GCTA rejects, so anything
# beyond MAF has to happen upstream.  The all-ancestry GRM is therefore built
# with --maf 0.01 and nothing else, which is fine if the 7.0 array data is
# well-called and quietly wrong if it is not.  That is a question with an
# answer, so this measures it rather than leaving the h2 resting on an
# assumption about someone else's QC.
#
# Reported, not acted on: if the numbers are clean there is nothing to do, and
# if they are not, the decision to filter should be a visible edit rather than a
# default nobody chose.  chr1 alone, because per-SNP call rate is a property of
# the assay and the batch, not of the chromosome, and one chromosome is enough
# to see a problem that would be worth acting on.
set -euo pipefail

# Locate the repository -- see the note in split_allanc_chroms.sbatch: under
# sbatch, ${BASH_SOURCE[0]} names a copy in /var/spool/slurm/slurmd/jobNNN/.
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

CHR="${QC_CHR:-1}"
D="$OUT/qc_allanc"
ensure_dirs "$D"

"$PLINK" --bfile "$(geno_allanc_prefix "$CHR")" --missing --out "$D/chr$CHR" --threads 1

echo
echo "=== per-SNP missingness (chr$CHR) ==="
awk 'NR>1{n++; s+=$5; if($5>0.02)a++; if($5>0.05)b++} END{
  printf "  %d variants, mean F_MISS %.5f, >2%%: %d (%.2f%%), >5%%: %d (%.2f%%)\n",
         n, s/n, a+0, 100*(a+0)/n, b+0, 100*(b+0)/n}' "$D/chr$CHR.lmiss"

echo "=== per-sample missingness (chr$CHR) ==="
awk 'NR>1{n++; s+=$6; if($6>0.02)a++; if($6>0.05)b++} END{
  printf "  %d samples, mean F_MISS %.5f, >2%%: %d, >5%%: %d\n",
         n, s/n, a+0, b+0}' "$D/chr$CHR.imiss"
