#!/bin/bash
#SBATCH --job-name=abcd_grm_pcair
#SBATCH --account=VERTES-SL3-CPU   # override at submit: sbatch -A <acct>
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.log
#
# Cut an "unrelated" GRM using GENESIS PC-AiR's sample list instead of GCTA's
# --grm-cutoff, so 02_reml can be run both ways without modification.
#
# WHY BOTH
# --------
# --grm-cutoff 0.05 thresholds relatedness computed against POOLED allele
# frequencies.  In one ancestry group that is relatedness.  Across several it is
# relatedness PLUS shared ancestry: two unrelated members of the same group are
# more alike than the pooled average and score positive, so a single global
# threshold prunes the groups unequally and quietly changes the ancestry
# composition of the "unrelated" subset.  h2 does not show this happening.
#
# PC-AiR (Conomos et al.) exists for exactly this: it partitions a structured
# sample into a maximally-unrelated subset using kinship estimates that are
# adjusted for ancestry, so shared ancestry is not mistaken for kinship.  ABCD
# 7.0 ships the result at genesis/unrelateds_individuals.txt (8,181 subjects),
# computed by the release rather than by us.
#
# Neither is "the" right answer, so both are run and reported:
#   --grm-cutoff  is comparable with the published EUR h2, which used it;
#   PC-AiR        is the estimator appropriate to a multi-ancestry sample,
#                 and retains 5,649 phenotyped subjects against 3,329 in EUR.
# A material gap between the two IS the finding, not a nuisance to resolve.
#
# GRM_ALLANC_PCAIR_KEEP defaults to the file work/ writes from the release list
# in FID<TAB>IID form; the 7.0 .fam has FID == IID on every row (verified), so
# the doubled column is correct rather than a placeholder.
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
[[ -n "$REPO_ROOT" ]] || { echo "FATAL: cannot find hpc/config.sh" >&2; exit 2; }
# shellcheck source=/dev/null
source "$REPO_ROOT/hpc/config.sh"

require_paths GRM_ALLANC GCTA
KEEP="${GRM_ALLANC_PCAIR_KEEP:-$ABCD_HPC_ROOT/results/ancestry/pcair_unrelated.keep}"
OUTP="${GRM_ALLANC_PCAIR:-$GRM_ALLANC.pcair}"
T="${SLURM_CPUS_PER_TASK:-4}"

[[ -f "$KEEP" ]] || { echo "FATAL: no PC-AiR keep list at $KEEP" >&2; exit 2; }
echo "keep list: $KEEP ($(wc -l < "$KEEP") subjects)"
echo "full GRM : $GRM_ALLANC ($(wc -l < "$GRM_ALLANC.grm.id") subjects)"

"$GCTA" --grm "$GRM_ALLANC" --keep "$KEEP" --make-grm --thread-num "$T" --out "$OUTP"

n_out=$(wc -l < "$OUTP.grm.id")
echo "PC-AiR unrelated GRM: $n_out subjects"

# A keep-list naming subjects the GRM does not hold is ignored silently by GCTA,
# so the intersection is stated rather than left to be inferred from a number
# that looks plausible either way.
echo "  (keep list had $(wc -l < "$KEEP"); the difference is subjects not in the GRM)"

echo "=== done  $(date) ==="
