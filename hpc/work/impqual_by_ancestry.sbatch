#!/bin/bash
#SBATCH --job-name=abcd_impqual
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.log
#
# Phase 1: is TOPMed imputation quality WORSE for the non-European children,
# and would a single R2 threshold therefore discard their variants
# preferentially?
#
# WHY THIS CANNOT USE THE R2 IN THE FILE
# --------------------------------------
# Minimac writes ONE R2 per variant, estimated across the whole cohort.  It has
# no per-group value, so "filter on R2" cannot be audited by reading R2 -- the
# number that would reveal the bias is the one that does not exist.  Recomputing
# Minimac's R2 per group is not possible after the fact either: it is an
# internal quantity of the imputation, not a function of the output.
#
# What IS recoverable is the thing R2 summarises.  An imputed dosage sits in
# [0,2]; a confident call sits near 0, 1 or 2, and an uncertain one sits between.
# So per variant and per group,
#
#     ambiguity = mean over individuals of  min(|DS-0|, |DS-1|, |DS-2|)
#
# is a direct, per-group measure of imputation uncertainty: 0 = every child
# called confidently, 0.5 = maximally uncertain.  It needs only the DS field and
# a group label, both of which we have.
#
# This matters because it is the SAME failure shape as README section 8.6,
# arriving through a different door: a single genome-wide threshold applied to a
# quantity that means different things in different ancestry groups, quietly
# re-Europeanising the analysis while reporting nothing.  Section 8.6 was caught
# because someone looked; this is looking.
#
# Per-group MAF is reported alongside, because a variant polymorphic in one
# group and monomorphic in another is a second way a pooled MAF filter can drop
# a group's variation without saying so.
#
# Sampled every IMPQUAL_STRIDE-th variant (default 500): parsing 11,670 dosage
# fields for every variant on chr22 would cost hours to answer a question that a
# few thousand variants settle. Sampling is systematic, not random, so it is
# reproducible without a seed.
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
CHR="${IMPQUAL_CHR:-22}"
STRIDE="${IMPQUAL_STRIDE:-500}"
STRATA="${STRATA_TSV:-$ABCD_HPC_ROOT/results/ancestry/strata.tsv}"
D="$OUT/impqual"; ensure_dirs "$D"

[[ -f "$STRATA" ]] || { echo "FATAL: no strata at $STRATA (run work/assign_ancestry.py)" >&2; exit 2; }
echo "chr$CHR, every ${STRIDE}th variant, strata from $STRATA"

zcat "$IMP_DIR/chr${CHR}_dose.vcf.gz" | awk -v STRATA="$STRATA" -v STRIDE="$STRIDE" -F'\t' '
BEGIN {
  while ((getline line < STRATA) > 0) {
    n_s = split(line, f, "\t")
    if (f[1] == "IID") continue            # header
    grp[f[1]] = f[3]
  }
}
# The #CHROM line carries the sample order; map each sample column to a stratum.
/^#CHROM/ {
  for (i = 10; i <= NF; i++) { col_grp[i] = ($i in grp) ? grp[$i] : "unassigned" }
  for (i = 10; i <= NF; i++) seen[col_grp[i]]++
  printf "sample columns: %d\n", NF - 9
  for (g in seen) printf "  %-14s %d\n", g, seen[g]
  next
}
/^#/ {next}
{
  nv++
  if (nv % STRIDE != 0) next
  # DS index within the colon-delimited sample field, read from FORMAT ($9).
  nf = split($9, fmt, ":"); ds_i = 0
  for (i = 1; i <= nf; i++) if (fmt[i] == "DS") ds_i = i
  if (!ds_i) next
  used++
  delete sum_amb; delete sum_ds; delete cnt
  for (i = 10; i <= NF; i++) {
    g = col_grp[i]
    split($i, s, ":")
    ds = s[ds_i] + 0
    a = ds; if (a > 1) a = 2 - a          # min(|ds-0|,|ds-1|,|ds-2|) on [0,2]
    if (a > 0.5) a = 1 - a
    sum_amb[g] += a; sum_ds[g] += ds; cnt[g]++
  }
  for (g in cnt) {
    if (cnt[g] < 20) continue
    amb = sum_amb[g] / cnt[g]
    af  = sum_ds[g] / (2 * cnt[g]); if (af > 0.5) af = 1 - af
    tot_amb[g] += amb; tot_n[g]++
    if (af >= 0.01) poly[g]++
    if (amb <= 0.02) conf[g]++            # near-perfectly called
    if (amb >= 0.15) bad[g]++             # substantially uncertain
  }
}
END {
  printf "\nvariants scanned %d, sampled %d\n\n", nv, used
  printf "%-14s %10s %10s %10s %10s\n", "stratum", "mean_amb", "conf<=.02", "unc>=.15", "MAF>=.01"
  for (g in tot_n)
    printf "%-14s %10.5f %9.1f%% %9.1f%% %9.1f%%\n", g, tot_amb[g]/tot_n[g],
           100*conf[g]/tot_n[g], 100*bad[g]/tot_n[g], 100*poly[g]/tot_n[g]
  printf "\nmean_amb: 0 = every child called confidently, 0.5 = maximally uncertain.\n"
  printf "A materially higher value in a non-EUR stratum means a single R2 threshold\n"
  printf "would discard that stratum'"'"'s variants preferentially.\n"
}'
echo "=== done  $(date) ==="
