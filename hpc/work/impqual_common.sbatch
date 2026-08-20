#!/bin/bash
#SBATCH --job-name=abcd_impcommon
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.log
#
# Phase 1b.  Phase 1 (job 33839205) returned mean dosage ambiguity of
# 0.0003-0.0007 in every ancestry stratum and concluded imputation quality is
# not ancestry-biased.  That conclusion is not yet earned, and the reason is a
# selection effect in the statistic itself:
#
#   97.5% of TOPMed variants are near-monomorphic in this cohort (Phase 0).
#   For those, essentially every child has dosage ~0, so ambiguity is ~0 BY
#   CONSTRUCTION.  Averaging over all variants therefore measures how rare the
#   variants are, not how well they were imputed, and it would return ~0.0004
#   even if every common variant were imputed badly.
#
# So this re-measures ambiguity over the variants that will ACTUALLY BE KEPT --
# those common in at least one stratum -- which is the only set whose imputation
# quality can affect the GRM.
#
# It also answers the question Phase 0 identified as the real risk.  The MAF in
# the INFO field is the POOLED cohort frequency, so a variant common in one
# ancestry group and absent elsewhere can fail a pooled threshold despite being
# informative for the children it is common in.  Phase 1 already hinted at the
# scale: 3.5% of variants clear MAF 0.01 within cluster1 against 2.0% within
# EURlike.  This counts the variants a pooled filter would discard that a
# per-group filter would keep, and names which stratum loses them.
#
# That is the same shape as README section 8.6 and section 8.12: one pooled
# number standing in for a quantity that differs by ancestry.
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
STRIDE="${IMPQUAL_STRIDE:-200}"
STRATA="${STRATA_TSV:-$ABCD_HPC_ROOT/results/ancestry/strata.tsv}"
D="$OUT/impqual"; ensure_dirs "$D"

echo "chr$CHR, every ${STRIDE}th variant"

zcat "$IMP_DIR/chr${CHR}_dose.vcf.gz" | awk -v STRATA="$STRATA" -v STRIDE="$STRIDE" -F'\t' '
BEGIN {
  while ((getline line < STRATA) > 0) { split(line,f,"\t"); if (f[1]!="IID") grp[f[1]]=f[3] }
  NG = split("EURlike cluster1 cluster2 cluster3", G, " ")
}
/^#CHROM/ { for (i=10;i<=NF;i++) col_grp[i] = ($i in grp) ? grp[$i] : "unassigned"; next }
/^#/ {next}
{
  nv++; if (nv % STRIDE != 0) next
  nf = split($9, fmt, ":"); ds_i = 0
  for (i=1;i<=nf;i++) if (fmt[i]=="DS") ds_i=i
  if (!ds_i) next
  # pooled R2 from INFO
  r2 = 0
  if (match($8, /(^|;)R2=[0-9.e+-]+/)) { r2s=substr($8,RSTART,RLENGTH); sub(/^[;]?R2=/,"",r2s); r2=r2s+0 }
  used++
  delete sum_ds; delete sum_amb; delete cnt
  for (i=10;i<=NF;i++) {
    g = col_grp[i]; split($i, s, ":"); ds = s[ds_i]+0
    a = ds; if (a>1) a = 2-a; if (a>0.5) a = 1-a
    sum_ds[g]+=ds; sum_amb[g]+=a; cnt[g]++
  }
  # pooled AF over assigned strata
  tds=0; tn=0
  for (j=1;j<=NG;j++) { g=G[j]; if (cnt[g]<1) continue; tds+=sum_ds[g]; tn+=cnt[g] }
  if (tn<1) next
  pooled = tds/(2*tn); if (pooled>0.5) pooled = 1-pooled

  anycommon = 0; whichmax=""; maxaf=0
  for (j=1;j<=NG;j++) {
    g=G[j]; if (cnt[g]<1) { gaf[g]=-1; continue }
    af = sum_ds[g]/(2*cnt[g]); if (af>0.5) af=1-af
    gaf[g]=af
    if (af>=0.01) { anycommon=1; ncommon_in[g]++ }
    if (af>maxaf) { maxaf=af; whichmax=g }
  }
  pooledpass = (pooled>=0.01)

  if (anycommon)  n_any++
  if (pooledpass) n_pool++
  if (anycommon && !pooledpass) { n_lost++; lost_by[whichmax]++ }
  if (!anycommon && pooledpass) n_gain++
  if (anycommon && pooledpass && r2>=0.8) n_both_r2++
  if (anycommon && r2>=0.8) n_any_r2++

  # ambiguity restricted to variants that are COMMON IN THAT STRATUM -- the only
  # variants whose imputation quality can influence that stratum'"'"'s GRM entries.
  if (anycommon) for (j=1;j<=NG;j++) {
    g=G[j]
    if (cnt[g] >= 1 && gaf[g] >= 0.01) { amb_c[g] += sum_amb[g]/cnt[g]; amb_n[g]++
                          if (sum_amb[g]/cnt[g] >= 0.05) amb_bad[g]++ }
  }
}
END {
  printf "\nsampled %d of %d variants\n", used, nv
  printf "\n--- what a POOLED MAF>=0.01 filter does ---\n"
  printf "  common in >=1 stratum            %8d\n", n_any
  printf "  passing pooled MAF>=0.01         %8d\n", n_pool
  printf "  LOST: common in a stratum but    %8d   (%.1f%% of stratum-common)\n", n_lost, (n_any?100*n_lost/n_any:0)
  printf "        failing the pooled filter\n"
  printf "  gained: pooled-common but not    %8d\n", n_gain
  printf "        common in any one stratum\n"
  printf "\n  which stratum the lost variants are commonest in:\n"
  for (g in lost_by) printf "    %-12s %6d  (%.1f%% of lost)\n", g, lost_by[g], (n_lost?100*lost_by[g]/n_lost:0)
  printf "\n--- imputation quality, RESTRICTED to variants common in that stratum ---\n"
  printf "  %-12s %10s %12s %10s\n", "stratum", "mean_amb", "frac_amb>=.05", "n_var"
  for (j=1;j<=NG;j++) { g=G[j]
    if (amb_n[g]>0) printf "  %-12s %10.5f %11.2f%% %10d\n", g, amb_c[g]/amb_n[g], 100*amb_bad[g]/amb_n[g], amb_n[g] }
  printf "\n  (Phase 1 got ~0.0004 for every stratum by averaging over all variants,\n"
  printf "   97.5%% of which are near-monomorphic and trivially unambiguous.)\n"
  printf "\n--- R2>=0.8 among stratum-common variants ---\n"
  printf "  common in >=1 stratum            %8d\n", n_any
  printf "  ... and R2>=0.8                  %8d  (%.1f%%)\n", n_any_r2, (n_any?100*n_any_r2/n_any:0)
}'
echo "=== done  $(date) ==="
