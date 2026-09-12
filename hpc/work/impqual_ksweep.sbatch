#!/bin/bash
#SBATCH --job-name=abcd_ksweep
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.log
#
# Are the pre-GRM ancestry diagnostics sensitive to k, the number of k-means
# clusters?  k=4 was chosen from a preview and never tested, which is exactly
# the kind of unexamined choice this project keeps getting caught by.
#
# All k are evaluated in ONE pass over the VCF.  A separate pass per k would
# cost an hour each and answer the same question; the per-variant dosage parse
# is the expensive part and it is shared.
#
# Emits machine-readable TSV (one row per k x stratum x metric) rather than a
# formatted table, because the report plots these and a human-readable table
# would have to be re-parsed to do it.
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
KLIST="${KLIST:-3 4 5}"
A="$ABCD_HPC_ROOT/results/ancestry"
D="$OUT/impqual"; ensure_dirs "$D"
OUTTSV="$D/ksweep_chr${CHR}.tsv"

# One file per k, passed as k:path pairs so awk can load them all up front.
SPEC=""
for k in $KLIST; do
  f="$A/strata_k${k}.tsv"
  [[ -f "$f" ]] || { echo "FATAL: missing $f" >&2; exit 2; }
  SPEC="${SPEC}${SPEC:+,}${k}:${f}"
done
echo "chr$CHR stride=$STRIDE  strata: $SPEC"

zcat "$IMP_DIR/chr${CHR}_dose.vcf.gz" | awk -v SPEC="$SPEC" -v STRIDE="$STRIDE" -F'\t' '
BEGIN {
  nk = split(SPEC, kk, ",")
  for (a = 1; a <= nk; a++) {
    split(kk[a], kv, ":"); K[a] = kv[1]
    while ((getline line < kv[2]) > 0) {
      split(line, f, "\t"); if (f[1] == "IID") continue
      lab[a, f[1]] = f[3]
      if (!( (a SUBSEP f[3]) in strseen )) { strseen[a, f[3]] = 1; nstr[a]++; STR[a, nstr[a]] = f[3] }
    }
    close(kv[2])
  }
}
/^#CHROM/ {
  for (i = 10; i <= NF; i++) for (a = 1; a <= nk; a++)
    cg[a, i] = ((a SUBSEP $i) in lab) ? lab[a, $i] : "unassigned"
  next
}
/^#/ { next }
{
  nv++; if (nv % STRIDE) next
  nf = split($9, fmt, ":"); ds_i = 0
  for (i = 1; i <= nf; i++) if (fmt[i] == "DS") ds_i = i
  if (!ds_i) next
  r2 = 0
  if (match($8, /(^|;)R2=[0-9.e+-]+/)) { r2s = substr($8, RSTART, RLENGTH); sub(/^[;]?R2=/, "", r2s); r2 = r2s + 0 }
  used++
  delete sds; delete samb; delete cnt
  for (i = 10; i <= NF; i++) {
    split($i, sf, ":"); ds = sf[ds_i] + 0
    amb = ds; if (amb > 1) amb = 2 - amb; if (amb > 0.5) amb = 1 - amb
    for (a = 1; a <= nk; a++) { g = cg[a, i]; sds[a, g] += ds; samb[a, g] += amb; cnt[a, g]++ }
  }
  for (a = 1; a <= nk; a++) {
    tds = 0; tn = 0
    for (j = 1; j <= nstr[a]; j++) { g = STR[a, j]; if (cnt[a, g] < 1) continue; tds += sds[a, g]; tn += cnt[a, g] }
    if (tn < 1) continue
    pooled = tds / (2 * tn); if (pooled > 0.5) pooled = 1 - pooled
    anyc = 0; maxaf = 0; whichmax = ""
    for (j = 1; j <= nstr[a]; j++) {
      g = STR[a, j]; if (cnt[a, g] < 1) { af[a, g] = -1; continue }
      x = sds[a, g] / (2 * cnt[a, g]); if (x > 0.5) x = 1 - x
      af[a, g] = x
      if (x >= 0.01) anyc = 1
      if (x > maxaf) { maxaf = x; whichmax = g }
    }
    pp = (pooled >= 0.01)
    if (anyc) { n_any[a]++; if (r2 >= 0.8) n_any_r2[a]++ }
    if (pp) n_pool[a]++
    if (anyc && !pp) { n_lost[a]++; lost_by[a, whichmax]++ }
    if (anyc) for (j = 1; j <= nstr[a]; j++) {
      g = STR[a, j]
      if (cnt[a, g] >= 1 && af[a, g] >= 0.01) {
        v = samb[a, g] / cnt[a, g]
        amb_sum[a, g] += v; amb_n[a, g]++
        if (v >= 0.05) amb_bad[a, g]++
      }
    }
  }
}
END {
  print "k\tstratum\tmetric\tvalue"
  for (a = 1; a <= nk; a++) {
    printf "%s\t%s\t%s\t%d\n", K[a], "ALL", "n_stratum_common", n_any[a]
    printf "%s\t%s\t%s\t%d\n", K[a], "ALL", "n_pooled_pass",    n_pool[a]
    printf "%s\t%s\t%s\t%d\n", K[a], "ALL", "n_lost_by_pooled", n_lost[a]
    printf "%s\t%s\t%s\t%d\n", K[a], "ALL", "n_stratum_common_r2pass", n_any_r2[a]
    printf "%s\t%s\t%s\t%d\n", K[a], "ALL", "n_sampled", used
    for (j = 1; j <= nstr[a]; j++) {
      g = STR[a, j]
      printf "%s\t%s\t%s\t%d\n", K[a], g, "lost_attributed", lost_by[a, g] + 0
      printf "%s\t%s\t%s\t%d\n", K[a], g, "n_common_variants", amb_n[a, g] + 0
      if (amb_n[a, g] > 0) {
        printf "%s\t%s\t%s\t%.6f\n", K[a], g, "mean_ambiguity", amb_sum[a, g] / amb_n[a, g]
        printf "%s\t%s\t%s\t%.4f\n", K[a], g, "pct_uncertain_ge05", 100 * amb_bad[a, g] / amb_n[a, g]
      }
    }
  }
}' > "$OUTTSV"

echo "wrote $OUTTSV ($(wc -l < "$OUTTSV") rows)"
column -t -s $'\t' "$OUTTSV" | head -40
echo "=== done  $(date) ==="
