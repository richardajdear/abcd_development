#!/bin/bash
#SBATCH --job-name=abcd_rsmap
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH --array=1-22%11
#SBATCH --output=%x_%A_%a.log
#
# Recover rsIDs for the imputed filesets.
#
# WHY THIS IS NEEDED: convert_imputed.sbatch passed
# `--set-all-var-ids '@:#:$r:$a'`, which overwrites EVERY variant ID with
# chr:pos:ref:alt.  That was wrong -- the intent was to name only the variants
# whose ID is missing, which is `--set-MISSING-var-ids`.  As written it discarded
# the rsIDs, and MAGMA, LDSC and the PRS all match on rsID, so all three of the
# downstream analyses were blocked by an ID column.
#
# WHY NOT JUST RE-CONVERT: the fix is one flag, but re-running the conversion
# means re-reading ~890 GB of VCF and re-parsing 11,670 dosage fields per
# variant -- the expensive part -- to change a string.  This job reads the same
# VCFs but takes only columns 1-5, so it never parses a genotype.  Same bytes off
# disk, a small fraction of the CPU.
#
# WHY rsIDs SOLVE THE BUILD PROBLEM TOO: the imputed data is GRCh38 and every
# reference this pipeline uses (g1000_eur, w_hm3.snplist, the sumstats) is hg19.
# An rsID names the same variant in both builds, so matching on rsID sidesteps
# liftover entirely for MAGMA, LDSC and the PRS.  Liftover is only needed for
# position-based work, and none of these three are.
#
# The map key is chr:pos:REF:ALT with the 'chr' prefix stripped, which is exactly
# how plink2's '@:#:$r:$a' rendered it into the .bim -- verified against a .bim
# row rather than assumed.
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
D="$OUT/rsid_map"; ensure_dirs "$D"
CHR=${SLURM_ARRAY_TASK_ID:?run as an array job}
BIM="$IMP_OUT/imp_union_chr${CHR}.bim"
[[ -f "$BIM" ]] || { echo "FATAL: no $BIM" >&2; exit 2; }

echo "=== chr$CHR  $(date) ==="
# Only variants present in our .bim are worth mapping; loading the .bim first
# keeps the output small and makes the coverage number meaningful.
zcat "$IMP_DIR/chr${CHR}_dose.vcf.gz" | cut -f1-5 | awk -v BIM="$BIM" -F'\t' '
BEGIN { while ((getline l < BIM) > 0) { split(l, b, "\t"); want[b[2]] = 1; nb++ } }
/^#/ { next }
{
  c = $1; sub(/^chr/, "", c)
  key = c ":" $2 ":" $4 ":" $5
  if (!(key in want)) next
  seen++
  if ($3 ~ /^rs/) { print key "\t" $3; got++ }
}
END {
  printf "  .bim variants %d, matched in VCF %d, with rsID %d (%.1f%%)\n",
         nb, seen, got, (seen ? 100*got/seen : 0) > "/dev/stderr"
}' > "$D/chr${CHR}.rsid.map"

echo "  wrote $D/chr${CHR}.rsid.map ($(wc -l < "$D/chr${CHR}.rsid.map") entries)"
echo "=== chr$CHR done  $(date) ==="
