#!/bin/bash
#SBATCH --job-name=v3_collect
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#
# Re-run the GENESIS collection for the four new phenotypes, both arms.
#
# WHY THIS EXISTS: run_all.sh's own final collect (jobs 35089437 pooled and
# 35089472 EUR) both died after 5 seconds with
#
#     gzip: stdout: No space left on device
#     Error in fread(cmd = sprintf("gzip -dc %s", shQuote(f))) :
#       ... the disk is full in the temporary directory ('/tmp/RtmpaNR80v')
#
# data.table's fread(cmd=) decompresses each chromosome file through R's
# tempdir(), which on these compute nodes is a small node-local /tmp shared
# with every other job on the node.  Nothing was wrong with the association
# scan: all 176 tasks COMPLETED and all 22 chromosome files are present for
# each of the four phenotypes in both arms.  Only the collector needs rerunning
# --- it is idempotent and rebuilds from whatever chromosome files exist.
#
# The fix is environmental, not a code change: point TMPDIR at /rds, which has
# hundreds of GB free, and give the job enough memory to hold one decompressed
# chromosome comfortably.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"

export TMPDIR="$OUT_V2/tmp_collect"
mkdir -p "$TMPDIR"
echo "TMPDIR=$TMPDIR  ($(df -h "$TMPDIR" | awk 'NR==2{print $4}') free)"

for arm in assoc assoc_eur; do
  D="$OUT_V2/$arm"
  for pheno in slope_topDelta slope_topC3 slope_projDelta slope_projC3; do
    n=$(ls "$D/${pheno}_chr"*.tsv.gz 2>/dev/null | wc -l)
    if (( n != 22 )); then
      echo "FATAL: $arm/$pheno has $n of 22 chromosome files" >&2; exit 3
    fi
    echo "=== collect $arm/$pheno  $(date) ==="
    "$RSCRIPT" "$RSCRIPT_DIR/05_collect_assoc.R" \
      --assoc-dir "$D" --pheno "$pheno" --out-dir "$D" \
      --maf "$ASSOC_MAF" --mac "$ASSOC_MAC"
    [[ -s "$D/$pheno.sumstats.tsv.gz" ]] \
      || { echo "FATAL: $D/$pheno.sumstats.tsv.gz missing/empty" >&2; exit 4; }
  done
  echo
  column -t -s $'\t' "$D/gwas_summary.tsv"
  echo
done
rmdir "$TMPDIR" 2>/dev/null || true
echo "=== done  $(date) ==="
