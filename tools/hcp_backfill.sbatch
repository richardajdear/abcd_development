#!/bin/bash
#SBATCH --job-name=hcp_backfill
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm/%x_%A_%a.log
#
# Backfill the HCP-MMP1.0 parcellation for sessions rr480's array never
# completed.  Each array task handles CHUNK consecutive lines of the session
# list (default 50 -> ~35 min at ~40 s/session, measured).
#
#   LIST=abcd-data-release-7.0/processed/hcp/sessions_to_reparcellate.txt   # 8,874 lines
#   N=$(wc -l < $LIST); CHUNK=50
#   sbatch --array=0-$(( (N + CHUNK - 1) / CHUNK - 1 ))%100 tools/hcp_backfill.sbatch $LIST work/parcellations_backfill/T1
#
# Then:  sbatch tools/hcp_extract.sbatch --out-dir abcd-data-release-7.0/processed/hcp \
#            --extra-parc-root work/parcellations_backfill/T1
#
# Idempotent: sessions whose output already has 180 ROI rows per hemisphere
# are skipped, so a partial array can simply be resubmitted.

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

LIST="${1:?session list}"
OUT_ROOT="${2:?output root}"
CHUNK="${CHUNK:-50}"
export TMPDIR="${TMPDIR:-/tmp}"

start=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
end=$(( start + CHUNK - 1 ))
echo "task $SLURM_ARRAY_TASK_ID: lines $start-$end of $LIST on $(hostname), $(date)"

ok=0; fail=0; skip=0
while read -r id; do
    [[ -z "$id" ]] && continue
    d="$OUT_ROOT/$id/HCP.fsaverage.aparc"
    if [[ -s "$d/lh.HCP.fsaverage.aparc.log" && -s "$d/rh.HCP.fsaverage.aparc.log" ]] \
       && [[ $(grep -c "_ROI$" "$d/lh.HCP.fsaverage.aparc.log") -eq 180 ]] \
       && [[ $(grep -c "_ROI$" "$d/rh.HCP.fsaverage.aparc.log") -eq 180 ]]; then
        echo "DONE $id (already complete)"; skip=$((skip+1)); continue
    fi
    t0=$(date +%s)
    if bash hpc/hcp_backfill_session.sh "$id" "$OUT_ROOT"; then
        ok=$((ok+1))
    else
        rc=$?; echo "FAILED rc=$rc $id"; fail=$((fail+1))
    fi
    echo "  $(( $(date +%s) - t0 ))s"
done < <(sed -n "${start},${end}p" "$LIST")

echo "task $SLURM_ARRAY_TASK_ID finished: ok=$ok fail=$fail skipped=$skip  $(date)"
