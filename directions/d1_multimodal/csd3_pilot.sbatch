#!/bin/bash
#SBATCH --job-name=d1_pilot
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/%x_%j.log
#
# D1 step 2 pilot on the 20 sessions in processed/hcp/surf_metrics_pilot_sessions.tsv
# (picked by directions/d1_multimodal/csd3_discovery.py):
#   1. abcd.surf_metrics on the pilot sessions (part "pilot");
#   2. gate C inputs: FreeSurfer mri_vol2surf of nu.mgz at projfrac 0.5 and projdist -1 (lh, 3 sessions);
#   3. T2 registration pilot: ABCC fMRIPrep desc-preproc_T1w -> DAIRC FreeSurfer subject
#      (mri_coreg, then bbregister --t1), timed per session.  fMRIPrep's T2w is already
#      resampled into that T1w space, so this one transform puts T2w on the DAIRC surfaces;
#   4. directions/d1_multimodal/csd3_pilot_gates.py: gates A-E and the T2 pilot summary.
#
#   sbatch directions/d1_multimodal/csd3_pilot.sbatch     (from the repo root)

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p slurm
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ABCD_CONDA_ENV:-legacy/hpc/work/envs/abcd}"
set -u
export PYTHONPATH=src:${PYTHONPATH:-}
export OMP_NUM_THREADS=1

REL=abcd-data-release-7.0/processed/hcp
PILOT=$REL/surf_metrics_pilot_sessions.tsv
WORK=$REL/surf_metrics_pilot
FS_ROOT=/rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer
FMRIPREP=/rds/project/rds-CeXlNYOYMxw/derivatives/abcc-fmriprep_v25.1.4
mkdir -p "$WORK"

echo "== 1. surf_metrics on pilot sessions"
rm -f "$REL/surf_metrics_parts/pilot.parquet"
python -m abcd.surf_metrics --sessions "$PILOT" --part-name pilot --workers 4

source /etc/profile.d/modules.sh
module load freesurfer/7.1.0
echo "FreeSurfer: $(cat "$FREESURFER_HOME/build-stamp.txt")"
export SUBJECTS_DIR="$WORK/subjects"
mkdir -p "$SUBJECTS_DIR"
tail -n +2 "$PILOT" | while IFS=$'\t' read -r sub ses _rest; do
    ln -sfn "$FS_ROOT/$sub/$ses" "$SUBJECTS_DIR/${sub}_${ses}"
done

echo "== 2. gate C: mri_vol2surf"
tail -n +2 "$PILOT" | head -3 | while IFS=$'\t' read -r sub ses _rest; do
    s=${sub}_${ses}
    mri_vol2surf --mov "$SUBJECTS_DIR/$s/mri/nu.mgz" --regheader "$s" --hemi lh \
        --projfrac 0.5 --interp trilinear --o "$WORK/$s.lh.nu.projfrac0.5.mgh" > "$WORK/$s.vol2surf.log" 2>&1
    mri_vol2surf --mov "$SUBJECTS_DIR/$s/mri/nu.mgz" --regheader "$s" --hemi lh \
        --projdist -1 --interp trilinear --o "$WORK/$s.lh.nu.projdist-1.mgh" >> "$WORK/$s.vol2surf.log" 2>&1
    # projfrac steps along the white normal; gm50 is the white-pial vertex midpoint.  For an
    # exact-definition check, write that midpoint surface into a scratch copy of the subject
    # (symlinks; the derivatives tree is read-only) and sample it at projdist 0.
    m="$WORK/midsubj/$s"; mkdir -p "$m/surf"
    for d in mri label; do ln -sfn "$FS_ROOT/$sub/$ses/$d" "$m/$d"; done
    for f in "$FS_ROOT/$sub/$ses/surf/"*; do ln -sfn "$f" "$m/surf/"; done
    python -c "
import sys, nibabel.freesurfer as fs
w, f, meta = fs.read_geometry(sys.argv[1] + '/lh.white', read_metadata=True)
p, _ = fs.read_geometry(sys.argv[1] + '/lh.pial')
fs.write_geometry(sys.argv[2] + '/lh.mid', w + 0.5 * (p - w), f, volume_info=meta)
" "$FS_ROOT/$sub/$ses/surf" "$m/surf"
    SUBJECTS_DIR="$WORK/midsubj" mri_vol2surf --mov "$m/mri/nu.mgz" --regheader "$s" --hemi lh \
        --surf mid --projdist 0 --interp trilinear --o "$WORK/$s.lh.nu.mid.mgh" >> "$WORK/$s.vol2surf.log" 2>&1
    echo "  $s done"
done

echo "== 3. T2 registration pilot (fMRIPrep T1w -> DAIRC FreeSurfer)"
reg_one() {
    sub=$1; ses=$2; s=${sub}_${ses}
    t1=$(ls "$FMRIPREP/$sub/$ses/anat/"*_desc-preproc_T1w.nii.gz 2>/dev/null | head -1)
    t2=$(ls "$FMRIPREP/$sub/$ses/anat/"*_desc-preproc_T2w.nii.gz 2>/dev/null | head -1)
    if [[ -z "$t1" || -z "$t2" ]]; then echo -e "$sub\t$ses\tno_fmriprep_T1w_or_T2w\tNA\tNA\tNA"; return; fi
    d="$WORK/reg/$s"; mkdir -p "$d"
    t0=$(date +%s.%N)
    mri_coreg --s "$s" --mov "$t1" --reg "$d/coreg.lta" --threads 1 > "$d/coreg.log" 2>&1 || { echo -e "$sub\t$ses\tmri_coreg_failed\tNA\tNA\tNA"; return; }
    t1s=$(date +%s.%N)
    bbregister --s "$s" --mov "$t1" --init-reg "$d/coreg.lta" --t1 --reg "$d/bbr.lta" > "$d/bbr.log" 2>&1 || { echo -e "$sub\t$ses\tbbregister_failed\tNA\tNA\tNA"; return; }
    t2s=$(date +%s.%N)
    cost=$(awk '{print $1}' "$d/bbr.dat.mincost" 2>/dev/null || echo NA)   # bbregister names it after the .dat
    echo -e "$sub\t$ses\tok\t$(echo "$t1s - $t0" | bc)\t$(echo "$t2s - $t1s" | bc)\t$cost"
}
export -f reg_one; export WORK FMRIPREP SUBJECTS_DIR
echo -e "participant_id\tsession_id\tstatus\tsec_coreg\tsec_bbregister\tbbr_mincost" > "$WORK/t2reg_pilot.tsv"
tail -n +2 "$PILOT" | cut -f1,2 | xargs -P "${SLURM_CPUS_PER_TASK:-8}" -n 2 bash -c 'reg_one "$0" "$1"' >> "$WORK/t2reg_pilot.tsv"
cat "$WORK/t2reg_pilot.tsv"

echo "== 4. gates"
python directions/d1_multimodal/csd3_pilot_gates.py
