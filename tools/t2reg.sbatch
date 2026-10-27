#!/bin/bash
#SBATCH --job-name=t2reg
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm/%x_%A_%a.log
#
# D1: register each session's ABCC fMRIPrep desc-preproc_T1w to its DAIRC FreeSurfer subject
# (orig.mgz), so fMRIPrep's T2w (already in that T1w space) can be sampled on the DAIRC
# surfaces by abcd.surf_metrics (measures t2__*, t1t2__*).  Same-session T1w -> T1w:
# mri_coreg (rigid, --sep 4, looser tolerances) to initialise, then bbregister --t1
# (boundary-based, 6 dof).  Settings chosen in directions/d1_multimodal/csd3_t2reg_speed.sbatch:
# within 0.006 mm of the default mri_coreg on the pilot, at about half the time.
#
# fMRIPrep's desc-preproc_T2w is NOT bias-corrected (its T1w is, with N4), so the T2w gets
# the same N4 here (ANTs N4BiasFieldCorrection, fMRIPrep's T1w settings, fMRIPrep brain
# mask).  Only the bias field is kept, resampled to 8 mm (corrected = T2w / field);
# surf_metrics divides the T2w samples by it.
#
# Writes <REG_ROOT>/<sub>/<ses>/{bbr.lta, bbr.dat.mincost, coreg.lta, t2_biasfield_8mm.nii.gz}
# (logs kept only on failure) and one status line per session to <REG_ROOT>/status/task_<i>.tsv.
# A session with all of these is skipped, so the array can be resubmitted.  Nothing is written to derivatives/
# (read-only): each session gets a private SUBJECTS_DIR of symlinks.
#
# Submit from the repo root, after writing the session list (sessions with fMRIPrep T1w+T2w):
#   N=$(( ($(wc -l < legacy/hpc/work/d1_t2reg/sessions.tsv) + 199) / 200 ))
#   sbatch --array=0-$((N-1)) tools/t2reg.sbatch

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p slurm
export REG_ROOT="${REG_ROOT:-legacy/hpc/work/d1_t2reg}"
export FS_ROOT=/rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer
export FMRIPREP=/rds/project/rds-CeXlNYOYMxw/derivatives/abcc-fmriprep_v25.1.4
export COREG_OPTS="${COREG_OPTS:---sep 4 --ftol 0.0001 --linmintol 0.01}"
CHUNK=${CHUNK:-200}
source /etc/profile.d/modules.sh
module load freesurfer/7.1.0
# ANTs 2.3.4 (spack build of module ants-2.3.4-gcc-5-lj6vm7c; the module does not reliably
# set PATH in batch shells, so point at it directly)
ANTS=/usr/local/software/spack/spack-git/opt/spack/linux-rhel7-broadwell/gcc-5/ants-2.3.4-lj6vm7ccmij4yywr4ynrttcko4aoflqt
export PATH="$ANTS/bin:$PATH" LD_LIBRARY_PATH="$ANTS/lib:${LD_LIBRARY_PATH:-}"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
mkdir -p "$REG_ROOT/status"

reg_one() {
    sub=$1; ses=$2; s=${sub}_${ses}
    d="$REG_ROOT/$sub/$ses"
    if [[ -s "$d/bbr.lta" && -s "$d/bbr.dat.mincost" && -s "$d/t2_biasfield_8mm.nii.gz" ]]; then
        echo -e "$sub\t$ses\tskipped_done\tNA\tNA\tNA\t$(awk '{print $1}' "$d/bbr.dat.mincost")"; return
    fi
    t1=$(ls "$FMRIPREP/$sub/$ses/anat/"*_desc-preproc_T1w.nii.gz 2>/dev/null | head -1)
    t2=$(ls "$FMRIPREP/$sub/$ses/anat/"*_desc-preproc_T2w.nii.gz 2>/dev/null | head -1)
    mask=$(ls "$FMRIPREP/$sub/$ses/anat/"*_desc-brain_mask.nii.gz 2>/dev/null | head -1)
    if [[ -z "$t1" || -z "$t2" || -z "$mask" || ! -s "$FS_ROOT/$sub/$ses/surf/lh.white" ]]; then
        echo -e "$sub\t$ses\tmissing_input\tNA\tNA\tNA\tNA"; return
    fi
    mkdir -p "$d"
    export SUBJECTS_DIR; SUBJECTS_DIR=$(mktemp -d "${TMPDIR:-/tmp}/t2reg.${s}.XXXX")
    ln -s "$FS_ROOT/$sub/$ses" "$SUBJECTS_DIR/$s"
    t0=$(date +%s.%N); t1s=$t0; t2s=$t0
    if [[ ! -s "$d/bbr.lta" || ! -s "$d/bbr.dat.mincost" ]]; then     # registration not yet done
        # shellcheck disable=SC2086
        if ! mri_coreg --s "$s" --mov "$t1" --reg "$d/coreg.lta" --threads 1 $COREG_OPTS > "$d/coreg.log" 2>&1; then
            echo -e "$sub\t$ses\tmri_coreg_failed\tNA\tNA\tNA\tNA"; rm -rf "$SUBJECTS_DIR"; return
        fi
        t1s=$(date +%s.%N)
        if ! bbregister --s "$s" --mov "$t1" --init-reg "$d/coreg.lta" --t1 --reg "$d/bbr.lta" > "$d/bbr.log" 2>&1; then
            echo -e "$sub\t$ses\tbbregister_failed\tNA\tNA\tNA\tNA"; rm -rf "$SUBJECTS_DIR"; return
        fi
        t2s=$(date +%s.%N)
    fi
    n4=$(mktemp -d "${TMPDIR:-/tmp}/n4.${s}.XXXX")
    # fMRIPrep's T2w has small negative values (resampling undershoot); N4 log-transforms and,
    # given them, silently returns a flat field.  Clip to >= 1 and let N4 rescale (-r 1).
    # (Checked with a synthetic +-30 % field: recovered with log-slope 0.99, r = 1.000.)
    if ! { ImageMath 3 "$n4/t2pos.nii.gz" ReplaceVoxelValue "$t2" -1e9 1 1 \
           && N4BiasFieldCorrection -d 3 -i "$n4/t2pos.nii.gz" -x "$mask" -r 1 -s 4 \
                  -c "[50x50x50x50,1e-7]" -b "[200]" -o "[$n4/t2n4.nii.gz,$n4/field.nii.gz]" \
           && ResampleImage 3 "$n4/field.nii.gz" "$d/t2_biasfield_8mm.nii.gz" 8x8x8 0 0; } > "$d/n4.log" 2>&1; then
        echo -e "$sub\t$ses\tn4_failed\tNA\tNA\tNA\tNA"; rm -rf "$SUBJECTS_DIR" "$n4" "$d/t2_biasfield_8mm.nii.gz"; return
    fi
    # a field that is exactly constant means N4 did not run properly
    rng=$(PrintHeader "$d/t2_biasfield_8mm.nii.gz" 2>/dev/null | awk -F'[][,]' '/Intensity Range/ {print $3 - $2}')
    if [[ -z "$rng" ]] || awk -v r="$rng" 'BEGIN {exit !(r < 1e-4)}'; then
        echo -e "$sub\t$ses\tn4_flat_field\tNA\tNA\tNA\tNA"; rm -rf "$SUBJECTS_DIR" "$n4" "$d/t2_biasfield_8mm.nii.gz"; return
    fi
    t3s=$(date +%s.%N)
    rm -rf "$SUBJECTS_DIR" "$n4" "$d"/bbr.dat~ "$d"/bbr.log.old "$d"/*.log "$d"/bbr.dat.{param,sum} "$d"/bbr.dat
    echo -e "$sub\t$ses\tok\t$(echo "$t1s - $t0" | bc)\t$(echo "$t2s - $t1s" | bc)\t$(echo "$t3s - $t2s" | bc)\t$(awk '{print $1}' "$d/bbr.dat.mincost")"
}
export -f reg_one

i=${SLURM_ARRAY_TASK_ID}
out="$REG_ROOT/status/task_$(printf %04d "$i").tsv"
echo -e "participant_id\tsession_id\tstatus\tsec_coreg\tsec_bbregister\tsec_n4\tbbr_mincost" > "$out"
tail -n +2 "$REG_ROOT/sessions.tsv" | sed -n "$((i * CHUNK + 1)),$(((i + 1) * CHUNK))p" | cut -f1,2 \
    | xargs -P "${SLURM_CPUS_PER_TASK:-8}" -n 2 bash -c 'reg_one "$0" "$1"' >> "$out"
echo "task $i: $(tail -n +2 "$out" | cut -f3 | sort | uniq -c | tr '\n' ' ')"
