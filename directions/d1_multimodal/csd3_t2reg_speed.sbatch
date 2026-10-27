#!/bin/bash
#SBATCH --job-name=d1_t2reg_speed
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm/%x_%j.log
#
# Faster fMRIPrep-T1w -> DAIRC FreeSurfer registration, tested against the pilot reference
# (mri_coreg default seps 4,2 + bbregister --t1; surf_metrics_pilot/reg/<s>/bbr.dat):
#   sep4    mri_coreg --sep 4                                       + bbregister --t1
#   sep4tol mri_coreg --sep 4 --ftol 1e-4 --linmintol 1e-2          + bbregister --t1
#   header  bbregister --t1 --init-header (no mri_coreg)
# Writes surf_metrics_pilot/t2reg_speed.tsv; compare with csd3_t2reg_compare.py.

set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
WORK=abcd-data-release-7.0/processed/hcp/surf_metrics_pilot
FMRIPREP=/rds/project/rds-CeXlNYOYMxw/derivatives/abcc-fmriprep_v25.1.4
source /etc/profile.d/modules.sh
module load freesurfer/7.1.0
export SUBJECTS_DIR="$WORK/subjects" WORK FMRIPREP

run_one() {
    sub=$1; ses=$2; v=$3; s=${sub}_${ses}
    t1=$(ls "$FMRIPREP/$sub/$ses/anat/"*_desc-preproc_T1w.nii.gz 2>/dev/null | head -1)
    [[ -z "$t1" ]] && return
    d="$WORK/reg_$v/$s"; mkdir -p "$d"
    t0=$(date +%s.%N)
    case $v in
        sep4)    mri_coreg --s "$s" --mov "$t1" --reg "$d/coreg.lta" --threads 1 --sep 4 > "$d/coreg.log" 2>&1
                 init=(--init-reg "$d/coreg.lta") ;;
        sep4tol) mri_coreg --s "$s" --mov "$t1" --reg "$d/coreg.lta" --threads 1 --sep 4 \
                     --ftol 0.0001 --linmintol 0.01 > "$d/coreg.log" 2>&1
                 init=(--init-reg "$d/coreg.lta") ;;
        header)  init=(--init-header) ;;
    esac
    t1s=$(date +%s.%N)
    bbregister --s "$s" --mov "$t1" "${init[@]}" --t1 --reg "$d/bbr.dat" > "$d/bbr.log" 2>&1 || { echo -e "$sub\t$ses\t$v\tfailed\tNA\tNA\tNA"; return; }
    t2s=$(date +%s.%N)
    echo -e "$sub\t$ses\t$v\tok\t$(echo "$t1s - $t0" | bc)\t$(echo "$t2s - $t1s" | bc)\t$(awk '{print $1}' "$d/bbr.dat.mincost")"
}
export -f run_one
out="$WORK/t2reg_speed.tsv"
echo -e "participant_id\tsession_id\tvariant\tstatus\tsec_coreg\tsec_bbregister\tbbr_mincost" > "$out"
for v in sep4 sep4tol header; do
    tail -n +2 "$WORK/../surf_metrics_pilot_sessions.tsv" | cut -f1,2 | sed "s/$/\t$v/"
done | xargs -P "${SLURM_CPUS_PER_TASK:-8}" -n 3 bash -c 'run_one "$0" "$1" "$2"' >> "$out"
cat "$out"
