#!/bin/bash
# Parcellate ONE FreeSurfer session into HCP-MMP1.0 without writing into the
# (read-only, managers-owned) derivatives tree.
#
#   hcp_backfill_session.sh sub-XXXXXXXX/ses-00A /path/to/output_root
#
# Writes  <output_root>/<sub>/<ses>/HCP.fsaverage.aparc/{lh,rh}.HCP.fsaverage.aparc.{annot,log}
# in exactly the layout of derivatives/parcellations/T1, so abcd.hcp_stats can
# read both trees.  Method = Code/parcellation_T1/parcellate.sh (R. Romero-
# Garcia), minus the volume parcellation and the other atlases, and with a
# PRIVATE SUBJECTS_DIR of symlinks per task -- the original script's shared
# per-subject fsaverageSubP symlink is what raced between concurrent sessions
# of one subject and left 5,439 sessions with empty stats (see
# docs/PLAN_HCP_thickness.md sec 1.4).
#
# Exit status 0 only if both hemispheres have 180 *_ROI rows.

set -eo pipefail

ID="$1"                         # sub-XXXX/ses-YYA
OUT_ROOT="$2"
FS_ROOT="${FS_ROOT:-/rds/project/rds-CeXlNYOYMxw/derivatives/freesurfer}"
FSAVG="${FSAVG:-/rds/project/rds-CeXlNYOYMxw/userdata/rr480/fsaverageSubP}"
ATLAS=HCP.fsaverage.aparc

sub="${ID%/*}"; ses="${ID#*/}"
SRC="$FS_ROOT/$sub/$ses"
OUT="$OUT_ROOT/$sub/$ses/$ATLAS"

if [[ ! -s "$SRC/surf/lh.sphere.reg" || ! -s "$SRC/surf/rh.sphere.reg" || ! -s "$SRC/surf/lh.thickness" ]]; then
    echo "SKIP $ID: surfaces incomplete" >&2
    exit 3
fi

source /etc/profile.d/modules.sh
module load freesurfer/7.1.0 >/dev/null 2>&1

# private SUBJECTS_DIR: two symlinks, nothing else
export SUBJECTS_DIR
SUBJECTS_DIR="$(mktemp -d "${TMPDIR:-/tmp}/fsdir.${sub}_${ses}.XXXX")"
trap 'rm -rf "$SUBJECTS_DIR"' EXIT
ln -s "$FSAVG" "$SUBJECTS_DIR/fsaverageSubP"
ln -s "$SRC"   "$SUBJECTS_DIR/${sub}_${ses}"

mkdir -p "$OUT"
for hemi in lh rh; do
    annot="$OUT/$hemi.$ATLAS.annot"
    log="$OUT/$hemi.$ATLAS.log"
    if [[ ! -s "$annot" ]]; then
        mri_surf2surf --srcsubject fsaverageSubP \
                      --sval-annot "$FSAVG/label/$hemi.$ATLAS.annot" \
                      --trgsubject "${sub}_${ses}" \
                      --trgsurfval "$annot" \
                      --hemi "$hemi" > "$OUT/$hemi.surf2surf.log" 2>&1
    fi
    mris_anatomical_stats -a "$annot" -b "${sub}_${ses}" "$hemi" > "$log" 2> "$OUT/$hemi.stats.err"
    n=$(grep -c "_ROI$" "$log" || true)
    if [[ "$n" -ne 180 ]]; then
        echo "FAIL $ID $hemi: $n ROI rows" >&2
        exit 4
    fi
done
echo "OK $ID"
