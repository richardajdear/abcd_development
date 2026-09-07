#!/bin/bash
#SBATCH --job-name=v2_bkpcrel
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
# Build the Zaitlen second GRM from PC-RELATE relatedness instead of from
# GCTA --make-bK on the pooled GRM.  See export_pcrelate_pairs.R for why.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/config.sh
# GCTA GRM off-diagonals are ~2*kinship, so the pipeline's BK_THRESH (on GRM
# values) becomes BK_THRESH/2 as a raw PC-Relate kinship floor.
BK_HALF=$(awk -v t="$BK_THRESH" 'BEGIN{printf "%.4f", t/2}')
OUTDIR="$OUT_V2/reml_zaitlen_pcrel"
ensure_dirs "$OUTDIR"
PAIRS="$OUTDIR/pcrelate_pairs_kin$BK_HALF.tsv"
BK="$OUTDIR/$(basename "$GRM_FULL").bKpcrel$BK_HALF"

echo "=== export PC-Relate pairs (kin >= $BK_HALF)  $(date) ==="
"$RSCRIPT" "$ABCD_HPC_ROOT/setup/export_pcrelate_pairs.R" \
  --pcrelate "$KIN_DIR/pcrelate.rds" --out "$PAIRS" --kin-thresh "$BK_HALF"

echo "=== rewrite GRM  $(date) ==="
python3 "$ABCD_HPC_ROOT/setup/make_bk_pcrelate.py" "$GRM_FULL" "$PAIRS" "$BK"
ln -sf "$GRM_FULL.grm.N.bin" "$BK.grm.N.bin"
cp "$GRM_FULL.grm.id" "$BK.grm.id"

# The two .bin files must be the same size (same triangle, only zeros differ).
a=$(wc -c < "$GRM_FULL.grm.bin"); b=$(wc -c < "$BK.grm.bin")
[[ "$a" == "$b" ]] || { echo "FATAL: size mismatch $a vs $b" >&2; exit 3; }
printf '%s\n%s\n' "$GRM_FULL" "$BK" > "$OUTDIR/mgrm.txt"
echo "=== done  $(date) ==="
cat "$OUTDIR/mgrm.txt"
