#!/bin/bash
# SCTWAS driver (idempotent; rerun as downloads finish). Submits, per cell type whose files are downloaded and not
# yet submitted: prep -> 22-chr SuSiE array -> collect+assoc, chained with afterok. Run from anywhere.
# Stage 0 (00_snp_filter.sbatch) must have finished (ld/panel_keep.bed present).
#   CTS="Ext1 Ext" bash run_sctwas.sh        # subset
#   EXTRA: Ext1 is additionally fitted at n = 983 (tag Ext1_n983) as the sample-size sensitivity.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; D=$REPO/genetic_analysis/sctwas; W=$REPO/genetic_analysis/work/sctwas
cd $W; mkdir -p slurm stage1 submitted
[[ -s ld/panel_keep.bed ]] || { echo "stage 0 not done"; exit 1; }
CTS="${CTS:-Ext1 Ext OPC IN OD Ast MG Ext2 Ext6 Ext3 Ext4 Ext5 Ext7 Ext8}"
jid() { sed -E 's/[^0-9]*([0-9]+).*/\1/'; }
for ct in $CTS; do
  [[ -e submitted/$ct ]] && continue
  [[ -e raw/${ct}_eqtl_full_assoc.tsv.gz.done && -e raw/${ct}_eqtl_top_assoc.tsv.gz.done ]] || { echo "$ct: not downloaded yet"; continue; }
  if [[ -s stage1/$ct/prep_summary.tsv ]]; then dep=""; else
    p=$(sbatch --parsable --export=ALL,CT=$ct $D/01_prep_eqtl.sbatch); dep="--dependency=afterok:$p"; fi
  f=$(sbatch --parsable $dep --export=ALL,CT=$ct,NSAMP=757,TAG=$ct $D/02_susie_fit.sbatch)
  a=$(sbatch --parsable --dependency=afterok:$f --export=ALL,TAG=$ct $D/03_collect_assoc.sbatch)
  echo "$ct prep=${p:-skip} fit=$f assoc=$a" | tee submitted/$ct
  if [[ $ct == Ext1 ]]; then
    f2=$(sbatch --parsable $dep --export=ALL,CT=Ext1,NSAMP=983,TAG=Ext1_n983 $D/02_susie_fit.sbatch)
    a2=$(sbatch --parsable --dependency=afterok:$f2 --export=ALL,TAG=Ext1_n983 $D/03_collect_assoc.sbatch)
    echo "Ext1_n983 fit=$f2 assoc=$a2" | tee -a submitted/$ct
  fi
  unset p
done
