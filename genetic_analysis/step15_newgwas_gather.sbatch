#!/bin/bash
#SBATCH --job-name=newgwas_gather
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=slurm/%x_%A.log
# Step 15 gather: concatenate + score PRS-CS per arm; within-ancestry-standardised
# (_zanc) copies of every score directory; manifest.
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/setup/paths.sh
ROOT="${NEWGWAS_ROOT:-$V2ROOT/work/scores_newgwas}"; export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR"
ZANC_PY="$V2ROOT/setup/standardise_within_ancestry.py"
for TRAIT in BIP_eur BIP_pooled ADHD INT; do
  TAG="${TRAIT//_/}"; OUTD="$ROOT/PRSCS/$TRAIT"; n=$(ls "$OUTD"/_chr/${TRAIT}_pst_eff_a1_b0.5_phiauto_chr*.txt 2>/dev/null | wc -l)
  (( n == 22 )) || { echo "FATAL: $TRAIT PRS-CS has $n of 22 chromosomes" >&2; exit 4; }
  cat "$OUTD"/_chr/${TRAIT}_pst_eff_a1_b0.5_phiauto_chr{1..22}.txt | awk 'BEGIN{OFS="\t"}{print $2,$4,$6}' > "$OUTD/${TAG}_prscs.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_prscs.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_prscs" > /dev/null
  echo "PRS-CS $TRAIT: $(wc -l < "$OUTD/${TAG}_prscs.weights") SNPs"
done
MAN="$ROOT/scores_manifest.tsv"; printf 'method\ttrait_arm\tn_profiles\n' > "$MAN"
for d in "$ROOT"/{CT,PRSCS,SBayesR,SBayesRC}/*; do
  [[ -d "$d" ]] || continue; n=$(ls "$d"/score_*.profile 2>/dev/null | wc -l); (( n > 0 )) || continue
  ZD="$d/_zanc"; mkdir -p "$ZD"; "$PYTHON" "$ZANC_PY" "$STRATA" "$d" "$ZD" | tail -n +2
  for wf in "$d"/*.weights; do [[ -e "$wf" ]] && ln -sf "$wf" "$ZD/"; done
  printf '%s\t%s\t%s\n' "$(basename "$(dirname "$d")")" "$(basename "$d")" "$n" >> "$MAN"
done
column -t "$MAN"; echo "=== gather done $(date)"
