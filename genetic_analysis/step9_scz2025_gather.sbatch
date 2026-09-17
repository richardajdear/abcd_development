#!/bin/bash
#SBATCH --job-name=scz25_gather
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=slurm/%x_%A.log
#
# Step 9, gather: concatenate the per-chromosome PRS-CS / PRS-CSx posteriors
# into weight files and score them; build the ancestry-matched composites;
# write within-ancestry-standardised (_zanc) copies of every score directory.
# After this, $SCZ25_ROOT holds one directory per (method, trait-arm) in the
# same shape as legacy prs_final/, ready for step9_scz2025_assoc.sbatch.
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/setup/paths.sh
ROOT="${SCZ25_ROOT:-$V2ROOT/work/scores_scz2025}"
ZANC_PY="$V2ROOT/setup/standardise_within_ancestry.py"
MATCH_PY="$V2ROOT/setup/build_matched_scores.py"
[[ -s "$ZANC_PY" && -s "$MATCH_PY" && -s "$STRATA" ]] || { echo "FATAL: helper or strata file missing" >&2; exit 2; }

gather22() {   # <chr dir> <file glob prefix> <weights out>
  local dir="$1" pref="$2" out="$3" n
  n=$(ls "$dir"/${pref}_pst_eff_a1_b0.5_phiauto_chr*.txt 2>/dev/null | wc -l)
  (( n == 22 )) || { echo "FATAL: $dir/$pref has $n of 22 chromosome files" >&2; return 4; }
  cat "$dir"/${pref}_pst_eff_a1_b0.5_phiauto_chr{1..22}.txt | awk 'BEGIN{OFS="\t"}{print $2,$4,$6}' > "$out"
  echo "  $(basename "$out"): $(wc -l < "$out") SNPs"
}

echo "=== PRS-CS: EUR, META  $(date)"
for TRAIT in SCZ25_EUR SCZ25_META; do
  TAG="${TRAIT//_/}"; OUTD="$ROOT/PRSCS/$TRAIT"
  gather22 "$OUTD/_chr" "$TRAIT" "$OUTD/${TAG}_prscs.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_prscs.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_prscs" > /dev/null
  [[ -s "$OUTD/score_${TAG}_prscs.profile" ]] || { echo "FATAL: no profile for $TRAIT" >&2; exit 4; }
done

echo "=== PRS-CSx: per-population posteriors + meta  $(date)"
for POP in EUR AFR EAS META; do
  OUTD="$ROOT/PRSCSX/SCZ25_$POP"; mkdir -p "$OUTD"
  gather22 "$ROOT/PRSCSX/_chr" "SCZ25_$POP" "$OUTD/SCZ25${POP}_csx.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/SCZ25${POP}_csx.weights" 1 2 3 sum --out "$OUTD/score_SCZ25${POP}_csx" > /dev/null
  [[ -s "$OUTD/score_SCZ25${POP}_csx.profile" ]] || { echo "FATAL: no profile for PRSCSX $POP" >&2; exit 4; }
done

echo "=== ancestry-matched composites  $(date)"
"$PYTHON" "$MATCH_PY" "$STRATA" "$ROOT"

echo "=== within-ancestry standardisation  $(date)"
for d in "$ROOT"/{CT,PRSCS,SBayesR,SBayesRC,PRSCSX}/SCZ25_*; do
  [[ -d "$d" ]] || continue
  [[ "$(basename "$d")" == SCZ25_MATCHED ]] && continue   # already z within cluster
  n=$(ls "$d"/score_*.profile 2>/dev/null | wc -l); (( n > 0 )) || continue
  ZD="$d/_zanc"; mkdir -p "$ZD"
  "$PYTHON" "$ZANC_PY" "$STRATA" "$d" "$ZD" | tail -n +2
  for wf in "$d"/*.weights; do [[ -e "$wf" ]] && ln -sf "$wf" "$ZD/"; done
done

echo "=== manifest"
MAN="$ROOT/scores_manifest.tsv"; printf 'method\ttrait_arm\tn_profiles\tzanc\n' > "$MAN"
for d in "$ROOT"/{CT,PRSCS,SBayesR,SBayesRC,PRSCSX}/SCZ25_*; do
  [[ -d "$d" ]] || continue
  n=$(ls "$d"/score_*.profile 2>/dev/null | wc -l); (( n > 0 )) || continue
  z=no; [[ -d "$d/_zanc" ]] && z=yes
  printf '%s\t%s\t%s\t%s\n' "$(basename "$(dirname "$d")")" "$(basename "$d")" "$n" "$z" >> "$MAN"
done
column -t "$MAN"
echo "=== done $(date)"
