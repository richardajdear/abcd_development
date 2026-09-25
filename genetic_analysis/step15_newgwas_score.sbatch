#!/bin/bash
#SBATCH --job-name=newgwas_score
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 15 scoring: C+T, PRS-CS, SBayesR, SBayesRC for the new discovery GWAS
# (BIP_eur, BIP_pooled, ADHD, INT), with exactly the settings of
# setup/prs_final.sbatch / step9_scz2025_score.sbatch so a difference between
# traits is the GWAS, not the method.  Scoring is phenotype-independent: runs
# once, both atlases and both phenotype constructions read the same profiles.
#
#   sbatch --array=1-12 step15_newgwas_score.sbatch heavy     # CT, SBayesR, SBayesRC x 4 arms
#   sbatch --array=1-88 --cpus-per-task=4 --mem=24G step15_newgwas_score.sbatch chr   # PRS-CS x 4 arms x 22 chr
# Optional 2nd argument: comma-separated cell indices to re-run.
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/setup/paths.sh
ROOT="${NEWGWAS_ROOT:-$V2ROOT/work/scores_newgwas}"; GWAS="$ROOT/gwas"
export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR"
ARMS=(BIP_eur BIP_pooled ADHD INT)
GROUP="${1:?usage: step15_newgwas_score.sbatch heavy|chr [cell,cell,...]}"
CELLS=()
case "$GROUP" in
  heavy) for m in CT SBayesR SBayesRC; do for a in "${ARMS[@]}"; do CELLS+=("$m:$a:"); done; done ;;
  chr)   for a in "${ARMS[@]}"; do for c in $(seq 1 22); do CELLS+=("PRSCS:$a:$c"); done; done ;;
  *) echo "FATAL: group must be heavy or chr" >&2; exit 2 ;;
esac
if [[ -n "${2:-}" ]]; then
  IFS=, read -r -a PICK < <(printf '%s\n' "$2")
  I=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 )); (( I < ${#PICK[@]} )) || { echo "FATAL: no entry $((I+1)) in '$2'" >&2; exit 2; }
  T=$(( PICK[I] - 1 ))
else
  T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
fi
(( T >= 0 && T < ${#CELLS[@]} )) || { echo "FATAL: cell $((T+1)) outside 1-${#CELLS[@]}" >&2; exit 2; }
SPEC="${CELLS[$T]}"; METH="${SPEC%%:*}"; REST="${SPEC#*:}"; TRAIT="${REST%%:*}"; CHR="${REST#*:}"
TAG="${TRAIT//_/}"
echo "=== $METH / $TRAIT ${CHR:+chr$CHR }$(date)  cpus=${SLURM_CPUS_PER_TASK:-?}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
MA="$GWAS/$TRAIT.ma"; [[ -s "$MA" ]] || { echo "FATAL: $MA missing (run step15_newgwas_normalise.py)" >&2; exit 3; }

case "$METH" in
CT)
  OUTD="$ROOT/CT/$TRAIT"; mkdir -p "$OUTD"
  awk 'BEGIN{FS=OFS=" "} NR==1{print "SNP","P"; next} {print $1,$7}' "$MA" > "$OUTD/$TRAIT.assoc"
  "$PLINK" --bfile "$GENO" --clump "$OUTD/$TRAIT.assoc" --clump-p1 "$CLUMP_P1" --clump-p2 "$CLUMP_P2" \
    --clump-r2 "$CLUMP_R2" --clump-kb "$CLUMP_KB" --out "$OUTD/$TRAIT" > /dev/null
  awk 'NR>1 && $3!="" {print $3}' "$OUTD/$TRAIT.clumped" > "$OUTD/$TRAIT.index"
  echo "index SNPs: $(wc -l < "$OUTD/$TRAIT.index")"
  for thr in $PRS_THRESHOLDS; do
    tag="${thr//./p}"
    awk -v t="$thr" 'NR==FNR{k[$1]=1; next} FNR==1{next} ($1 in k) && ($7+0)<(t+0) {print $1"\t"$2"\t"$5}' \
      "$OUTD/$TRAIT.index" "$MA" > "$OUTD/${TAG}_$tag.weights"
    if [[ -s "$OUTD/${TAG}_$tag.weights" ]]; then
      "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_$tag.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_$tag" >/dev/null 2>&1 \
        || echo "WARN: score at $thr failed" >&2
    else
      rm -f "$OUTD/${TAG}_$tag.weights"; echo "  thr $thr: no SNPs, skipped"; continue
    fi
    echo "  thr $thr: $(wc -l < "$OUTD/${TAG}_$tag.weights") SNPs"
  done
  rm -f "$OUTD/$TRAIT.assoc" ;;
SBayesR)
  OUTD="$ROOT/SBayesR/$TRAIT"; mkdir -p "$OUTD"
  for c in $(seq 1 22); do m="$SBR_LDM/ukbEURu_hm3_chr${c}_v3_50k.ldm.sparse"; [[ -s "$m.info" ]] && echo "$m"; done > "$OUTD/ldm.txt"
  "$GCTB" --sbayes R --mldm "$OUTD/ldm.txt" --gwas-summary "$MA" --robust --chain-length 10000 --burn-in 4000 \
    --seed 2026 --thread "${SLURM_CPUS_PER_TASK:-8}" --out "$OUTD/$TRAIT" 2>&1 | tail -25
  [[ -s "$OUTD/$TRAIT.snpRes" ]] || { echo "FATAL: no snpRes" >&2; exit 4; }
  awk 'NR>1{print $2"\t"$5"\t"$8}' "$OUTD/$TRAIT.snpRes" > "$OUTD/${TAG}_sbr.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_sbr.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_sbr" > /dev/null ;;
SBayesRC)
  OUTD="$ROOT/SBayesRC/$TRAIT"; mkdir -p "$OUTD"; IMP="$OUTD/${TRAIT}_imp"
  [[ -s "$IMP.imputed.ma" ]] || "$GCTB" --ldm-eigen "$SBRC_EIGEN" --gwas-summary "$MA" --impute-summary \
      --thread "${SLURM_CPUS_PER_TASK:-8}" --out "$IMP" 2>&1 | tail -6
  [[ -s "$IMP.imputed.ma" ]] || { echo "FATAL: imputation failed" >&2; exit 3; }
  "$GCTB" --sbayes RC --ldm-eigen "$SBRC_EIGEN" --gwas-summary "$IMP.imputed.ma" --annot "$ANNOT" \
    --chain-length 3000 --burn-in 1000 --seed 2026 --thread "${SLURM_CPUS_PER_TASK:-8}" --out "$OUTD/$TRAIT" 2>&1 | tail -25
  [[ -s "$OUTD/$TRAIT.snpRes" ]] || { echo "FATAL: no snpRes" >&2; exit 4; }
  awk 'NR>1{print $2"\t"$5"\t"$8}' "$OUTD/$TRAIT.snpRes" > "$OUTD/${TAG}_sbrc.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_sbrc.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_sbrc" > /dev/null ;;
PRSCS)
  SST="$GWAS/$TRAIT.sst"; NG=$(cat "$GWAS/$TRAIT.n_gwas")
  OUTD="$ROOT/PRSCS/$TRAIT/_chr"; mkdir -p "$OUTD"
  out="$OUTD/${TRAIT}_pst_eff_a1_b0.5_phiauto_chr${CHR}.txt"
  if [[ -s "$out" ]]; then echo "exists: $out"; exit 0; fi
  echo "n_gwas=$NG"
  "$PYTHON" "$PRSCS" --ref_dir="$LDBLK_UKBB" --bim_prefix="$GENO" --sst_file="$SST" --n_gwas="$NG" \
    --chrom="$CHR" --out_dir="$OUTD/$TRAIT" --seed=2026
  [[ -s "$out" ]] || { echo "FATAL: PRS-CS wrote nothing for chr$CHR" >&2; exit 4; } ;;
esac
echo "=== $METH / $TRAIT ${CHR:+chr$CHR }done $(date)"
