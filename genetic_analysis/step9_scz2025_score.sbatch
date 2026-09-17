#!/bin/bash
#SBATCH --job-name=scz25_score
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/%x_%A_%a.log
#
# Step 9, scoring: polygenic weights and PLINK scores for the Nature 2025 SCZ
# GWAS.  Scoring is phenotype-independent (README_HPC.md 3), so this runs
# ONCE and both parcellation arms (DK, HCP-MMP) read the same profiles.
#
# Every method's settings are the ones setup/prs_final.sbatch used for PGC3
# (clump p1=p2=1, r2 0.1, 250 kb, 8 thresholds; PRS-CS phi auto on UKB EUR
# LD; SBayesR robust 10k/4k on the shrunk-sparse UKB HM3 LDM; SBayesRC 3k/1k
# with baseline 2.2 annotations on the UKB eigen LDM), so a change in beta is
# a change in the discovery GWAS and nothing else.
#
# ANCESTRY DESIGN (the point of this step; README_HPC.md 8, 2026-09-17):
#   SCZ25_EUR   -> EUR target (rule 4 matched; replaces PGC3 european)     CT, PRS-CS, SBayesR, SBayesRC
#   SCZ25_META  -> pooled target (rule 4 matched; replaces PGC3 primary)   CT, PRS-CS, SBayesR, SBayesRC
#                  (the Bayesian trio use EUR LD on a multi-ancestry file --
#                  the same discovery-side LD caveat as before; C+T clumps
#                  on our own sample and carries none)
#   SCZ25_AFR, SCZ25_EAS -> C+T only, for the per-stratum transferability
#                  tests and the ancestry-matched composite (gather step).
#                  SBayesR/RC have only a EUR LD reference here, so they are
#                  not run on the AFR/EAS files.
#   PRS-CSx     -> EUR + AFR + EAS GWAS jointly, each with its own 1000G LD
#                  panel; writes per-population posteriors and a meta weight.
#                  This is the principled multi-ancestry method and the new
#                  thing the 2025 release makes possible.
#
# Two task groups, because they want different resources:
#   sbatch --array=1-8  step9_scz2025_score.sbatch heavy
#       CT x {EUR,AFR,EAS,META}, SBayesR x {EUR,META}, SBayesRC x {EUR,META}
#   sbatch --array=1-66 --cpus-per-task=4 --mem=24G step9_scz2025_score.sbatch chr
#       PRS-CS {EUR,META} x chr1-22 (44) and PRS-CSx x chr1-22 (22); the
#       per-chromosome split turns PGC3's 7-hour serial PRS-CS run into ~1 h.
#       step9_scz2025_gather.sbatch concatenates and scores them.
set -euo pipefail
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/setup/paths.sh
ROOT="${SCZ25_ROOT:-$V2ROOT/work/scores_scz2025}"
GWAS="$ROOT/gwas"
LD1KG=/rds/user/rajd2/hpc-work/prscs            # ldblk_1kg_{eur,afr,eas,amr} + snpinfo_mult_1kg_hm3
PRSCSX=/rds/user/rajd2/hpc-work/prscs/PRScsx/PRScsx.py

GROUP="${1:?usage: step9_scz2025_score.sbatch heavy|chr}"
CELLS=()
case "$GROUP" in
  heavy) CELLS=(CT:SCZ25_EUR CT:SCZ25_AFR CT:SCZ25_EAS CT:SCZ25_META
                SBayesR:SCZ25_EUR SBayesR:SCZ25_META SBayesRC:SCZ25_EUR SBayesRC:SCZ25_META) ;;
  chr)   for t in SCZ25_EUR SCZ25_META; do for c in $(seq 1 22); do CELLS+=("PRSCS:$t:$c"); done; done
         for c in $(seq 1 22); do CELLS+=("PRSCSX:ALL:$c"); done ;;
  *) echo "FATAL: group must be heavy or chr" >&2; exit 2 ;;
esac
# Optional 2nd argument: comma-separated CELL INDICES (1-based, into the group
# list) so a subset can be re-run: array task i runs the i-th index given.
#   sbatch --array=1-3 step9_scz2025_score.sbatch chr 1,5,40
if [[ -n "${2:-}" ]]; then
  IFS=, read -r -a PICK < <(printf '%s\n' "$2")
  I=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
  (( I >= 0 && I < ${#PICK[@]} )) || { echo "FATAL: task $((I+1)) has no entry in '$2'" >&2; exit 2; }
  T=$(( PICK[I] - 1 ))
else
  T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
fi
(( T >= 0 && T < ${#CELLS[@]} )) || { echo "FATAL: cell index $((T+1)) outside 1-${#CELLS[@]} for group $GROUP" >&2; exit 2; }
# No here-strings: bash backs them with a temp file in $TMPDIR, and a compute
# node with a full /tmp killed 30 tasks at this line ("cannot create temp
# file for here-document").  Parameter expansion needs no file.
SPEC="${CELLS[$T]}"
METH="${SPEC%%:*}"; REST="${SPEC#*:}"
TRAIT="${REST%%:*}"; CHR="${REST#*:}"; [[ "$CHR" == "$REST" ]] && CHR=""
# and point every tool's scratch at rds rather than the node (rule 10)
export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR"
TAG="${TRAIT//_/}"
echo "=== $METH / $TRAIT ${CHR:+chr$CHR }: $(date)  cpus=${SLURM_CPUS_PER_TASK:-?} ==="
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

case "$METH" in
CT|SBayesR|SBayesRC)
  MA="$GWAS/$TRAIT.ma"
  [[ -s "$MA" ]] || { echo "FATAL: $MA missing (run step9_scz2025_normalise.sbatch)" >&2; exit 3; }
  OUTD="$ROOT/$METH/$TRAIT"; mkdir -p "$OUTD"
  echo "$(( $(wc -l < "$MA") - 1 )) SNPs in $MA" ;;
esac

case "$METH" in
CT)
  awk 'BEGIN{FS=OFS=" "} NR==1{print "SNP","P"; next} {print $1,$7}' "$MA" > "$OUTD/$TRAIT.assoc"
  "$PLINK" --bfile "$GENO" --clump "$OUTD/$TRAIT.assoc" \
    --clump-p1 "$CLUMP_P1" --clump-p2 "$CLUMP_P2" \
    --clump-r2 "$CLUMP_R2" --clump-kb "$CLUMP_KB" --out "$OUTD/$TRAIT"
  awk 'NR>1 && $3!="" {print $3}' "$OUTD/$TRAIT.clumped" > "$OUTD/$TRAIT.index"
  echo "index SNPs: $(wc -l < "$OUTD/$TRAIT.index")"
  for thr in $PRS_THRESHOLDS; do
    tag="${thr//./p}"
    awk -v t="$thr" 'NR==FNR{k[$1]=1; next} FNR==1{next}
                     ($1 in k) && ($7+0)<(t+0) {print $1"\t"$2"\t"$5}' \
      "$OUTD/$TRAIT.index" "$MA" > "$OUTD/${TAG}_$tag.weights"
    "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_$tag.weights" 1 2 3 sum \
      --out "$OUTD/score_${TAG}_$tag" >/dev/null 2>&1 || echo "WARN: score at $thr failed" >&2
    echo "  thr $thr: $(wc -l < "$OUTD/${TAG}_$tag.weights") SNPs"
  done
  rm -f "$OUTD/$TRAIT.assoc" ;;
SBayesR)
  for c in $(seq 1 22); do
    m="$SBR_LDM/ukbEURu_hm3_chr${c}_v3_50k.ldm.sparse"
    [[ -s "$m.info" ]] && echo "$m"
  done > "$OUTD/ldm.txt"
  "$GCTB" --sbayes R --mldm "$OUTD/ldm.txt" --gwas-summary "$MA" --robust \
    --chain-length 10000 --burn-in 4000 --seed 2026 --thread "${SLURM_CPUS_PER_TASK:-8}" \
    --out "$OUTD/$TRAIT" 2>&1 | tail -25
  [[ -s "$OUTD/$TRAIT.snpRes" ]] || { echo "FATAL: no snpRes" >&2; exit 4; }
  awk 'NR>1{print $2"\t"$5"\t"$8}' "$OUTD/$TRAIT.snpRes" > "$OUTD/${TAG}_sbr.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_sbr.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_sbr" ;;
SBayesRC)
  IMP="$OUTD/${TRAIT}_imp"
  [[ -s "$IMP.imputed.ma" ]] || "$GCTB" --ldm-eigen "$SBRC_EIGEN" \
      --gwas-summary "$MA" --impute-summary --thread "${SLURM_CPUS_PER_TASK:-8}" --out "$IMP" 2>&1 | tail -6
  [[ -s "$IMP.imputed.ma" ]] || { echo "FATAL: imputation failed" >&2; exit 3; }
  "$GCTB" --sbayes RC --ldm-eigen "$SBRC_EIGEN" --gwas-summary "$IMP.imputed.ma" \
    --annot "$ANNOT" --chain-length 3000 --burn-in 1000 --seed 2026 --thread "${SLURM_CPUS_PER_TASK:-8}" \
    --out "$OUTD/$TRAIT" 2>&1 | tail -25
  [[ -s "$OUTD/$TRAIT.snpRes" ]] || { echo "FATAL: no snpRes" >&2; exit 4; }
  awk 'NR>1{print $2"\t"$5"\t"$8}' "$OUTD/$TRAIT.snpRes" > "$OUTD/${TAG}_sbrc.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_sbrc.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_sbrc" ;;
PRSCS)
  SST="$GWAS/$TRAIT.sst"; NG=$(cat "$GWAS/$TRAIT.n_gwas")
  [[ -s "$SST" && -n "$NG" ]] || { echo "FATAL: $SST / n_gwas missing" >&2; exit 3; }
  OUTD="$ROOT/PRSCS/$TRAIT/_chr"; mkdir -p "$OUTD"
  out="$OUTD/${TRAIT}_pst_eff_a1_b0.5_phiauto_chr${CHR}.txt"
  if [[ -s "$out" ]]; then echo "exists: $out"; exit 0; fi
  SEED=(); grep -q "'seed='" "$PRSCS" && SEED=(--seed=2026)
  echo "n_gwas=$NG  ref=$LDBLK_UKBB"
  "$PYTHON" "$PRSCS" --ref_dir="$LDBLK_UKBB" --bim_prefix="$GENO" --sst_file="$SST" \
    --n_gwas="$NG" --chrom="$CHR" --out_dir="$OUTD/$TRAIT" "${SEED[@]}"
  [[ -s "$out" ]] || { echo "FATAL: PRS-CS wrote nothing for chr$CHR" >&2; exit 4; } ;;
PRSCSX)
  OUTD="$ROOT/PRSCSX/_chr"; mkdir -p "$OUTD"
  out="$OUTD/SCZ25_META_pst_eff_a1_b0.5_phiauto_chr${CHR}.txt"
  if [[ -s "$out" ]]; then echo "exists: $out"; exit 0; fi
  for a in EUR AFR EAS; do [[ -s "$GWAS/SCZ25_$a.sst" ]] || { echo "FATAL: $GWAS/SCZ25_$a.sst missing" >&2; exit 3; }; done
  NG="$(cat "$GWAS/SCZ25_EUR.n_gwas"),$(cat "$GWAS/SCZ25_AFR.n_gwas"),$(cat "$GWAS/SCZ25_EAS.n_gwas")"
  echo "n_gwas(EUR,AFR,EAS)=$NG  ref=$LD1KG"
  "$PYTHON" "$PRSCSX" --ref_dir="$LD1KG" --bim_prefix="$GENO" \
    --sst_file="$GWAS/SCZ25_EUR.sst,$GWAS/SCZ25_AFR.sst,$GWAS/SCZ25_EAS.sst" \
    --n_gwas="$NG" --pop=EUR,AFR,EAS --chrom="$CHR" --out_dir="$OUTD" --out_name=SCZ25 \
    --meta=TRUE --seed=2026
  [[ -s "$out" ]] || { echo "FATAL: PRS-CSx wrote no META file for chr$CHR" >&2; exit 4; } ;;
*) echo "FATAL: unknown method $METH" >&2; exit 2 ;;
esac
echo "=== $METH / $TRAIT ${CHR:+chr$CHR }done  $(date) ==="
