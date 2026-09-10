#!/bin/bash
#SBATCH --job-name=v2_final
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-24%8
#
# FINAL PRS MATRIX: four methods x the discovery GWAS matched to the target arm.
#
#   C+T        clumping + thresholding, best of 8 thresholds
#   PRS-CS     continuous shrinkage, UKB LD reference (ldblk_ukbb_eur)
#   SBayesR    mixture of normals, shrunk-sparse UKB LD
#   SBayesRC   mixture of normals + 97 functional annotations, eigen UKB LD
#
# ARM MATCHING, which is the point of this run:
#   pooled arm (n=8,082, multi-ancestry target) uses the MULTI-ANCESTRY GWAS
#     SCZ primary (EUR+EAS+AFR+LAT), MDD div (trans-ancestry)
#   EUR arm (n=4,116) uses the EUROPEAN GWAS
#     SCZ european, MDD eur
# Earlier runs mixed these up: C+T inherited primary/div from v1 while every
# Bayesian method used european/eur, so "method" and "discovery GWAS" were
# confounded in the comparison.
#
# CAVEAT, stated because it is the mirror of the problem this project has been
# chasing: PRS-CS, SBayesR and SBayesRC all use a EUROPEAN LD reference (UKB
# EUR), because no multi-ancestry LD panel exists for them.  Applying them to
# the multi-ancestry primary/div sumstats is therefore LD-mismatched on the
# DISCOVERY side.  The pooled-arm Bayesian numbers carry that caveat; C+T does
# not, since it clumps on our own sample.
#
# ASD and ALZ have no ancestry-stratified release, so one file serves both arms.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
MA_DIR="$FINAL/gwas"
mkdir -p "$FINAL"

# task -> (method, trait-arm label)
T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
METHODS=(CT PRSCS SBayesR SBayesRC)
TRAITS=(SCZ_pooled SCZ_eur MDD_pooled MDD_eur ASD ALZ)
# The full 4 x 6 grid.  Every cell is recomputed from the unified .ma files so
# that no number in the final tables inherits an older, differently-prepared
# input -- ASD and ALZ in particular changed a lot once their allele
# frequencies came from our own target rather than the HapMap3 panel
# (ASD 378,844 -> 2,389,935 SNPs; ALZ 1,110,485 -> 5,994,563).
TODO=()
for m in "${METHODS[@]}"; do for t in "${TRAITS[@]}"; do TODO+=("$m:$t"); done; done
SPEC="${TODO[$T]}"
METH="${SPEC%%:*}"; TRAIT="${SPEC##*:}"
MA="$MA_DIR/${TRAIT}.ma"
[[ -s "$MA" ]] || { echo "FATAL: $MA missing (run normalise_gwas.py)" >&2; exit 3; }
# prs_assoc.R splits the score filename on "_" to get disorder/threshold, so
# the trait token must not contain one: SCZ_pooled -> SCZpooled.
TAG="${TRAIT//_/}"
OUTD="$FINAL/${METH}/${TRAIT}"; mkdir -p "$OUTD"
echo "=== $METH / $TRAIT : $(( $(wc -l < "$MA") - 1 )) SNPs  $(date) ==="

case "$METH" in
CT)
  # PLINK --clump wants columns literally named SNP and P.
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
      --out "$OUTD/score_${TAG}_$tag" >/dev/null 2>&1 || true
  done ;;
PRSCS)
  awk 'BEGIN{FS=" "; OFS="\t"} NR==1{print "SNP","A1","A2","BETA","P"; next}
       {print $1,$2,$3,$5,$7}' "$MA" > "$OUTD/sst.txt"
  # per-SNP N varies; PRS-CS takes one number, so use the max (the full
  # discovery sample) rather than whichever row happens to come first.
  NG=$(awk 'NR>1 && $8+0>m{m=$8+0} END{printf "%d", m}' "$MA")
  echo "n_gwas=$NG"
  export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
  for chr in $(seq 1 22); do
    out="$OUTD/${TRAIT}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
    [[ -s "$out" ]] && continue
    "$PYTHON" "$PRSCS" --ref_dir="$LDBLK_UKBB" --bim_prefix="$GENO" \
      --sst_file="$OUTD/sst.txt" --n_gwas="$NG" --chrom="$chr" \
      --out_dir="$OUTD/$TRAIT" || echo "chr$chr FAILED" >&2
  done
  cat "$OUTD"/${TRAIT}_pst_eff_*chr*.txt | awk 'BEGIN{OFS="\t"}{print $2,$4,$6}' > "$OUTD/${TAG}_prscs.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_prscs.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_prscs" ;;
SBayesR)
  for c in $(seq 1 22); do
    m="$SBR_LDM/ukbEURu_hm3_chr${c}_v3_50k.ldm.sparse"
    [[ -s "$m.info" ]] && echo "$m"
  done > "$OUTD/ldm.txt"
  "$GCTB" --sbayes R --mldm "$OUTD/ldm.txt" --gwas-summary "$MA" --robust \
    --chain-length 10000 --burn-in 4000 --seed 2026 --thread 8 \
    --out "$OUTD/$TRAIT" 2>&1 | tail -25
  [[ -s "$OUTD/$TRAIT.snpRes" ]] || { echo "FATAL: no snpRes" >&2; exit 4; }
  awk 'NR>1{print $2"\t"$5"\t"$8}' "$OUTD/$TRAIT.snpRes" > "$OUTD/${TAG}_sbr.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_sbr.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_sbr" ;;
SBayesRC)
  IMP="$OUTD/${TRAIT}_imp"
  [[ -s "$IMP.imputed.ma" ]] || "$GCTB" --ldm-eigen "$SBRC_EIGEN" \
      --gwas-summary "$MA" --impute-summary --thread 8 --out "$IMP" 2>&1 | tail -6
  [[ -s "$IMP.imputed.ma" ]] || { echo "FATAL: imputation failed" >&2; exit 3; }
  "$GCTB" --sbayes RC --ldm-eigen "$SBRC_EIGEN" --gwas-summary "$IMP.imputed.ma" \
    --annot "$ANNOT" --chain-length 3000 --burn-in 1000 --seed 2026 --thread 8 \
    --out "$OUTD/$TRAIT" 2>&1 | tail -25
  [[ -s "$OUTD/$TRAIT.snpRes" ]] || { echo "FATAL: no snpRes" >&2; exit 4; }
  awk 'NR>1{print $2"\t"$5"\t"$8}' "$OUTD/$TRAIT.snpRes" > "$OUTD/${TAG}_sbrc.weights"
  "$PLINK" --bfile "$GENO" --score "$OUTD/${TAG}_sbrc.weights" 1 2 3 sum --out "$OUTD/score_${TAG}_sbrc" ;;
esac

# --- associate: raw, and with the score standardised within ancestry cluster --
"$RSCRIPT" "$ASSOC_R" --prs-dir "$OUTD" --pheno "$PHENO" --covar-quant "$COVQ" \
  --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EURKEEP" \
  --out "$FINAL/assoc_${METH}_${TRAIT}.tsv" || true
"$RSCRIPT" "$FAM_R" --prs-dir "$OUTD" --pheno "$PHENO" --covar-quant "$COVQ" \
  --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EURKEEP" \
  --out "$FINAL/fam_${METH}_${TRAIT}.tsv" || true
ZD="$OUTD/_zanc"; mkdir -p "$ZD"
"$PYTHON" "$V2ROOT/work/setup/standardise_within_ancestry.py" "$STRATA" "$OUTD" "$ZD"
"$RSCRIPT" "$ASSOC_R" --prs-dir "$ZD" --pheno "$PHENO" --covar-quant "$COVQ" \
  --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EURKEEP" \
  --out "$FINAL/assoc_${METH}_${TRAIT}_zanc.tsv" || true
"$RSCRIPT" "$FAM_R" --prs-dir "$ZD" --pheno "$PHENO" --covar-quant "$COVQ" \
  --covar-cat "$COVC" --manifest "$MANIF" --eur-ids "$EURKEEP" \
  --out "$FINAL/fam_${METH}_${TRAIT}_zanc.tsv" || true
echo "=== $METH / $TRAIT done  $(date) ==="
