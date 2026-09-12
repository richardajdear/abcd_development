#!/bin/bash
#SBATCH --job-name=v2_sbrc
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#
# SBayesRC (Zheng et al. 2024, Nat Genet) -- SBayesR plus functional annotations.
#
# WHY IT MIGHT BEAT SBayesR HERE.  SBayesRC extends the mixture prior with
# per-SNP functional annotations (S-LDSC BaselineLDv2.2, ~96 annotations), so
# variants in functionally important regions are shrunk less.  That helps most
# exactly when discovery power is limited -- the annotations supply information
# the GWAS alone cannot -- which is our situation.  It also runs a low-rank
# model on eigen-decomposed LD, which handles LD better than the sparse
# approximation SBayesR uses.
#
# Same EUR-only .ma inputs as SBayesR (13.8), so the two differ only in the
# prior and the LD representation.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
# GCTB 2.5.5, not 2.5.2: 2.5.2 segfaulted immediately after loading the eigen
# LD and annotations (all inputs matched -- 1,153,913 SNPs, 97 annotations,
# 591 blocks -- then died), and 2.5.5 is the version the SBayesRC tutorial
# now ships.
GCTB=/rds/user/rajd2/hpc-work/gctb/gctb_2.5.5_Linux/gctb
BASE=/rds/user/rajd2/hpc-work/gctb
OUTD="$OUT_V2/sbayesrc"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"
mkdir -p "$OUTD"

case "${SLURM_ARRAY_TASK_ID:-1}" in
  1) D=SCZ ;; 2) D=MDD ;; 3) D=ASD ;; 4) D=ALZ ;;
esac
MA="$OUT_V2/sbayesr/$D.ma"
[[ -s "$MA" ]] || { echo "FATAL: $MA missing (run make_ma_eur.py)" >&2; exit 3; }

# The eigen archive and the annotation archive each unpack to their own
# directory; resolve them rather than assuming a layout.
LDDIR="$BASE/ukbEUR_HM3"
# Name the annotation file explicitly.  A glob for "*annot*baseline*" also
# matched the download's own .status marker, and GCTB read that as the
# annotation: "0 matched SNPs in the annotation file", after which the
# included SNP set was empty and the GWAS matched 0 of 7,659,767 SNPs.
ANNOT="$BASE/annot_baseline2.2.txt"
echo "LD eigen dir: $LDDIR"
echo "annotations : $ANNOT"
[[ -n "$LDDIR" && -d "$LDDIR" ]] || { echo "FATAL: eigen LD dir not found" >&2; exit 3; }
[[ -s "$ANNOT" ]] || { echo "FATAL: annotation file not found at $ANNOT" >&2; exit 3; }
echo "annotation rows: $(wc -l < "$ANNOT"), columns: $(head -1 "$ANNOT" | awk '{print NF}')"
echo "=== SBayesRC $D: $(( $(wc -l < "$MA") - 1 )) SNPs  $(date) ==="

# SBayesRC is a TWO-STEP workflow and skipping step 1 is what produced the
# segfault.  The method requires a summary statistic for every SNP in the LD
# reference; --impute-summary fills in the ones the GWAS lacks from the LD
# structure.  Only then can --sbayes RC run.
# NOTE the output name: GCTB writes "<--out>.imputed.ma", so with
# --out X_imputed the file is X_imputed.imputed.ma.  Checking for X_imputed.ma
# silently left MA_USE pointing at the UN-imputed file, which is exactly the
# input that segfaults step 2.
IMP="$OUTD/${D}_imputed"
IMPMA="$IMP.imputed.ma"
if [[ ! -s "$IMPMA" ]]; then
  echo "--- step 1: impute summary statistics  $(date) ---"
  "$GCTB" --ldm-eigen "$LDDIR" --gwas-summary "$MA" --impute-summary \
    --thread 8 --out "$IMP" 2>&1 | tail -15
fi
[[ -s "$IMPMA" ]] || { echo "FATAL: imputation produced no $IMPMA" >&2; exit 3; }
MA_USE="$IMPMA"
echo "using imputed sumstats: $(( $(wc -l < "$IMPMA") - 1 )) SNPs"

echo "--- step 2: SBayesRC  $(date) ---"
"$GCTB" --sbayes RC \
  --ldm-eigen "$LDDIR" \
  --gwas-summary "$MA_USE" \
  --annot "$ANNOT" \
  --chain-length 3000 --burn-in 1000 \
  --seed 2026 --thread 8 \
  --out "$OUTD/$D" 2>&1 | tail -40

SNPRES="$OUTD/$D.snpRes"
[[ -s "$SNPRES" ]] || { echo "FATAL: $D produced no .snpRes" >&2; exit 4; }
echo "posterior effects: $(( $(wc -l < "$SNPRES") - 1 )) SNPs"

mkdir -p "$OUTD/$D"
awk 'NR>1{print $2"\t"$5"\t"$8}' "$SNPRES" > "$OUTD/$D/score_${D}_sbayesrc.txt"
"$W/bin/plink" --bfile "$PRS_GENO" --score "$OUTD/$D/score_${D}_sbayesrc.txt" 1 2 3 sum \
  --out "$OUTD/$D/score_${D}_sbrc"

"$RSCRIPT" "$REPO/tools/prs_assoc.R" --prs-dir "$OUTD/$D" --pheno "$PHENO" \
  --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
  --eur-ids "$EUR_IDS" --out "$OUTD/prs_association_${D}_sbayesrc.tsv" || true
"$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" --prs-dir "$OUTD/$D" --pheno "$PHENO" \
  --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
  --eur-ids "$EUR_IDS" --out "$OUTD/prs_withinfamily_${D}_sbayesrc.tsv" || true
echo "=== $D done  $(date) ==="
