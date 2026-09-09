#!/bin/bash
#SBATCH --job-name=v2_sbr
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#
# SBayesR (GCTB) on EUROPEAN-ONLY discovery GWAS, all four traits.
#
# WHY A THIRD SCORING METHOD.  C+T (best of 8 thresholds, selected on the
# outcome) and PRS-CS (one untuned Bayesian fit) disagree about SCZ ->
# global_slope: 4.3e-03 against 0.195.  SBayesR is a genuinely different prior
# -- a finite mixture of normals rather than continuous shrinkage -- fitted
# against a shrunk sparse UKB LD matrix rather than 1000G blocks.  If it lands
# with PRS-CS, the C+T result is threshold selection; if it lands with C+T, the
# PRS-CS attenuation is over-shrinkage.  Either way it adjudicates.
#
# ANCESTRY.  Every input is the European-ancestry release, verified rather than
# assumed -- see make_ma_eur.py's header.  Notably SCZ uses `.european.`, NOT
# `.core.`: the PGC3 README defines core as "cohorts primarily of east asian
# and european ancestry", so the earlier arms that used core were not the
# EUR-only comparison they were labelled as.
#
# SBayesR is known to be sensitive to summary-statistic QC and to LD mismatch.
# --robust is used (the recommended parameterisation for public sumstats).
# --num-chains is NOT used: GCTB 2.5.2 answers it with "R is not currently
# available in the multi-chain Bayesian analysis" and writes no output.  A
# failure to converge is reported rather than silently accepted.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
GCTB=/rds/user/rajd2/hpc-work/gctb/gctb_2.5.2_Linux/gctb
# The archive unpacks to ukbEURu_hm3_shrunk_sparse (not ..._sparse), with files
# named ukbEURu_hm3_chr<N>_v3_50k.ldm.sparse.{info,bin}.
LDM=/rds/user/rajd2/hpc-work/gctb/ukbEURu_hm3_shrunk_sparse
OUTD="$OUT_V2/sbayesr"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"

case "${SLURM_ARRAY_TASK_ID:-1}" in
  1) D=SCZ ;; 2) D=MDD ;; 3) D=ASD ;; 4) D=ALZ ;;
esac
MA="$OUTD/$D.ma"
[[ -s "$MA" ]] || { echo "FATAL: $MA missing" >&2; exit 3; }
echo "=== SBayesR $D: $(( $(wc -l < "$MA") - 1 )) SNPs  $(date) ==="

# The shrunk sparse reference ships one .ldm.sparse per chromosome; GCTB takes
# a multi-chromosome LD matrix list.
LDLIST="$OUTD/${D}_ldm.txt"
# Autosomes only -- our genotypes carry no chrX.
for c in $(seq 1 22); do
  m="$LDM/ukbEURu_hm3_chr${c}_v3_50k.ldm.sparse"
  [[ -s "$m.info" && -s "$m.bin" ]] && echo "$m"
done > "$LDLIST"
echo "LD matrices: $(wc -l < "$LDLIST")"
head -2 "$LDLIST"

"$GCTB" --sbayes R \
  --mldm "$LDLIST" \
  --gwas-summary "$MA" \
  --robust \
  --chain-length 10000 --burn-in 4000 \
  --seed 2026 \
  --thread 8 \
  --out "$OUTD/$D" 2>&1 | tail -40

SNPRES="$OUTD/$D.snpRes"
if [[ ! -s "$SNPRES" ]]; then
  echo "FATAL: $D produced no .snpRes -- SBayesR did not converge or the input was rejected" >&2
  exit 4
fi
echo "posterior effects: $(( $(wc -l < "$SNPRES") - 1 )) SNPs"

# .snpRes columns: Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP
mkdir -p "$OUTD/$D"
awk 'NR>1{print $2"\t"$5"\t"$8}' "$SNPRES" > "$OUTD/$D/score_${D}_sbayesr.txt"
"$W/bin/plink" --bfile "$PRS_GENO" \
  --score "$OUTD/$D/score_${D}_sbayesr.txt" 1 2 3 sum \
  --out "$OUTD/$D/score_${D}_sbr"

"$RSCRIPT" "$REPO/tools/prs_assoc.R" --prs-dir "$OUTD/$D" \
  --pheno "$PHENO" --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" \
  --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_association_${D}_sbayesr.tsv" || true
"$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" --prs-dir "$OUTD/$D" \
  --pheno "$PHENO" --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" \
  --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_withinfamily_${D}_sbayesr.tsv" || true
echo "=== $D done  $(date) ==="
