#!/bin/bash
#SBATCH --job-name=v2_csukbb
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-88
#
# PRS-CS re-run against the UK BIOBANK LD reference (ldblk_ukbb_eur).
#
# The 1000G-reference run in prscs/ is kept.  Only the LD panel changes
# (503 individuals -> ~375k), plus SCZ moving from .core. to .european. so all
# three methods share one discovery GWAS.  See prscs_ukbb_prep.sbatch.
#
# PRS-CS: Bayesian polygenic scores with a continuous-shrinkage prior.
#
# WHY REPLACE CLUMPING+THRESHOLDING.  v1's scores use P+T: clump to
# approximately independent index SNPs, keep those below a p-value threshold,
# sum their effect sizes.  That discards almost every variant (at p<0.05,
# ~74k of 7M) and treats each survivor's marginal effect as if it were causal,
# even though it is inflated by the winner's curse and by LD with neighbours.
#
# PRS-CS instead keeps ALL HapMap3 SNPs and re-estimates every effect jointly,
# placing a continuous-shrinkage prior on the effects and modelling local LD
# from an external reference.  Effects with weak evidence shrink toward zero;
# effects supported across an LD block survive.  The global shrinkage parameter
# phi is learned from the data (PRS-CS-auto), so there is no threshold to tune
# and no multiple-testing over eight thresholds.  Published gains are typically
# 1.5-2x the variance explained on the same discovery data, which at our effect
# size (SCZ PRS r2 ~ 0.17%) is the cheapest real power available.
#
# The comparison this enables is like-for-like: same discovery GWAS, same target
# genotypes, same association model downstream -- only the weighting changes.
#
# Array: 1 = SCZ, 2 = MDD, 3 = ALZ (late-onset control), 4 = ASD (childhood
# control).  Running the controls through the SAME method matters; a specificity
# claim built on a better-powered score for SCZ and a weaker one for the
# controls would be an artefact of the method, not of biology.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
PRSCS_DIR=/rds/user/rajd2/hpc-work/prscs
PY="$W/envs/abcd/bin/python"
OUTD="$OUT_V2/prscs_ukbb"
mkdir -p "$OUTD"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"

# One task per (disorder, chromosome): 4 x 22 = 88.  PRS-CS's MCMC is serial
# and takes ~1-3 h for a large chromosome, so chromosome-level parallelism is
# what keeps this inside the 12 h QOS wall limit.
T=$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))
case $(( T / 22 )) in
  0) D=SCZ ;;
  1) D=MDD ;;
  2) D=ALZ ;;
  3) D=ASD ;;
esac
CHR=$(( T % 22 + 1 ))
echo "=== PRS-CS $D  $(date) ==="

# --- 1. sumstats in PRS-CS format: SNP A1 A2 BETA P ---------------------------
# A1 is the EFFECT allele throughout.  PRS-CS also needs the discovery N as a
# scalar (--n_gwas), so the per-SNP N column is summarised to its median.
SS="$OUTD/${D}_prscs_input.txt"
[[ -s "$SS" ]] || { echo "FATAL: $SS missing -- run prscs_ukbb_prep.sbatch first" >&2; exit 3; }
NGWAS=$(cat "$OUTD/${D}_n_gwas.txt")
echo "$D: $(( $(wc -l < "$SS") - 1 )) SNPs, n_gwas=$NGWAS"

# --- 2. posterior effect sizes, one job per chromosome ------------------------
# Limiting BLAS threads is PRS-CS's own recommendation: the MCMC is serial and
# oversubscribed BLAS makes it slower, not faster.
export MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$OUTD/$D"
out="$OUTD/$D/${D}_pst_eff_a1_b0.5_phiauto_chr${CHR}.txt"
if [[ ! -s "$out" ]]; then
  "$PY" "$PRSCS_DIR/PRScs/PRScs.py" \
    --ref_dir="$PRSCS_DIR/ldblk_ukbb_eur" \
    --bim_prefix="$PRS_GENO" \
    --sst_file="$SS" \
    --n_gwas="$NGWAS" \
    --chrom="$CHR" \
    --out_dir="$OUTD/$D/$D"
else
  echo "  chr$CHR already done"
fi

# --- 3. score + associate, but ONLY from the task holding the last chromosome,
# and only once every chromosome for this disorder has landed.  Any task may be
# last, so each one checks rather than assuming.
n_done=$(ls "$OUTD/$D"/${D}_pst_eff_a1_b0.5_phiauto_chr*.txt 2>/dev/null | wc -l)
if [[ "$n_done" -lt 22 ]]; then
  echo "$D: $n_done/22 chromosomes done; scoring deferred to the last task"
  exit 0
fi
echo "$D: all 22 chromosomes present, scoring"
cat "$OUTD/$D"/${D}_pst_eff_a1_b0.5_phiauto_chr*.txt > "$OUTD/$D/${D}_all.weights"
echo "$D posterior weights: $(wc -l < "$OUTD/$D/${D}_all.weights") SNPs"
# PRS-CS output columns: CHR SNP BP A1 A2 BETA -> --score needs SNP A1 BETA
awk 'BEGIN{OFS="\t"} {print $2,$4,$6}' "$OUTD/$D/${D}_all.weights" > "$OUTD/$D/score_${D}_prscs.txt"
"$W/bin/plink" --bfile "$PRS_GENO" \
  --score "$OUTD/$D/score_${D}_prscs.txt" 1 2 3 sum \
  --out "$OUTD/$D/score_${D}_auto"

# prs_assoc.R globs score_*.profile in one directory, so each disorder's
# PRS-CS score lives in its own subdirectory and is associated separately.
"$RSCRIPT" "$REPO/tools/prs_assoc.R" \
  --prs-dir "$OUTD/$D" --pheno "$PHENO" --covar-quant "$COVAR_QUANT" \
  --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_association_${D}_prscs_ukbb.tsv" || true
"$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" \
  --prs-dir "$OUTD/$D" --pheno "$PHENO" --covar-quant "$COVAR_QUANT" \
  --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_withinfamily_${D}_prscs_ukbb.tsv" || true
echo "=== $D done  $(date) ==="
