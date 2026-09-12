#!/bin/bash
#SBATCH --job-name=v2_wancrc
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#
# Within-ancestry standardisation for SBayesRC, which finished AFTER
# within_ancestry_all.sbatch had already run and so was not covered by it.
# Same procedure, same strata file, so the SBayesRC rows are directly
# comparable with the other methods' corrected numbers.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
STRATA="$W/results/ancestry/strata_k4.tsv"
OUTD="$OUT_V2/within_ancestry"
SETUP="$REPO/hpc_v2/work/setup"
mkdir -p "$OUTD"
for trait in SCZ MDD ASD ALZ ALZ_noAPOE; do
  dir="$OUT_V2/sbayesrc/$trait"
  [[ -d "$dir" ]] || { echo "SKIP $trait"; continue; }
  ls "$dir"/score_*.profile >/dev/null 2>&1 || { echo "SKIP $trait (no profile)"; continue; }
  lbl="${trait/_/}"
  wdir="$OUTD/z_SBayesRC_${lbl}"
  echo "=== SBayesRC / $trait  $(date) ==="
  python3 "$SETUP/standardise_within_ancestry.py" "$STRATA" "$dir" "$wdir"
  "$RSCRIPT" "$REPO/tools/prs_assoc.R" --prs-dir "$wdir" --pheno "$PHENO" \
    --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
    --eur-ids "$EUR_IDS" --out "$OUTD/prs_association_SBayesRC_${lbl}_wanc.tsv" || true
  "$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" --prs-dir "$wdir" --pheno "$PHENO" \
    --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
    --eur-ids "$EUR_IDS" --out "$OUTD/prs_withinfamily_SBayesRC_${lbl}_wanc.tsv" || true
done
echo "=== done  $(date) ==="
