#!/bin/bash
#SBATCH --job-name=v2_wanc
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#
# Re-associate EVERY pooled polygenic score after standardising it WITHIN
# ancestry cluster.  See standardise_within_ancestry.py for the rationale.
#
# Covers all five scoring methods so the correction is applied uniformly and
# the methods stay comparable:
#   C+T            SCZ, MDD (v1's scores) and ASD, ALZ, ALZnoAPOE (v2 controls)
#   PRS-CS 1000G   SCZ, MDD, ASD, ALZ
#   PRS-CS UKB     SCZ, MDD, ASD, ALZ, ALZnoAPOE
#   SBayesR        SCZ, MDD, ASD, ALZ, ALZnoAPOE
#   PRS-CSx        EUR, EAS, EUR+EAS sum
#
# Output tables land in results_v2/within_ancestry/ as
#   prs_association_<method>_<trait>_wanc.tsv
#   prs_withinfamily_<method>_<trait>_wanc.tsv
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
STRATA="$W/results/ancestry/strata_k4.tsv"
OUTD="$OUT_V2/within_ancestry"
SETUP="$REPO/hpc_v2/work/setup"
mkdir -p "$OUTD"

run_one () {   # run_one <method> <trait> <profile_dir>
  local meth="$1" trait="$2" dir="$3"
  [[ -d "$dir" ]] || { echo "SKIP $meth/$trait: $dir absent"; return 0; }
  ls "$dir"/score_*.profile >/dev/null 2>&1 || { echo "SKIP $meth/$trait: no profiles in $dir"; return 0; }
  local wdir="$OUTD/z_${meth}_${trait}"
  echo "=== $meth / $trait  $(date) ==="
  python3 "$SETUP/standardise_within_ancestry.py" "$STRATA" "$dir" "$wdir"
  "$RSCRIPT" "$REPO/tools/prs_assoc.R" --prs-dir "$wdir" --pheno "$PHENO" \
    --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
    --eur-ids "$EUR_IDS" --out "$OUTD/prs_association_${meth}_${trait}_wanc.tsv" || true
  "$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" --prs-dir "$wdir" --pheno "$PHENO" \
    --covar-quant "$COVAR_QUANT" --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
    --eur-ids "$EUR_IDS" --out "$OUTD/prs_withinfamily_${meth}_${trait}_wanc.tsv" || true
}

# C+T -- v1's SCZ/MDD scores live together in one directory
run_one CT      SCZMDD    "$W/results/prs_imp"
run_one CT      ASD       "$OUT_V2/control_ASD/prs"
run_one CT      ALZ       "$OUT_V2/control_ALZ/prs"
# PRS-CS on the 1000G reference
for D in SCZ MDD ASD ALZ; do run_one PRSCS1kG "$D" "$OUT_V2/prscs/$D"; done
# PRS-CS on the UKB reference
for D in SCZ MDD ASD ALZ; do run_one PRSCSukbb "$D" "$OUT_V2/prscs_ukbb/$D"; done
run_one PRSCSukbb ALZnoAPOE "$OUT_V2/prscs_ukbb/ALZ_noAPOE"
# SBayesR
for D in SCZ MDD ASD ALZ; do run_one SBayesR "$D" "$OUT_V2/sbayesr/$D"; done
run_one SBayesR ALZnoAPOE "$OUT_V2/sbayesr/ALZ_noAPOE"
# PRS-CSx
for v in EUR EAS SUM; do run_one PRSCSx "$v" "$OUT_V2/prscsx/csx_$v"; done

echo "=== all done  $(date) ==="
ls "$OUTD"/*.tsv | wc -l
