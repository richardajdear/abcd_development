#!/bin/bash
#SBATCH --job-name=prs_1lmm_controls
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-2
#SBATCH --output=slurm/%x_%A_%a.log
#
# Order-of-operations arm, control panel (rule 7): the ASD, ALZ (Wightman,
# +/- APOE), ALZ_IGAP (Kunkle, +/- APOE) and EA scores -- the legacy PGC3-era
# score profiles step 6 reuses -- against the single-LMM phenotypes
# (global_slope_1lmm, baseline_thickness_1lmm) built by the 2026-09-15
# order-of-operations check.  The SCZ/MDD cells were run by
# work/pc_compare/prs_1lmm.sbatch; this completes the panel so the single-LMM
# construction is read against the same controls as the per-region one.
# Task 1 = DK, task 2 = HCP-MMP.  Output: results_70tab*/prs_final_1lmm/assoc_*.tsv,
# then step6_orderops_collect.py writes prs_final_1lmm/table_order_of_operations_all.tsv.
set -eo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development; cd "$REPO"
export TMPDIR="$REPO/genetic_analysis/work/tmp"
R=/home/rajd2/rds/hpc-work/envs/genesis/bin/Rscript
PY="$REPO/legacy/hpc/work/envs/abcd/bin/python"
SUF=""; [[ "${SLURM_ARRAY_TASK_ID:-1}" == 2 ]] && SUF="_hcp"
OUT="$REPO/genetic_analysis/work/results_70tab$SUF/prs_final_1lmm"; PH="$OUT/pheno"
[[ -s "$PH/phenotypes_gcta.txt" ]] || { echo "FATAL: $PH/phenotypes_gcta.txt missing" >&2; exit 2; }
EUR="$REPO/legacy/hpc/work/results/ancestry/eur_anchor.keep"; SC="$REPO/legacy/hpc_v2/work/results_v2/prs_final"
n=0
for m in CT PRSCS SBayesR SBayesRC; do for t in ASD ALZ ALZ_noAPOE ALZ_IGAP ALZ_IGAP_noAPOE EA; do
  for z in "" _zanc; do
    d="$SC/$m/$t${z:+/_zanc}"; [[ -d "$d" ]] || { echo "  (no dir: $d)"; continue; }
    echo "=== $m/$t$z  $(date +%T)"
    "$R" tools/prs_assoc.R --prs-dir "$d" --pheno "$PH/phenotypes_gcta.txt" --covar-quant "$PH/covar_quant.txt" \
      --covar-cat "$PH/covar_categorical.txt" --manifest "$PH/phenotype_manifest.tsv" --eur-ids "$EUR" \
      --out "$OUT/assoc_${m}_${t}${z}.tsv" >/dev/null 2>&1 && n=$((n+1)) || echo "  (failed: $m/$t$z)"
  done
done; done
echo "=== $n new association tables in $OUT  $(date +%T)"
"$PY" genetic_analysis/step6_orderops_collect.py
