#!/bin/bash
# 03 -- C3-C PRS screen on whatever per-child score files are present locally (README_HPC.md 8.5).
# Same model and inputs as the cluster: tools/prs_assoc.R on the c3axis export, with and without
# the global slope as a covariate. Then 03_prs_screen_collect.py writes the matched-cell tables.
# Run from the repo root after build_c3axis_pheno.py has written $C3/pheno{,_adjG}:
#   python genetic_analysis/c3axis/build_c3axis_pheno.py out/thickness_hcp_70_aa6e91efba82_c3axis \
#       genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno genetic_analysis/work/results_70tab_hcp/c3axis/pheno
#   bash genetic_analysis/c3axis/03_prs_screen_local.sh && python genetic_analysis/c3axis/03_prs_screen_collect.py
set -eo pipefail
C3=genetic_analysis/work/results_70tab_hcp/c3axis
EUR=legacy/hpc/work/results/ancestry/eur_anchor.keep
RSCRIPT=${RSCRIPT:-Rscript}
mkdir -p "$C3/prs_local" "$C3/prs_local_adjG"
for d in $(find genetic_analysis/work/scores_scz2025 legacy/hpc_v2/work/results_v2/prs_final -name 'score_*.profile' -exec dirname {} \; | sort -u); do
  tag=$(echo "$d" | sed 's|.*scores_scz2025/||; s|.*prs_final/||; s|/|__|g')
  for v in "" _adjG; do
    extra=(); [[ $v == _adjG ]] && extra=(--extra-covar gcov_global_slope)
    "$RSCRIPT" tools/prs_assoc.R --prs-dir "$d" --pheno "$C3/pheno$v/phenotypes_gcta.txt" \
      --covar-quant "$C3/pheno$v/covar_quant.txt" --covar-cat "$C3/pheno$v/covar_categorical.txt" \
      --manifest "$C3/pheno$v/phenotype_manifest.tsv" --eur-ids "$EUR" "${extra[@]}" \
      --out "$C3/prs_local$v/assoc_$tag.tsv" > "$C3/prs_local$v/log_$tag.txt" 2>&1 || echo "FAIL $tag$v"
  done
done
echo "$(ls $C3/prs_local/assoc_*.tsv | wc -l) score cells"
