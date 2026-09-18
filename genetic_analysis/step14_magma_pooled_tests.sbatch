#!/bin/bash
#SBATCH --job-name=magma_pooled_tests
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/%x_%A.log
# Step 14: the step-13 panel tests on the pooled-arm phenotype gene results -> results_70tab*/magma_pooled/
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"; PY="$REPO/legacy/hpc/work/envs/abcd/bin/python"
ROOT="$REPO/genetic_analysis/work/magma_scz2025"; G="$ROOT/genes"; export TMPDIR="$ROOT/tmp"
DIS=(SCZ25_EUR SCZ25_META SCZ25_METAnaive PGC3_EUR PGC3_primary MDD_EUR MDD_div ASD ALZ ALZ_noAPOE EA)
SETS_SCZ="$ROOT/genesets_scz2025.txt"; SETS_MDD="$REPO/legacy/hpc/work/genesets/mdd_highconf.txt"
genez() { awk -v g="$2" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; print "GENE\t" g; next} {print $c["GENE"] "\t" $c["ZSTAT"]}' "$1" > "$3"; }
for PARC in dk hcp; do
  S=""; [[ "$PARC" == hcp ]] && S="_hcp"; OUT="$REPO/genetic_analysis/work/results_70tab$S/magma_pooled"
  for p in global_slope baseline_thickness global_slope_1lmm baseline_thickness_1lmm; do
    raw="$OUT/genes/$p.genes.raw"; go="$OUT/genes/$p.genes.out"; [[ -s "$raw" ]] || { echo "FATAL: $raw" >&2; exit 3; }
    genez "$go" "$p" "$OUT/${p}_genez.txt"
    for D in "${DIS[@]}"; do
      [[ -s "$OUT/${D}_vs_${p}.gsa.out" ]] || "$MAGMA" --gene-results "$G/$D.genes.raw" --gene-covar "$OUT/${p}_genez.txt" --model direction=both --out "$OUT/${D}_vs_${p}" > /dev/null
      [[ -s "$OUT/${p}_on_${D}.gsa.out" ]] || "$MAGMA" --gene-results "$raw" --gene-covar "$G/${D}_genez.txt" --model direction=both --out "$OUT/${p}_on_${D}" > /dev/null
    done
    "$MAGMA" --gene-results "$raw" --set-annot "$SETS_SCZ" col=1,2 --model direction=greater --out "$OUT/${p}_sets_scz" > /dev/null
    "$MAGMA" --gene-results "$raw" --set-annot "$SETS_MDD" col=1,2 --model direction=greater --out "$OUT/${p}_sets_mdd" > /dev/null
    echo "  $PARC $p done $(date +%T)"
  done
done
MAGMA_PANEL_DIR=magma_pooled "$PY" "$REPO/genetic_analysis/step13_magma_panel_collect.py"
