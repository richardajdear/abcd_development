#!/bin/bash
#SBATCH --job-name=magma_ahba_c123
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=slurm/%x_%A.log
#
# Step 15: MAGMA gene-property of the AHBA C1-C3 loadings on the SINGLE-LMM
# phenotype gene Z, both atlases.
#
# The gap this closes.  Step 7 ran the C1-C3 gene-property test
# (magma_eur/magma_summary.tsv, rows *_ahba_components_posneg) but only on the
# PER-REGION constructions -- its loop is over magma_eur/, which holds
# per-region .genes.raw only.  Step 13 computed the single-LMM gene analyses
# (magma_panel/genes/*_1lmm.genes.raw) but spent them on the disorder panel,
# never on C1-C3.  Both halves exist; the cell does not.  This script is the
# cell, and it costs seconds: the expensive gene analysis is already done and
# is reused, exactly as step 8 reuses step 7's.
#
# Arms, all reusing committed .genes.raw:
#   eur/1lmm       magma_panel/genes/{baseline_thickness,global_slope}_1lmm   <- the target
#   eur/perregion  magma_eur/{baseline_thickness,global_slope,slope_PC1..3}   <- control, SAME covar file
#   pooled/*       magma_pooled/genes/* (step 14, ABCD's own LD reference)    <- if present
#
# The per-region arm is re-run here rather than quoted from step 7 because
# step 7's covariate file was built from the sibling AHBA repo's Entrez-keyed
# loadings, which are not in this repo; this one is built from committed data
# (make_ahba_c123_gene_covar.py) and maps 7,289 of 7,973 genes.  A single-LMM
# number is only interpretable against a per-region number computed on the
# same gene set, so both are produced in one run.  Expect the per-region rows
# to reproduce step 7 closely but not exactly.
#
# Models per .genes.raw:
#   posneg            C1+ C2+ C3+ C1- C2- C3- jointly -- replicates step 7's
#                     `ahba_components_posneg`.  The split is not cosmetic: a
#                     component's poles are different biology (C1 is sensory
#                     vs association cortex) and one signed covariate averages
#                     an enrichment at one pole against the opposite at the other.
#   marginal_C1/C2/C3 each signed component alone; the C3 row is directly
#                     comparable with step 8's AHBA_C3 (magma_ahba_pls_h3).
# --model direction=both throughout, as in legacy/hpc/04_magma.sbatch:73 -- a
# maturation axis can be enriched at either tail and a one-sided test would
# report half of a real result as null.
#
# STATED PRIOR (carried from step 8 / ahba_pls/FOLLOWUP_GENETICS.md): expect
# null.  C1-C3 is null on every per-region construction (best p 0.068, C3- on
# baseline thickness) and the 1lmm and per-region gene Z track each other very
# closely across the whole step-13 panel.  Report null as null.
#
#   sbatch genetic_analysis/step15_magma_ahba_c123.sbatch
#   bash   genetic_analysis/step15_magma_ahba_c123.sbatch   # runs fine on a login node
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"
COVAR="$REPO/genetic_analysis/hpc/ahba_c123_gene_covar_entrez.txt"
export TMPDIR="${TMPDIR:-$REPO/genetic_analysis/work/magma_scz2025/tmp}"; mkdir -p "$TMPDIR"
[[ -x "$MAGMA" ]] || { echo "FATAL: $MAGMA missing" >&2; exit 2; }
[[ -s "$COVAR" ]] || { echo "FATAL: $COVAR missing -- run make_ahba_c123_gene_covar.py" >&2; exit 2; }

for PARC in dk hcp; do
  SUF=""; [[ "$PARC" == hcp ]] && SUF="_hcp"
  R="$REPO/genetic_analysis/work/results_70tab$SUF"
  OUT="$R/magma_ahba_c123"; CV="$OUT/covar"; mkdir -p "$CV"

  # One file per model.  Passing the 9-column covar file whole and selecting
  # with use= is avoided deliberately: the column names carry + and -, and
  # step 8 (2026-09-15) was already bitten once by MAGMA silently dropping
  # covariate columns rather than failing.  Splitting is unambiguous.
  awk -F'\t' 'NR==1{print "GENE\tC1+\tC2+\tC3+\tC1-\tC2-\tC3-"; next}
              {print $1"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' "$COVAR" > "$CV/posneg.txt"
  for i in 1 2 3; do
    awk -F'\t' -v c=$((i+1)) -v n="C$i" 'NR==1{print "GENE\t" n; next}{print $1"\t"$c}' \
      "$COVAR" > "$CV/C$i.txt"
  done

  # arm:construction:phenotype:genes.raw
  CELLS=()
  for p in baseline_thickness global_slope; do
    CELLS+=("eur:1lmm:$p:$R/magma_panel/genes/${p}_1lmm.genes.raw")
  done
  for p in baseline_thickness global_slope slope_PC1 slope_PC2 slope_PC3; do
    CELLS+=("eur:perregion:$p:$R/magma_eur/$p.genes.raw")
  done
  for p in baseline_thickness global_slope; do
    CELLS+=("pooled:1lmm:$p:$R/magma_pooled/genes/${p}_1lmm.genes.raw")
    CELLS+=("pooled:perregion:$p:$R/magma_pooled/genes/$p.genes.raw")
  done

  echo "=== $PARC  $(date)"
  for SPEC in "${CELLS[@]}"; do
    ARM="${SPEC%%:*}"; S1="${SPEC#*:}"; CON="${S1%%:*}"; S2="${S1#*:}"; P="${S2%%:*}"; RAW="${S2#*:}"
    # The pooled arm (step 14) may not have landed yet; the EUR arm is required.
    if [[ ! -s "$RAW" ]]; then
      [[ "$ARM" == pooled ]] && { echo "  skip $ARM $CON $P (no $(basename "$RAW"))"; continue; }
      echo "FATAL: $RAW missing" >&2; exit 3
    fi
    for m in posneg C1 C2 C3; do
      tag="$OUT/${ARM}_${CON}_${P}_${m}"
      [[ -s "$tag.gsa.out" ]] || "$MAGMA" --gene-results "$RAW" \
        --gene-covar "$CV/$m.txt" --model direction=both --out "$tag" > /dev/null
    done
    echo "  $ARM $CON $P: 4 models done"
  done
done

echo "=== collect  $(date)"
for PARC in dk hcp; do
  SUF=""; [[ "$PARC" == hcp ]] && SUF="_hcp"
  OUT="$REPO/genetic_analysis/work/results_70tab$SUF/magma_ahba_c123"
  { printf 'parcellation\tarm\tconstruction\tphenotype\tmodel\tvariable\tn_genes\tbeta\tbeta_std\tse\tp\n'
    for f in "$OUT"/*.gsa.out; do
      [[ -s "$f" ]] || continue
      b=$(basename "$f" .gsa.out)
      m="${b##*_}"; rest="${b%_*}"               # model is the last field
      arm="${rest%%_*}"; r1="${rest#*_}"
      con="${r1%%_*}"; ph="${r1#*_}"
      [[ "$m" == posneg ]] && model=posneg || model="marginal_$m"
      awk -v pc="$PARC" -v a="$arm" -v c="$con" -v ph="$ph" -v m="$model" '
        /^#/ {next}
        $1=="VARIABLE" {for(i=1;i<=NF;i++) col[$i]=i; next}
        NF>=7 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               pc, a, c, ph, m, $col["VARIABLE"], $col["NGENES"],
               $col["BETA"], $col["BETA_STD"], $col["SE"], $col["P"]}' "$f"
    done; } > "$OUT/table_ahba_c123.tsv"
  echo "wrote $OUT/table_ahba_c123.tsv ($(( $(wc -l < "$OUT/table_ahba_c123.tsv") - 1 )) rows)"
done
echo "=== step 15 done  $(date)"
