#!/bin/bash
#SBATCH --job-name=magma_panel
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/%x_%A.log
#
# Step 13: MAGMA, full disorder panel x both phenotype constructions x both
# atlases.  Phenotype side: EUR-arm GWAS (LD-matched to the 1000G EUR panel);
# the per-region gene results are step 7's (magma_eur/), the single-LMM ones
# are computed here from scan_1lmm/assoc_eur on the step-10 union annotation.
# Disorder side: step 10/11 gene results in work/magma_scz2025/genes/.
# Tests per (phenotype, disorder):
#   reverse gene-property  disorder genes ~ phenotype gene Z   (legacy 04_magma step 3)
#   forward gene-property  phenotype genes ~ disorder gene Z
# and per phenotype the competitive gene-set test on the SCZ sets (step 10 file)
# and the MDD sets (legacy mdd_highconf.txt), direction=greater.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"; PY="$REPO/legacy/hpc/work/envs/abcd/bin/python"
REF=/rds/user/rajd2/hpc-work/magma/reference_data
ROOT="$REPO/genetic_analysis/work/magma_scz2025"; G="$ROOT/genes"
ANNOT="$ROOT/ref/ncbi37.window35-10.union3.genes.annot"
export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR"
DIS=(SCZ25_EUR SCZ25_META SCZ25_METAnaive PGC3_EUR PGC3_primary MDD_EUR MDD_div ASD ALZ ALZ_noAPOE EA)
SETS_SCZ="$ROOT/genesets_scz2025.txt"; SETS_MDD="$REPO/legacy/hpc/work/genesets/mdd_highconf.txt"
genez() { awk -v g="$2" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; print "GENE\t" g; next} {print $c["GENE"] "\t" $c["ZSTAT"]}' "$1" > "$3"; }
for D in "${DIS[@]}"; do [[ -s "$G/$D.genes.raw" ]] || { echo "FATAL: $G/$D.genes.raw missing" >&2; exit 2; }; genez "$G/$D.genes.out" "$D" "$G/${D}_genez.txt"; done

for PARC in dk hcp; do
  SUF=""; [[ "$PARC" == hcp ]] && SUF="_hcp"; R="$REPO/genetic_analysis/work/results_70tab$SUF"
  OUT="$R/magma_panel"; mkdir -p "$OUT/genes"
  echo "=== $PARC: single-LMM phenotype gene analyses  $(date)"
  for p in global_slope_1lmm baseline_thickness_1lmm; do
    if [[ ! -s "$OUT/genes/$p.genes.raw" ]]; then
      ss="$R/scan_1lmm/assoc_eur/$p.sumstats.tsv.gz"; [[ -s "$ss" ]] || { echo "FATAL: $ss missing" >&2; exit 3; }
      zcat "$ss" | awk 'BEGIN{OFS="\t"} NR==1{print "SNP","P","N"; next}{print $1,$10,$7}' > "$OUT/genes/$p.snp_p.tsv"
      "$MAGMA" --bfile "$REF/g1000_eur" --gene-annot "$ANNOT" --pval "$OUT/genes/$p.snp_p.tsv" use=SNP,P ncol=N --out "$OUT/genes/$p" | tail -3
      rm -f "$OUT/genes/$p.snp_p.tsv"
    fi
    ln -sf "$OUT/genes/$p.genes.raw" "$OUT/genes/${p}.link.raw" 2>/dev/null || true
  done
  # phenotype gene results: per-region from step 7, single LMM from above
  declare -A PRAW=( [global_slope]="$R/magma_eur/global_slope.genes.raw" [baseline_thickness]="$R/magma_eur/baseline_thickness.genes.raw"
                    [global_slope_1lmm]="$OUT/genes/global_slope_1lmm.genes.raw" [baseline_thickness_1lmm]="$OUT/genes/baseline_thickness_1lmm.genes.raw" )
  for p in global_slope baseline_thickness global_slope_1lmm baseline_thickness_1lmm; do
    raw="${PRAW[$p]}"; go="${raw%.raw}.out"; [[ -s "$raw" && -s "$go" ]] || { echo "FATAL: $raw missing" >&2; exit 3; }
    genez "$go" "$p" "$OUT/${p}_genez.txt"
    for D in "${DIS[@]}"; do
      [[ -s "$OUT/${D}_vs_${p}.gsa.out" ]] || "$MAGMA" --gene-results "$G/$D.genes.raw" --gene-covar "$OUT/${p}_genez.txt" --model direction=both --out "$OUT/${D}_vs_${p}" > /dev/null
      [[ -s "$OUT/${p}_on_${D}.gsa.out" ]] || "$MAGMA" --gene-results "$raw" --gene-covar "$G/${D}_genez.txt" --model direction=both --out "$OUT/${p}_on_${D}" > /dev/null
    done
    "$MAGMA" --gene-results "$raw" --set-annot "$SETS_SCZ" col=1,2 --model direction=greater --out "$OUT/${p}_sets_scz" > /dev/null
    "$MAGMA" --gene-results "$raw" --set-annot "$SETS_MDD" col=1,2 --model direction=greater --out "$OUT/${p}_sets_mdd" > /dev/null
    echo "  $PARC $p done  $(date +%T)"
  done
done
echo "=== collect  $(date)"
"$PY" "$REPO/genetic_analysis/step13_magma_panel_collect.py"
echo "=== step 13 done  $(date)"
