#!/bin/bash
#SBATCH --job-name=scz25_magma_finish
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/%x_%A.log
#
# Step 10, finish: the gene meta-analysis, the gene sets, and the tests
# against the phenotype gene results on both parcellations.
#
#   A. --meta raw=SCZ25_EUR,SCZ25_AFR,SCZ25_EAS -> SCZ25_META.genes.{out,raw}
#      (the principled multi-ancestry gene result)
#   B. gene sets (setup/build_scz2025_genesets.py)
#   C. per atlas (DK, HCP), phenotype gene results = step 7's EUR-arm MAGMA
#      (results_70tab*/magma_eur/<pheno>.genes.{raw,out}); for each disorder
#      gene result D in {SCZ25_EUR, SCZ25_META, SCZ25_METAnaive, PGC3_EUR,
#      PGC3_primary (= step 7's SCZ, naive)}:
#        reverse gene-property  D.genes.raw    ~ phenotype gene Z   (D_vs_<pheno>,   as legacy 04_magma step 3)
#        forward gene-property  pheno.genes.raw ~ D gene Z          (<pheno>_on_D)
#      and for each phenotype the competitive gene-set test on all sets
#      (direction=greater, as prioritised_gsa.sbatch), plus the 2025 sets
#      conditioned on their PGC3 analogue (does the new GWAS add beyond the
#      old?), and the circular positive control on SCZ25_META itself.
#   D. collector -> results_70tab*/magma_scz2025/table_magma_scz2025.tsv and
#      work/magma_scz2025/table_gene_level_comparison.tsv
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"
PY="$REPO/legacy/hpc/work/envs/abcd/bin/python"
ROOT="${SCZ25_MAGMA_ROOT:-$REPO/genetic_analysis/work/magma_scz2025}"
G="$ROOT/genes"; export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR"
ANNOT="$ROOT/ref/ncbi37.window35-10.union3.genes.annot"
GW25="$REPO/genetic_analysis/work/scores_scz2025/gwas"; GW3="$REPO/legacy/hpc_v2/work/results_v2/prs_final/gwas"
LEGACY_SETS="$REPO/legacy/hpc/work/genesets/scz_prioritised.txt"
PHENOS=(baseline_thickness global_slope slope_PC1 slope_PC2 slope_PC3)

echo "=== A. gene meta-analysis EUR+AFR+EAS  $(date)"
for n in SCZ25_EUR SCZ25_AFR SCZ25_EAS SCZ25_METAnaive PGC3_EUR; do [[ -s "$G/$n.genes.raw" ]] || { echo "FATAL: $G/$n.genes.raw missing" >&2; exit 2; }; done
# raw= writes the .genes.raw (Z, N and merged gene correlations) for the
# downstream gene-set / gene-property analyses; genes= writes the .genes.out
# table.  MAGMA does one or the other per call, so both are run.
"$MAGMA" --meta raw="$G/SCZ25_EUR,$G/SCZ25_AFR,$G/SCZ25_EAS" prefix --out "$G/SCZ25_META" | tail -2
"$MAGMA" --meta genes="$G/SCZ25_EUR,$G/SCZ25_AFR,$G/SCZ25_EAS" prefix --out "$G/SCZ25_META" | tail -2
[[ -s "$G/SCZ25_META.genes.raw" && -s "$G/SCZ25_META.genes.out" ]] || { echo "FATAL: meta produced no genes.raw/.genes.out" >&2; exit 4; }
echo "  meta genes: $(( $(wc -l < "$G/SCZ25_META.genes.out") - 1 ))"

echo "=== B. gene sets  $(date)"
SETS="$ROOT/genesets_scz2025.txt"
"$PY" "$REPO/genetic_analysis/setup/build_scz2025_genesets.py" "$ANNOT" "$GW25/SCZ25_META.ma" "$G/SCZ25_META.genes.out" \
      "$GW3/SCZ_pooled.ma" "$G/PGC3_EUR.genes.out" "$LEGACY_SETS" "$SETS"
cut -f2 "$SETS" | sort | uniq -c | awk '{printf "  %-24s %s\n", $2, $1}'

genez() {  # <genes.out> <label> <out.txt>: GENE + ZSTAT with a header MAGMA accepts
  awk -v g="$2" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; print "GENE\t" g; next} {print $c["GENE"] "\t" $c["ZSTAT"]}' "$1" > "$3"
}
for D in SCZ25_EUR SCZ25_META SCZ25_METAnaive PGC3_EUR; do genez "$G/$D.genes.out" "$D" "$G/${D}_genez.txt"; done

for PARC in dk hcp; do
  SUF=""; [[ "$PARC" == hcp ]] && SUF="_hcp"
  PH="$REPO/genetic_analysis/work/results_70tab$SUF/magma_eur"
  OUT="$REPO/genetic_analysis/work/results_70tab$SUF/magma_scz2025"; mkdir -p "$OUT"
  echo "=== C. $PARC  phenotype gene results: $PH  $(date)"
  # step 7's SCZ (PGC3 primary x g1000_eur) under a uniform name
  ln -sf "$PH/SCZ.genes.raw" "$G/PGC3_primary.genes.raw"; ln -sf "$PH/SCZ.genes.out" "$G/PGC3_primary.genes.out"
  genez "$G/PGC3_primary.genes.out" PGC3_primary "$G/PGC3_primary_genez.txt"
  for p in "${PHENOS[@]}"; do
    [[ -s "$PH/$p.genes.raw" ]] || { echo "FATAL: $PH/$p.genes.raw missing (step 7)" >&2; exit 3; }
    genez "$PH/$p.genes.out" "$p" "$OUT/${p}_genez.txt"
    for D in SCZ25_EUR SCZ25_META SCZ25_METAnaive PGC3_EUR PGC3_primary; do
      "$MAGMA" --gene-results "$G/$D.genes.raw" --gene-covar "$OUT/${p}_genez.txt" --model direction=both --out "$OUT/${D}_vs_${p}" > /dev/null
      "$MAGMA" --gene-results "$PH/$p.genes.raw" --gene-covar "$G/${D}_genez.txt" --model direction=both --out "$OUT/${p}_on_${D}" > /dev/null
    done
    "$MAGMA" --gene-results "$PH/$p.genes.raw" --set-annot "$SETS" col=1,2 --model direction=greater --out "$OUT/${p}_sets" > /dev/null
    "$MAGMA" --gene-results "$PH/$p.genes.raw" --set-annot "$SETS" col=1,2 --model direction=greater condition=PGC3_locus_pool --out "$OUT/${p}_sets_condPGC3pool" > /dev/null || echo "  WARN: conditional (PGC3_locus_pool) failed for $p"
    "$MAGMA" --gene-results "$PH/$p.genes.raw" --set-annot "$SETS" col=1,2 --model direction=greater condition=SCZ_locus_pool --out "$OUT/${p}_sets_condST12pool" > /dev/null || echo "  WARN: conditional (SCZ_locus_pool) failed for $p"
    echo "  $p done"
  done
  # circular positive control: the 2025 sets on the 2025 gene result itself (plumbing)
  "$MAGMA" --gene-results "$G/SCZ25_META.genes.raw" --set-annot "$SETS" col=1,2 --model direction=greater --out "$OUT/SCZ25_META_sets" > /dev/null
done

echo "=== D. collect  $(date)"
"$PY" "$REPO/genetic_analysis/step10_scz2025_magma_collect.py"
echo "=== step 10 done  $(date)"
