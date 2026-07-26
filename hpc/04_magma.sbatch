#!/bin/bash
#SBATCH -A VERTES-SL2-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --job-name=abcd-magma
#SBATCH -o slurm/abcd-magma-%j.out
#SBATCH -e slurm/abcd-magma-%j.err
source "$(dirname "$0")/config.sh"

ANNOT="$MAGMA_REF/snp_gene_annotation/ncbi37.window35-10.genes.annot"
REF="$MAGMA_REF/g1000_eur/g1000_eur"

# Step 1: aggregate SNP p-values to genes, for each ABCD phenotype GWAS.
for PHENO in "${PHENOTYPES[@]}"; do
  N=$(awk 'NR==2{print $6}' "$OUT/gwas/${PHENO}.fastGWA")   # N column
  awk 'NR==1{print "SNP\tP"} NR>1{print $2"\t"$10}' \
      "$OUT/gwas/${PHENO}.fastGWA" > "$OUT/magma/${PHENO}.pval"
  "$MAGMA" --bfile "$REF" --gene-annot "$ANNOT" \
           --pval "$OUT/magma/${PHENO}.pval" N="$N" \
           --out "$OUT/magma/${PHENO}"
done

# Step 2: competitive gene-set tests. The gene-set files come from
# `python -m abcd.magma_export` and contain, as covariates rather than
# hard sets: AHBA C1-C3 loadings (signed, split into +/- halves) and the
# snRNA-seq PC1 maturation loadings.
#
# --model direction=both because a maturation axis can be enriched at either
# tail and forcing one-sided would hide half the hypothesis.
for PHENO in "${PHENOTYPES[@]}"; do
  "$MAGMA" --gene-results "$OUT/magma/${PHENO}.genes.raw" \
           --gene-covar "$MAGMA_REF/genesets/ahba_components_posneg.txt" \
           --model direction=both \
           --out "$OUT/magma/${PHENO}_ahba"
  "$MAGMA" --gene-results "$OUT/magma/${PHENO}.genes.raw" \
           --gene-covar "$MAGMA_REF/genesets/snrnaseq_pc1.txt" \
           --model direction=both \
           --out "$OUT/magma/${PHENO}_pc1"
done

# Step 3: the reverse direction, which is the actual hypothesis -- are SCZ/MDD
# risk genes enriched in the ABCD developmental phenotype's association signal?
for D in SCZ MDD; do
  case $D in
    SCZ) F="$GWAS_DIR/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.fixed.tsv"; S=ID;   P=PVAL ;;
    MDD) F="$GWAS_DIR/pgc-mdd2025_no23andMe_div_v3-49-46-01_formatted.tsv";     S=rsid; P=p_value ;;
  esac
  "$MAGMA" --bfile "$REF" --gene-annot "$ANNOT" \
           --pval "$F" use="$S,$P" ncol=N \
           --out "$OUT/magma/${D}"
  for PHENO in "${PHENOTYPES[@]}"; do
    # ABCD gene Z-scores as a continuous covariate on disorder gene Z-scores.
    awk 'NR>1{print $1"\t"$8}' "$OUT/magma/${PHENO}.genes.out" \
        > "$OUT/magma/${PHENO}_genez.txt"
    "$MAGMA" --gene-results "$OUT/magma/${D}.genes.raw" \
             --gene-covar "$OUT/magma/${PHENO}_genez.txt" \
             --model direction=both \
             --out "$OUT/magma/${D}_vs_${PHENO}"
  done
done

tail -n +4 "$OUT"/magma/*_vs_*.gsa.out
