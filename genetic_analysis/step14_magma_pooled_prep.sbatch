#!/bin/bash
#SBATCH --job-name=magma_pooled_prep
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=slurm/%x_%A.log
#
# Step 14 prep: MAGMA LD reference = the ABCD analysis sample itself.
# The manual's requirement is that the reference's ancestry match the sample
# the p-values came from; for the pooled (admixed) ABCD scan no continental
# panel does, but the GWAS sample's own genotypes match by construction.
# Restrict the imputed PRS fileset (7.07 M SNPs, 11,670 children) to the 8,596
# phenotyped-and-genotyped analysis set, and build the 35/10 kb NCBI37.3
# annotation on its SNP positions (hg19; same gene coordinates as the union
# annotation used for the disorder side, so gene results join).
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
PLINK="$REPO/legacy/hpc/work/bin/plink"; MAGMA="$REPO/legacy/hpc/work/bin/magma"
GENO="$REPO/genetic_analysis/work/inputs/geno/abcd_imp_prs"
GENELOC=/rds/user/rajd2/hpc-work/magma/gene_locations/NCBI37.3.gene.loc
ROOT="$REPO/genetic_analysis/work/magma_scz2025"; OUT="$ROOT/ref_abcd"; mkdir -p "$OUT"
export TMPDIR="$ROOT/tmp"
awk 'NR>1{print $1, $2}' "$REPO/genetic_analysis/work/pheno_70tab_fidiid/phenotypes_gcta.txt" > "$OUT/analysis.keep"
echo "keep: $(wc -l < "$OUT/analysis.keep") children"
if [[ ! -s "$OUT/abcd_analysis.bed" ]]; then
  "$PLINK" --bfile "$GENO" --keep "$OUT/analysis.keep" --make-bed --threads 8 --memory 40000 --out "$OUT/abcd_analysis" | grep -E "people|variants|pass" | head -5
fi
echo "reference: $(wc -l < "$OUT/abcd_analysis.fam") individuals, $(wc -l < "$OUT/abcd_analysis.bim") SNPs"
if [[ ! -s "$OUT/ncbi37.window35-10.abcd.genes.annot" ]]; then
  awk '{print $2, $1, $4}' "$OUT/abcd_analysis.bim" > "$OUT/abcd.snploc"
  "$MAGMA" --annotate window=35,10 --snp-loc "$OUT/abcd.snploc" --gene-loc "$GENELOC" --out "$OUT/ncbi37.window35-10.abcd" | tail -3
fi
echo "annotated genes: $(grep -vc '^#' "$OUT/ncbi37.window35-10.abcd.genes.annot")"
echo "=== prep done $(date)"
