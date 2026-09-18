#!/bin/bash
#SBATCH --job-name=scz25_magma_prep
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/%x_%A.log
#
# Step 10 (2026-09-17), MAGMA on the Nature 2025 SCZ GWAS -- preparation.
#
# WHY MAGMA NEEDS SPECIAL HANDLING FOR A MULTI-ANCESTRY GWAS.  MAGMA's gene
# analysis turns SNP p-values into a gene statistic whose null distribution
# depends on the LD among the gene's SNPs, which it takes from a reference
# panel.  Step 7 used 1000G EUR for everything, including the PGC3 *primary*
# (multi-ancestry) file -- an LD mismatch inherited from v1.  For the 2025
# release there is no single panel that matches the AFR+EUR+EAS meta, and the
# MAGMA manual's answer is the one used here: run each ancestry's GWAS against
# its own 1000G panel, then combine the gene results with `--meta`
# (weighted Stouffer's Z, weights sqrt(N)), which requires the SAME SNP-gene
# annotation across cohorts.  This job therefore
#   1. unpacks the 1000G AFR and EAS panels (downloaded from the MAGMA site,
#      Build 37, dbSNP151 synonyms) next to the existing EUR one;
#   2. builds one annotation (NCBI37.3 genes, 35 kb up / 10 kb down -- the
#      convention every gene-set file in this project was prepared against)
#      over the UNION of the three panels' SNP positions, so AFR- and
#      EAS-specific SNPs are annotated too.  Gene coordinates are identical
#      to the EUR-only annotation step 7 used for the phenotype GWAS, so
#      gene-level results join and meta-analyse without complaint.
set -euo pipefail
MG=/rds/user/rajd2/hpc-work/magma
REF="$MG/reference_data"; GENELOC="$MG/gene_locations/NCBI37.3.gene.loc"
MAGMA=/home/rajd2/rds/hpc-work/abcd_development/legacy/hpc/work/bin/magma
ROOT="${SCZ25_MAGMA_ROOT:-/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/magma_scz2025}"
mkdir -p "$ROOT/ref"; export TMPDIR="$ROOT/tmp"; mkdir -p "$TMPDIR"

echo "=== 1. reference panels  $(date)"
for pop in afr eas; do
  if [[ -s "$REF/g1000_$pop.bed" && -s "$REF/g1000_$pop.bim" && -s "$REF/g1000_$pop.synonyms" ]]; then
    echo "  g1000_$pop present"; continue
  fi
  [[ -s "$REF/g1000_$pop.zip" ]] || { echo "FATAL: $REF/g1000_$pop.zip missing" >&2; exit 2; }
  ( cd "$REF" && unzip -o -q "g1000_$pop.zip" && mv -f README "README_g1000_$pop" 2>/dev/null || true )
  echo "  g1000_$pop unpacked: $(wc -l < "$REF/g1000_$pop.bim") SNPs, $(wc -l < "$REF/g1000_$pop.fam") individuals"
done
for pop in eur afr eas; do [[ -s "$REF/g1000_$pop.bed" ]] || { echo "FATAL: g1000_$pop.bed missing" >&2; exit 2; }; done

echo "=== 2. union SNP locations  $(date)"
UNION="$ROOT/ref/g1000_union3.snploc"
if [[ ! -s "$UNION" ]]; then
  # SNP CHR BP; first occurrence wins (EUR first, so positions agree with the
  # step-7 annotation wherever a SNP is shared)
  awk '!($2 in s){s[$2]=1; print $2, $1, $4}' "$REF/g1000_eur.bim" "$REF/g1000_afr.bim" "$REF/g1000_eas.bim" > "$UNION"
fi
echo "  $(wc -l < "$UNION") SNPs (EUR $(wc -l < "$REF/g1000_eur.bim"), AFR $(wc -l < "$REF/g1000_afr.bim"), EAS $(wc -l < "$REF/g1000_eas.bim"))"

echo "=== 3. annotation, window 35/10 kb  $(date)"
ANNOT="$ROOT/ref/ncbi37.window35-10.union3"
if [[ ! -s "$ANNOT.genes.annot" ]]; then
  "$MAGMA" --annotate window=35,10 --snp-loc "$UNION" --gene-loc "$GENELOC" --out "$ANNOT" | tail -5
fi
echo "  $(grep -vc '^#' "$ANNOT.genes.annot") genes annotated"
echo "=== prep done  $(date)"
