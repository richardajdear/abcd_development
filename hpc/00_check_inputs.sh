#!/bin/bash
# Preflight. Verifies every input exists and that IDs actually intersect
# before any expensive job is submitted. Run interactively on a login node.
source "$(dirname "$0")/config.sh"

fail=0
check_file () { [ -f "$1" ] && echo "  ok   $1" || { echo "  MISSING $1"; fail=1; }; }
check_exec () { [ -x "$1" ] && echo "  ok   $1" || { echo "  MISSING/NOT EXECUTABLE $1"; fail=1; }; }

echo "genotypes:"; for e in bed bim fam; do check_file "$GENO.$e"; done
echo "binaries:";  check_exec "$GCTA"; check_exec "$MAGMA"
echo "phenotypes:"; check_file "$PHENO_DIR/phenotypes_gcta.txt"
echo "sumstats:"
check_file "$GWAS_DIR/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.fixed.tsv"
check_file "$GWAS_DIR/pgc-mdd2025_no23andMe_div_v3-49-46-01_formatted.tsv"

# The failure mode that costs the most time: phenotype IDs that do not match
# the .fam file. ABCD writes NDAR_INV... in genetics tables and
# sub-NDARINV... in imaging tables; abcd.gcta_export normalises to the
# genetics form, but verify rather than trust.
if [ -f "$PHENO_DIR/phenotypes_gcta.txt" ] && [ -f "$GENO.fam" ]; then
  n=$(awk 'NR>1{print $2}' "$PHENO_DIR/phenotypes_gcta.txt" | sort -u \
      | comm -12 - <(awk '{print $2}' "$GENO.fam" | sort -u) | wc -l)
  echo "ID intersection with .fam: $n"
  [ "$n" -lt 1000 ] && { echo "  TOO FEW -- check ID format (NDAR_INV vs sub-NDARINV)"; fail=1; }
fi

exit $fail
