#!/bin/bash
#SBATCH --job-name=abcd_imp_probe
#SBATCH --account=VERTES-SL3-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log
#
# Size the imputed-genotype project before committing to it.
#
# Everything about that project -- disk, wall-clock, whether it fits in the
# remaining ~180 GB at all -- follows from ONE number nobody has measured: how
# many TOPMed r3 variants survive an imputation-quality and frequency filter.
# The raw count is irrelevant (most imputed variants are near-monomorphic
# guesses); the surviving count is what becomes a PLINK fileset and a GRM.
#
# chr22 is the smallest autosome (12.6 GB) and its variant count scales to the
# genome closely enough to plan from.  Measuring it costs one job; guessing it
# wrong costs a week and possibly a full filesystem.
#
# R2 is Minimac4's imputation r-squared: the correlation between the imputed
# dosage and the (unobserved) true genotype, estimated internally.  0.8 is the
# conventional floor; 0.3 is the permissive one.  Both are reported so the
# threshold is a visible choice.
set -euo pipefail
D=/rds/project/rds-CeXlNYOYMxw/Data_Genetics/genotype_microarray/imputed
V=$D/chr22_dose.vcf.gz

echo "=== chr22 variant census  $(date) ==="
# One decompression pass; parse INFO for R2 and MAF rather than calling bcftools,
# which is not installed on this account.
zcat "$V" | awk -F'\t' '
  /^#/ {next}
  {
    n++
    r2 = 0; maf = 0
    if (match($8, /[;\t]?R2=[0-9.e+-]+/))  { s=substr($8, RSTART, RLENGTH); sub(/^[;]?R2=/,"",s); r2=s+0 }
    if (match($8, /[;\t]?MAF=[0-9.e+-]+/)) { s=substr($8, RSTART, RLENGTH); sub(/^[;]?MAF=/,"",s); maf=s+0 }
    if (r2>=0.3) a++
    if (r2>=0.8) b++
    if (maf>=0.01) c++
    if (r2>=0.3 && maf>=0.01) d++
    if (r2>=0.8 && maf>=0.01) e++
    if ($8 ~ /TYPED/) t++
  }
  END {
    printf "  total variants          %12d\n", n
    printf "  R2>=0.3                 %12d  (%5.1f%%)\n", a, 100*a/n
    printf "  R2>=0.8                 %12d  (%5.1f%%)\n", b, 100*b/n
    printf "  MAF>=0.01               %12d  (%5.1f%%)\n", c, 100*c/n
    printf "  R2>=0.3 & MAF>=0.01     %12d  (%5.1f%%)   <- permissive\n", d, 100*d/n
    printf "  R2>=0.8 & MAF>=0.01     %12d  (%5.1f%%)   <- conventional\n", e, 100*e/n
    printf "  genotyped (TYPED flag)  %12d\n", t
  }'
echo "=== done  $(date) ==="
