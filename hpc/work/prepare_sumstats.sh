#!/bin/bash
# Normalise the two disorder GWAS to the layout hpc/ assumes.
#
# Site-specific setup, not part of the pipeline: it exists because the copies of
# these files on CSD3 do not have the shape hpc/config.sh's defaults describe.
# Two concrete mismatches, both of which would have surfaced as a MAGMA/LDSC
# parse error rather than as a wrong number:
#
#   SCZ  PGC3_SCZ_wave3.primary...vcf.fixed.tsv is a PGC "sumstats VCF": 73
#        `##` metadata lines precede the real header.  Every consumer in hpc/
#        (MAGMA --pval, LDSC munge, 06_prs's awk) reads the header from line 1,
#        so all three would read `##fileFormat=PGCsumstatsVCFv1.0` as the header.
#        It also has NO N column -- it carries NCAS, NCON and NEFFDIV2 -- while
#        04_magma passes `ncol=N` and config's SCZ_N_SPEC says `--N-col NEFF`.
#
#   MDD  pgc-mdd2025... has the right columns but names the sample size `n`,
#        lowercase, against `ncol=N` / `--N-col N`.
#
# Rather than special-case each consumer, write one normalised copy of each with
# a real header on line 1 and an `N` column, and point $SCZ_SUMSTATS /
# $MDD_SUMSTATS at those.  N is the EFFECTIVE sample size in both cases
# (SCZ: 2 x NEFFDIV2), which is what a case/control GWAS should contribute to a
# gene-level or LD-score analysis; using NCAS+NCON would overstate it.
#
# Idempotent: skips a file that already exists.  ~15 min, ~1.5 GB.
set -euo pipefail

SRC="${SRC:-/rds/user/rajd2/hpc-work/magma/gwas}"
DST="${DST:-$(cd "$(dirname "$0")" && pwd)/sumstats}"
mkdir -p "$DST"

scz_in="$SRC/PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.fixed.tsv"
mdd_in="$SRC/pgc-mdd2025_no23andMe_div_v3-49-46-01_formatted.tsv"

if [[ ! -f "$DST/SCZ.tsv" ]]; then
  echo "SCZ: stripping ## metadata and adding N = 2 x NEFFDIV2"
  awk -F'\t' -v OFS='\t' '
    /^##/ {next}
    !hdr { for (i=1;i<=NF;i++) c[$i]=i
           if (!c["NEFFDIV2"]) { print "no NEFFDIV2 column" > "/dev/stderr"; exit 1 }
           ne=c["NEFFDIV2"]; print $0, "N"; hdr=1; next }
    { print $0, 2*$ne }' "$scz_in" > "$DST/SCZ.tsv"
fi
head -1 "$DST/SCZ.tsv"
echo "SCZ rows: $(( $(wc -l < "$DST/SCZ.tsv") - 1 ))"

if [[ ! -f "$DST/MDD.tsv" ]]; then
  echo "MDD: renaming n -> N"
  awk -F'\t' -v OFS='\t' '
    NR==1 { for (i=1;i<=NF;i++) if ($i=="n") $i="N"; print; next }
    { print }' "$mdd_in" > "$DST/MDD.tsv"
fi
head -1 "$DST/MDD.tsv"
echo "MDD rows: $(( $(wc -l < "$DST/MDD.tsv") - 1 ))"
