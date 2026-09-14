#!/bin/bash
# Collect the two-GRM .hsq files, parsing BY ROW NAME (the two-GRM layout is
# not the single-GRM one; a single-GRM .hsq must yield no row rather than be
# mislabelled as Zaitlen).  Same columns as step 05's own collector.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/config.sh
OUTDIR="${1:-$OUT_V2/reml_zaitlen_pcrel}"
summary="$OUTDIR/reml_zaitlen_summary.tsv"
{
  echo -e "phenotype\th2_snp\tse_snp\th2_ped\tse_ped\tpval\tn"
  while IFS= read -r pheno; do
    f="$OUTDIR/$pheno.hsq"
    [[ -f "$f" ]] || continue
    awk -v p="$pheno" -F'\t' '
      $1=="V(G1)/Vp"       {h1=$2; s1=$3}
      $1=="Sum of V(G)/Vp" {hs=$2; ss=$3}
      $1=="Pval"           {pv=$2}
      $1=="n"              {n=$2}
      END {if (h1!="" && hs!="") printf "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\n", p, h1, s1, hs, ss, pv, n}
    ' "$f"
  done < <(phenotype_names)
} > "$summary"
column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
