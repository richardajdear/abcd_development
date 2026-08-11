#!/bin/bash
#SBATCH --job-name=abcd_gwas
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=%x_%A_%a.log
#SBATCH --array=1-5
#
# Step 3: GWAS by fastGWA-MLM, one phenotype per array task.
#
# Uses the SPARSE GRM and therefore the FULL sample -- fastGWA models the
# relatedness rather than discarding it, so unlike 02_reml the twin and sibling
# pairs are retained.  The two steps' N differ by design.
#
# Power caveat that governs how these results may be read: the phenotype is a
# shrunken per-subject slope whose reliability is ~0.21 (2 visits).  With
# N~8,000, genome-wide significance is not a realistic expectation; the
# defensible readouts are lambda_GC (is the model calibrated?) and the
# gene-level enrichment in 04_magma, which aggregates over SNPs.  Treat a null
# genome-wide scan as the expected outcome, not as a failed run.
# Locate config.sh.  Under sbatch, $0 is a COPY of this script in
# /var/spool/slurm/slurmd/jobNNN/, so `dirname "$0"` does not contain
# config.sh -- run_all.sh's SLURM mode therefore sourced nothing, and because
# `set -euo pipefail` lives INSIDE config.sh it was never enabled either, so
# every step ran to completion with undefined functions and exited 0.  Search
# the submit directory too ($PWD, which run_all.sh sets via --chdir).
set -e
_cfg=""
for _d in "$(dirname "$0")" "$PWD" "${SLURM_SUBMIT_DIR:-}"; do
  if [[ -n "$_d" && -f "$_d/config.sh" ]]; then _cfg="$_d/config.sh"; break; fi
done
[[ -n "$_cfg" ]] || { echo "FATAL: cannot locate hpc/config.sh (looked in $(dirname "$0"), $PWD, ${SLURM_SUBMIT_DIR:-unset})" >&2; exit 2; }
source "$_cfg"
require_paths GENO GRM_SPARSE PHENO COVAR_QUANT COVAR_CAT MANIFEST GCTA
ensure_dirs "$GWAS_DIR"

names=()   # not mapfile: bash 3.2 compatibility, see 02_reml.sbatch
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=$((SLURM_ARRAY_TASK_ID - 1))
  if (( idx >= ${#names[@]} )); then
    echo "array index $SLURM_ARRAY_TASK_ID exceeds ${#names[@]} phenotypes; nothing to do"
    exit 0
  fi
  todo=("${names[$idx]}")
else
  todo=("${names[@]}")
fi

for pheno in "${todo[@]}"; do
  mpheno=$(mpheno_for "$pheno") || { echo "ERROR: '$pheno' not in $MANIFEST" >&2; exit 3; }
  echo "=== 03_gwas  $pheno  (--mpheno $mpheno)  $(date) ==="

  run "$GCTA" \
    --bfile "$GENO" \
    --grm-sparse "$GRM_SPARSE" \
    --fastGWA-mlm \
    --pheno "$PHENO" \
    --mpheno "$mpheno" \
    --qcovar "$COVAR_QUANT" \
    --covar "$COVAR_CAT" \
    --maf "$MAF" \
    --thread-num "$THREADS" \
    --out "$GWAS_DIR/$pheno"
done

# ---------------------------------------------------------------------------
# Collect: lambda_GC, N, and hit counts per phenotype.
#
# lambda_GC is computed from the median chi-square, obtained by sorting the
# column and taking the middle value(s) with `sort -g`.  The previous version
# took a "median" with a single awk pass that indexed the middle of the *input*
# order rather than of sorted values -- which returns an arbitrary element, not
# the median, and additionally mishandled even-length input.  Sorting is not the
# clever way to do it but it is the correct one, and at ~1M SNPs it costs
# seconds.
#
# chi2 = (BETA/SE)^2 rather than an inverse-normal transform of P: it avoids
# depending on how many digits fastGWA printed for very small p-values.
# 0.4549364 is the median of a chi-square with 1 df.
# ---------------------------------------------------------------------------
if [[ "$DRY_RUN" != "1" ]]; then
  summary="$GWAS_DIR/gwas_summary.tsv"
  {
    echo -e "phenotype\tn_snps\tn_mean\tlambda_gc\tn_p5e8\tn_p1e5\tmin_p"
    for pheno in "${names[@]}"; do
      f="$GWAS_DIR/$pheno.fastGWA"
      [[ -f "$f" ]] || continue

      # Resolve columns by NAME from the header.  fastGWA's column order has
      # changed between GCTA versions; positional indexing is a silent hazard.
      read -r cN cBETA cSE cP < <(awk 'NR==1{
        for(i=1;i<=NF;i++) c[$i]=i
        print c["N"], c["BETA"], c["SE"], c["P"]; exit}' "$f")
      if [[ -z "$cP" || "$cP" == "0" ]]; then
        echo -e "$pheno\tERROR: could not locate P/BETA/SE columns in header" >&2
        continue
      fi

      n_snps=$(awk 'NR>1' "$f" | wc -l | tr -d ' ')
      median_chi2=$(awk -v b="$cBETA" -v s="$cSE" \
                      'NR>1 && $s>0 {print ($b/$s)^2}' "$f" \
                    | sort -g \
                    | awk '{v[NR]=$1} END{
                        if(NR==0){print "NA"; exit}
                        if(NR%2) print v[(NR+1)/2]
                        else print (v[NR/2]+v[NR/2+1])/2 }')
      lambda=$(awk -v m="$median_chi2" 'BEGIN{ if(m=="NA") print "NA"; else printf "%.4f", m/0.4549364 }')

      awk -v p="$pheno" -v cn="$cN" -v cp="$cP" -v ns="$n_snps" -v lam="$lambda" '
        NR>1 {
          nsum += $cn; nn++
          if ($cp < 5e-8) g++
          if ($cp < 1e-5) s++
          if (minp=="" || $cp < minp) minp = $cp
        }
        END {
          printf "%s\t%s\t%d\t%s\t%d\t%d\t%s\n", p, ns, (nn?nsum/nn:0), lam, g+0, s+0, (minp==""?"NA":minp)
        }' "$f"
    done
  } > "$summary"
  echo
  column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
  echo
  echo "lambda_gc near 1.00 means the mixed model is calibrated.  Above ~1.10"
  echo "suggests residual structure -- check that the ancestry PCs are in"
  echo "covar_quant.txt and that the sparse GRM was actually used."
fi

echo "=== 03_gwas done  $(date) ==="
