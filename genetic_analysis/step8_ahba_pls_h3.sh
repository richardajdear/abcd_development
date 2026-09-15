#!/bin/bash
# README_HPC.md 3 step 8 / ahba_pls/FOLLOWUP_GENETICS.md H3: MAGMA gene-property
# of each ABCD phenotype's gene-level Z on the NSPN-PLS2 / thinning-signature
# gene weights (three downsampling columns) and AHBA C3 -- marginal, and the
# DS columns conditioned on AHBA_C3 (the test that matters).  Seconds per run;
# the expensive gene analysis is step 7's .genes.raw, reused here.
# Stated prior (FOLLOWUP_GENETICS.md): expect null; report null as null.
#   PARC=hcp bash genetic_analysis/step8_ahba_pls_h3.sh   # HCP arm, once its step 7 exists
set -euo pipefail
cd "$(dirname "$0")/.."
source genetic_analysis/config.sh
COVAR="$REPO_ROOT/ahba_pls/hpc/lead_pls2_gene_covar_entrez.txt"
MAG="$OUT_V2/magma_eur"; OUT="$OUT_V2/magma_ahba_pls_h3"; mkdir -p "$OUT"
[[ -f "$COVAR" ]] || { echo "FATAL: $COVAR missing" >&2; exit 2; }
# MAGMA drops any covariate column that contains missing values ("variable
# ... was removed during preprocessing"): ds25/ds50/AHBA_C3 have 3,447 /
# 6,846 / 6,862 NA genes of 13,789, so a single --gene-covar run silently tests
# ds0 only (first attempt, 2026-09-15).  Split the file per variable with NA
# genes removed -- which is also literally what H3 asks for: "the three DS
# columns as separate marginal models", plus each conditioned on AHBA_C3.
CV="$OUT/covar"; mkdir -p "$CV"
for v in thinning_Z_ds0 thinning_Z_ds25 thinning_Z_ds50 AHBA_C3; do
  awk -F'\t' -v v="$v" 'NR==1{for(i=1;i<=NF;i++) if($i==v) c=i; print "GENE\t" v; next} $c!="NA"{print $1"\t"$c}' "$COVAR" > "$CV/$v.txt"
done
for v in thinning_Z_ds0 thinning_Z_ds25 thinning_Z_ds50; do
  awk -F'\t' -v v="$v" 'NR==1{for(i=1;i<=NF;i++){if($i==v) c=i; if($i=="AHBA_C3") k=i}; print "GENE\t" v "\tAHBA_C3"; next} $c!="NA" && $k!="NA"{print $1"\t"$c"\t"$k}' "$COVAR" > "$CV/${v}_with_C3.txt"
done
for pheno in $(phenotype_names); do
  raw="$MAG/$pheno.genes.raw"; [[ -f "$raw" ]] || { echo "SKIP $pheno (no $raw)"; continue; }
  for v in thinning_Z_ds0 thinning_Z_ds25 thinning_Z_ds50 AHBA_C3; do
    "$MAGMA" --gene-results "$raw" --gene-covar "$CV/$v.txt" --out "$OUT/${pheno}_marginal_$v" >/dev/null
  done
  for v in thinning_Z_ds0 thinning_Z_ds25 thinning_Z_ds50; do
    "$MAGMA" --gene-results "$raw" --gene-covar "$CV/${v}_with_C3.txt" --model condition=AHBA_C3 --out "$OUT/${pheno}_condC3_$v" >/dev/null
  done
  echo "  $pheno: 4 marginal + 3 conditional done"
done
# collect: MAGMA .gsa.out -> one table
{ printf 'phenotype\tmodel\tvariable\tn_genes\tbeta\tbeta_std\tse\tp\n'
  for f in "$OUT"/*_marginal_*.gsa.out "$OUT"/*_condC3_*.gsa.out; do
    b=$(basename "$f" .gsa.out)
    case "$b" in *_marginal_*) ph=${b%%_marginal_*}; model=marginal ;; *) ph=${b%%_condC3_*}; model=condC3 ;; esac
    # header-driven: the conditional output carries an extra MODEL column
    awk -v ph="$ph" -v m="$model" '/^#/ {next}
      $1=="VARIABLE" {for(i=1;i<=NF;i++) c[$i]=i; next}
      NF>=7 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", ph, m, $c["VARIABLE"], $c["NGENES"], $c["BETA"], $c["BETA_STD"], $c["SE"], $c["P"]}' "$f"
  done; } > "$OUT/table_h3.tsv"
echo "wrote $OUT/table_h3.tsv ($(( $(wc -l < "$OUT/table_h3.tsv") - 1 )) rows)"
