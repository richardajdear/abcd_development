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
require_paths MAGMA_REF MAGMA_GENE_LOC GWAS_DIR MANIFEST MAGMA
ensure_dirs "$MAGMA_DIR"

names=()   # not mapfile: bash 3.2 compatibility, see 02_reml.sbatch
while IFS= read -r line; do names+=("$line"); done < <(phenotype_names)

# Step 0: SNP->gene annotation.  Built here if absent rather than assumed to
# exist, because it depends on the reference panel's .bim and so must be rebuilt
# whenever that panel changes.  window=35,10 kb is the MAGMA default convention
# (35 kb upstream, 10 kb downstream) and matches what the AHBA/snRNA-seq
# gene-set files were prepared against.
if [[ ! -f "$MAGMA_ANNOT" ]]; then
  echo "=== 04_magma  building SNP-gene annotation  $(date) ==="
  run "$MAGMA" --annotate window=35,10 \
      --snp-loc "$MAGMA_REF.bim" \
      --gene-loc "$MAGMA_GENE_LOC" \
      --out "${MAGMA_ANNOT%.genes.annot}"
fi

# Step 1: aggregate SNP p-values to genes, for each ABCD phenotype GWAS.
#
# Columns are resolved by NAME from the fastGWA header, as in 03_gwas: the
# previous version read SNP from field 2, P from field 10 and N from field 6 of
# row 2, which silently produces a garbage gene analysis if GCTA's column order
# changes.  MAGMA is given the .fastGWA file directly with use=/ncol= rather
# than an awk-extracted two-column file, so there is one less intermediate to
# get wrong.
for pheno in "${names[@]}"; do
  f="$GWAS_DIR/$pheno.fastGWA"
  if [[ ! -f "$f" ]]; then
    echo "SKIP $pheno: $f not found (run 03_gwas first)" >&2
    continue
  fi
  # Field-wise check, not a substring match on the header line: fastGWA is
  # tab-separated, so a `case " $hdr " in *" SNP "*)` test never matches and
  # would reject every valid file.
  awk 'NR==1{for(i=1;i<=NF;i++) c[$i]=1; exit !(c["SNP"] && c["P"] && c["N"])}' "$f" \
    || { echo "ERROR: $f header lacks SNP, P or N: $(head -1 "$f")" >&2; exit 3; }
  echo "=== 04_magma  gene analysis  $pheno  $(date) ==="
  run "$MAGMA" --bfile "$MAGMA_REF" --gene-annot "$MAGMA_ANNOT" \
      --pval "$f" use=SNP,P ncol=N \
      --out "$MAGMA_DIR/$pheno"
done

# Step 2: competitive gene-set tests. The gene-set files come from
# `python -m abcd.magma_export` and contain, as covariates rather than
# hard sets: AHBA C1-C3 loadings (signed, split into +/- halves) and the
# snRNA-seq PC1 maturation loadings.
#
# --model direction=both because a maturation axis can be enriched at either
# tail and forcing one-sided would hide half the hypothesis.
#
# Gene-covariate files come from $GENESET_DIR (written by abcd.magma_export),
# not from under $MAGMA_REF -- they are analysis outputs, not reference data,
# and the previous version looked for them inside the LD panel directory.
for pheno in "${names[@]}"; do
  raw="$MAGMA_DIR/$pheno.genes.raw"
  [[ -f "$raw" || "$DRY_RUN" == "1" ]] || { echo "SKIP $pheno: no $raw" >&2; continue; }
  for gs in ahba_components_posneg snrnaseq_pc1; do
    gsf="$GENESET_DIR/$gs.txt"
    if [[ ! -f "$gsf" ]]; then
      echo "SKIP $pheno x $gs: $gsf not found (run python -m abcd.magma_export)" >&2
      continue
    fi
    run "$MAGMA" --gene-results "$raw" \
        --gene-covar "$gsf" \
        --model direction=both \
        --out "$MAGMA_DIR/${pheno}_${gs}"
  done
done

# Step 3: the reverse direction, which is the actual hypothesis -- are SCZ/MDD
# risk genes enriched in the ABCD developmental phenotype's association signal?
for D in SCZ MDD; do
  case $D in
    SCZ) F="$SCZ_SUMSTATS"; S="$SCZ_SNP_COL"; P="$SCZ_P_COL" ;;
    MDD) F="$MDD_SUMSTATS"; S="$MDD_SNP_COL"; P="$MDD_P_COL" ;;
  esac
  if [[ ! -f "$F" && "$DRY_RUN" != "1" ]]; then
    echo "SKIP $D: $F not found (set ${D}_SUMSTATS)" >&2
    continue
  fi
  echo "=== 04_magma  $D gene analysis  $(date) ==="
  run "$MAGMA" --bfile "$MAGMA_REF" --gene-annot "$MAGMA_ANNOT" \
      --pval "$F" use="$S,$P" ncol=N \
      --out "$MAGMA_DIR/$D"

  for pheno in "${names[@]}"; do
    go="$MAGMA_DIR/$pheno.genes.out"
    [[ -f "$go" || "$DRY_RUN" == "1" ]] || { echo "SKIP $D x $pheno: no $go" >&2; continue; }
    # ABCD gene Z-scores as a continuous covariate on disorder gene Z-scores.
    # ZSTAT is located by name; it was field 8 positionally, which is a silent
    # hazard if MAGMA's .genes.out layout changes.  MAGMA wants a header.
    if [[ -f "$go" ]]; then
      awk 'NR==1{for(i=1;i<=NF;i++) c[$i]=i;
                 if(!c["GENE"]||!c["ZSTAT"]){print "ERR" > "/dev/stderr"; exit 1}
                 print "GENE\t" g; next}
           {print $c["GENE"] "\t" $c["ZSTAT"]}' g="$pheno" \
        OFS='\t' "$go" > "$MAGMA_DIR/${pheno}_genez.txt" \
        || { echo "ERROR: GENE/ZSTAT not found in $go" >&2; exit 3; }
    fi
    run "$MAGMA" --gene-results "$MAGMA_DIR/$D.genes.raw" \
        --gene-covar "$MAGMA_DIR/${pheno}_genez.txt" \
        --model direction=both \
        --out "$MAGMA_DIR/${D}_vs_${pheno}"
  done
done

# ---------------------------------------------------------------------------
# Collect every gene-set / gene-covariate test into one table.
#
# MAGMA .gsa.out has 3 comment lines before the header, hence `tail -n +4` in
# the previous version -- but that printed raw blocks per file with no
# indication of which file each came from.  This tags every row with its
# comparison so the enrichment results are readable as a single table.
# ---------------------------------------------------------------------------
if [[ "$DRY_RUN" != "1" ]]; then
  summary="$MAGMA_DIR/magma_summary.tsv"
  {
    echo -e "comparison\tvariable\ttype\tn_genes\tbeta\tse\tp"
    for f in "$MAGMA_DIR"/*.gsa.out; do
      [[ -f "$f" ]] || continue
      cmp=$(basename "$f" .gsa.out)
      awk -v c="$cmp" '
        /^#/ {next}
        $1=="VARIABLE" {for(i=1;i<=NF;i++) col[$i]=i; have=1; next}
        have && NF>3 {
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", c, $col["VARIABLE"],
            (col["TYPE"]?$col["TYPE"]:"NA"), (col["NGENES"]?$col["NGENES"]:"NA"),
            (col["BETA"]?$col["BETA"]:"NA"), (col["SE"]?$col["SE"]:"NA"),
            (col["P"]?$col["P"]:"NA")
        }' "$f"
    done
  } > "$summary"
  echo
  column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
  echo
  echo "Rows named <disorder>_vs_<phenotype> are the reverse-direction test:"
  echo "positive beta = the ABCD phenotype's gene-level signal is enriched in"
  echo "that disorder's risk genes."
fi

echo "=== 04_magma done  $(date) ==="
