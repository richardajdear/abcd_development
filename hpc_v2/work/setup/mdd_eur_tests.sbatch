#!/bin/bash
#SBATCH --job-name=v2_mddeur
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#
# Re-run the MDD tests with the EUROPEAN-ONLY MDD GWAS.
#
# WHY: everything MDD in this project used pgc-mdd2025_no23andMe_DIV -- the
# trans-ancestry release (105 cohorts, 547,355 cases, 75.3% European by case
# count).  Paired with our EUR-stratified ABCD arm that leaves an ancestry
# mismatch on the disorder side.  The European release of the SAME meta-analysis
# (pgc-mdd2025_no23andMe_EUR: 75 cohorts, 412,305 cases / 1,588,397 controls) is
# the exact counterpart -- same consortium, same version family, same
# no-23andMe cohort exclusion, same GWAS-SSF column layout -- so ancestry is the
# ONLY thing that differs between the two analyses.  Downloaded from
# figshare 10.6084/m9.figshare.27061255.v5, md5 9dc86e3e4ff6f81fb46b4b5577335d9e
# verified against the release manifest.
#
# WHAT ACTUALLY CHANGES.  Only tests that consume MDD's own summary statistics:
#   - LDSC rg (MDD x each phenotype)
#   - MAGMA MDD_vs_<pheno>: our phenotype gene Z as a covariate on MDD gene Z
# The MDD_hc_finemap / MDD_highconf / MDD_pool -> phenotype tests do NOT change:
# those are fixed GENE LISTS scored against OUR gene results, with no MDD
# sumstats involved.  Reported unchanged rather than re-run.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
source "$REPO/hpc/config.local.sh"      # MAGMA_REF, MAGMA_ANNOT, HPCWORK
MAGMA_REF="${MAGMA_REF:-$HPCWORK/magma/reference_data/g1000_eur}"
MAGMA_ANNOT="${MAGMA_ANNOT:-$HPCWORK/magma/snp_gene_annotation/ncbi37.window35-10.genes.annot}"
W="$REPO/hpc/work"
SRC=/rds/user/rajd2/hpc-work/magma/gwas
OUTD="$OUT_V2/mdd_eur"
mkdir -p "$OUTD"

MDD_EUR="$OUTD/MDD_EUR.tsv"
if [[ ! -s "$MDD_EUR" ]]; then
  echo "=== prepare MDD-EUR (strip ## pre-header, rename n -> N)  $(date) ==="
  zcat "$SRC/pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz" \
    | grep -v '^##' \
    | awk -F'\t' -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="n") $i="N"; print; next} {print}' \
    > "$MDD_EUR"
  echo "rows: $(( $(wc -l < "$MDD_EUR") - 1 ))"
fi

# --- LDSC -------------------------------------------------------------------
LD="$W/ldsc_ref/eur_w_ld_chr/"
if [[ ! -f "$OUTD/MDD_EUR.sumstats.gz" ]]; then
  echo "=== munge MDD-EUR  $(date) ==="
  "$W/bin/munge_sumstats.py" --sumstats "$MDD_EUR" --snp rsid --p p_value \
    --a1 effect_allele --a2 other_allele --signed-sumstats beta,0 --N-col N \
    --merge-alleles "$W/ldsc_ref/eur_w_ld_chr/w_hm3.snplist" \
    --chunksize 500000 --out "$OUTD/MDD_EUR"
fi
for arm in v2 v1; do
  case $arm in
    v2) SS="$OUT_V2/ldsc_eur" ;;
    v1) SS="$W/results/ldsc_imp_eur" ;;
  esac
  list=""
  for p in $(phenotype_names); do
    [[ -f "$SS/$p.sumstats.gz" ]] && list="$list,$SS/$p.sumstats.gz"
  done
  echo "=== rg MDD-EUR x $arm EUR  $(date) ==="
  "$W/bin/ldsc.py" --rg "$OUTD/MDD_EUR.sumstats.gz${list}" \
    --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUTD/rg_MDDeur_${arm}" || true
done

# --- MAGMA ------------------------------------------------------------------
MAGMA="$W/bin/magma"
if [[ ! -f "$OUTD/MDD_EUR.genes.raw" ]]; then
  echo "=== MAGMA gene analysis, MDD-EUR  $(date) ==="
  # Paths come from v1's config (MAGMA_REF / MAGMA_ANNOT), not from a relative
  # walk out of hpc/work -- the first version of this script guessed the latter
  # and MAGMA died on a missing .bed.
  "$MAGMA" --bfile "$MAGMA_REF" --gene-annot "$MAGMA_ANNOT" \
    --pval "$MDD_EUR" use=rsid,p_value ncol=N --out "$OUTD/MDD_EUR"
fi
# Positive control: the MDD high-confidence sets must be enriched in MDD's own
# gene results.  If this is null, nothing below is interpretable.
"$MAGMA" --gene-results "$OUTD/MDD_EUR.genes.raw" \
  --set-annot "$W/genesets/mdd_highconf.txt" col=1,2 \
  --model direction=greater --out "$OUTD/MDD_EUR_mddhc" || true
# The direction that matters: our phenotype gene Z as covariate on MDD gene Z.
for arm in v2 v1; do
  case $arm in
    v2) GD="$OUT_V2/magma_eur" ;;
    v1) GD="$W/results/magma_imp_eur" ;;
  esac
  for p in $(phenotype_names); do
    gz="$GD/${p}_genez.txt"
    [[ -f "$gz" ]] || { echo "SKIP $arm $p: $gz absent"; continue; }
    "$MAGMA" --gene-results "$OUTD/MDD_EUR.genes.raw" --gene-covar "$gz" \
      --model direction=both --out "$OUTD/MDDeur_${arm}_vs_${p}" || true
  done
done
echo "=== done  $(date) ==="
