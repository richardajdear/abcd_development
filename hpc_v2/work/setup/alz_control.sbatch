#!/bin/bash
#SBATCH --job-name=v2_alz
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#
# ALZHEIMER'S AS A NEGATIVE CONTROL.
#
# The project's hypothesis is specific: the genetics of ADOLESCENT cortical
# development should relate to risk for ADOLESCENT-ONSET psychiatric illness.
# That prediction is only informative if a LATE-onset brain disorder does NOT
# show the same association.  If Alzheimer's PRS predicts adolescent cortical
# thinning as well as SCZ PRS does, the SCZ result is much more likely to be a
# generic confound (stratification, or any brain-related polygenicity) than
# anything developmentally specific.
#
# GWAS: PGCALZ2sumstatsExcluding23andMe.txt -- Wightman et al. 2021, the PGC-ALZ
# release excluding 23andMe, already on disk (12,688,340 SNPs).  Chosen over
# Bellenguez 2022 (larger, but EADB/GR@ACE rather than PGC) to keep the control
# inside the same consortium as SCZ and MDD.
#
# CAVEATS, both stated rather than hidden:
#  1. The file has Z but no BETA and no allele frequency, so no per-SNP beta can
#     be reconstructed.  Z is used as the PRS weight.  Since prs_assoc.R
#     standardises the score before fitting, only the RELATIVE weighting across
#     SNPs matters, and Z differs from beta by the per-SNP SE only.
#  2. APOE dominates Alzheimer's genetics.  Both an all-SNP score and an
#     APOE-excluded score (chr19:44.4-46.5 Mb, GRCh37) are built, because a
#     control that is really "is APOE associated with cortical thinning" is a
#     different question from "is polygenic AD risk associated with it".
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
SETUP="$REPO/hpc_v2/work/setup"          # resolve BEFORE sourcing v1's config:
source "$REPO/hpc/config.local.sh"       # it sets ABCD_HPC_ROOT by plain
                                         # assignment and would clobber v2's.
MAGMA_REF="${MAGMA_REF:-$HPCWORK/magma/reference_data/g1000_eur}"
MAGMA_ANNOT="${MAGMA_ANNOT:-$HPCWORK/magma/snp_gene_annotation/ncbi37.window35-10.genes.annot}"
W="$REPO/hpc/work"
SRC=/rds/user/rajd2/hpc-work/magma/gwas
OUTD="$OUT_V2/alz"
PRSD="$OUTD/prs"
mkdir -p "$OUTD" "$PRSD"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"

# --- 1. rsID mapping ---------------------------------------------------------
if [[ ! -s "$OUTD/ALZ.tsv" ]]; then
  echo "=== map PGC-ALZ2 to rsIDs  $(date) ==="
  python3 "$SETUP/prepare_alz.py" \
    "$PRS_GENO.bim" "$SRC/PGCALZ2sumstatsExcluding23andMe.txt" "$OUTD/ALZ.tsv"
fi
echo "ALZ.tsv: $(( $(wc -l < "$OUTD/ALZ.tsv") - 1 )) SNPs"

# --- 2. LDSC ------------------------------------------------------------------
LD="$W/ldsc_ref/eur_w_ld_chr/"
if [[ ! -f "$OUTD/ALZ.sumstats.gz" ]]; then
  "$W/bin/munge_sumstats.py" --sumstats "$OUTD/ALZ.tsv" --snp SNP --p P \
    --a1 A1 --a2 A2 --signed-sumstats Z,0 --N-col N \
    --merge-alleles "$W/ldsc_ref/eur_w_ld_chr/w_hm3.snplist" \
    --chunksize 500000 --out "$OUTD/ALZ"
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
  echo "=== rg ALZ x $arm EUR  $(date) ==="
  "$W/bin/ldsc.py" --rg "$OUTD/ALZ.sumstats.gz${list}" \
    --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUTD/rg_ALZ_${arm}" || true
done

# --- 3. MAGMA -----------------------------------------------------------------
MAGMA="$W/bin/magma"
if [[ ! -f "$OUTD/ALZ.genes.raw" ]]; then
  echo "=== MAGMA gene analysis, ALZ  $(date) ==="
  "$MAGMA" --bfile "$MAGMA_REF" --gene-annot "$MAGMA_ANNOT" \
    --pval "$OUTD/ALZ.tsv" use=SNP,P ncol=N --out "$OUTD/ALZ"
fi
for arm in v2 v1; do
  case $arm in
    v2) GD="$OUT_V2/magma_eur" ;;
    v1) GD="$W/results/magma_imp_eur" ;;
  esac
  for p in $(phenotype_names); do
    [[ -f "$GD/${p}_genez.txt" ]] || continue
    "$MAGMA" --gene-results "$OUTD/ALZ.genes.raw" --gene-covar "$GD/${p}_genez.txt" \
      --model direction=both --out "$OUTD/ALZ_${arm}_vs_${p}" || true
  done
done

# --- 4. PRS, same P+T settings as v1's SCZ/MDD scores -------------------------
# APOE-excluded copy first, so both variants clump independently.
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {print}' "$OUTD/ALZ.tsv" > "$PRSD/ALZ_norm.tsv"
awk 'NR==FNR{if(FNR>1 && $1=="19" && $4>=44400000 && $4<=46500000) apoe[$2]=1; next}
     FNR==1{print; next} !($1 in apoe)' \
     "$PRS_GENO.bim" "$OUTD/ALZ.tsv" > "$PRSD/ALZnoAPOE_norm.tsv"
echo "ALZ SNPs: $(( $(wc -l < "$PRSD/ALZ_norm.tsv") - 1 ));  minus APOE: $(( $(wc -l < "$PRSD/ALZnoAPOE_norm.tsv") - 1 ))"

for D in ALZ ALZnoAPOE; do
  norm="$PRSD/${D}_norm.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR==1{print "SNP","P"; next} {print $1,$5}' "$norm" > "$PRSD/$D.assoc"
  echo "=== clump $D  $(date) ==="
  "$PLINK" --bfile "$PRS_GENO" --clump "$PRSD/$D.assoc" \
    --clump-p1 "$CLUMP_P1" --clump-p2 "$CLUMP_P2" \
    --clump-r2 "$CLUMP_R2" --clump-kb "$CLUMP_KB" \
    --out "$PRSD/$D" || true
  [[ -f "$PRSD/$D.clumped" ]] || { echo "no clumps for $D" >&2; continue; }
  awk 'NR>1 && $3!="" {print $3}' "$PRSD/$D.clumped" > "$PRSD/$D.index_snps"
  echo "  $D: $(wc -l < "$PRSD/$D.index_snps") index SNPs"
  for thr in $PRS_THRESHOLDS; do
    tag="${thr//./p}"
    awk -v t="$thr" 'NR==FNR{keep[$1]=1; next} FNR==1{next}
                     ($1 in keep) && ($5+0) < (t+0) {print $1"\t"$2"\t"$4}' \
      "$PRSD/$D.index_snps" "$norm" > "$PRSD/${D}_$tag.weights"
    "$PLINK" --bfile "$PRS_GENO" --score "$PRSD/${D}_$tag.weights" 1 2 3 sum \
      --out "$PRSD/score_${D}_$tag" || true
  done
done

# --- 5. Associate: population and within-family, same scripts as SCZ/MDD ------
"$RSCRIPT" "$REPO/tools/prs_assoc.R" \
  --prs-dir "$PRSD" --pheno "$PHENO" --covar-quant "$COVAR_QUANT" \
  --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_association_ALZ.tsv" || true
"$RSCRIPT" "$RSCRIPT_DIR/06_prs_family.R" \
  --prs-dir "$PRSD" --pheno "$PHENO" --covar-quant "$COVAR_QUANT" \
  --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_withinfamily_ALZ.tsv" || true
echo "=== done  $(date) ==="
