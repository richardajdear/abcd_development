#!/bin/bash
#SBATCH --job-name=v2_ctrl
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --array=1-2
#
# CONTROL DISORDERS: an onset-age gradient, which is the point.
#
# The hypothesis is specific -- the genetics of ADOLESCENT cortical development
# should relate to risk for ADOLESCENT-ONSET illness.  A single association with
# SCZ cannot demonstrate that specificity; it is equally consistent with generic
# brain-related polygenicity or residual stratification.  What discriminates is
# a GRADIENT across onset age:
#
#   ASD    childhood onset   (SPARK + iPSYCH + PGC, N = 58,794)
#   SCZ    adolescent/young adult   -- the hypothesis, already run
#   MDD    adolescent/adult         -- already run
#   ALZ    late life         (PGC-ALZ2, Wightman 2021 excl. 23andMe)
#
# A result that holds for SCZ/MDD and fails for BOTH flanking controls is a far
# stronger claim than any single association.  A result that holds for all four
# is evidence the finding is not developmental at all.
#
# Array task 1 = ALZ, 2 = ASD.  ALZ needs chr:pos -> rsID mapping first
# (see prepare_alz.py; it is GRCh37 and our own .bim is GRCh38); ASD ships
# rsIDs, so it needs none.
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
SETUP="$REPO/hpc_v2/work/setup"       # resolve BEFORE sourcing v1's config,
source "$REPO/hpc_v2/config.sh"       # which sets ABCD_HPC_ROOT by plain
source "$REPO/hpc/config.sh"          # assignment.  hpc/config.sh brings
                                      # CLUMP_*/PRS_THRESHOLDS and MAGMA_REF.
W="$REPO/hpc/work"
SRC=/rds/user/rajd2/hpc-work/magma/gwas
# hpc/config.sh redefines OUT_V2's siblings; re-derive what we need from v2.
OUT_V2="$REPO/hpc_v2/work/results_v2"
PHENO="$W/pheno_allanc_prs/phenotypes_gcta.txt"
COVAR_QUANT="$W/pheno_allanc_prs/covar_quant.txt"
COVAR_CAT="$W/pheno_allanc_prs/covar_categorical.txt"
MANIFEST="$W/pheno_allanc_prs/phenotype_manifest.tsv"
EUR_IDS="$W/results/ancestry/eur_anchor.keep"
RSCRIPT="/rds/user/rajd2/hpc-work/envs/genesis/bin/Rscript"
PLINK="$W/bin/plink"
MAGMA="$W/bin/magma"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"

case "${SLURM_ARRAY_TASK_ID:-1}" in
  1) D=ALZ ;;
  2) D=ASD ;;
  *) echo "bad array id" >&2; exit 2 ;;
esac
OUTD="$OUT_V2/control_$D"; PRSD="$OUTD/prs"
mkdir -p "$OUTD" "$PRSD"
echo "=== control disorder $D  $(date) ==="

# --- 1. normalise to SNP / A1 / A2 / EFFECT / P / N ---------------------------
NORM="$OUTD/$D.tsv"
if [[ ! -s "$NORM" ]]; then
  case $D in
    ALZ)
      # No marker ID column; map chr:pos -> rsID against the GRCh37 MAGMA panel.
      python3 "$SETUP/prepare_alz.py" "$MAGMA_REF.bim" \
        "$SRC/PGCALZ2sumstatsExcluding23andMe.txt" "$OUTD/ALZ_raw.tsv"
      # Z is the only signed statistic available (no BETA, no frequency).
      awk 'BEGIN{FS=OFS="\t"} NR==1{print "SNP","A1","A2","EFFECT","P","N"; next}
           {print $1,$2,$3,$4,$5,$6}' "$OUTD/ALZ_raw.tsv" > "$NORM" ;;
    ASD)
      # SPARK+iPSYCH+PGC: rsIDs present, Effect is log-OR, N per SNP.
      awk 'BEGIN{FS=OFS="\t"} NR==1{
             for(i=1;i<=NF;i++) c[$i]=i
             print "SNP","A1","A2","EFFECT","P","N"; next}
           $c["rsID"]!="" && $c["rsID"]!="NA" {
             print $c["rsID"], toupper($c["Allele1"]), toupper($c["Allele2"]),
                   $c["Effect"], $c["P-value"], $c["TotalSampleSize"]}' \
        "$SRC/ASD_SPARK_iPSYCH_PGC.with_rsID.tsv" > "$NORM" ;;
  esac
fi
echo "$D: $(( $(wc -l < "$NORM") - 1 )) SNPs"

# --- 2. LDSC ------------------------------------------------------------------
LD="$W/ldsc_ref/eur_w_ld_chr/"
if [[ ! -f "$OUTD/$D.sumstats.gz" ]]; then
  "$W/bin/munge_sumstats.py" --sumstats "$NORM" --snp SNP --p P \
    --a1 A1 --a2 A2 --signed-sumstats EFFECT,0 --N-col N \
    --merge-alleles "$W/ldsc_ref/eur_w_ld_chr/w_hm3.snplist" \
    --chunksize 500000 --out "$OUTD/$D"
fi
for arm in v2 v1; do
  case $arm in
    v2) SS="$OUT_V2/ldsc_eur" ;;
    v1) SS="$W/results/ldsc_imp_eur" ;;
  esac
  list=""
  for p in baseline_thickness global_slope slope_PC3 slope_PC2 slope_PC1; do
    [[ -f "$SS/$p.sumstats.gz" ]] && list="$list,$SS/$p.sumstats.gz"
  done
  "$W/bin/ldsc.py" --rg "$OUTD/$D.sumstats.gz${list}" \
    --ref-ld-chr "$LD" --w-ld-chr "$LD" --out "$OUTD/rg_${D}_${arm}" || true
done

# --- 3. MAGMA -----------------------------------------------------------------
if [[ ! -f "$OUTD/$D.genes.raw" ]]; then
  "$MAGMA" --bfile "$MAGMA_REF" --gene-annot "$MAGMA_ANNOT" \
    --pval "$NORM" use=SNP,P ncol=N --out "$OUTD/$D"
fi
for arm in v2 v1; do
  case $arm in
    v2) GD="$OUT_V2/magma_eur" ;;
    v1) GD="$W/results/magma_imp_eur" ;;
  esac
  for p in baseline_thickness global_slope slope_PC3 slope_PC2 slope_PC1; do
    [[ -f "$GD/${p}_genez.txt" ]] || continue
    "$MAGMA" --gene-results "$OUTD/$D.genes.raw" --gene-covar "$GD/${p}_genez.txt" \
      --model direction=both --out "$OUTD/${D}_${arm}_vs_${p}" || true
  done
done

# --- 4. PRS, identical P+T settings to v1's SCZ/MDD scores --------------------
variants=("$D")
if [[ "$D" == "ALZ" ]]; then
  # APOE dominates AD genetics.  "Is APOE associated with cortical thinning" is
  # a different question from "is polygenic AD risk associated with it", so the
  # control is built both ways.  chr19:44.4-46.5 Mb, GRCh37.
  awk 'NR==FNR{if($1=="19" && $4>=44400000 && $4<=46500000) apoe[$2]=1; next}
       FNR==1{print; next} !($1 in apoe)' \
       "$MAGMA_REF.bim" "$NORM" > "$OUTD/ALZnoAPOE.tsv"
  variants+=("ALZnoAPOE")
fi
for V in "${variants[@]}"; do
  src="$OUTD/$V.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR==1{print "SNP","P"; next} {print $1,$5}' "$src" > "$PRSD/$V.assoc"
  "$PLINK" --bfile "$PRS_GENO" --clump "$PRSD/$V.assoc" \
    --clump-p1 "$CLUMP_P1" --clump-p2 "$CLUMP_P2" \
    --clump-r2 "$CLUMP_R2" --clump-kb "$CLUMP_KB" --out "$PRSD/$V" || true
  [[ -f "$PRSD/$V.clumped" ]] || { echo "no clumps for $V" >&2; continue; }
  awk 'NR>1 && $3!="" {print $3}' "$PRSD/$V.clumped" > "$PRSD/$V.index_snps"
  echo "  $V: $(wc -l < "$PRSD/$V.index_snps") index SNPs"
  for thr in $PRS_THRESHOLDS; do
    tag="${thr//./p}"
    awk -v t="$thr" 'NR==FNR{keep[$1]=1; next} FNR==1{next}
                     ($1 in keep) && ($5+0) < (t+0) {print $1"\t"$2"\t"$4}' \
      "$PRSD/$V.index_snps" "$src" > "$PRSD/${V}_$tag.weights"
    "$PLINK" --bfile "$PRS_GENO" --score "$PRSD/${V}_$tag.weights" 1 2 3 sum \
      --out "$PRSD/score_${V}_$tag" || true
  done
done

# --- 5. Associate, population and within-family -------------------------------
"$RSCRIPT" "$REPO/tools/prs_assoc.R" \
  --prs-dir "$PRSD" --pheno "$PHENO" --covar-quant "$COVAR_QUANT" \
  --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_association_$D.tsv" || true
"$RSCRIPT" "$REPO/hpc_v2/R/06_prs_family.R" \
  --prs-dir "$PRSD" --pheno "$PHENO" --covar-quant "$COVAR_QUANT" \
  --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" --eur-ids "$EUR_IDS" \
  --out "$OUTD/prs_withinfamily_$D.tsv" || true
echo "=== $D done  $(date) ==="
