#!/bin/bash
#SBATCH --job-name=v2_csx
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-22
#
# PRS-CSx: multi-ancestry polygenic scores for SCHIZOPHRENIA.
#
# WHY ONLY SCZ.  PRS-CSx needs ancestry-STRATIFIED discovery GWAS, and only
# schizophrenia has them.  The PGC3 SCZ release (figshare 19426775) ships four
# ancestry-specific meta-analyses -- european, asian, afram, latino -- and its
# README is explicit that `core` is "cohorts primarily of east asian and
# european ancestry" and `primary` adds african-american and latino.  The PGC
# MDD 2025 release (figshare 27061255) ships only European-ancestry strata plus
# the pooled `div` file, so PRS-CSx is not possible for MDD.  ASD and ALZ have
# no per-ancestry release either.
#
# WHAT IT BUYS.  Our pooled target is 32% non-European while every score so far
# was built on a European LD reference.  PRS-CSx couples effects across
# populations, sharing information through a common global shrinkage parameter
# while modelling each population's own LD, and returns one posterior per
# population.  That fixes a concrete, identified defect in the current scores
# (reference mismatch for non-European subjects); it does NOT raise the
# phenotype's heritability, which is the actual limiting factor.
#
# One task per chromosome (the MCMC is serial).
set -e
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/hpc_v2/config.sh"
W="$REPO/hpc/work"
CSX=/rds/user/rajd2/hpc-work/prscs
PY="$W/envs/abcd/bin/python"
OUTD="$OUT_V2/prscsx"; mkdir -p "$OUTD"
PRS_GENO="$W/genotype_imputed/prs/abcd_imp_prs"
CHR="${SLURM_ARRAY_TASK_ID:-22}"

# --- inputs, one per ancestry, in PRS-CS format: SNP A1 A2 BETA P ------------
# N per ancestry = 2 x NEFF (the PGC files carry Neff/2), taken from the file.
# AFR and AMR are EXCLUDED, and the reason is a property of the release, not a
# choice: PGC3_SCZ_wave3.afram and .latino ship only 11 columns
# (CHROM ID POS A1 A2 FCAS FCON IMPINFO BETA SE PVAL) with no NCAS/NCON/NEFF,
# and no matching daner file exists.  PRS-CSx requires a per-population sample
# size, and inventing one from memory is not acceptable here.  The cost is
# real and worth stating: our target's non-European subjects are mostly of
# African and Hispanic ancestry, so dropping those two strata removes exactly
# the populations PRS-CSx would most have helped.  EUR + EAS is what the
# release actually supports.
POPS="EUR EAS"
G=/rds/user/rajd2/hpc-work/magma/gwas

sst=""; nn=""; pp=""
for pop in $POPS; do
  case $pop in
    EUR) f="$G/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz" ;;
    EAS) f="$G/PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.gz" ;;
    AFR) f="$G/PGC3_SCZ_wave3.afram.autosome.public.v3.vcf.tsv.gz" ;;
    AMR) f="$G/PGC3_SCZ_wave3.latino.autosome.public.v3.vcf.tsv.gz" ;;
  esac
  out="$OUTD/SCZ_${pop}.txt"
  if [[ ! -s "$out" ]]; then
    zcat "$f" | grep -v '^##' | awk 'BEGIN{FS=OFS="\t"}
      NR==1{for(i=1;i<=NF;i++)c[$i]=i; print "SNP","A1","A2","BETA","P"; next}
      { snp=$c["ID"]; a1=toupper($c["A1"]); a2=toupper($c["A2"])
        b=$c["BETA"]; p=$c["PVAL"]
        if(snp==""||b==""||p=="") next
        if(a1!~/^[ACGT]$/||a2!~/^[ACGT]$/) next
        if(b+0!=b+0||b=="inf"||b=="-inf") next
        if(p+0<=0||p+0>1) next
        print snp,a1,a2,b,p }' > "$out"
  fi
  # config.sh sets `pipefail`, and an awk that exits at row 2 closes the pipe
  # under zcat -- which dies on SIGPIPE and takes the whole job with it (every
  # one of the 22 tasks failed with exit 141 and an empty log).  Read the second
  # data row without terminating the pipeline early.
  n=$(zcat "$f" | grep -v '^##' | awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} NR==2{v=2*$c["NEFF"]} END{printf "%d", v}')
  echo "  $pop: $(( $(wc -l < "$out") - 1 )) SNPs, n_gwas=$n"
  sst="${sst},$out"; nn="${nn},$n"; pp="${pp},$pop"
done

export MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "=== PRS-CSx chr$CHR  $(date) ==="
"$PY" "$CSX/PRScsx/PRScsx.py" \
  --ref_dir="$CSX" \
  --bim_prefix="$PRS_GENO" \
  --sst_file="${sst#,}" \
  --n_gwas="${nn#,}" \
  --pop="${pp#,}" \
  --chrom="$CHR" \
  --out_dir="$OUTD" \
  --out_name=SCZ
echo "=== chr$CHR done  $(date) ==="
