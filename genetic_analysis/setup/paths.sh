# Canonical paths for the final PRS analysis.  Everything resolves inside
# legacy/hpc_v2/work/inputs/, which is symlinks to the real files -- so no script
# below reaches into hpc/work/ or /rds/user/.../magma/ directly.
#
# DESIGN: the discovery GWAS is matched to the target arm.
#   pooled arm (n=8,082, multi-ancestry)  -> SCZ primary, MDD div
#   EUR arm    (n=4,116, European)        -> SCZ european, MDD eur
# ASD and ALZ have no ancestry-stratified release, so the same file serves both
# arms and that is flagged wherever they appear.
V2ROOT=/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis
IN="$V2ROOT/work/inputs"
RES="${OUT_V2:-$V2ROOT/work/results_70tab}"   # fresh dirs for the 7.0 re-run (README_HPC.md §4 rule 9)
FINAL="$RES/prs_final"

GWAS_DIR="$IN/gwas"
PHENO="$IN/pheno/phenotypes_gcta.txt"
COVQ="$IN/pheno/covar_quant.txt"
COVC="$IN/pheno/covar_categorical.txt"
MANIF="$IN/pheno/phenotype_manifest.tsv"
EURKEEP="$IN/ancestry/eur_anchor.keep"
STRATA="$IN/ancestry/strata_k4.tsv"
GENO="$IN/geno/abcd_imp_prs"

PLINK="$IN/bin/plink"
GCTB="$IN/bin/gctb"
RSCRIPT="$IN/bin/Rscript"
PYTHON="$IN/bin/python"
PRSCS="$IN/bin/PRScs/PRScs.py"
LDBLK_UKBB="$IN/ref/ldblk_ukbb_eur"
SBR_LDM="$IN/ref/sbayesr_ldm"
SBRC_EIGEN="$IN/ref/sbayesrc_eigen"
ANNOT="$IN/ref/annot_baseline2.2.txt"

ASSOC_R=/home/rajd2/rds/hpc-work/abcd_development/tools/prs_assoc.R
FAM_R=/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/R/06_prs_family.R

# (trait, arm) -> normalised sumstats built by normalise_gwas.py
CLUMP_P1=1; CLUMP_P2=1; CLUMP_R2=0.1; CLUMP_KB=250
PRS_THRESHOLDS="5e-8 1e-5 0.001 0.01 0.05 0.1 0.5 1"
