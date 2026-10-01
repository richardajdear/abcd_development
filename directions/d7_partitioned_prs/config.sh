# D7 partitioned SBayesRC scores: CSD3 paths.  Sourced by every d7 sbatch.
# Everything that is not D7-specific comes from the genetics pipeline's
# canonical path file, so the genotypes, strata and binaries are the ones every
# Figure-1 score was built with.
source /home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/setup/paths.sh

REPO=/home/rajd2/rds/hpc-work/abcd_development
D7="$REPO/directions/d7_partitioned_prs"
WORK="${D7_WORK:-$D7/work}"                 # gitignored: SNP-level and per-subject files
RESULTS="$D7/results"                       # committed: summary tables only
export TMPDIR="$WORK/tmp"; mkdir -p "$TMPDIR" "$WORK" "$RESULTS"

# python env with bed-reader (see env_d7.yml)
D7PY="${D7PY:-/home/rajd2/rds/hpc-work/envs/d7/bin/python}"

# MAGMA gene locations, the file steps 10 and 14 annotate with (NCBI37.3, GRCh37)
GENELOC="${GENELOC:-/rds/user/rajd2/hpc-work/magma/gene_locations/NCBI37.3.gene.loc}"

# EUR arm definition used by every step-9 association (step9_scz2025_assoc_1lmm.sbatch)
EUR_KEEP="$REPO/legacy/hpc/work/results/ancestry/eur_anchor.keep"

# Score arms: NAME|weights file|raw .profile next to it|cell
# cell = pooled (multi-ancestry GWAS -> full sample, within-cluster standardised; rule 4)
#      | EUR    (European GWAS -> EUR arm, raw score)
SCZ25="$REPO/genetic_analysis/work/scores_scz2025/SBayesRC"
# MDD/EA: the v2 scores Figure 1 used (step6_prs_assoc.sbatch SCORE_ROOT=$V2_LEGACY/prs_final);
# $RES/prs_final holds only the 7.0 association tables, not the weights.
PRSF="${D7_PRSF:-$REPO/legacy/hpc_v2/work/results_v2/prs_final/SBayesRC}"
ARMS=(
  "SCZ25_META|$SCZ25/SCZ25_META/SCZ25META_sbrc.weights|$SCZ25/SCZ25_META/score_SCZ25META_sbrc.profile|pooled"
  "SCZ25_EUR|$SCZ25/SCZ25_EUR/SCZ25EUR_sbrc.weights|$SCZ25/SCZ25_EUR/score_SCZ25EUR_sbrc.profile|EUR"
  "MDD_pooled|$PRSF/MDD_pooled/MDDpooled_sbrc.weights|$PRSF/MDD_pooled/score_MDDpooled_sbrc.profile|pooled"
  "MDD_eur|$PRSF/MDD_eur/MDDeur_sbrc.weights|$PRSF/MDD_eur/score_MDDeur_sbrc.profile|EUR"
  "EA|$PRSF/EA/EA_sbrc.weights|$PRSF/EA/score_EA_sbrc.profile|EUR"
)

# Phenotype exports (single LMM; FID = family id).  HCP-MMP is primary, DK secondary.
PHENO_HCP="$REPO/genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno"
PHENO_DK="$REPO/genetic_analysis/work/results_70tab/prs_final_1lmm/pheno"
SYMPTOMS="$WORK/d7_symptom_outcomes.tsv"     # laptop step L1, copied by scp
C3PHENO="$REPO/genetic_analysis/work/results_70tab_hcp/c3axis/pheno"   # optional (README_HPC 8.3 C3-A)
