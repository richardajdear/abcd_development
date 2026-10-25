# Shared paths for the MAGMA window-robustness runs (gwas_cells/csd3/MAGMA_WINDOWS_INSTRUCTIONS.md).
REPO=/home/rajd2/rds/hpc-work/abcd_development
MAGMA="$REPO/legacy/hpc/work/bin/magma"; PY="$REPO/legacy/hpc/work/envs/abcd/bin/python"
MG=/rds/user/rajd2/hpc-work/magma; REF="$MG/reference_data"
GENELOC="$MG/gene_locations/NCBI37.3.gene.loc"
UNION="$REPO/genetic_analysis/work/magma_scz2025/ref/g1000_union3.snploc"
W="$REPO/genetic_analysis/work/magma_windows"
TAGS=(w0_0 w10_10 w35_10 w50_50 w100_100)
TRAITS=(ABCD_thin_EUR ABCD_base_EUR SCZ25_EUR MDD_EUR)
A70="$REPO/genetic_analysis/work/results_70tab_hcp/scan_1lmm/assoc_eur"
declare -A SRC=( [ABCD_thin_EUR]="$A70/global_slope_1lmm.sumstats.tsv.gz"
                 [ABCD_base_EUR]="$A70/baseline_thickness_1lmm.sumstats.tsv.gz"
                 [SCZ25_EUR]="$REPO/genetic_analysis/work/scores_scz2025/gwas/SCZ25_EUR.ma"
                 [MDD_EUR]="$REPO/legacy/hpc_v2/work/results_v2/prs_final/gwas/MDD_eur.ma" )
mkdir -p "$W/ref" "$W/tmp"; export TMPDIR="$W/tmp"
win() { local t=${1#w}; echo "${t%%_multi}" | tr _ ,; }   # w35_10 -> 35,10
