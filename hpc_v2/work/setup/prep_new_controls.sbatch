#!/bin/bash
#SBATCH --job-name=v2_newctrl
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#
# Normalise the two GWAS added to test whether the ALZ_noAPOE association (and
# therefore the SCZ one) is an SES/education signal rather than disorder
# biology, then build the APOE-excluded versions.
#
#   ALZ_IGAP  Kunkle 2019 stage 1, 21,982 CLINICALLY DIAGNOSED cases / 41,944
#             controls, European.  No UK Biobank by-proxy cases -- which is the
#             whole point, since Wightman 2021 has them and a proxy phenotype
#             (parental dementia, self-reported) carries parental longevity,
#             SES and education with it.
#   EA        Okbay 2016 years of education, N = 405,072.  A positive control
#             for the confound itself, not a disorder control.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
"$PYTHON" "$V2ROOT/work/setup/normalise_gwas.py" "$GWAS_DIR" "$FINAL/gwas" \
  "$IN/ref/maf_union.tsv" ALZ_IGAP,EA
bash "$V2ROOT/work/setup/make_alz_noapoe.sbatch" || true
wc -l "$FINAL"/gwas/ALZ_IGAP*.ma "$FINAL"/gwas/EA.ma
