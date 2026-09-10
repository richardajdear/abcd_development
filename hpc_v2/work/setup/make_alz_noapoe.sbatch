#!/bin/bash
#SBATCH --job-name=v2_noapoe
#SBATCH --partition=icelake
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#
# ALZ with the APOE region removed, chr19:44.4-46.5 Mb GRCh37 -- the same
# window and the same reference .bim as the earlier control run, so the two are
# comparable.  GRCh37 is correct here even though our target is GRCh38: the ALZ
# rsIDs were mapped through this GRCh37 reference in the first place (the build
# mismatch that cost us 99.7% of the file once already, README §13.x).
#
# WHY this cell matters.  With APOE in, Alzheimer's looks as strongly associated
# with faster thinning as schizophrenia does -- C+T -0.041 p=0.0087, SBayesR
# -0.025, SBayesRC -0.028.  That is not evidence that polygenic AD risk affects
# adolescent cortex; APOE is one locus of very large effect and "is APOE
# associated with cortical thinning" is a different question.  The specificity
# claim (SCZ > MDD > ALZ ~ ASD) rests on the APOE-excluded score, so it has to
# be a cell in the grid rather than a footnote.
set -e
source /home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/setup/paths.sh
awk 'NR==FNR{if($1=="19" && $4>=44400000 && $4<=46500000) a[$2]=1; next}
     FNR==1{print; next} !($1 in a)' \
     "$IN/ref/g1000_eur.bim" "$FINAL/gwas/ALZ.ma" > "$FINAL/gwas/ALZ_noAPOE.ma"
echo "ALZ        : $(( $(wc -l < "$FINAL/gwas/ALZ.ma") - 1 )) SNPs"
echo "ALZ_noAPOE : $(( $(wc -l < "$FINAL/gwas/ALZ_noAPOE.ma") - 1 )) SNPs"
