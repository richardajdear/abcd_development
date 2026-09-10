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
for SRC in ALZ ALZ_IGAP; do
  [[ -s "$FINAL/gwas/$SRC.ma" ]] || { echo "skip $SRC (not built yet)"; continue; }
  awk 'NR==FNR{if($1=="19" && $4>=44400000 && $4<=46500000) a[$2]=1; next}
       FNR==1{print; next} !($1 in a)' \
       "$IN/ref/g1000_eur.bim" "$FINAL/gwas/$SRC.ma" > "$FINAL/gwas/${SRC}_noAPOE.ma"
  echo "$SRC        : $(( $(wc -l < "$FINAL/gwas/$SRC.ma") - 1 )) SNPs"
  echo "${SRC}_noAPOE : $(( $(wc -l < "$FINAL/gwas/${SRC}_noAPOE.ma") - 1 )) SNPs"
  # independent check in the TARGET's build (GRCh38), since the window above is
  # GRCh37: nothing genome-wide significant may survive in the APOE region
  awk 'NR==FNR{if($1=="19" && $4>=43500000 && $4<=46500000) w[$2]=1; next}
       FNR>1 && ($1 in w) && $7+0<1e-8 {n++} END{printf "  residual APOE p<1e-8 in %s_noAPOE: %d\n", "'"$SRC"'", n+0}' \
       "$GENO.bim" "$FINAL/gwas/${SRC}_noAPOE.ma"
done
