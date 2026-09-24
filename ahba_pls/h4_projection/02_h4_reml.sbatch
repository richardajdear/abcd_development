#!/bin/bash
#SBATCH --job-name=h4_reml
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-6
#SBATCH --output=slurm/%x_%A_%a.log
#
# H4 on CSD3: GCTA GREML of the projection phenotypes, IDENTICAL to
# genetic_analysis/step5_reml.sbatch (the project's canonical h2: dense imputed
# GRM restricted to the PC-AiR unrelated set, sex + site discrete, baseline_age
# + n_visits + 10 PCs quantitative, from the FID = IID HCP export).  Only the
# phenotype file differs: work/pheno_h4.tsv (laptop, IID = NDAR_INV...) joined
# onto the aligned export by the 8-char NDAR token and written in the export's
# own FID/IID spelling, so GCTA's literal FID+IID match (rule 1) holds.
# global_slope is carried through the SAME join as the positive control: its
# h2 here must equal results_70tab_hcp/reml_imp_pooled/global_slope.hsq.
#
#   cd ~/rds/hpc-work/abcd_development
#   PARC=hcp sbatch --export=ALL,PARC=hcp ahba_pls/h4_projection/02_h4_reml.sbatch
set -euo pipefail
export PARC=hcp
REPO=/home/rajd2/rds/hpc-work/abcd_development
source "$REPO/genetic_analysis/config.sh"
HERE="$REPO/ahba_pls/h4_projection"
PHENO_FID="${PHENO_FID:-${PHENO_DIR}_fidiid}"
GRM_UNREL="${GRM_UNREL:-$V1_LEGACY/results/grm_imp_pooled/abcd_imp.pcair}"
OUT="$HERE/results/reml"
[[ "$(basename "$PHENO_DIR")" == *hcp* ]] || { echo "PHENO_DIR=$PHENO_DIR is not the HCP export" >&2; exit 2; }
COVAR_QUANT="$PHENO_FID/covar_quant.txt"; COVAR_CAT="$PHENO_FID/covar_categorical.txt"
EXPORT="$PHENO_FID/phenotypes_gcta.txt"; H4="$HERE/work/pheno_h4.tsv"
require_paths EXPORT COVAR_QUANT COVAR_CAT H4 GCTA
[[ -f "$GRM_UNREL.grm.bin" ]] || { echo "FATAL: $GRM_UNREL incomplete" >&2; exit 2; }
mkdir -p "$OUT" slurm

# ---- join once (task 1 writes, the others wait for the file) ----------------
PHENO="$HERE/work/pheno_h4_fidiid.txt"
if [[ "${SLURM_ARRAY_TASK_ID:-1}" == "1" ]]; then
  awk 'BEGIN{OFS=" "}
    function tok(s){ return substr(s, length(s)-7) }
    FNR==NR { if (FNR>1) { row[tok($1)]=$0 } ; if (FNR==1) { hdr=$0 } ; next }
    FNR==1 { split(hdr, h, "\t"); printf "FID IID"; for (i=2;i<=length(h);i++) printf " %s", h[i]; print ""; next }
    { t=tok($2); if (t in row) { split(row[t], v, "\t"); printf "%s %s", $1, $2; for (i=2;i<=length(v);i++) printf " %s", v[i]; print "" } }
  ' "$H4" "$EXPORT" > "$PHENO"
  echo "joined $(($(wc -l < "$PHENO")-1)) of $(($(wc -l < "$EXPORT")-1)) export rows onto $(($(wc -l < "$H4")-1)) laptop rows"
  head -2 "$PHENO"
else
  for i in $(seq 1 60); do [[ -s "$PHENO" && $(wc -l < "$PHENO") -gt 8000 ]] && break; sleep 10; done
fi
names=($(head -1 "$PHENO" | cut -d' ' -f3-))
idx=$((SLURM_ARRAY_TASK_ID - 1))
(( idx < ${#names[@]} )) || { echo "no phenotype for task $SLURM_ARRAY_TASK_ID"; exit 0; }
pheno="${names[$idx]}"; mpheno=$((idx + 1))
echo "=== h4_reml  $pheno  (--mpheno $mpheno)  $(date) ==="
run "$GCTA" --reml --grm "$GRM_UNREL" --pheno "$PHENO" --mpheno "$mpheno" \
  --qcovar "$COVAR_QUANT" --covar "$COVAR_CAT" --thread-num "${THREADS:-8}" --out "$OUT/$pheno"
{
  printf 'phenotype\th2\tse\tpval\tn\n'
  for p in "${names[@]}"; do f="$OUT/$p.hsq"; [[ -f "$f" ]] || continue
    awk -v p="$p" '$1=="V(G)/Vp"{h2=$2;se=$3} $1=="Pval"{pv=$2} $1=="n"{n=$2} END{if(h2!="")printf "%s\t%.4f\t%.4f\t%s\t%s\n",p,h2,se,pv,n}' "$f"; done
} > "$OUT/reml_summary.tsv"
cat "$OUT/reml_summary.tsv"; echo "=== done $(date) ==="
