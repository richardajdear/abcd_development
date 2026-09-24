#!/bin/bash
# Local end-to-end test of steps 1-3 on the reference panel itself (no ABCD
# data), plus a planted-effect test of step 4.  Writes the aggregate accuracy
# table results/smoketest_accuracy.tsv (committed); everything else goes to
# smoketest/out/ (gitignored).
#
#   PLINK=<plink 1.9> BCFTOOLS=bcftools JAVA=java RSCRIPT=Rscript PY=python \
#     bash c4_imputation/smoketest/run_smoketest.sh
#
# Imputation leg: 5-fold leave-fold-out.  Each fold's ~22 held-out individuals
# are converted to an array-like PLINK fileset (C4 removed, SNPs thinned,
# unphased, 5 % strand-flipped, 50 decoy SNPs), then run through 01 -> 02 -> 03
# exactly as the ABCD data will be, and scored against the true C4 alleles.
# Two densities: every panel SNP (thin 1) and every third (thin 3, close to the
# ~2-3 k MHC SNPs an array carries).  The reference in each fold is the other
# ~89 individuals, so accuracy here is a LOWER bound on the full 111-person panel.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; C4="$(dirname "$HERE")"; O="$HERE/out"; mkdir -p "$O"
PLINK="${PLINK:-plink}"; BCFTOOLS="${BCFTOOLS:-bcftools}"; JAVA="${JAVA:-java}"; PY="${PY:-python}"; RSCRIPT="${RSCRIPT:-Rscript}"
export PLINK BCFTOOLS JAVA
export BEAGLE_JAR="$C4/resources/beagle.27Feb25.75f.jar" CONFORM_JAR="$C4/resources/conform-gt.24May16.cee.jar"
FULLPANEL="$C4/resources/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37.vcf.gz"
ACC="$C4/results/smoketest_accuracy.tsv"

for THIN in 1 3; do
  prefixes=()
  for FOLD in 1 2 3 4 5; do
    P="$O/fold${FOLD}_thin${THIN}"
    "$PY" "$HERE/make_test_data.py" --panel "$FULLPANEL" --out "$P" --fold "$FOLD" --nfold 5 --thin "$THIN" --bcftools "$BCFTOOLS"
    "$PLINK" --vcf "$P.target.vcf" --double-id --keep-allele-order --make-bed --out "$P.target" > /dev/null
    bash "$C4/01_extract_mhc.sh" "$P.target" "$P.mhc"
    PANEL="$P.ref.vcf.gz" bash "$C4/02_impute_c4.sh" "$P.mhc.vcf.gz" "$P.imputed" 4 4 | tail -n 1
    "$PY" "$C4/03_c4_grex.py" --vcf "$P.imputed.vcf.gz" --out-calls "$P.calls.tsv" \
      --out-summary "$P.summary.tsv" --panel "$P.ref.vcf.gz" --bcftools "$BCFTOOLS"
    prefixes+=("$P")
  done
  "$PY" "$HERE/evaluate.py" --pairs "${prefixes[@]}" --label "leave_fold_out_thin${THIN}" --out "$ACC"
done

# Association leg: planted effect on synthetic children built from panel haplotypes.
"$PY" "$HERE/make_assoc_fixture.py" --panel "$FULLPANEL" --out "$O/assoc_fixture" --n 3000 --beta -0.10
"$RSCRIPT" "$C4/04_c4_assoc.R" --c4 "$O/assoc_fixture/c4_calls.tsv" --pheno-dir "$O/assoc_fixture/pheno" \
  --eur-ids "$O/assoc_fixture/eur.keep" --prs-file "$O/assoc_fixture/prs.tsv" \
  --out "$O/assoc_fixture/c4_assoc.tsv"
"$PY" - "$O/assoc_fixture/c4_assoc.tsv" <<'EOF'
import sys, pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t")
r = t[(t.model == "M1_primary") & (t.term == "C4A_GREx_z")].iloc[0]
lo, hi = r.beta - 1.96 * r.se, r.beta + 1.96 * r.se
assert lo <= -0.10 <= hi, f"planted -0.10 not recovered: {r.beta:.3f} +/- {r.se:.3f}"
b = t[(t.model == "M5_control_baseline") & (t.term == "C4A_GREx_z")].iloc[0]
assert abs(b.beta / b.se) < 3.5, f"null control phenotype shows an effect: {b.beta:.3f}"
print(f"planted -0.100 recovered as {r.beta:.3f} (95% CI {lo:.3f}, {hi:.3f}), n={int(r.n)}; control null ok")
EOF
# Staging leg: the two CSD3 sbatch scripts, run as plain bash with every path
# overridden (fold-1 target as "the array", the fixture as "the phenotypes").
SB="$O/sbatch"; mkdir -p "$SB"
awk '{print $2}' "$O/fold1_thin3.target.fam" > "$SB/eur.keep"
REPO="$(dirname "$C4")" GENO="$O/fold1_thin3.target" C4_WORK="$SB/work" C4_RESULTS="$SB/results" \
  PY="$PY" EUR_KEEP="$SB/eur.keep" SLURM_CPUS_PER_TASK=2 bash "$C4/c4_impute.sbatch" | tail -n 4
[[ -s "$SB/results/c4_imputation_summary.tsv" && -s "$SB/work/c4_calls.tsv" ]]
F="$O/assoc_fixture"
REPO="$(dirname "$C4")" C4_WORK="$SB/work" C4_RESULTS="$SB/results" RSCRIPT="$RSCRIPT" PY="$PY" \
  CALLS="$F/c4_calls.tsv" EUR_KEEP="$F/eur.keep" PRS_DIR="$F/scores/SBayesRC/SCZ25_EUR" \
  PH_HCP_1LMM="$F/pheno" PH_DK_1LMM="$F/pheno" PH_HCP_PERREGION="$F/pheno" \
  SCORE_ROOT="$F/scores" MHC_BIM="$F/chr6.bim" bash "$C4/c4_assoc.sbatch" | grep -E '^===|^M1|n_mhc|SBayesRC|done'
for c in hcp_1lmm dk_1lmm hcp_perregion; do [[ -s "$SB/results/c4_assoc_$c.tsv" ]]; done
got=$(awk -F'\t' 'NR==2{print $5}' "$SB/results/mhc_in_scz_scores.tsv"); [[ "$got" == "$(cat "$F/expected_mhc.txt")" ]] \
  || { echo "MHC count $got != expected $(cat "$F/expected_mhc.txt")" >&2; exit 4; }
echo "smoketest PASSED"
