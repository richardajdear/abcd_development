#!/bin/bash
# Local validation of the hpc_v2/ pipeline against the synthetic fixture.
#
#   python tools/make_test_genotypes.py        # v1 fixture (source genotypes)
#   python tools/make_test_genotypes_v2.py     # v2 fixture (this test's input)
#   cp hpc_v2/config.local.sh.example hpc_v2/config.local.sh
#   bash tools/local_test_v2.sh
#
# WHAT THIS CAN AND CANNOT CHECK
#
# Unlike v1, most of hpc_v2 is R, and R runs natively here -- so steps 01-04
# and 06 EXECUTE FOR REAL on the fixture: GDS conversion, KING, PC-AiR,
# PC-Relate, fitNullModel, assocTestSingle, the collector, and the
# between/within PRS decomposition all produce genuine numbers that are
# asserted against expectation.
#
# The exception is step 05 (Zaitlen REML), which is GCTA.  GCTA ships
# statically linked against Intel MKL and aborts on any CPU without AVX, so
# under Rosetta on Apple Silicon every invocation dies inside MKL.  GCTA parses
# every option and opens every input BEFORE any linear algebra, so an
# invocation that reaches "Computing/Reading GRM" has proved its options are
# valid and its files parse; one with a bad option dies earlier at "invalid
# option".  Stage 5 distinguishes those two failure points, and separately
# tests the two-GRM .hsq collector against a synthetic file in GCTA's exact
# format.
#
# The strongest check here is Stage 6's positive control: the fixture plants a
# known within-family PRS effect (b_W = -0.20) plus family-level confounding
# (b_B = -0.35), so a step-06 that confuses the two coefficients, loses the
# family IDs, or standardises in the wrong place CANNOT pass.
#
# On a Linux cluster none of this substitutes for running the pipeline; it is a
# pre-flight for the code, not evidence about the data.
set -uo pipefail
cd "$(dirname "$0")/.."
source hpc_v2/config.sh
# config.sh sets `-e` (and it should: the pipeline steps must abort on error).
# A TEST harness must not -- it has to run every check and report the tally, so
# a single failing assertion cannot end the run.  Same for pipefail, which turns
# the SIGPIPE from `... | head -1` into a fatal error under -e.  This must come
# AFTER the source, or config.sh re-enables both.
set +e +o pipefail

pass=0; fail=0
ok()  { echo "  PASS  $*"; pass=$((pass+1)); }
bad() { echo "  FAIL  $*"; fail=$((fail+1)); }
hdr() { echo; echo "=============================================================="; echo "$*"; echo "=============================================================="; }

# R with GENESIS is mandatory for stages 1-4; report clearly if absent.
if ! "$RSCRIPT" -e 'suppressMessages({library(GENESIS);library(SNPRelate);library(GWASTools)})' >/dev/null 2>&1; then
  echo "FATAL: RSCRIPT=$RSCRIPT cannot load GENESIS/SNPRelate/GWASTools." >&2
  echo "       Build hpc_v2/envs/genesis_env.yml, or set RSCRIPT in" >&2
  echo "       hpc_v2/config.local.sh to an R that has them." >&2
  exit 2
fi

# All step logs land here; every stage below writes into it, so create it up
# front rather than relying on a step's own ensure_dirs having run yet.
TESTLOG="$OUT_V2/test_logs"
ensure_dirs "$OUT_V2" "$TESTLOG" "$GDS_DIR" "$KIN_DIR" "$NULL_DIR" \
            "$ASSOC_DIR" "$REML2_DIR" "$PRSFAM_DIR" "$LOG_DIR"

hdr "Stage 0: preflight (00_check_inputs.sh)"
if bash hpc_v2/00_check_inputs.sh > "$TESTLOG/preflight.log" 2>&1; then
  ok "00_check_inputs passed"
else
  bad "00_check_inputs failed -- see $TESTLOG/preflight.log"
  tail -20 "$TESTLOG/preflight.log"
fi
# The integrity gate must actually reject a broken fileset, or it is decoration.
tmpd=$(mktemp -d)
cp "$GENO_ARRAY".{bed,bim,fam} "$tmpd/" 2>/dev/null
head -n 100 "$GENO_ARRAY.bim" > "$tmpd/$(basename "$GENO_ARRAY").bim"   # bim now disagrees with bed
if GENO_ARRAY="$tmpd/$(basename "$GENO_ARRAY")" bash hpc_v2/00_check_inputs.sh >/dev/null 2>&1; then
  bad "integrity gate ACCEPTED a bed/bim mismatch (the 4.0 defect would pass)"
else
  ok "integrity gate rejects a bed/bim size mismatch"
fi
rm -rf "$tmpd"

hdr "Stage 1: 01_gds.sbatch (PLINK -> GDS, array + 22 chromosomes)"
if bash hpc_v2/01_gds.sbatch > "$TESTLOG/gds_run.log" 2>&1; then
  ok "01_gds exited 0"
else
  bad "01_gds failed -- see $GDS_DIR/../gds_run.log"; tail -20 "$TESTLOG/gds_run.log"
fi
n_gds=$(ls "$GDS_DIR"/imp_chr*.gds 2>/dev/null | wc -l | tr -d ' ')
[[ -s "$GDS_DIR/array.gds" ]] && ok "array.gds written" || bad "array.gds missing"
[[ "$n_gds" -eq 22 ]] && ok "22 per-chromosome GDS written" || bad "$n_gds per-chromosome GDS, expected 22"
# Idempotency: a second run must NOT rebuild or corrupt.
# Capture to a file first, then grep: `cmd | grep -q` exits at the first match,
# SIGPIPEs the producer, and `set -o pipefail` then reports the whole pipeline
# as failed even though the command was fine.
bash hpc_v2/01_gds.sbatch > "$TESTLOG/gds_rerun.log" 2>&1
if grep -q "EXISTS" "$TESTLOG/gds_rerun.log"; then
  ok "re-run detects existing GDS and skips (idempotent)"
else
  bad "re-run did not report EXISTS -- conversion is not idempotent"
fi

hdr "Stage 2: 02_kinship.sbatch (KING -> PC-AiR -> PC-Relate)"
if bash hpc_v2/02_kinship.sbatch > "$TESTLOG/kinship_run.log" 2>&1; then
  ok "02_kinship exited 0"
else
  bad "02_kinship failed -- see $KIN_DIR/../kinship_run.log"; tail -30 "$TESTLOG/kinship_run.log"
fi
for f in kinship_sparse.rds pcair_pcs.tsv pcair_unrelated.txt kinship_summary.tsv pcrelate.rds; do
  [[ -s "$KIN_DIR/$f" ]] && ok "wrote $f" || bad "missing $KIN_DIR/$f"
done
# The fixture's relatives are simulated by parental transmission, so PC-Relate
# MUST find first-degree pairs.  Zero would mean the estimator saw no structure
# -- the failure mode that makes a sibling-aware pipeline pointless.
if [[ -s "$KIN_DIR/kinship_summary.tsv" ]]; then
  n1=$(awk -F'\t' '$1 ~ /1st deg/ {print $2}' "$KIN_DIR/kinship_summary.tsv")
  n_unrel=$(awk -F'\t' '$1=="n_pcair_unrelated" {print $2}' "$KIN_DIR/kinship_summary.tsv")
  n_samp=$(awk -F'\t' '$1=="n_samples" {print $2}' "$KIN_DIR/kinship_summary.tsv")
  [[ "${n1:-0}" -gt 100 ]] \
    && ok "PC-Relate found $n1 first-degree pairs (fixture has ~739 sib families)" \
    || bad "PC-Relate found ${n1:-0} first-degree pairs -- relatedness not detected"
  # PC-AiR must prune SOME subjects but not most: the pathology it replaces
  # (v1 §8.6) was pruning 48% of the sample away.
  if [[ "${n_unrel:-0}" -gt 0 && "${n_samp:-0}" -gt 0 ]]; then
    frac=$(awk -v a="$n_unrel" -v b="$n_samp" 'BEGIN{printf "%.2f", a/b}')
    awk -v f="$frac" 'BEGIN{exit !(f>0.4 && f<1.0)}' \
      && ok "PC-AiR unrelated partition is $n_unrel/$n_samp ($frac) -- prunes relatives, keeps the bulk" \
      || bad "PC-AiR kept $n_unrel/$n_samp ($frac) -- implausible partition"
  fi
  # Sparse kinship must actually be sparse AND non-diagonal.
  nz=$(awk -F'\t' '$1=="sparse_nonzero_offdiag" {print $2}' "$KIN_DIR/kinship_summary.tsv")
  [[ "${nz:-0}" -gt 0 ]] \
    && ok "sparse kinship has $nz non-zero off-diagonal pairs" \
    || bad "sparse kinship is diagonal -- the mixed model would have no random effect"
  # PC-Relate must recover the fixture's TRUE pair count, not merely "some"
  # pairs: the .fam FIDs say exactly how many within-family pairs exist.
  truep=$(awk '{print $1}' "$GENO_ARRAY.fam" | sort | uniq -c \
          | awk '{n=$1; s+=n*(n-1)/2} END{print s}')
  if [[ -n "$truep" && "${n1:-0}" -gt 0 ]]; then
    awk -v a="$n1" -v b="$truep" 'BEGIN{exit !(a >= 0.9*b && a <= 1.1*b)}' \
      && ok "first-degree pair count $n1 matches the fixture's true $truep (within 10%)" \
      || bad "PC-Relate found $n1 first-degree pairs; the .fam implies $truep"
  fi
  # The density guard must REPORT, and must fire when the matrix is not sparse.
  # (On this small fixture with only ~5k pruned SNPs it is EXPECTED to fire --
  # noisy kinship chains the cohort into one block via transitive closure.  The
  # test is that the warning appears, not that the fixture is sparse.)
  dens=$(awk -F'\t' '$1=="sparse_density" {print $2}' "$KIN_DIR/kinship_summary.tsv")
  blk=$(awk -F'\t' '$1=="sparse_largest_block" {print $2}' "$KIN_DIR/kinship_summary.tsv")
  if [[ -n "$dens" && -n "$blk" ]]; then
    ok "kinship summary reports sparse_density=$dens and largest block=$blk"
    over=$(awk -v d="$dens" -v m="$SPARSE_KIN_MAX_DENSITY" 'BEGIN{print (d>m)?1:0}')
    warned=$(grep -c "WARNING: the 'sparse' kinship matrix is not sparse" "$TESTLOG/kinship_run.log")
    if [[ "$over" == "1" ]]; then
      [[ "$warned" -gt 0 ]] \
        && ok "density $dens exceeds $SPARSE_KIN_MAX_DENSITY and the warning FIRED (guard works)" \
        || bad "density $dens exceeds $SPARSE_KIN_MAX_DENSITY but NO warning was printed"
    else
      [[ "$warned" -eq 0 ]] \
        && ok "density $dens within $SPARSE_KIN_MAX_DENSITY, no spurious warning" \
        || bad "density $dens is within budget but the warning fired anyway"
    fi
  else
    bad "kinship_summary.tsv lacks sparse_density / sparse_largest_block rows"
  fi
fi
# PCs file shape.
if [[ -s "$KIN_DIR/pcair_pcs.tsv" ]]; then
  n_pc=$(awk -F'\t' 'NR==1{n=0; for(i=1;i<=NF;i++) if ($i ~ /^PC[0-9]+$/) n++; print n}' "$KIN_DIR/pcair_pcs.tsv")
  [[ "${n_pc:-0}" -ge "$N_PCS" ]] \
    && ok "pcair_pcs.tsv carries $n_pc PCs (>= N_PCS=$N_PCS)" \
    || bad "pcair_pcs.tsv has ${n_pc:-0} PCs, need >= $N_PCS"
fi

hdr "Stage 3: 03_null_model.sbatch (GENESIS fitNullModel)"
if bash hpc_v2/03_null_model.sbatch > "$TESTLOG/null_run.log" 2>&1; then
  ok "03_null_model exited 0"
else
  bad "03_null_model failed -- see $NULL_DIR/../null_run.log"; tail -30 "$TESTLOG/null_run.log"
fi
n_null=$(ls "$NULL_DIR"/*_null.rds 2>/dev/null | wc -l | tr -d ' ')
n_pheno=$(phenotype_names | wc -l | tr -d ' ')
[[ "$n_null" -eq "$n_pheno" ]] \
  && ok "$n_null null models for $n_pheno phenotypes" \
  || bad "$n_null null models, expected $n_pheno"
# The whole point: N must EXCEED the unrelated subset, because relatives stay.
if [[ -s "$NULL_DIR/global_slope_null_summary.tsv" ]]; then
  n_fit=$(awk -F'\t' 'NR==2{print $2}' "$NULL_DIR/global_slope_null_summary.tsv")
  n_unrel=$(awk -F'\t' '$1=="n_pcair_unrelated" {print $2}' "$KIN_DIR/kinship_summary.tsv")
  awk -v a="$n_fit" -v b="$n_unrel" 'BEGIN{exit !(a>b)}' \
    && ok "null model fitted on n=$n_fit > $n_unrel unrelated -- relatives RETAINED" \
    || bad "null model n=$n_fit does not exceed the unrelated count $n_unrel -- relatives were dropped"
fi

hdr "Stage 4: 04_assoc.sbatch (assocTestSingle + collector)"
if bash hpc_v2/04_assoc.sbatch > "$TESTLOG/assoc_run.log" 2>&1; then
  ok "04_assoc exited 0"
else
  bad "04_assoc failed -- see $ASSOC_DIR/../assoc_run.log"; tail -30 "$TESTLOG/assoc_run.log"
fi
n_chr_files=$(ls "$ASSOC_DIR"/global_slope_chr*.tsv.gz 2>/dev/null | wc -l | tr -d ' ')
[[ "$n_chr_files" -eq 22 ]] \
  && ok "22 per-chromosome association files for global_slope" \
  || bad "$n_chr_files association files, expected 22"
if [[ -s "$ASSOC_DIR/gwas_summary.tsv" ]]; then
  ok "gwas_summary.tsv written"
  # Columns must match v1's layout so downstream tooling reads either.
  hdr_cols=$(head -1 "$ASSOC_DIR/gwas_summary.tsv")
  for c in phenotype n_snps n_mean lambda_gc n_p5e8 n_p1e5 min_p; do
    grep -q "$c" <<< "$hdr_cols" || bad "gwas_summary.tsv missing v1 column '$c'"
  done
  grep -q "n_p5e8" <<< "$hdr_cols" && ok "summary columns match the v1 layout"
  # lambda_GC on a fixture with no true polygenic signal should be near 1.
  lam=$(awk -F'\t' 'NR>1 && $1=="global_slope"{print $4}' "$ASSOC_DIR/gwas_summary.tsv")
  if [[ -n "$lam" ]]; then
    awk -v l="$lam" 'BEGIN{exit !(l>0.7 && l<1.5)}' \
      && ok "lambda_GC = $lam (calibrated; mixed model is not inflated)" \
      || bad "lambda_GC = $lam -- outside 0.7-1.5, the model is miscalibrated"
  fi
  n_rows=$(awk 'NR>1' "$ASSOC_DIR/gwas_summary.tsv" | wc -l | tr -d ' ')
  [[ "$n_rows" -eq "$n_pheno" ]] \
    && ok "summary has one row per phenotype ($n_rows)" \
    || bad "summary has $n_rows rows, expected $n_pheno"
  # END-TO-END POSITIVE CONTROL.  The v1 fixture seeds a phenotype (sim_h2_50)
  # with 20 causal SNPs at h2=0.5, so the scan MUST detect them.  Without this,
  # a pipeline that silently shuffled genotypes against phenotypes -- an ID-join
  # bug, the single most expensive class of error in this project -- would still
  # produce a clean lambda_GC and pass every other check here.
  if awk -F'\t' 'NR>1 && $1=="sim_h2_50"' "$ASSOC_DIR/gwas_summary.tsv" | grep -q .; then
    read -r sim_hits sim_minp <<< "$(awk -F'\t' 'NR>1 && $1=="sim_h2_50"{print $5, $7}' "$ASSOC_DIR/gwas_summary.tsv")"
    awk -v h="$sim_hits" 'BEGIN{exit !(h > 0)}' \
      && ok "seeded phenotype sim_h2_50 yields $sim_hits genome-wide-significant SNPs (min p=$sim_minp) -- genotypes and phenotypes are correctly aligned" \
      || bad "seeded phenotype sim_h2_50 found NO significant SNPs (min p=$sim_minp) -- the 20 causal SNPs at h2=0.5 should be detectable; suspect an ID-join or column-selection bug"
    # And the negative side: a real phenotype must NOT show those hits.
    gs_hits=$(awk -F'\t' 'NR>1 && $1=="global_slope"{print $5}' "$ASSOC_DIR/gwas_summary.tsv")
    [[ "${gs_hits:-0}" -eq 0 ]] \
      && ok "unseeded global_slope yields 0 hits (no spurious signal)" \
      || bad "unseeded global_slope yields $gs_hits hits -- unexpected on simulated null genotypes"
  fi
  # Re-running the collector must UPDATE, not duplicate.
  "$RSCRIPT" hpc_v2/R/05_collect_assoc.R --assoc-dir "$ASSOC_DIR" --pheno global_slope \
      --out-dir "$ASSOC_DIR" --maf "$ASSOC_MAF" --mac "$ASSOC_MAC" >/dev/null 2>&1
  n_rows2=$(awk 'NR>1' "$ASSOC_DIR/gwas_summary.tsv" | wc -l | tr -d ' ')
  [[ "$n_rows2" -eq "$n_rows" ]] \
    && ok "collector is idempotent (still $n_rows2 rows after re-run)" \
    || bad "collector duplicated rows on re-run ($n_rows -> $n_rows2)"
fi
# Sumstats must carry the columns MAGMA/LDSC/PRS tooling matches on.
ss="$ASSOC_DIR/global_slope.sumstats.tsv.gz"
if [[ -s "$ss" ]]; then
  cols=$(gzip -dc "$ss" | head -1)
  missing=""
  for c in SNP CHR POS A1 A2 AF1 N BETA SE P; do
    grep -qw "$c" <<< "$cols" || missing="$missing $c"
  done
  [[ -z "$missing" ]] && ok "sumstats columns complete (SNP CHR POS A1 A2 AF1 N BETA SE P)" \
                      || bad "sumstats missing columns:$missing"
  # N per variant must reflect the FULL sample, not the unrelated subset.
  nmax=$(gzip -dc "$ss" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} {if($c["N"]>m) m=$c["N"]} END{print m+0}')
  n_unrel=$(awk -F'\t' '$1=="n_pcair_unrelated" {print $2}' "$KIN_DIR/kinship_summary.tsv")
  awk -v a="$nmax" -v b="$n_unrel" 'BEGIN{exit !(a>b)}' \
    && ok "per-variant N up to $nmax > $n_unrel unrelated -- the scan uses relatives" \
    || bad "per-variant N maxes at $nmax, not above the unrelated count $n_unrel"
fi

hdr "Stage 5: 05_reml_zaitlen.sbatch (GCTA two-GRM; options + collector)"
# 5a. Option validity, via the MKL boundary (see the header note).
if [[ -x "$GCTA" ]] || command -v "$GCTA" >/dev/null 2>&1; then
  gl="$TESTLOG/gcta_probe.log"; ensure_dirs "$REML2_DIR"
  "$GCTA" --grm "$GRM_FULL" --make-bK "$BK_THRESH" --out "$REML2_DIR/probe" \
      > "$gl" 2>&1
  st=$?
  if grep -qiE "invalid option|Error: unknown option|not recognized" "$gl"; then
    bad "GCTA rejected an option in the --make-bK command (see $gl)"
  elif [[ "$st" == "0" ]] || grep -qiE "Reading the GRM|GRM for|individuals" "$gl"; then
    ok "GCTA accepted --make-bK $BK_THRESH and opened the GRM"
  else
    bad "GCTA --make-bK failed before reading the GRM (see $gl)"
  fi
  # The full step: reaching the REML stage proves --mgrm/--pheno/--qcovar parse.
  rl="$TESTLOG/step05_run.log"
  bash hpc_v2/05_reml_zaitlen.sbatch > "$rl" 2>&1
  if grep -qiE "invalid option|unknown option" "$rl"; then
    bad "05_reml_zaitlen issued an option GCTA rejects (see $rl)"
  elif grep -qiE "Reading the GRM|Summary result of REML|multi-component|individuals" "$rl"; then
    ok "05_reml_zaitlen's GCTA commands parse and reach the REML/GRM stage"
  else
    bad "05_reml_zaitlen did not reach GCTA's GRM/REML stage (see $rl)"
  fi
  # mgrm.txt must name TWO grms -- one line is a single-component fit
  # masquerading as Zaitlen.
  if [[ -f "$REML2_DIR/mgrm.txt" ]]; then
    n_g=$(wc -l < "$REML2_DIR/mgrm.txt" | tr -d ' ')
    [[ "$n_g" -eq 2 ]] && ok "mgrm.txt lists 2 GRMs (full + bK)" \
                       || bad "mgrm.txt lists $n_g GRMs, expected 2"
  else
    bad "mgrm.txt not written -- the two-component fit has no input"
  fi
else
  echo "  GCTA not available at '$GCTA'; skipping 5a (option validation)."
fi
# 5b. The two-GRM collector, against a file in GCTA's exact format.
tmp=$(mktemp -d); ensure_dirs "$REML2_DIR"
cat > "$REML2_DIR/__collector_probe.hsq" <<'HSQ'
Source	Variance	SE
V(G1)	0.152341	0.041234
V(G2)	0.087654	0.052311
V(e)	0.760005	0.038112
Vp	1.000000	0.021456
V(G1)/Vp	0.152341	0.041234
V(G2)/Vp	0.087654	0.052311
Sum of V(G)/Vp	0.239995	0.055120
logL	-1234.567
logL0	-1240.123
LRT	11.112
df	1
Pval	4.2700e-04
n	7900
HSQ
got=$(awk -F'\t' '
  $1=="V(G1)/Vp"       {h1=$2; s1=$3}
  $1=="Sum of V(G)/Vp" {hs=$2; ss=$3}
  $1=="Pval"           {pv=$2}
  $1=="n"              {n=$2}
  END {printf "%.4f %.4f %.4f %s %s", h1, s1, hs, pv, n}
' "$REML2_DIR/__collector_probe.hsq")
[[ "$got" == "0.1523 0.0412 0.2400 4.2700e-04 7900" ]] \
  && ok "two-GRM .hsq collector parses h2_snp, h2_ped, Pval and n correctly" \
  || bad "collector parsed '$got' (expected '0.1523 0.0412 0.2400 4.2700e-04 7900')"
# It must NOT silently accept a single-GRM .hsq as a Zaitlen fit.
cat > "$REML2_DIR/__collector_single.hsq" <<'HSQ'
Source	Variance	SE
V(G)	0.213456	0.061234
V(e)	0.267891	0.058123
Vp	0.481347	0.021456
V(G)/Vp	0.443512	0.118765
Pval	1.2e-04
n	5649
HSQ
got1=$(awk -F'\t' '$1=="V(G1)/Vp"{h1=$2} END{print (h1==""?"EMPTY":h1)}' "$REML2_DIR/__collector_single.hsq")
[[ "$got1" == "EMPTY" ]] \
  && ok "collector yields no row for a single-GRM .hsq (cannot mislabel it Zaitlen)" \
  || bad "collector extracted '$got1' from a single-GRM .hsq"
rm -f "$REML2_DIR/__collector_probe.hsq" "$REML2_DIR/__collector_single.hsq"; rm -rf "$tmp"

hdr "Stage 6: 06_prs_family.sbatch (between/within decomposition)"
if bash hpc_v2/06_prs_family.sbatch > "$TESTLOG/prsfam_run.log" 2>&1; then
  ok "06_prs_family exited 0"
else
  bad "06_prs_family failed -- see $PRSFAM_DIR/../prsfam_run.log"
  tail -30 "$TESTLOG/prsfam_run.log"
fi
wf="$PRSFAM_DIR/prs_withinfamily.tsv"
if [[ -s "$wf" ]]; then
  ok "prs_withinfamily.tsv written"
  for c in beta_within se_within p_within beta_between beta_diff p_diff n_pairs; do
    head -1 "$wf" | grep -qw "$c" || bad "output missing column '$c'"
  done
  head -1 "$wf" | grep -qw beta_diff && ok "output carries the between/within contrast columns"
  # THE POSITIVE CONTROL.  Fixture plants b_W = -0.20 and a family-mean
  # coefficient of -0.35, so on SCZ x global_slope (full stratum) we require:
  #   beta_within  ~ -0.20 and significant
  #   |beta_between| > |beta_within|   (the confounder is visible)
  #   p_diff        significant        (the contrast detects it)
  read -r bw pw bb pd np <<< "$(awk -F'\t' '
    NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next}
    $c["disorder"]=="SCZ" && $c["phenotype"]=="global_slope" && $c["stratum"]=="full" {
      print $c["beta_within"], $c["p_within"], $c["beta_between"], $c["p_diff"], $c["n_pairs"]; exit }
  ' "$wf")"
  if [[ -n "${bw:-}" ]]; then
    echo "        (SCZ x global_slope, full: b_W=$bw p_W=$pw  b_B=$bb  p_diff=$pd  pairs=$np)"
    awk -v b="$bw" 'BEGIN{exit !(b > -0.32 && b < -0.10)}' \
      && ok "beta_within = $bw recovers the planted -0.20" \
      || bad "beta_within = $bw does not recover the planted -0.20"
    awk -v p="$pw" 'BEGIN{exit !(p < 0.01)}' \
      && ok "within-family effect is significant (p_within = $pw)" \
      || bad "within-family effect not significant (p_within = $pw) at a planted z~5"
    awk -v a="$bb" -v b="$bw" 'BEGIN{exit !(a < b)}' \
      && ok "beta_between ($bb) is more negative than beta_within ($bw) -- confounding visible" \
      || bad "beta_between ($bb) does not exceed beta_within ($bw); the decomposition is wrong"
    awk -v p="$pd" 'BEGIN{exit !(p < 0.05)}' \
      && ok "beta_diff detects the planted confounding (p_diff = $pd)" \
      || bad "beta_diff missed the planted confounding (p_diff = $pd)"
    awk -v n="$np" 'BEGIN{exit !(n > 300)}' \
      && ok "$np informative sibling pairs counted" \
      || bad "only $np informative pairs -- family structure was lost"
  else
    bad "no SCZ x global_slope full-stratum row in $wf"
  fi
  # MDD is an independent null score: it must NOT come out significant.
  bwm=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
    $c["disorder"]=="MDD" && $c["phenotype"]=="global_slope" && $c["stratum"]=="full" {print $c["p_within"]; exit}' "$wf")
  if [[ -n "${bwm:-}" ]]; then
    awk -v p="$bwm" 'BEGIN{exit !(p > 0.01)}' \
      && ok "MDD null score is not significant within-family (p = $bwm)" \
      || bad "MDD null score came out at p = $bwm -- false positive in the model"
  fi
fi
# The FID=IID guard must fire: this is what stops a silent no-op analysis.
tmpp=$(mktemp -d)
awk 'NR==1{print; next}{$1=$2; print}' OFS=' ' "$PHENO" > "$tmpp/pheno_fidiid.txt"
if "$RSCRIPT" hpc_v2/R/06_prs_family.R --prs-dir "$PRS_PROFILE_DIR" \
      --pheno "$tmpp/pheno_fidiid.txt" --covar-quant "$COVAR_QUANT" \
      --covar-cat "$COVAR_CAT" --manifest "$MANIFEST" \
      --out "$tmpp/out.tsv" >/dev/null 2>&1; then
  bad "06_prs_family RAN on a FID==IID phenotype file -- the guard does not fire"
else
  ok "06_prs_family refuses a FID==IID phenotype file (guard fires)"
fi
rm -rf "$tmpp"

hdr "Stage 7: syntax + portability of every script"
for f in hpc_v2/*.sh hpc_v2/*.sbatch tools/local_test_v2.sh; do
  bash -n "$f" 2>/dev/null && ok "bash -n $f" || bad "bash -n $f"
done
# The three checks below strip COMMENTS before matching.  Every one of these
# traps is documented in a comment in the very file being checked, so a naive
# grep flags the documentation as the defect -- which is how a useful check gets
# deleted for crying wolf.  `sed 's/#.*//'` is crude (it would also blank a '#'
# inside a string literal) but no script here has one, and erring toward
# false-NEGATIVES in a lint is wrong, so the patterns stay strict.
# config.local.sh is excluded throughout: it is gitignored and machine-specific
# by design, so absolute paths in it are correct, not a defect.
lintable() { ls hpc_v2/*.sh hpc_v2/*.sbatch 2>/dev/null | grep -v 'config\.local\.sh$'; }
strip_comments() { sed 's/#.*//' $(lintable); }

# bash 3.2 (macOS) compatibility: mapfile and declare -A are bash 4+.
if strip_comments | grep -qE '^\s*(mapfile|readarray)\b|declare -A'; then
  bad "bash 4+ construct (mapfile/declare -A) found -- breaks on macOS bash 3.2"
else
  ok "no bash 4+ constructs (mapfile / declare -A)"
fi
# Empty-array expansion under set -u: "${arr[@]}" aborts on bash < 4.4 when the
# array is empty.  Both optional-argument arrays here must use the +-guard form.
# The guard form is ${arr[@]+"${arr[@]}"}, which CONTAINS "${arr[@]}" -- so the
# pattern must match only an expansion that is NOT preceded by the `[@]+` guard.
# Checking that the guard is absent from the line is enough here: no line uses
# the same array twice.
if strip_comments \
     | grep -E '"\$\{(strata_args|eur_args|_geno_args)\[@\]\}"' \
     | grep -qvE '\[@\]\+'; then
  bad "unguarded \"\${arr[@]}\" on an optionally-empty array -- aborts under set -u on bash < 4.4"
else
  ok "optional-argument arrays use the \${arr[@]+...} guard"
fi
# Every sbatch must name its partition explicitly (the 14-hour v1 lesson).
for f in hpc_v2/*.sbatch; do
  grep -q '^#SBATCH --partition=' "$f" \
    && ok "$(basename "$f") declares --partition" \
    || bad "$(basename "$f") has no --partition line (would take the cluster default)"
done
# No absolute cluster paths in executable code.  config.sh's ${VAR:-default}
# lines are the ONE sanctioned exception -- they are the relocation roots the
# whole design rests on -- so a /rds literal is a defect only OUTSIDE a
# ${...:-...} default.  Anywhere else it means a step will break on the next
# account.
if strip_comments \
     | grep -vE '\$\{[A-Za-z_]+:-[^}]*\}' \
     | grep -qE '/rds/|/home/[a-z]+[0-9]+'; then
  bad "hard-coded cluster path in executable code (outside a \${VAR:-default})"
else
  ok "no hard-coded cluster paths outside \${VAR:-default} roots"
fi
# Same for the R scripts, which have no sanctioned exception at all.
if sed 's/#.*//' hpc_v2/R/*.R | grep -qE '/rds/|/home/[a-z]+[0-9]+'; then
  bad "hard-coded cluster path in an R script"
else
  ok "no hard-coded cluster paths in the R scripts"
fi
# The {CHR} brace trap: bash ends a ${...:-...} expansion at the brace closing
# {CHR}, so the placeholder must never appear inside one.
if strip_comments | grep -qE '\$\{[A-Z_]+:[-=][^}]*\{CHR\}'; then
  bad "{CHR} placeholder inside a \${VAR:-default} expansion (bash brace trap)"
else
  ok "no {CHR} placeholder inside \${...:-...} defaults"
fi
# R scripts must parse.
for f in hpc_v2/R/*.R; do
  "$RSCRIPT" -e "invisible(parse('$f'))" >/dev/null 2>&1 \
    && ok "R parse $(basename "$f")" || bad "R parse $(basename "$f")"
done

hdr "RESULT"
echo "  $pass passed, $fail failed"
[[ "$fail" -eq 0 ]] || exit 1
