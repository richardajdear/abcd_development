#!/bin/bash
# Local validation of the hpc/ pipeline against the synthetic fixture.
#
#   python tools/make_test_genotypes.py     # build the fixture first
#   bash tools/local_test.sh
#
# WHAT THIS CAN AND CANNOT CHECK
#
# GCTA ships statically linked against Intel MKL, which aborts on any CPU
# without AVX.  Under Rosetta on Apple Silicon that is every invocation, so the
# numerical steps (GRM, REML, fastGWA) cannot be executed on this machine at all.
# The failure is in MKL, not in GCTA's argument handling, and that distinction is
# what makes a useful test possible:
#
#   GCTA parses every option and opens every input file BEFORE it does any
#   linear algebra.
#
# So an invocation that reaches "Computing GRM..." and then dies in MKL has
# already proved that its options are valid, its files exist and parse, its IDs
# intersect, and --mpheno is in range.  An invocation with a bad option dies
# earlier, at "invalid option", without touching MKL.  Distinguishing those two
# failure points is a real test of every command the pipeline issues -- it is how
# the --geno/--hwe options (PLINK's, not GCTA's) were caught here rather than on
# the cluster.
#
# Stage 1 below runs that check on every GCTA command in the pipeline.
# Stage 2 computes a GRM with PLINK 1.9, which runs natively, to confirm the
# fixture's relatedness structure is what the GRM step will see.
# Stage 2b is the exception to the MKL problem: --grm-cutoff and --make-bK-sparse
# read an existing GRM and threshold it without entering MKL, so GCTA's own
# relatedness handling DOES execute here and is checked against expectation.
# Stage 3 tests the summary collectors against output files in GCTA's format.
#
# So the coverage is better than "options only": the relatedness logic runs for
# real, and only the three MKL-bound numerical steps (GRM element computation,
# REML, fastGWA) remain unexecuted.
#
# On a Linux cluster none of this substitutes for actually running the pipeline;
# it is a pre-flight for the code, not evidence about the data.
set -uo pipefail
cd "$(dirname "$0")/.."
source hpc/config.sh

pass=0; fail=0
ok()   { echo "  PASS  $*"; pass=$((pass+1)); }
bad()  { echo "  FAIL  $*"; fail=$((fail+1)); }

echo "=============================================================="
echo "Stage 0: build a GRM with PLINK so 02/03 have their input"
echo "=============================================================="
# 02_reml and 03_gwas require_paths GRM, and GCTA cannot build one here, so
# without this they abort before issuing a single command and stages 1-3 would
# silently test only 01.  PLINK writes GCTA's own binary GRM format, so the
# files it produces are the files GCTA would read.
if ! command -v plink >/dev/null 2>&1; then
  echo "  plink not on PATH -- cannot build a GRM, so 02/03 cannot be validated." >&2
  echo "  conda install -c bioconda plink" >&2
  exit 2
fi
ensure_dirs "$GRM_DIR"
for spec in "$GRM:" "$GRM_UNREL:--rel-cutoff $GRM_CUTOFF" "$GRM_SPARSE:"; do
  prefix="${spec%%:*}"; extra="${spec#*:}"
  # shellcheck disable=SC2086
  plink --bfile "$GENO" --autosome --maf "$MAF" $extra \
        --make-grm-bin --out "$prefix" >/dev/null 2>&1
  [[ -f "$prefix.grm.bin" ]] || { echo "  FAILED to build $prefix.grm.bin" >&2; exit 2; }
done
echo "  built $(ls "$GRM_DIR"/*.grm.bin | wc -l | tr -d ' ') GRMs"
echo "  full  $(wc -l < "$GRM.grm.id" | tr -d ' ') subjects"
echo "  unrel $(wc -l < "$GRM_UNREL.grm.id" | tr -d ' ') subjects after --rel-cutoff $GRM_CUTOFF"
echo
echo "NOTE: PLINK's --rel-cutoff and GCTA's --grm-cutoff use different pruning"
echo "      algorithms, so the retained set here is not identical to what the"
echo "      cluster will produce.  This GRM exists to unblock stages 1-3, not to"
echo "      stand in for the real one."

echo
echo "=============================================================="
echo "Stage 1: GCTA option and input validation"
echo "=============================================================="
# Collect every GCTA command the pipeline would run, via DRY_RUN.
cmds=$(mktemp)
for s in hpc/01_grm.sbatch hpc/02_reml.sbatch hpc/03_gwas.sbatch; do
  DRY_RUN=1 bash "$s" 2>/dev/null | sed -n 's/^DRY_RUN: //p'
done | grep -F "$GCTA" > "$cmds"
n_cmds=$(wc -l < "$cmds" | tr -d ' ')
echo "$n_cmds GCTA invocations to validate"
# 3 from 01, then one per phenotype from each of 02 and 03.
n_expected=$(( 3 + 2 * $(phenotype_names | wc -l | tr -d ' ') ))
if [[ "$n_cmds" -ne "$n_expected" ]]; then
  bad "collected $n_cmds commands but expected $n_expected -- a step aborted before issuing any"
fi
echo

n=0
while IFS= read -r cmd; do
  n=$((n+1))
  label=$(echo "$cmd" | grep -oE '\-\-(make-grm|make-bK-sparse|reml|fastGWA-mlm)' | head -1 || true)
  mph=$(echo "$cmd" | grep -oE '\-\-mpheno [0-9]+' | head -1 || true)
  # `|| true`: every invocation here is EXPECTED to exit non-zero (MKL aborts on
  # this CPU), and under `pipefail` an unguarded failure would end the loop after
  # the first command -- silently testing 1 of 15.
  out=$(eval "$cmd" 2>&1 || true)

  if echo "$out" | grep -qi "invalid option"; then
    bad "$n ${label} ${mph}: $(echo "$out" | grep -i 'invalid option' | head -1 | tr -s ' ')"
  elif echo "$out" | grep -qiE "Error:|error occurs" && ! echo "$out" | grep -qi "MKL"; then
    bad "$n ${label} ${mph}: $(echo "$out" | grep -iE 'Error:' | head -1 | tr -s ' ')"
  elif echo "$out" | grep -qi "MKL FATAL"; then
    # Reached the linear algebra: options, files and IDs all accepted.
    ok "$n ${label} ${mph} -- options and inputs accepted (stopped at MKL)"
  else
    # No MKL error and no error at all: on a machine with AVX this is a real run.
    ok "$n ${label} ${mph} -- ran to completion"
  fi
done < "$cmds"
rm -f "$cmds"

echo
echo "=============================================================="
echo "Stage 2: relatedness structure of the fixture"
echo "=============================================================="
# Read the Stage 0 GRM in list format and compare within- against
# between-family pairs.  If this fails, the fixture has no relatives and
# --grm-cutoff, the REML/fastGWA N difference, and the sparse GRM are all
# untested -- so it fails loudly rather than being treated as cosmetic.
plink --bfile "$GENO" --autosome --maf "$MAF" \
      --make-grm-gz no-gz --out "$GRM_DIR/plink_list" >/dev/null 2>&1
if [[ -f "$GRM_DIR/plink_list.grm" ]]; then
  # GCTA .grm list format: i j n_snps relatedness (1-based indices into .fam).
  read -r within between npairs <<< "$(awk '
    NR==FNR {fid[FNR]=$1; next}
    $1!=$2 { if (fid[$1]==fid[$2]) {w+=$4; nw++} else {b+=$4; nb++} }
    END {printf "%.3f %.4f %d", (nw?w/nw:0), (nb?b/nb:0), nw}
  ' "$GENO.fam" "$GRM_DIR/plink_list.grm")"
  echo "  within-family  $within  ($npairs pairs)"
  echo "  between-family $between"
  awk -v w="$within" 'BEGIN{exit !(w>0.35 && w<0.65)}' \
    && ok "within-family relatedness is sibling-like" \
    || bad "within-family relatedness $within is not sibling-like -- fixture has no relatives, so --grm-cutoff is untested"
  awk -v b="$between" 'BEGIN{exit !(b>-0.05 && b<0.05)}' \
    && ok "between-family relatedness is ~0" \
    || bad "between-family relatedness $between is not ~0"
  n_full=$(wc -l < "$GRM.grm.id" | tr -d ' ')
  n_unrel=$(wc -l < "$GRM_UNREL.grm.id" | tr -d ' ')
  [[ "$n_unrel" -lt "$n_full" ]] \
    && ok "relatedness pruning removes subjects ($n_full -> $n_unrel)" \
    || bad "pruning removed nobody ($n_full -> $n_unrel): the REML/fastGWA N difference is untested"
else
  bad "plink --make-grm-gz produced no .grm"
fi

echo
echo "=============================================================="
echo "Stage 2b: GCTA's own relatedness handling (runs natively)"
echo "=============================================================="
# --grm-cutoff and --make-bK-sparse read an existing GRM and threshold it
# without calling MKL, so unlike every other GCTA step these DO execute here.
# That makes the relatedness logic -- the part most likely to be silently wrong
# -- genuinely testable on a laptop.
if [[ -f "$GRM_UNREL.grm.id" ]]; then
  n_gcta=$(wc -l < "$GRM_UNREL.grm.id" | tr -d ' ')
  echo "  GCTA --grm-cutoff $GRM_CUTOFF retained $n_gcta of 1500"
  [[ "$n_gcta" -lt 1500 && "$n_gcta" -gt 300 ]] \
    && ok "GCTA pruning retained a plausible unrelated subset ($n_gcta)" \
    || bad "GCTA pruning retained $n_gcta -- implausible"
fi
if [[ -f "$GRM_SPARSE.grm.sp" ]]; then
  npair=$(wc -l < "$GRM_SPARSE.grm.sp" | tr -d ' ')
  # .grm.sp is i j value, 0-based.  Off-diagonal entries are the retained
  # relative pairs; they should be near 0.5, not near the sparse cutoff.
  offdiag=$(awk '$1!=$2 {s+=$3; n++} END {if(n) printf "%.3f", s/n; else print "NA"}' \
            "$GRM_SPARSE.grm.sp")
  echo "  sparse GRM: $npair retained pairs, mean off-diagonal $offdiag"
  awk -v v="$offdiag" 'BEGIN{exit !(v>0.35 && v<0.65)}' \
    && ok "sparse GRM off-diagonal is sibling-like ($offdiag)" \
    || bad "sparse GRM off-diagonal $offdiag -- fastGWA would model the wrong relatedness"
fi

echo
echo "=============================================================="
echo "Stage 3: summary collectors"
echo "=============================================================="
# The collectors parse GCTA's output formats with awk.  GCTA cannot run here, so
# feed them files in those formats and check the parsed numbers.  This tests the
# parsing, not GCTA's output -- if a future GCTA changes its layout, only the
# cluster run will reveal it, which is why 02/03 resolve columns by name.
tmp=$(mktemp -d)
cat > "$tmp/x.hsq" <<'HSQ'
Source	Variance	SE
V(G)	0.213456	0.061234
V(e)	0.267891	0.058123
Vp	0.481347	0.021456
V(G)/Vp	0.443512	0.118765
logL	-1234.567
logL0	-1240.123
LRT	11.112
df	1
Pval	4.2e-04
n	892
HSQ
read -r h2 se pv nn <<< "$(awk '
  $1=="V(G)/Vp" {h2=$2; se=$3} $1=="Pval" {pv=$2} $1=="n" {n=$2}
  END {printf "%.4f %.4f %s %s", h2, se, pv, n}' "$tmp/x.hsq")"
[[ "$h2" == "0.4435" && "$se" == "0.1188" && "$pv" == "4.2e-04" && "$nn" == "892" ]] \
  && ok ".hsq parser: h2=$h2 se=$se p=$pv n=$nn" \
  || bad ".hsq parser got h2=$h2 se=$se p=$pv n=$nn"

# fastGWA output with a known median chi-square.  Columns deliberately in a
# different order than GCTA's default, to prove name-based resolution works.
{
  echo -e "CHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tSE\tBETA\tP"
  # BETA/SE ratios 1,2,3,4,5 -> chi2 1,4,9,16,25 -> median 9
  for i in 1 2 3 4 5; do
    echo -e "1\trs$i\t${i}000\tA\tG\t1000\t0.3\t1\t$i\t0.5"
  done
} > "$tmp/y.fastGWA"
read -r cN cBETA cSE cP < <(awk 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; print c["N"], c["BETA"], c["SE"], c["P"]; exit}' "$tmp/y.fastGWA")
[[ "$cN" == "6" && "$cBETA" == "9" && "$cSE" == "8" && "$cP" == "10" ]] \
  && ok "column resolution by name: N=$cN BETA=$cBETA SE=$cSE P=$cP" \
  || bad "column resolution gave N=$cN BETA=$cBETA SE=$cSE P=$cP"

med=$(awk -v b="$cBETA" -v s="$cSE" 'NR>1 && $s>0 {print ($b/$s)^2}' "$tmp/y.fastGWA" \
      | sort -g | awk '{v[NR]=$1} END{if(NR%2) print v[(NR+1)/2]; else print (v[NR/2]+v[NR/2+1])/2}')
[[ "$med" == "9" ]] \
  && ok "median chi-square = $med (expected 9)" \
  || bad "median chi-square = $med, expected 9"

# The bug this replaces: taking the middle of INPUT order rather than of sorted
# values.  Reverse the rows; a correct median is unchanged.
{ head -1 "$tmp/y.fastGWA"; tail -n +2 "$tmp/y.fastGWA" | tail -r 2>/dev/null || tail -n +2 "$tmp/y.fastGWA" | tac; } > "$tmp/y_rev.fastGWA"
med_rev=$(awk -v b="$cBETA" -v s="$cSE" 'NR>1 && $s>0 {print ($b/$s)^2}' "$tmp/y_rev.fastGWA" \
      | sort -g | awk '{v[NR]=$1} END{if(NR%2) print v[(NR+1)/2]; else print (v[NR/2]+v[NR/2+1])/2}')
[[ "$med_rev" == "$med" ]] \
  && ok "median is invariant to input order" \
  || bad "median changed with input order: $med vs $med_rev (not sorting)"
rm -rf "$tmp"

echo
echo "=============================================================="
echo "Stage 4: 04_magma against a MAGMA stub"
echo "=============================================================="
# Real MAGMA has no macOS build, so $MAGMA here is scratch/bin/magma: a stub
# that records its arguments and emits files in MAGMA's output formats.  This
# tests that 04_magma builds valid invocations, locates GENE/ZSTAT by name, and
# parses .gsa.out -- NOT any statistical result.  Every number below is the
# stub's, so only the structure is meaningful.
if [[ -x "${MAGMA:-}" ]] && [[ -f "$MAGMA_REF.bed" ]]; then
  ensure_dirs "$GWAS_DIR"
  # 04_magma needs fastGWA input; synthesise it in GCTA's format.
  # BETA must be signed and centred near zero: LDSC's munge_sumstats refuses a
  # signed statistic whose median is far from the stated null (it reads that as
  # a mislabelled column), and an all-positive BETA is not a plausible GWAS.
  # Every SNP also needs both alleles and a p-value consistent with BETA/SE.
  for p in $(phenotype_names); do
    { printf 'CHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP\n'
      awk -v OFS='\t' 'NR<=300{
        b = ((NR % 2) ? 1 : -1) * (0.002 + 0.004 * (NR % 23))
        se = 0.02
        z = b / se; if (z < 0) z = -z
        # Two-sided normal tail, good enough for a fixture.
        p = 2 * exp(-0.717 * z - 0.416 * z * z); if (p > 1) p = 1
        print $1, $2, $4, "A", "G", 7800, 0.31, b, se, p
      }' "$GENO.bim"
    } > "$GWAS_DIR/$p.fastGWA"
  done
  if bash hpc/04_magma.sbatch >/dev/null 2>&1; then
    ok "04_magma ran to completion"
  else
    bad "04_magma exited non-zero"
  fi
  n_pheno=$(phenotype_names | wc -l | tr -d ' ')
  # One gene analysis per phenotype, 2 gene-set tests each, 2 disorders x
  # n_pheno reverse tests.
  n_gsa=$(ls "$MAGMA_DIR"/*.gsa.out 2>/dev/null | wc -l | tr -d ' ')
  n_exp=$(( 2 * n_pheno + 2 * n_pheno ))
  [[ "$n_gsa" -eq "$n_exp" ]] \
    && ok "gene-set tests: $n_gsa .gsa.out files (expected $n_exp)" \
    || bad "gene-set tests: $n_gsa .gsa.out files, expected $n_exp"
  # The genez file is the one real transform in the step: GENE + ZSTAT by name,
  # with the phenotype as the covariate column name.
  gz="$MAGMA_DIR/global_slope_genez.txt"
  if [[ -f "$gz" ]]; then
    [[ "$(head -1 "$gz")" == "GENE	global_slope" ]] \
      && ok "genez header names the phenotype as the covariate" \
      || bad "genez header is '$(head -1 "$gz")'"
    # ZSTAT is column 8 of the stub's .genes.out; check the value tracks it.
    z_src=$(awk 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} NR==2{print $c["ZSTAT"]}' \
            "$MAGMA_DIR/global_slope.genes.out")
    z_out=$(awk 'NR==2{print $2}' "$gz")
    [[ "$z_src" == "$z_out" ]] \
      && ok "genez ZSTAT resolved by name ($z_out)" \
      || bad "genez ZSTAT is $z_out but .genes.out has $z_src"
  else
    bad "no $gz written"
  fi
  n_rows=$(( $(wc -l < "$MAGMA_DIR/magma_summary.tsv" | tr -d ' ') - 1 ))
  [[ "$n_rows" -eq $(( 2 * n_gsa )) ]] \
    && ok "magma_summary.tsv collected $n_rows rows from $n_gsa files" \
    || bad "magma_summary.tsv has $n_rows rows from $n_gsa files (expected $(( 2 * n_gsa )))"
else
  echo "  SKIP  no MAGMA stub or reference panel; see hpc/config.local.sh.example"
fi

echo
echo "=============================================================="
echo "Stage 5: 05_ldsc_rg (LDSC runs natively)"
echo "=============================================================="
# LDSC is pure Python, so unlike GCTA it runs here in full.  What the fixture
# CANNOT support is a defined rg: 5,000 SNPs with synthetic LD scores give LDSC
# far too little information, so rg comes back NA by design.  The checks below
# therefore test the machinery -- munge accepts the fastGWA and disorder files,
# the sign check passes, the log is parsed by name -- and deliberately do NOT
# assert anything about rg's value.
if [[ -x "${LDSC_MUNGE:-}" && -f "$HM3_SNPLIST" ]]; then
  ensure_dirs "$LDSC_DIR"
  if bash hpc/05_ldsc_rg.sbatch >"$LDSC_DIR/ldsc_run.log" 2>&1; then
    ok "05_ldsc ran to completion"
  else
    bad "05_ldsc exited non-zero (see $LDSC_DIR/ldsc_run.log)"
  fi
  n_ss=$(ls "$LDSC_DIR"/*.sumstats.gz 2>/dev/null | wc -l | tr -d ' ')
  n_exp=$(( $(phenotype_names | wc -l | tr -d ' ') + 2 ))
  [[ "$n_ss" -eq "$n_exp" ]] \
    && ok "munged $n_ss sumstats files (phenotypes + SCZ + MDD)" \
    || bad "munged $n_ss sumstats files, expected $n_exp"
  # The sign check catches a mislabelled effect column, which would silently
  # reverse rg.  Assert it against the OUTPUT file rather than against log text:
  # munge_one reuses an existing .sumstats.gz (correct on the cluster, where
  # re-munging a 10M-SNP file is expensive), so on a second run the log contains
  # no munge output at all and a grep-based check passes or fails depending on
  # whether the tree happened to be clean.
  ss="$LDSC_DIR/global_slope.sumstats.gz"
  if [[ -f "$ss" ]]; then
    # LDSC sumstats format: SNP A1 A2 N Z.  Z must be signed and centred.
    read -r zmed nz <<< "$(gzip -dc "$ss" | awk '
      NR==1 {for(i=1;i<=NF;i++) c[$i]=i; next}
      c["Z"] && $c["Z"] != "NA" {v[++n]=$c["Z"]+0}
      END {
        if (!n) {print "NA 0"; exit}
        # median of sorted values
        for (i=1;i<=n;i++) for (j=i+1;j<=n;j++) if (v[j]<v[i]) {t=v[i];v[i]=v[j];v[j]=t}
        printf "%.4f %d", (n%2 ? v[(n+1)/2] : (v[n/2]+v[n/2+1])/2), n
      }')"
    if [[ "$nz" -gt 0 ]] && awk -v m="$zmed" 'BEGIN{exit !(m>-0.5 && m<0.5)}'; then
      ok "munged Z is signed and centred (median $zmed over $nz SNPs)"
    else
      bad "munged Z median is $zmed over $nz SNPs -- effect column may be mislabelled"
    fi
  else
    bad "no $ss to check the signed statistic"
  fi
  s="$LDSC_DIR/ldsc_rg_summary.tsv"
  if [[ -f "$s" ]]; then
    n_rows=$(( $(wc -l < "$s" | tr -d ' ') - 1 ))
    n_want=$(( 2 * $(phenotype_names | wc -l | tr -d ' ') ))
    [[ "$n_rows" -eq "$n_want" ]] \
      && ok "rg summary has $n_rows rows (2 disorders x phenotypes)" \
      || bad "rg summary has $n_rows rows, expected $n_want"
    # Phenotype names must be recovered from the .sumstats.gz paths.
    if awk -F'\t' 'NR>1 && $2 ~ /\// {found=1} END{exit !found}' "$s"; then
      bad "phenotype column still contains file paths"
    else
      ok "phenotype names resolved from sumstats paths"
    fi
    # The underpowered flag is the guard against over-reading a null rg.
    awk -F'\t' 'NR>1 && ($11=="yes" || $11=="no") {n++} END{exit !(n>0)}' "$s" \
      && ok "underpowered flag computed from h2 z-score" \
      || bad "underpowered flag missing or malformed"
  else
    bad "no ldsc_rg_summary.tsv written"
  fi
else
  echo "  SKIP  no LDSC binary or HapMap3 snplist; see hpc/config.local.sh.example"
fi

echo
echo "=============================================================="
echo "Stage 6: 06_prs (PLINK clump+score, then lmer in R)"
echo "=============================================================="
if [[ -x "${PLINK:-}" ]] && [[ -f "$MAGMA_REF.bed" ]]; then
  ensure_dirs "$PRS_DIR"
  if bash hpc/06_prs.sbatch >"$PRS_DIR/prs_run.log" 2>&1; then
    ok "06_prs ran to completion"
  else
    bad "06_prs exited non-zero (see $PRS_DIR/prs_run.log)"
  fi
  n_thr=$(echo $PRS_THRESHOLDS | wc -w | tr -d ' ')
  n_prof=$(ls "$PRS_DIR"/score_*.profile 2>/dev/null | wc -l | tr -d ' ')
  [[ "$n_prof" -eq $(( 2 * n_thr )) ]] \
    && ok "scored $n_prof profiles (2 disorders x $n_thr thresholds)" \
    || bad "scored $n_prof profiles, expected $(( 2 * n_thr ))"
  # A score with no variance would make every association exactly null while
  # looking superficially fine, so check the scores actually differ by subject.
  f=$(ls "$PRS_DIR"/score_*.profile 2>/dev/null | head -1)
  if [[ -n "$f" ]]; then
    nuniq=$(awk 'NR>1{print $NF}' "$f" | sort -u | wc -l | tr -d ' ')
    [[ "$nuniq" -gt 100 ]] \
      && ok "scores vary across subjects ($nuniq distinct values)" \
      || bad "only $nuniq distinct score values -- score is near-constant"
  fi
  # The association step needs R; skip rather than fail if it is absent.
  a="$PRS_DIR/prs_association.tsv"
  if [[ -f "$a" ]]; then
    n_rows=$(( $(wc -l < "$a" | tr -d ' ') - 1 ))
    ok "association table has $n_rows rows"
    # Both strata must be present -- the EUR/full contrast is the design.
    for st in EUR full; do
      awk -F'\t' -v s="$st" 'NR>1 && $4==s {n++} END{exit !(n>0)}' "$a" \
        && ok "stratum '$st' present" \
        || bad "stratum '$st' missing from association table"
    done
    # EUR must be a strict subset, or the ancestry cut did nothing.
    read -r n_eur n_full <<< "$(awk -F'\t' 'NR>1 {if($4=="EUR") e=$5; else f=$5} END{print e, f}' "$a")"
    [[ -n "$n_eur" && -n "$n_full" && "$n_eur" -lt "$n_full" ]] \
      && ok "EUR subset is smaller than full sample ($n_eur < $n_full)" \
      || bad "EUR n=$n_eur vs full n=$n_full -- ancestry cut did not restrict"
    # Random weights against a simulated phenotype: essentially nothing should
    # survive correction.  A pile of hits here means the model is miscalibrated
    # (most likely the family random effect is not being fitted).
    n_sig=$(awk -F'\t' 'NR>1 && $12+0 < 0.05 {n++} END{print n+0}' "$a")
    [[ "$n_sig" -le 2 ]] \
      && ok "null calibration: $n_sig significant after correction (expect ~0)" \
      || bad "null calibration: $n_sig significant after correction -- model likely miscalibrated"
  else
    echo "  SKIP  prs_association.tsv absent (Rscript unavailable in this env)"
    echo "        run in the 'r' env: Rscript tools/prs_assoc.R --prs-dir $PRS_DIR ..."
  fi
else
  echo "  SKIP  no PLINK or reference panel"
fi

echo
echo "=============================================================="
echo "Stage 7: run_all.sh reaches every step"
echo "=============================================================="
# The entry point is where bash-4-isms hide: `declare -A` for the step table
# aborted on bash 3.2, and a missing $LOG_DIR made `tee` kill the chain at step
# 01.  Both failures LOOK like a clean early exit, so assert that a dry run
# actually announces all six steps.
seen=$(DRY_RUN=1 bash hpc/run_all.sh 2>&1 | grep -cE '^### 0[1-6] ')
[[ "$seen" -eq 6 ]] \
  && ok "run_all.sh dry run reached all 6 steps" \
  || bad "run_all.sh dry run reached $seen of 6 steps"

echo
echo "=============================================================="
echo "$pass passed, $fail failed"
echo "=============================================================="
[[ "$fail" == "0" ]] || exit 1
