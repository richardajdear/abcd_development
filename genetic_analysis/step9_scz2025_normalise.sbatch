#!/bin/bash
#SBATCH --job-name=scz25_norm
#SBATCH --account=vertes-sl3-cpu
#SBATCH --partition=icelake
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/%x_%A.log
#
# Step 9 (2026-09-17) -- the Nature 2025 multi-ancestry schizophrenia GWAS
# (https://www.nature.com/articles/s41586-025-10000-6; SCZ + schizoaffective,
# MVP + PGC3 + All of Us + FinnGen for EUR; MVP + AoU + GPC + MGS for AFR;
# PGC3 for EAS) -> the pipeline's `.ma` format (SNP A1 A2 freq b se p N), one
# file per discovery arm, exactly as setup/normalise_gwas.py did for PGC3.
#
# Raw files: $SCZ25_SRC (GWAS-VCF BCFs, hg38, rsID in ID; README.txt there).
# Output:    $SCZ25_ROOT/gwas/SCZ25_{EUR,AFR,EAS,META}.ma  (+ .sst for PRS-CS,
#            .n_gwas = max per-variant effective N, gwas_summary.tsv)
#
# Conventions, matched to the PGC3 files so the comparison is like for like:
#   * A1 = ALT = the allele ES refers to; A2 = REF.
#   * freq = discovery ALT frequency (AF), as the PGC3 SCZ files used their own
#     FCAS/FCON-weighted frequency.  (Rule 5's target-MAF borrowing was for
#     files with no frequency column; these have one.)
#   * N = NE (bcftools +metal's per-variant effective N) PLUS the missing
#     half of the PGC3 cohorts' Neff -- see EFFECTIVE N CORRECTION below.
#   * p = 10^-LP.
#   * Keep: rsID-named biallelic SNPs (single-base REF and ALT), FILTER PASS
#     or unset (drop IFFY / REF_MISMATCH), finite ES/SE/LP/AF/NE, se > 0,
#     0 < AF < 1.  rsIDs occurring more than once (multi-allelic splits) are
#     dropped entirely -- PLINK --score and GCTB both stop on a duplicate ID.
#   * Autosomes only: the target fileset (abcd_imp_prs) has no chrX.
#   * Build does not matter: every downstream join is on rsID.
set -euo pipefail
REPO=/home/rajd2/rds/hpc-work/abcd_development
SRC="${SCZ25_SRC:-/rds/user/rajd2/hpc-work/magma/gwas/SCZ_2025_multiancestry}"
ROOT="${SCZ25_ROOT:-$REPO/genetic_analysis/work/scores_scz2025}"
GWAS="$ROOT/gwas"
BIM="$REPO/genetic_analysis/work/inputs/geno/abcd_imp_prs.bim"
BCFTOOLS="${BCFTOOLS:-/usr/local/software/spack/spack-views/rocky8-icelake-20220710/bcftools-1.14/gcc-11.3.0/intel-oneapi-mkl-2022.1.0/xrcih6dysuehfmkkuvr77tfgyx75sf52/bin/bcftools}"
if [[ ! -x "$BCFTOOLS" ]]; then
  source /usr/share/Modules/init/bash 2>/dev/null || true
  module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
  BCFTOOLS=$(command -v bcftools)
fi
"$BCFTOOLS" --version | head -1
[[ -s "$BIM" ]] || { echo "FATAL: target bim $BIM missing" >&2; exit 2; }
mkdir -p "$GWAS"

declare -A FILE=( [EUR]=SCZ_EUR_autosomes.bcf [AFR]=SCZ_AFR_autosome.bcf
                  [EAS]=SCZ_EAS_autosome.bcf  [META]=SCZ_AFR_EU_EAS_autosome.bcf )
# --- EFFECTIVE N CORRECTION -------------------------------------------------
# bcftools +metal summed each study's own N column into NE.  For the three
# non-PGC cohorts that column is a proper Neff (4*NC*(NS-NC)/NS -- verified on
# single-study rows: AoU EUR 161,682/6,579 -> 25,245).  For the two PGC3 files
# it is PGC's NEFF = Neff/2 (single-study rows: PGC3 european 130,644/53,386
# -> 58,749.1, which is exactly half the 117,498 our PGC3 .ma carries; PGC3
# asian 30,761/14,004 -> 14,424.8).  So NE under-counts the PGC3 half of the
# discovery sample.  N is what PRS-CS (n_gwas) and SBayesR/RC read, and the
# PGC3 run used full Neff, so add the missing half back per variant using the
# per-study direction string ED (position = study order in the metal command,
# verified on single-study rows):
#   EUR  : ED[2]  = PGC3 european  -> + N_pgc3eur(rsID)/2   (from legacy SCZ_eur.ma)
#   META : ED[9]  = PGC3 european  -> + N_pgc3eur(rsID)/2
#          ED[12] = PGC3 asian     -> + NE_eas(rsID)         (from SCZ25_EAS.ma, = Neff/2)
#   EAS  : PGC3 asian only         -> N = 2 x NE
#   AFR  : no PGC3 cohort          -> N = NE
# Constant fallbacks (58,749 / 14,425) when the rsID is not in the lookup.
# The raw NE is kept in gwas_summary.tsv so the size of the correction is on
# record.  EAS is processed first because META needs its table.
PGC3_EUR_MA="$REPO/legacy/hpc_v2/work/results_v2/prs_final/gwas/SCZ_eur.ma"
[[ -s "$PGC3_EUR_MA" ]] || { echo "FATAL: $PGC3_EUR_MA missing (needed for the N correction)" >&2; exit 2; }
declare -A EDLEN=( [EUR]=4 [META]=12 [AFR]=7 [EAS]=0 )
SUMMARY="$GWAS/gwas_summary.tsv"
printf 'arm\tsource\tn_records\tn_kept\tn_dup_rsid_dropped\tn_in_target_bim\tmedian_NE_raw\tmedian_N_corrected\tn_gwas_max_N\tn_p5e8\n' > "$SUMMARY"

for ARM in EAS EUR AFR META; do
  src="$SRC/${FILE[$ARM]}"; out="$GWAS/SCZ25_$ARM.ma"
  tmp="$GWAS/.SCZ25_$ARM.tmp"; dup="$GWAS/.SCZ25_$ARM.dups"; cnt="$GWAS/.SCZ25_$ARM.cnt"
  [[ -s "$src" ]] || { echo "FATAL: $src missing" >&2; exit 3; }
  echo "=== $ARM  $(basename "$src")  $(date)"
  : > "$dup"
  FMT='%ID\t%REF\t%ALT\t%FILTER[\t%ES\t%SE\t%LP\t%AF\t%NE\t%ED]\n'
  [[ "$ARM" == EAS ]] && FMT='%ID\t%REF\t%ALT\t%FILTER[\t%ES\t%SE\t%LP\t%AF\t%NE]\n'
  "$BCFTOOLS" query -f "$FMT" "$src" \
  | awk -v dup="$dup" -v cnt="$cnt" -v mode="$ARM" -v edlen="${EDLEN[$ARM]}" \
        -v pgc3f="$PGC3_EUR_MA" -v easf="$GWAS/SCZ25_EAS.ma" '
      BEGIN {
        if (mode=="EUR" || mode=="META") { while ((getline l < pgc3f) > 0) { split(l, a, " "); if (a[1] != "SNP") h[a[1]] = a[8]/2 } close(pgc3f) }
        if (mode=="META") { while ((getline l < easf) > 0) { split(l, a, " "); if (a[1] != "SNP") e[a[1]] = a[8]/2 } close(easf) }
      }
      { nrec++ }
      $1 ~ /^rs[0-9]+$/ && ($4=="." || $4=="PASS") && $2 ~ /^[ACGT]$/ && $3 ~ /^[ACGT]$/ \
        && $5!="." && $6!="." && $7!="." && $8!="." && $9!="." {
        se=$6+0; af=$8+0; ne=$9+0; lp=$7+0; ed=$10
        if (se<=0 || af<=0 || af>=1 || ne<=0) next
        if (edlen > 0 && length(ed) != edlen) { bad++; next }
        p = (lp > 299) ? 1e-300 : 10^(-lp)
        N = ne
        if (mode=="EAS") N = 2*ne
        else if (mode=="EUR")  { if (substr(ed,2,1) != "?") N += ($1 in h ? h[$1] : 58749) }
        else if (mode=="META") { if (substr(ed,9,1) != "?") N += ($1 in h ? h[$1] : 58749)
                                 if (substr(ed,12,1) != "?") N += ($1 in e ? e[$1] : 14425) }
        if ($1 in seen) d[$1]=1; seen[$1]=1
        printf "%s %s %s %.6f %.6g %.6g %.6g %d %d\n", $1, $3, $2, af, $5, se, p, int(N+0.5), int(ne+0.5)
      }
      END { for (k in d) print k > dup; print nrec, bad+0 > cnt }' > "$tmp"
  read -r nrec nbad < "$cnt"
  (( nbad == 0 )) || { echo "FATAL: $nbad rows with an ED string not of length ${EDLEN[$ARM]} -- study order assumption broken" >&2; exit 5; }
  ndup=$(wc -l < "$dup")
  # column 9 (raw NE) is only for the summary; the .ma keeps the 8 standard columns
  { echo "SNP A1 A2 freq b se p N"
    awk -v dupf="$dup" 'FILENAME==dupf{d[$1]; next} !($1 in d){print $1,$2,$3,$4,$5,$6,$7,$8}' "$dup" "$tmp"; } > "$out"
  read -r nmed_raw nmed_cor nmax n58 < <(awk -v dupf="$dup" 'FILENAME==dupf{d[$1]; next} !($1 in d){i++; r[i]=$9; c[i]=$8; if($8>m)m=$8; if($7<5e-8)h++}
        END{asort(r); asort(c); printf "%d %d %d %d\n", r[int(i/2)+1], c[int(i/2)+1], m, h+0}' "$dup" "$tmp")
  rm -f "$tmp"
  nkept=$(( $(wc -l < "$out") - 1 ))
  nbim=$(awk 'NR==FNR{b[$2]; next} FNR>1 && ($1 in b)' "$BIM" "$out" | wc -l)
  echo "$nmax" > "$GWAS/SCZ25_$ARM.n_gwas"
  # PRS-CS / PRS-CSx summary-statistic file (SNP A1 A2 BETA P); p floored so it
  # survives the p -> z round trip (setup/prs_final.sbatch PRSCS note).
  awk 'BEGIN{OFS="\t"} NR==1{print "SNP","A1","A2","BETA","P"; next}
       {p=$7+0; if (p<=0) p=1e-300; print $1,$2,$3,$5,p}' "$out" > "$GWAS/SCZ25_$ARM.sst"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$ARM" "${FILE[$ARM]}" "$nrec" "$nkept" "$ndup" "$nbim" "$nmed_raw" "$nmed_cor" "$nmax" "$n58" >> "$SUMMARY"
  echo "  records $nrec -> kept $nkept (dup rsIDs dropped $ndup); in target bim $nbim; median NE raw $nmed_raw -> N corrected $nmed_cor; n_gwas (max N) $nmax; p<5e-8: $n58"
  rm -f "$dup" "$cnt"
done
echo; column -t "$SUMMARY"
echo "=== done $(date)"
