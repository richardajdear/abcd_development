"""Put all six discovery GWAS into one format: SNP A1 A2 freq b se p N.

One normaliser for every trait and arm, so every downstream method reads the
same columns and any difference between methods is the method.

  SCZ_pooled  PGC3_SCZ_wave3.primary   EUR+EAS+AFR+LAT (multi-ancestry)
  SCZ_eur     PGC3_SCZ_wave3.european  European only
  MDD_pooled  pgc-mdd2025 ..._div      trans-ancestry
  MDD_eur     pgc-mdd2025 ..._eur      European only
  ASD         SPARK+iPSYCH+PGC         no ancestry-stratified release exists
  ALZ_IGAP    Kunkle 2019 IGAP stage1  clinically diagnosed only, NO proxy cases
  EA          Okbay 2016 EduYears      positive control for the SES/education
                                       confound, not a disorder control
  ALZ         PGC-ALZ2 (Wightman)      no ancestry-stratified release exists

Missing pieces, filled the same way for every trait that needs it:
  freq  from the UKB HapMap3 panel when the file carries none (ASD, ALZ)
  se    from |b| and p when the file carries none (ASD)
  b,se  from Z when the file carries only Z (ALZ), via the standard
        b = z/sqrt(2p(1-p)(n+z^2)), se = 1/sqrt(2p(1-p)(n+z^2))
  N     EFFECTIVE sample size, which both case/control files already supply:
          SCZ  2 x NEFFDIV2 (the PGC convention stores Neff/2)
          MDD  the "n" column, which is non-integer and so already Neff
        ASD and ALZ supply no case/control split, only a total, so those two
        remain on raw N.  That over-states their precision, so the Bayesian
        methods shrink them less than they should and their scores are noisier
        than a like-for-like Neff version would be.  Both are negative controls
        and this works against them, so it is recorded rather than hidden --
        but it is a reason not to read their nulls as strongly as SCZ's signal.
"""
import csv, gzip, math, os, sys
from statistics import NormalDist

IN, OUT, MAFFILE = sys.argv[1], sys.argv[2], sys.argv[3]
# optional 4th arg: comma-separated subset of traits to emit (default all)
ONLY = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
os.makedirs(OUT, exist_ok=True)
MAF = {}
with open(MAFFILE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        MAF[r["SNP"]] = float(r["MAF"])
print(f"MAF panel: {len(MAF):,} SNPs", flush=True)

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def emit(name, gen):
    if ONLY is not None and name not in ONLY: return
    p = os.path.join(OUT, f"{name}.ma")
    n = 0
    with open(p, "w") as o:
        o.write("SNP A1 A2 freq b se p N\n")
        for snp, a1, a2, fr, b, se, pv, N in gen:
            o.write(f"{snp} {a1} {a2} {fr:.6f} {b:.6g} {se:.6g} {pv:.6g} {int(N)}\n")
            n += 1
    print(f"  {name}.ma: {n:,} SNPs", flush=True)

def pgc_vcf(path, has_n=True):
    """PGC sumstats-VCF style.  Both files carry NCAS/NCON and FCAS/FCON; they
    differ only in naming the effective size NEFF (european) or NEFFDIV2
    (primary).  Either way it is Neff/2 by PGC convention, so N = 2x it."""
    with opener(path) as f:
        for line in f:
            if not line.startswith("##"):
                hdr = line.rstrip("\n").split("\t"); break
        c = {k: i for i, k in enumerate(hdr)}
        for line in f:
            v = line.rstrip("\n").split("\t")
            try:
                a1, a2 = v[c["A1"]].upper(), v[c["A2"]].upper()
                b, se, pv = float(v[c["BETA"]]), float(v[c["SE"]]), float(v[c["PVAL"]])
                ncas, ncon = float(v[c["NCAS"]]), float(v[c["NCON"]])
                fr = (float(v[c["FCAS"]])*ncas + float(v[c["FCON"]])*ncon)/(ncas+ncon)
                nefcol = "NEFF" if "NEFF" in c else "NEFFDIV2"
                N = 2*float(v[c[nefcol]])
            except (ValueError, KeyError, ZeroDivisionError):
                continue
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or se <= 0:
                continue
            yield v[c["ID"]], a1, a2, fr, b, se, pv, N

def mdd(path):
    with opener(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                fr = float(r["effect_allele_frequency"]); b = float(r["beta"])
                se = float(r["standard_error"]); pv = float(r["p_value"])
                # The file's "n" is ALREADY an effective sample size, not a
                # headcount: it is non-integer (e.g. 1183604.88 where
                # ncases+ncontrols = 2288955).  It is the per-cohort sum of
                # 4/(1/ncas_i + 1/ncon_i), which is the right quantity and is
                # strictly SMALLER than the same formula applied to the pooled
                # totals (1.18M vs 1.36M here) because pooling hides each
                # cohort's case/control imbalance.  So use it as given.
                #
                # I briefly replaced this with the pooled-ratio formula on the
                # theory that MDD was on a raw-N convention while SCZ was on
                # Neff (2 x NEFFDIV2).  It is not -- both are Neff, and the
                # substitution made MDD's N too large, which would have made
                # PRS-CS/SBayesR shrink MDD too little.  Recorded because the
                # giveaway is cheap to check and easy to miss: an N column with
                # a fractional part is not a headcount.
                N = float(r.get("N") or r.get("n"))
            except (ValueError, TypeError):
                continue
            a1, a2 = r["effect_allele"].upper(), r["other_allele"].upper()
            snp = r["rsid"]
            if not snp or a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or se <= 0:
                continue
            yield snp, a1, a2, fr, b, se, pv, N

def asd(path):
    with opener(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            snp = r.get("rsID") or ""
            fr = MAF.get(snp)
            if not snp or fr is None: continue
            try:
                b = float(r["Effect"]); se = float(r["StdErr"])
                pv = float(r["P-value"]); N = float(r["TotalSampleSize"])
            except (ValueError, TypeError, KeyError):
                continue
            a1, a2 = r["Allele1"].upper(), r["Allele2"].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or se <= 0: continue
            yield snp, a1, a2, fr, b, se, pv, N

def alz(path):
    """PGC-ALZ2 carries Z only, and no rsID in the raw file -- so this reads the
    rsID-mapped copy prepare_alz.py already produced (SNP A1 A2 Z P N) and
    reconstructs b and se from Z using the UKB panel's MAF."""
    with opener(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            snp = r["SNP"]; fr = MAF.get(snp)
            if fr is None: continue
            try:
                z = float(r["EFFECT"]); pv = float(r["P"]); N = float(r["N"])
            except (ValueError, TypeError):
                continue
            if z != z or abs(z) == float("inf") or N <= 0: continue
            a1, a2 = r["A1"].upper(), r["A2"].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1): continue
            d = math.sqrt(2*fr*(1-fr)*(N + z*z))
            if d <= 0: continue
            yield snp, a1, a2, fr, z/d, 1.0/d, pv, N

def kunkle(path):
    """Kunkle et al. 2019, IGAP stage 1 (GWAS Catalog GCST007511): 21,982
    CLINICALLY DIAGNOSED late-onset AD cases vs 41,944 controls, European.

    The point of adding it is what it does NOT contain.  Wightman 2021 (our
    other ALZ file) includes UK Biobank by-proxy cases, where the phenotype is
    parental dementia reported by the participant -- contaminated by parental
    longevity, SES and education.  An AD-by-proxy score therefore partly indexes
    educational attainment, and EA associates with cortical structure in ABCD.
    Kunkle has no proxy cases, so if the ALZ_noAPOE association survives here it
    is not proxy contamination; if it vanishes, it is.

    No frequency column (borrow from the target, as for ASD/ALZ) and no N
    column: N is constant by design, Neff = 4/(1/21982 + 1/41944) = 57,706.
    That is a third of Wightman's Neff, so a null here is weaker evidence than
    a null in a comparably powered file -- stated so the comparison is not
    over-read."""
    NEFF = 4.0 / (1.0/21982 + 1.0/41944)
    with opener(path) as f:
        for r in csv.DictReader(f, delimiter=" ", skipinitialspace=True):
            snp = r.get("MarkerName") or ""
            fr = MAF.get(snp)
            if not snp or fr is None: continue
            try:
                b = float(r["Beta"]); se = float(r["SE"]); pv = float(r["Pvalue"])
            except (ValueError, TypeError, KeyError):
                continue
            a1, a2 = r["Effect_allele"].upper(), r["Non_Effect_allele"].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or se <= 0:
                continue
            yield snp, a1, a2, fr, b, se, pv, NEFF

def okbay(path, cols):
    """Okbay et al. 2016 educational attainment, years of education
    (GCST003676), N = 405,072 European.  Not a disorder control -- a POSITIVE
    control for the confound itself.  If polygenic EA associates with
    global_slope at the magnitudes we are reporting for SCZ, then the whole set
    of disorder associations has to be read as possibly an SES/education signal,
    and conditioning on it becomes the test that matters.

    Lee 2018 EA3 (N up to 1.1M) would be far better powered, but its full
    summary statistics are 23andMe-restricted behind an SSGAC data-use
    agreement, which is the user's to sign rather than mine to click through.
    Column names vary between Okbay releases, so they are passed in."""
    c_snp, c_a1, c_a2, c_frq, c_b, c_se, c_p = cols
    with opener(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            snp = r.get(c_snp) or ""
            if not snp: continue
            try:
                b = float(r[c_b]); se = float(r[c_se]); pv = float(r[c_p])
            except (ValueError, TypeError, KeyError):
                continue
            fr = None
            if c_frq and r.get(c_frq) not in (None, "", "NA"):
                try: fr = float(r[c_frq])
                except ValueError: fr = None
            if fr is None: fr = MAF.get(snp)
            if fr is None: continue
            a1, a2 = r[c_a1].upper(), r[c_a2].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or se <= 0:
                continue
            yield snp, a1, a2, fr, b, se, pv, 405072.0

emit("SCZ_pooled", pgc_vcf(os.path.join(IN, "SCZ_primary.tsv")))
emit("SCZ_eur",    pgc_vcf(os.path.join(IN, "SCZ_european.tsv.gz")))
emit("MDD_pooled", mdd(os.path.join(IN, "MDD_div_raw.tsv")))
emit("MDD_eur",    mdd(os.path.join(IN, "MDD_eur_raw.tsv.gz")))
emit("ASD",        asd(os.path.join(IN, "ASD_raw.tsv")))
emit("ALZ",        alz(os.path.join(IN, "ALZ_rsid.tsv")))
emit("ALZ_IGAP",   kunkle(os.path.join(IN, "Kunkle_etal_Stage1_results.txt")))
# Okbay column names are supplied by the caller via OKBAY_COLS so a different
# release can be swapped in without editing this file.
_oc = os.environ.get("OKBAY_COLS", "MarkerName,A1,A2,EAF,Beta,SE,Pval").split(",")
emit("EA",         okbay(os.path.join(IN, "Okbay_EduYears_Main.txt.gz"), _oc))
