"""Put all six discovery GWAS into one format: SNP A1 A2 freq b se p N.

One normaliser for every trait and arm, so every downstream method reads the
same columns and any difference between methods is the method.

  SCZ_pooled  PGC3_SCZ_wave3.primary   EUR+EAS+AFR+LAT (multi-ancestry)
  SCZ_eur     PGC3_SCZ_wave3.european  European only
  MDD_pooled  pgc-mdd2025 ..._div      trans-ancestry
  MDD_eur     pgc-mdd2025 ..._eur      European only
  ASD         SPARK+iPSYCH+PGC         no ancestry-stratified release exists
  ALZ         PGC-ALZ2 (Wightman)      no ancestry-stratified release exists

Missing pieces, filled the same way for every trait that needs it:
  freq  from the UKB HapMap3 panel when the file carries none (ASD, ALZ)
  se    from |b| and p when the file carries none (ASD)
  b,se  from Z when the file carries only Z (ALZ), via the standard
        b = z/sqrt(2p(1-p)(n+z^2)), se = 1/sqrt(2p(1-p)(n+z^2))
  N     2 x NEFF for the PGC case/control files (their NEFF is Neff/2)
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

emit("SCZ_pooled", pgc_vcf(os.path.join(IN, "SCZ_primary.tsv")))
emit("SCZ_eur",    pgc_vcf(os.path.join(IN, "SCZ_european.tsv.gz")))
emit("MDD_pooled", mdd(os.path.join(IN, "MDD_div_raw.tsv")))
emit("MDD_eur",    mdd(os.path.join(IN, "MDD_eur_raw.tsv.gz")))
emit("ASD",        asd(os.path.join(IN, "ASD_raw.tsv")))
emit("ALZ",        alz(os.path.join(IN, "ALZ_rsid.tsv")))
