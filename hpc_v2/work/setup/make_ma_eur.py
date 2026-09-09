"""Build GCTB/COJO .ma files for SBayesR, from EUROPEAN-ONLY discovery GWAS.

.ma format: SNP A1 A2 freq b se p N   (A1 = effect allele)

ANCESTRY AUDIT -- what "EUR-only" means for each trait, checked not assumed:

  SCZ  PGC3_SCZ_wave3.EUROPEAN (53,386 cases / 77,258 controls).  NOT `core`,
       which the release README defines as "cohorts primarily of east asian and
       european ancestry" -- core carries ~1,700 extra non-European cases and
       was used by mistake in the earlier matched-rg and PRS-CS arms.
  MDD  pgc-mdd2025_no23andMe_EUR (412,305 / 1,588,397), 75 EUR cohorts.
  ASD  SPARK+iPSYCH+PGC.  NO ancestry-stratified release exists for this file;
       it is treated as European-ancestry because its constituent cohorts are,
       but that is an ASSUMPTION and is flagged in the output.
  ALZ  PGC-ALZ2 (Wightman 2021, excl. 23andMe).  Same caveat -- no per-ancestry
       release; predominantly but not exclusively European.

MISSING COLUMNS, and how they are filled:
  ASD and ALZ carry no allele frequency, so MAF comes from the 1000G EUR
  HapMap3 panel (the same panel PRS-CS used).  That restricts them to HapMap3,
  which is what SBayesR's LD reference covers anyway.
  ALZ also carries no beta/se, only Z.  Both are reconstructed with the
  standard approximation (Zhu et al. 2016), which is exactly what GCTA-COJO
  and SBayesR expect when only a Z is available:
      b  = z / sqrt(2p(1-p)(n + z^2))
      se = 1 / sqrt(2p(1-p)(n + z^2))
"""
import csv, gzip, sys, math
from statistics import NormalDist

maf_path, out_dir = sys.argv[1], sys.argv[2]
MAF = {}
with open(maf_path) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        MAF[row["SNP"]] = float(row["MAF"])
print(f"MAF reference: {len(MAF):,} HapMap3 SNPs", flush=True)

def ok(*vals):
    for v in vals:
        if v is None or v == "" or v != v:
            return False
    return True

def write(name, rows):
    p = f"{out_dir}/{name}.ma"
    n = 0
    with open(p, "w") as o:
        o.write("SNP A1 A2 freq b se p N\n")
        for snp, a1, a2, fr, b, se, pv, N in rows:
            o.write(f"{snp} {a1} {a2} {fr:.6f} {b:.6g} {se:.6g} {pv:.6g} {int(N)}\n")
            n += 1
    print(f"  {name}.ma: {n:,} SNPs", flush=True)

# ---- SCZ: European-ancestry meta-analysis -----------------------------------
def scz():
    src = "/rds/user/rajd2/hpc-work/magma/gwas/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
    with gzip.open(src, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                hdr = line.rstrip("\n").split("\t"); break
        c = {k: i for i, k in enumerate(hdr)}
        for line in f:
            v = line.rstrip("\n").split("\t")
            try:
                ncas, ncon = float(v[c["NCAS"]]), float(v[c["NCON"]])
                # allele frequency across cases and controls, weighted by n
                fr = (float(v[c["FCAS"]])*ncas + float(v[c["FCON"]])*ncon)/(ncas+ncon)
                b, se, pv = float(v[c["BETA"]]), float(v[c["SE"]]), float(v[c["PVAL"]])
                # NEFF here is Neff/2 (PGC convention), so N = 2*NEFF
                N = 2*float(v[c["NEFF"]])
            except (ValueError, KeyError, ZeroDivisionError):
                continue
            a1, a2 = v[c["A1"]].upper(), v[c["A2"]].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or se <= 0:
                continue
            yield v[c["ID"]], a1, a2, fr, b, se, pv, N

# ---- MDD: European-ancestry release -----------------------------------------
def mdd():
    src = "hpc_v2/work/results_v2/mdd_eur/MDD_EUR.tsv"
    with open(src) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                fr = float(row["effect_allele_frequency"]); b = float(row["beta"])
                se = float(row["standard_error"]); pv = float(row["p_value"]); N = float(row["N"])
            except (ValueError, TypeError):
                continue
            a1, a2 = row["effect_allele"].upper(), row["other_allele"].upper()
            snp = row["rsid"]
            if not snp or a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or se <= 0:
                continue
            yield snp, a1, a2, fr, b, se, pv, N

# ---- ASD: beta/se present, frequency taken from the 1000G EUR panel ---------
def asd():
    with open("hpc_v2/work/results_v2/control_ASD/ASD.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            snp = row["SNP"]; fr = MAF.get(snp)
            if fr is None: continue
            try:
                b = float(row["EFFECT"]); pv = float(row["P"]); N = float(row["N"])
            except (ValueError, TypeError):
                continue
            a1, a2 = row["A1"].upper(), row["A2"].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1): continue
            # SE is not carried through ASD.tsv; recover it from beta and p.
            if pv <= 0 or pv >= 1 or b == 0: continue
            z = NormalDist().inv_cdf(1 - pv/2)
            if z <= 0: continue
            se = abs(b)/z
            if se <= 0: continue
            yield snp, a1, a2, fr, b, se, pv, N

# ---- ALZ: Z only, so beta and se are both reconstructed ---------------------
def alz():
    with open("hpc_v2/work/results_v2/control_ALZ/ALZ.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            snp = row["SNP"]; fr = MAF.get(snp)
            if fr is None: continue
            try:
                z = float(row["EFFECT"]); pv = float(row["P"]); N = float(row["N"])
            except (ValueError, TypeError):
                continue
            if z != z or abs(z) == float("inf"): continue
            a1, a2 = row["A1"].upper(), row["A2"].upper()
            if a1 not in "ACGT" or a2 not in "ACGT" or not (0 < fr < 1) or N <= 0: continue
            denom = math.sqrt(2*fr*(1-fr)*(N + z*z))
            if denom <= 0: continue
            yield snp, a1, a2, fr, z/denom, 1.0/denom, pv, N

for name, gen in (("SCZ", scz), ("MDD", mdd), ("ASD", asd), ("ALZ", alz)):
    print(f"building {name}...", flush=True)
    write(name, gen())
