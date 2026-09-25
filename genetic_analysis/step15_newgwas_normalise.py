#!/usr/bin/env python
"""Step 15: normalise the new discovery GWAS to the pipeline's .ma format
(SNP A1 A2 freq b se p N), one file per trait arm, as setup/normalise_gwas.py
did for the step-6 panel.

  BIP_eur     O'Connell et al. 2025, PGC bip2024_eur_no23andMe        EUR           daner, hg19
  BIP_pooled  O'Connell et al. 2025, PGC bip2024_multianc_no23andMe   EUR+AFR+EAS+LAT daner, hg19
  ADHD        Demontis et al. 2023, ADHD2022_iPSYCH_deCODE_PGC.meta   EUR           daner-like, hg19
  INT         Savage, Jansen et al. 2018, intelligence meta-analysis  EUR           METAL, GRCh37

Raw files: ~/rds/hpc-work/magma/gwas/{BIP_2024,ADHD_2023,INT_2018}/ (README there;
checksums verified 2026-09-25).  The ADHD and BIP data-use terms prohibit
reposting: nothing per-SNP from these files may enter the public repo.

Conventions, matched to the existing arms:
  * b = ln(OR) for case-control files; stdBeta for INT (continuous; the PRS
    methods only need b and se on a common scale within a file).
  * A1 = the allele b refers to (daner: A1; Savage: A1, lower-case in the file).
  * freq = discovery A1 frequency: BIP HRC_FRQ_A1 (EUR); BIP multi-ancestry
    HRC_EUR_FRQ_A1 -- the frequency consistent with the EUR LD references every
    Bayesian method uses (GCTB checks freq against the UKB EUR LD); ADHD the
    case/control-weighted FRQ_A/FRQ_U; INT EAF_HRC.
  * N = EFFECTIVE sample size per SNP:
      BIP   2 x Neff_half (the PGC convention, as for SCZ)
      ADHD  the file has only Nca/Nco.  The pooled-ratio formula 4/(1/Nca+1/Nco)
            OVERSTATES Neff when cohorts differ in case fraction (the note in
            normalise_gwas.py on MDD), so it is rescaled to the summary-statistic
            estimate Neff ~ 4 / (2 p (1-p) se^2) (median over well-behaved SNPs:
            0.05 <= freq <= 0.5 by MAF, INFO >= 0.9) and CALIBRATED by the ratio of
            BIP_eur's own Neff to the same estimator on BIP_eur (1.144 on
            2026-09-25; the estimator runs low on PGC case-control SEs).  Result:
            max per-SNP N 103,447, against the paper's neff_half = 51,568, i.e.
            Neff = 103,136 (Demontis 2023, Results para. 1): 0.3 % high.
            Intelligence checks the uncalibrated continuous form exactly.
      INT   N_analyzed (continuous trait: N is N).
  * Keep rsID-named biallelic SNPs with valid A/C/G/T alleles, finite b/se/p,
    se > 0, 0 < freq < 1; drop rsIDs that occur more than once.

    python genetic_analysis/step15_newgwas_normalise.py [ARM ...]
"""
from __future__ import annotations
import csv, gzip, math, os, statistics, sys

RAW = "/rds/user/rajd2/hpc-work/magma/gwas"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "genetic_analysis/work/scores_newgwas/gwas")
os.makedirs(OUT, exist_ok=True)
ACGT = set("ACGT")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def daner(path, freq_col=None, neff="half"):
    """yield (snp, a1, a2, freq, b, se, p, n_raw_pooled, neff_file_or_None, info)"""
    with opener(path) as f:
        hdr = None
        for line in f:
            if line.startswith("##"):
                continue
            hdr = line.split(); break
        c = {k: i for i, k in enumerate(hdr)}
        fa = [k for k in hdr if k.startswith("FRQ_A_")]
        fu = [k for k in hdr if k.startswith("FRQ_U_")]
        for line in f:
            v = line.split()
            try:
                a1, a2 = v[c["A1"]].upper(), v[c["A2"]].upper()
                OR, se, p = float(v[c["OR"]]), float(v[c["SE"]]), float(v[c["P"]])
                nca, nco = float(v[c["Nca"]]), float(v[c["Nco"]])
                if freq_col:
                    fr = float(v[c[freq_col]])
                else:
                    fr = (float(v[c[fa[0]]]) * nca + float(v[c[fu[0]]]) * nco) / (nca + nco)
                info = float(v[c["INFO"]]) if "INFO" in c else 1.0
                nf = 2 * float(v[c["Neff_half"]]) if neff == "half" else None
            except (ValueError, KeyError, IndexError, ZeroDivisionError):
                continue
            if OR <= 0:
                continue
            yield v[c["SNP"]], a1, a2, fr, math.log(OR), se, p, 4.0 / (1.0 / nca + 1.0 / nco), nf, info


def savage(path):
    with opener(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            try:
                fr, b, se, p, n = float(d["EAF_HRC"]), float(d["stdBeta"]), float(d["SE"]), float(d["P"]), float(d["N_analyzed"])
                info = float(d["minINFO"])
            except (ValueError, TypeError):
                continue
            yield d["SNP"], d["A1"].upper(), d["A2"].upper(), fr, b, se, p, n, n, info


def neff_estimate(rows, k=4.0):
    """median over well-behaved SNPs of 4 / (2 p (1-p) se^2), and the median of the
    file's pooled-ratio N over the same SNPs -> (estimate, pooled, file_neff_median)"""
    est, pool, fil = [], [], []
    for snp, a1, a2, fr, b, se, p, npool, nf, info in rows:
        maf = min(fr, 1 - fr)
        if 0.05 <= maf and info >= 0.9 and se > 0:
            est.append(k / (2 * fr * (1 - fr) * se * se)); pool.append(npool)
            if nf is not None:
                fil.append(nf)
    return statistics.median(est), statistics.median(pool), (statistics.median(fil) if fil else float("nan")), len(est)


ARMS = {
    "BIP_eur":    (os.path.join(RAW, "BIP_2024/bip2024_eur_no23andMe.gz"), lambda p: daner(p, "HRC_FRQ_A1", "half"), "file"),
    "BIP_pooled": (os.path.join(RAW, "BIP_2024/bip2024_multianc_no23andMe.gz"), lambda p: daner(p, "HRC_EUR_FRQ_A1", "half"), "file"),
    "ADHD":       (os.path.join(RAW, "ADHD_2023/ADHD2022_iPSYCH_deCODE_PGC.meta.gz"), lambda p: daner(p, None, None), "rescaled"),
    "INT":        (os.path.join(RAW, "INT_2018/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt"), savage, "file"),
}


def main():
    arms = sys.argv[1:] or list(ARMS)
    summ = os.path.join(OUT, "newgwas_summary.tsv")
    rows_out = []
    if os.path.exists(summ):
        with open(summ) as f:
            rows_out = [r for r in csv.DictReader(f, delimiter="\t") if r["arm"] not in arms]
    for arm in arms:
        path, reader, nmode = ARMS[arm]
        print(f"=== {arm}: {path}", flush=True)
        rows = [r for r in reader(path)
                if r[0].startswith("rs") and r[1] in ACGT and r[2] in ACGT and len(r[1]) == 1 and len(r[2]) == 1
                and 0 < r[3] < 1 and r[5] > 0 and 0 < r[6] <= 1 and math.isfinite(r[4])]
        n_read = len(rows)
        seen, dup = set(), set()
        for r in rows:
            (dup if r[0] in seen else seen).add(r[0])
        rows = [r for r in rows if r[0] not in dup]
        # case-control on the log-OR scale: k = 4; standardised continuous trait: k = 1
        est, pool, filem, n_est = neff_estimate(rows, 1.0 if arm == "INT" else 4.0)
        # The estimator runs ~13 % low for PGC case-control metas (SEs carry GC and
        # INFO inflation): for BIP_eur it gives 142,785 against the file's own
        # 2 x Neff_half = 163,367.  Calibrate the ADHD estimate by that ratio, read
        # from this run's BIP_eur row (so it is reproducible, not a typed constant).
        calib = 1.0
        if nmode == "rescaled":
            with open(summ) as f:
                b = [r for r in csv.DictReader(f, delimiter="\t") if r["arm"] == "BIP_eur"]
            if not b:
                sys.exit("ADHD calibration needs the BIP_eur row: normalise BIP_eur first")
            calib = float(b[0]["neff_file_median"]) / float(b[0]["neff_sumstat_estimate"])
        scale = est * calib / pool if nmode == "rescaled" else float("nan")
        with open(os.path.join(OUT, f"{arm}.ma"), "w") as o:
            o.write("SNP A1 A2 freq b se p N\n")
            nmax = 0
            for snp, a1, a2, fr, b, se, p, npool, nf, info in rows:
                N = nf if nmode == "file" else npool * scale
                nmax = max(nmax, N)
                o.write(f"{snp} {a1} {a2} {fr:.6f} {b:.6g} {se:.6g} {max(p, 1e-300):.6g} {int(round(N))}\n")
        with open(os.path.join(OUT, f"{arm}.n_gwas"), "w") as o:
            o.write(f"{int(round(nmax))}\n")
        with open(os.path.join(OUT, f"{arm}.ma")) as fi, open(os.path.join(OUT, f"{arm}.sst"), "w") as so:
            next(fi); so.write("SNP\tA1\tA2\tBETA\tP\n")
            for line in fi:
                v = line.split(); so.write(f"{v[0]}\t{v[1]}\t{v[2]}\t{v[4]}\t{v[6]}\n")
        n58 = sum(1 for r in rows if r[6] < 5e-8)
        rec = dict(arm=arm, source=os.path.basename(path), n_valid=n_read, n_dup_dropped=len(dup), n_kept=len(rows),
                   n_p5e8=n58, n_mode=nmode, neff_sumstat_estimate=round(est), n_pooled_ratio_median=round(pool),
                   neff_file_median=("" if filem != filem else round(filem)), rescale=("" if scale != scale else round(scale, 4)),
                   n_gwas_max=int(round(nmax)), n_snps_for_estimate=n_est,
                   calibration=(round(calib, 4) if nmode == "rescaled" else ""))
        print("  " + "  ".join(f"{k}={v}" for k, v in rec.items()), flush=True)
        rows_out.append(rec)
    keys = list(rows_out[-1].keys())
    with open(summ, "w") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t"); w.writeheader(); w.writerows(rows_out)
    print(f"wrote {summ}")


if __name__ == "__main__":
    main()
