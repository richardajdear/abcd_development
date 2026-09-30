#!/usr/bin/env python
"""Step 16.4: genome-wide MOSTest and minP for one measure family.

Reads zmat/<family>_chr*.npz (real) and zmat/<family>_perm_chr*.npz (one
genotype permutation), and writes

  sumstats/mostest_<family>.sumstats.tsv.gz       real scan, per SNP
      SNP CHR POS A1 A2 FREQ N STAT P P_CHI2 MINP P_MINP    (P = MOSTest, gamma null)
  sumstats/mostest_<family>_perm.sumstats.tsv.gz  permuted scan, same columns
  <family>_mostest.json                            fits, diagnostics, gate inputs
  qq_<family>.tsv                                  binned QQ points (real, perm; MOSTest and minP)
  univariate_<family>.tsv                          per-measure lambda_GC and hit counts

The per-SNP files are gitignored (SNP-level, large); 06_collect.py turns the
JSON / TSV files into the tracked table_mostest_*.tsv.

Method (mostest_core.py): R = corr(permuted z); T = z' R_reg^-1 z; p from a
gamma fitted to the permuted T.  Also reported: the analytic chi2(k) p (should
agree when no eigenvalue was floored) and a second estimate of R from the REAL
z (should agree with the permuted R to within polygenic inflation; gate G3).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))   # run from any cwd
import mostest_core as mc  # noqa: E402


def load(zm: Path, run: str):
    parts = []
    for c in range(1, 23):
        f = zm / f"{run}_chr{c}.npz"
        if not f.exists():
            raise SystemExit(f"FATAL: missing {f}")
        d = np.load(f)
        parts.append(d)
    traits = list(parts[0]["traits"])
    meta = pd.concat([pd.DataFrame({k: d[k] for k in ("snp", "chrom", "pos", "a0", "a1", "freq", "n")})
                      for d in parts], ignore_index=True)
    z = np.concatenate([d["z"] for d in parts]).astype(np.float64)
    return traits, meta, z


def lam_from_z(z):
    return float(np.median(z ** 2) / stats.chi2.ppf(0.5, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", required=True)
    ap.add_argument("--root", required=True, help="step-16 output root (M16)")
    ap.add_argument("--floor-frac", type=float, default=1e-4)
    a = ap.parse_args()
    root = Path(a.root); zm = root / "zmat"; ss = root / "sumstats"; ss.mkdir(exist_ok=True)
    fam = a.family

    tr, meta, z = load(zm, fam)
    tr_p, meta_p, zp = load(zm, f"{fam}_perm")
    if tr != tr_p or not np.array_equal(meta.snp.values, meta_p.snp.values):
        raise SystemExit("FATAL: real and permuted scans differ in traits or variants")
    k = len(tr)
    print(f"{fam}: {len(meta):,} SNPs x {k} measures")

    accp = mc.CorrAccumulator(k); accp.add(zp); R = accp.corr()
    accr = mc.CorrAccumulator(k); accr.add(z); R_real = accr.corr()
    Rinv, rinfo = mc.regularised_inverse(R, a.floor_frac)
    iu = np.triu_indices(k, 1)

    t_real = mc.most_stat(z, Rinv); t_perm = mc.most_stat(zp, Rinv)
    gfit = mc.fit_gamma_null(t_perm)
    m_real = mc.minp_stat(z); m_perm = mc.minp_stat(zp)
    mfit = mc.fit_minp_null(m_perm)

    res = {}
    for tag, t, m, md in (("real", t_real, m_real, meta), ("perm", t_perm, m_perm, meta_p)):
        p = mc.p_most(t, gfit); pc = mc.p_chi2(t, k); pm = mc.p_minp(m, mfit)
        out = pd.DataFrame(dict(SNP=md.snp, CHR=md.chrom, POS=md.pos, A1=md.a1, A2=md.a0,
                                FREQ=md.freq.round(4), N=md.n, STAT=t.round(4),
                                P=p, P_CHI2=pc, MINP=m, P_MINP=pm))
        name = f"mostest_{fam}{'' if tag == 'real' else '_perm'}.sumstats.tsv.gz"
        out.to_csv(ss / name, sep="\t", index=False, float_format="%.4g", compression="gzip")
        res[tag] = dict(lambda_mostest=mc.lambda_gc(p), lambda_chi2=mc.lambda_gc(pc),
                        lambda_minp=mc.lambda_gc(pm),
                        min_p_mostest=float(np.nanmin(p)), min_p_minp=float(np.nanmin(pm)),
                        **{f"mostest_{kk}": v for kk, v in mc.tail_counts(p).items()},
                        **{f"minp_{kk}": v for kk, v in mc.tail_counts(pm).items()})
        res[tag]["_p"], res[tag]["_pm"] = p, pm

    # QQ points
    rows = []
    for tag in ("real", "perm"):
        for stat, key in (("MOSTest", "_p"), ("minP", "_pm")):
            q = mc.qq_table(res[tag][key])
            rows.append(pd.DataFrame(dict(family=fam, scan=tag, statistic=stat,
                                          exp_log10p=q[:, 0].round(4), obs_log10p=q[:, 1].round(4))))
        res[tag].pop("_p"); res[tag].pop("_pm")
    pd.concat(rows).to_csv(root / f"qq_{fam}.tsv", sep="\t", index=False)

    # per-measure univariate summaries (the scans MOSTest combines)
    p_uni = 2 * stats.norm.sf(np.abs(z))
    pd.DataFrame(dict(family=fam, measure=tr,
                      lambda_gc=[lam_from_z(z[:, j]) for j in range(k)],
                      lambda_gc_perm=[lam_from_z(zp[:, j]) for j in range(k)],
                      n_p_lt_5e8=(p_uni < 5e-8).sum(axis=0),
                      min_p=p_uni.min(axis=0))).to_csv(root / f"univariate_{fam}.tsv",
                                                       sep="\t", index=False, float_format="%.4g")

    summ = dict(family=fam, k=k, n_snps=int(len(meta)), n_median=float(np.median(meta.n)),
                R=dict(**rinfo, median_offdiag=float(np.median(R[iu])),
                       max_abs_diff_real_vs_perm=float(np.abs(R - R_real)[iu].max()),
                       mean_abs_diff_real_vs_perm=float(np.abs(R - R_real)[iu].mean())),
                gamma_fit=gfit, minp_fit={kk: v for kk, v in mfit.items() if kk != "sorted_perm"},
                real=res["real"], perm=res["perm"])
    (root / f"{fam}_mostest.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps({kk: summ[kk] for kk in ("R", "gamma_fit", "minp_fit", "real", "perm")}, indent=1))


if __name__ == "__main__":
    main()
