"""Synthetic children with a planted C4A effect, for testing 04_c4_assoc.R.

  python make_assoc_fixture.py --panel ../resources/MHC_...vcf.gz --out out/assoc_fixture --n 3000 --beta -0.10

Each child gets two C4 alleles drawn from the panel's haplotypes, so C4A/C4B
GREx have realistic distributions and correlation.  Files mirror the layout of
genetic_analysis/work/results_70tab_hcp/prs_final_1lmm/pheno/ (space-delimited,
FID = family id, IID in the sub-xxxxxxxx spelling) plus a c4_calls.tsv in the
03_c4_grex.py format.  Planted: global_slope_1lmm = beta * z(C4A_GREx) + family
effect + noise; baseline_thickness_1lmm and the PRS are independent of C4.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import re
import string
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("c4grex", Path(__file__).resolve().parents[1] / "03_c4_grex.py")
c4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c4)  # noqa: E702


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--beta", type=float, default=-0.10)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out); (out / "pheno").mkdir(parents=True, exist_ok=True)

    with gzip.open(a.panel, "rt") as fh:
        row = next(l for l in fh if not l.startswith("#") and l.split("\t")[2] == "C4").rstrip("\n").split("\t")
    alleles = ["REF"] + [x.strip("<>") for x in row[4].split(",")]
    haps = [alleles[int(h)] for g in row[9:] for h in re.split(r"[|/]", g)]

    n = a.n
    alph = np.array(list(string.ascii_uppercase + string.digits))
    toks = ["".join(rng.choice(alph, 8)) for _ in range(n)]
    fam = np.arange(n); pair = rng.random(n) < 0.2                   # ~20 % in sibling pairs
    fam[1:][pair[1:]] = fam[:-1][pair[1:]]
    rows = []
    for t in toks:
        h1, h2 = rng.choice(haps, 2)
        e = {s: c4.seg_counts(h1)[s] + c4.seg_counts(h2)[s] for s in c4.SEGMENTS}
        r = dict(IID=f"sub-{t}", tok=t, hap1=c4.structure(h1), hap2=c4.structure(h2),
                 post_mean=rng.uniform(0.6, 1.0), **{f"E_{s}": v for s, v in e.items()})
        r["C4A_copies"] = e["AL"] + e["AS"]; r["C4B_copies"] = e["BL"] + e["BS"]; r["HERV_copies"] = e["AL"] + e["BL"]
        for g, w in c4.GREX_WEIGHTS.items():
            r[g] = sum(wt * e[s] for s, wt in w.items())
        r["common5"] = int(r["hap1"] in c4.COMMON5 and r["hap2"] in c4.COMMON5)
        rows.append(r)
    calls = pd.DataFrame(rows)
    calls.to_csv(out / "c4_calls.tsv", sep="\t", index=False)

    z = (calls.C4A_GREx - calls.C4A_GREx.mean()) / calls.C4A_GREx.std()
    fe = pd.Series(rng.normal(0, 0.4, n)).groupby(fam).transform("first").to_numpy()
    fid = [f"F{f:05d}" for f in fam]
    ph = pd.DataFrame(dict(FID=fid, IID=calls.IID, global_slope_1lmm=a.beta * z + fe + rng.normal(0, 1, n),
                           baseline_thickness_1lmm=fe + rng.normal(0, 1, n)))
    ph.to_csv(out / "pheno/phenotypes_gcta.txt", sep=" ", index=False)
    q = pd.DataFrame(dict(FID=fid, IID=calls.IID, baseline_age=rng.uniform(9, 11, n),
                          **{f"PC{k}": rng.normal(0, 1, n) for k in range(1, 11)}))
    q.to_csv(out / "pheno/covar_quant.txt", sep=" ", index=False)
    pd.DataFrame(dict(FID=fid, IID=calls.IID, sex=rng.integers(1, 3, n), site=rng.integers(1, 22, n))
                 ).to_csv(out / "pheno/covar_categorical.txt", sep=" ", index=False)
    pd.DataFrame(dict(FID=fid, IID=calls.IID, SCORE=rng.normal(0, 1, n))).to_csv(out / "prs.tsv", sep="\t", index=False)
    # a PLINK-style score directory (for --prs-dir) and a fake weight table + chr6 .bim (for 05_mhc_in_scores.py)
    sd = out / "scores/SBayesRC/SCZ25_EUR"; sd.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dict(FID=fid, IID=calls.IID, CNT=1000, SCORE=pd.read_csv(out / "prs.tsv", sep="\t").SCORE)
                 ).to_csv(sd / "score_SCZ25_EUR.profile", sep=" ", index=False)
    bp = np.sort(rng.integers(20_000_000, 40_000_000, 500))
    snp = [f"rs{k}" for k in range(500)]
    pd.DataFrame({0: 6, 1: snp, 2: 0, 3: bp, 4: "A", 5: "G"}).to_csv(out / "chr6.bim", sep="\t", header=False, index=False)
    pd.DataFrame(dict(SNP=snp[::2], A1="A", BETA=0.001)).to_csv(sd / "weights.snpRes", sep="\t", index=False)
    (out / "expected_mhc.txt").write_text(str(int(((bp[::2] >= 25_000_000) & (bp[::2] <= 34_000_000)).sum())))
    keep = calls.IID[rng.random(n) < 0.9]
    pd.DataFrame(dict(FID=keep, IID=keep)).to_csv(out / "eur.keep", sep=" ", index=False, header=False)
    print(f"fixture: {n} children, {len(set(fam))} families, planted beta {a.beta}, EUR {len(keep)}")


if __name__ == "__main__":
    main()
