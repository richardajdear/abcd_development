"""Step 5 (CSD3): how much of the extended MHC do our SCZ polygenic scores carry?

  python 05_mhc_in_scores.py --score-root genetic_analysis/work/scores_scz2025 \
      --bim <imputed chr6 .bim> --out results/mhc_in_scz_scores.tsv

M4 of 04_c4_assoc.R asks whether C4A is independent of the SCZ score.  That
reading depends on whether the score already contains MHC variants (C+T
clumping keeps them unless excluded; the Bayesian methods may or may not,
depending on their LD reference).  This script finds every per-SNP weight
table under <score-root>/<METHOD>/SCZ25_EUR/ (any delimited text file with a
SNP/rsID column, excluding .profile/.log), and counts the SNPs that fall in
the extended MHC (GRCh37 chr6:25,000,000-34,000,000) by rsID against the
chr6 .bim.  Output is summary-level (counts only).
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd

MHC = (25_000_000, 34_000_000)
ID_COLS = ("SNP", "snp", "rsid", "RSID", "rsID", "ID", "Name", "SNPID", "variant")
SKIP = (".profile", ".log", ".nosex", ".sscore", ".tbi", ".bed", ".bim", ".fam")


def open_any(p: Path):
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-root", required=True)
    ap.add_argument("--bim", required=True)
    ap.add_argument("--arm", default="SCZ25_EUR")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    bim = pd.read_csv(a.bim, sep=r"\s+", header=None, usecols=[1, 3], names=["SNP", "BP"])
    mhc_ids = set(bim.SNP[(bim.BP >= MHC[0]) & (bim.BP <= MHC[1])])
    rows = []
    for mdir in sorted(Path(a.score_root).glob(f"*/{a.arm}")):
        files = [f for f in mdir.rglob("*") if f.is_file() and not f.name.endswith(SKIP) and f.stat().st_size > 0]
        found = False
        for f in files:
            try:
                with open_any(f) as fh:
                    head = fh.readline()
                sep = "\t" if "\t" in head else (r"\s+" if " " in head.strip() else ",")
                cols = pd.read_csv(f, sep=sep, nrows=0).columns
                idc = next((c for c in ID_COLS if c in cols), None)
                if idc is None:
                    continue
                ids = pd.read_csv(f, sep=sep, usecols=[idc])[idc].astype(str)
            except Exception:
                continue
            found = True
            rows.append(dict(method=mdir.parent.name, arm=a.arm, file=str(f.relative_to(mdir)),
                             n_snps=len(ids), n_mhc=int(ids.isin(mhc_ids).sum())))
        if not found:
            rows.append(dict(method=mdir.parent.name, arm=a.arm, file="(no weight table found)",
                             n_snps=None, n_mhc=None))
    out = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
