"""Disorder gene sets defined by the discovery GWAS's OWN gene-level results.

Usage:  python build_topgene_sets.py <genes_dir> <out.txt>
  genes_dir holds MDD_EUR.genes.out and SCZ25_EUR.genes.out (LD-matched gene
  analyses of the European discovery GWAS; CSD3 work/magma_scz2025/genes/).
Writes MAGMA long-format set annotation (GENE <tab> SET; use col=1,2):
  <D>_genesig   Bonferroni-significant genes (p < 0.05 / n genes)
  <D>_top{100,250,500}   the N most significant genes
for D in MDD, SCZ25.  The MHC (chr6 25-34 Mb) is excluded from every set so
the SCZ sets are not dominated by one LD block; the same rule applies to MDD
so the two disorders' sets are built identically and sizes are comparable.
"""
import csv
import sys
from pathlib import Path

SOURCES = {"MDD": "MDD_EUR.genes.out", "SCZ25": "SCZ25_EUR.genes.out"}
TOPN = (100, 250, 500)


def genes(path):
    rows = list(csv.DictReader(open(path), delimiter=" ", skipinitialspace=True))
    rows = [r for r in rows
            if not (r["CHR"] == "6" and 25e6 < float(r["START"]) < 34e6)]
    return sorted(rows, key=lambda r: float(r["P"])), len(rows)


def main(gdir, out):
    lines = []
    for tag, f in SOURCES.items():
        g, n = genes(Path(gdir) / f)
        sig = [r["GENE"] for r in g if float(r["P"]) < 0.05 / n]
        lines += [f"{x}\t{tag}_genesig" for x in sig]
        for k in TOPN:
            lines += [f"{r['GENE']}\t{tag}_top{k}" for r in g[:k]]
        print(f"{tag}: {n} genes (MHC excluded), {len(sig)} Bonferroni")
    Path(out).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
