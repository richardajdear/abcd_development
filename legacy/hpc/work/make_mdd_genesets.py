"""Build MAGMA gene sets from the MDD high-confidence gene table.

Source: Table S21 of the 2024 Cell depression GWAS (doi 10.1016/j.cell.2024.12.002,
PMC11829167), supplementary file `NIHMS2048064-supplement-11.xlsx`, obtained via
the Europe PMC supplementaryFiles API -- cell.com and pmc.ncbi.nlm.nih.gov both
serve a bot challenge instead of the file.

Two sheets, mirroring the SCZ table's structure (see make_prioritised_genesets.py):

  `High-confidence Gene List`      308 genes -- the prioritised set, defined by
                                   the readme as those "identified from
                                   finemapping, expression, or protein".
  `Table S21 Gene Mapping Methods` 4,600 genes -- every gene flagged by at least
                                   one of seven mapping methods.

An asymmetry with the SCZ analysis worth stating rather than burying
--------------------------------------------------------------------
The SCZ pool is 685 genes at 239 genome-wide-significant loci: a candidate list
built by asking "which genes sit at these loci".  The MDD pool is 4,600 genes,
roughly a quarter of all genes MAGMA tests, because two of its seven methods
(fastBAT, H-MAGMA) are themselves genome-wide gene-based association tests
rather than locus annotations.  A competitive test against a set that large is
a weak contrast -- it approaches asking whether the signal is in a quarter of
the genome versus the other three quarters -- and it is not the same question
the SCZ pool answers.  The 308-gene high-confidence set is the like-for-like
comparator for SCZ's 120 prioritised genes; the pools are not comparable.

The three defining criteria are also kept as separate sets, since fine-mapping,
expression and protein evidence are different kinds of claim.

    python make_mdd_genesets.py <supplement-11.xlsx> <gene.loc> <out.txt>

Writes a two-column MAGMA set-annot file (geneID, set name) for use as
`--set-annot <out.txt> col=1,2`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HC_SHEET = "High-confidence Gene List"
POOL_SHEET = "Table S21 Gene Mapping Methods"

CRITERIA = {
    "MDD_hc_finemap": "Fine_mapping",
    "MDD_hc_expression": "Expression",
    "MDD_hc_protein": "Protein",
}


def main(xlsx: str, gene_loc: str, out: str) -> int:
    sym2id: dict[str, str] = {}
    for line in Path(gene_loc).read_text().splitlines():
        f = line.split()
        if len(f) >= 6:
            sym2id.setdefault(f[5], f[0])
    print(f"{gene_loc}: {len(sym2id)} symbols")

    x = pd.ExcelFile(xlsx)
    hc = x.parse(HC_SHEET)
    pool = x.parse(POOL_SHEET)
    print(f"{HC_SHEET}: {len(hc)} genes")
    print(f"{POOL_SHEET}: {len(pool)} genes")

    rows: list[tuple[str, str]] = []

    def add(name: str, symbols) -> None:
        symbols = [s for s in symbols if isinstance(s, str)]
        ids = sorted({sym2id[s] for s in symbols if s in sym2id})
        missing = [s for s in symbols if s not in sym2id]
        rows.extend((i, name) for i in ids)
        print(f"  {name:<26} {len(ids):>4} mapped of {len(symbols):>4}"
              f"   ({len(missing)} unmapped)")

    add("MDD_highconf", hc["Gene"])
    for name, col in CRITERIA.items():
        if col in hc.columns:
            add(name, hc.loc[hc[col].astype(str).str.lower() == "true", "Gene"])

    add("MDD_pool", pool["Gene"])
    hc_syms = set(hc["Gene"].dropna())
    add("MDD_pool_not_hc",
        [s for s in pool["Gene"] if isinstance(s, str) and s not in hc_syms])

    with open(out, "w") as fh:
        for gid, name in rows:
            fh.write(f"{gid}\t{name}\n")
    n_sets = len({n for _, n in rows})
    print(f"\nwrote {out}: {len(rows)} memberships across {n_sets} sets")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(*sys.argv[1:4]))
