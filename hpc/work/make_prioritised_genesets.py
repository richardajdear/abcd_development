"""Build MAGMA gene sets from the SCZ prioritised-gene table (Trubetskoy 2022).

Source: Extended Data Table 1 of Trubetskoy et al. 2022, Nature 604:502
("Mapping genomic loci implicates genes and synaptic biology in
schizophrenia"), figshare 19426775 file 35775617, md5
8fee2faee10ddf2c9f8a16a54f2c5b84.  Sheet `Extended.Data.Table.1` is the 120
prioritised genes; sheet `ST12 all criteria` is the 685-gene candidate pool
those 120 were selected from, with a `Prioritised` flag.

Why this test is worth running when the genome-wide rg was null
---------------------------------------------------------------
`05_ldsc` asks whether SCZ and a cortical phenotype share genetic architecture
*averaged over the whole genome*, and answers no (rg = 0.036 +/- 0.047).  A
competitive gene-set test asks something much narrower: within the ABCD gene
results, is the association signal *concentrated* in these 120 genes relative to
all other genes?  That is one degree of freedom over a set chosen by orthogonal
evidence (fine-mapping, SMR, rare variants), so it is far better powered at
N = 4,119 than either a genome-wide rg or a GWAS scan -- and a null rg does not
predict a null here, because rg averages over the ~18,000 genes that carry no
SCZ signal at all.

ID mapping
----------
MAGMA's gene results are keyed by Entrez ID (from NCBI37.3.gene.loc); the table
supplies Ensembl IDs and symbols.  We map on symbol, which is what gene.loc
carries.  Genes that fail to map are reported rather than silently dropped --
mostly lincRNAs and antisense transcripts (`AC068490.2`-style names) that have
no Entrez entry in the NCBI37.3 build, so the mapped set is biased toward
protein-coding genes.  That is a real limitation of the annotation, not of the
selection, and it is stated in the output.

    python make_prioritised_genesets.py <table.xlsx> <gene.loc> <out.txt>

Writes a two-column MAGMA set-annot file (geneID, set name) for use as
`--set-annot <out.txt> col=1,2`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# Sets built from Extended Data Table 1.  The three priority flags are not
# mutually exclusive and are kept separate as well as pooled, because they rest
# on different evidence: FINEMAP is statistical fine-mapping of the GWAS itself,
# SMR is eQTL-based causal inference, and Rare is exome/CNV burden -- so a
# signal appearing in only one of them means something different from a signal
# in all three.
SUBSETS = {
    "SCZ_prio_finemap": "FINEMAP.priority.gene",
    "SCZ_prio_smr": "SMR.priority.gene",
    "SCZ_prio_rare": "Rare.priority.gene",
}

# Symbols the 2022 table uses that NCBI37.3 spells differently.  Checked one by
# one against gene.loc rather than pulled from an alias database: of the 17
# unmapped prioritised genes this is the only one that is present under another
# name.  The rest are lincRNAs and clone-based names (`RP11-...`, `LINC01088`)
# that the NCBI37.3 build simply does not contain -- it has 3 `LINC*` entries in
# total -- plus MLXIP, which has no NCBI37.3 entry under any spelling.  Genes
# absent from gene.loc were never annotated, so MAGMA has no result for them and
# they could not enter a gene-set test however they were named.
ALIASES = {"GPR98": "ADGRV1"}


def main(xlsx: str, gene_loc: str, out: str) -> int:
    # symbol -> Entrez, from the same gene.loc MAGMA annotated the GWAS with,
    # so a gene present here is a gene MAGMA could have tested.
    sym2id: dict[str, str] = {}
    for line in Path(gene_loc).read_text().splitlines():
        f = line.split()
        if len(f) >= 6:
            sym2id.setdefault(f[5], f[0])
    print(f"{gene_loc}: {len(sym2id)} symbols")

    x = pd.ExcelFile(xlsx)
    t1 = x.parse("Extended.Data.Table.1")
    print(f"Extended.Data.Table.1: {len(t1)} prioritised genes")

    rows: list[tuple[str, str]] = []

    def add(name: str, symbols) -> None:
        symbols = [s for s in symbols if isinstance(s, str)]
        ids, missing = [], []
        for s in symbols:
            key = ALIASES.get(s, s)
            if key in sym2id:
                ids.append(sym2id[key])
            else:
                missing.append(s)
        ids = sorted(set(ids))
        rows.extend((i, name) for i in ids)
        print(f"  {name:<24} {len(ids):>4} mapped of {len(symbols):>4}"
              f"   ({len(missing)} unmapped)")
        if missing:
            print(f"      unmapped: {', '.join(sorted(missing)[:12])}"
                  f"{' ...' if len(missing) > 12 else ''}")

    add("SCZ_prioritised", t1["Symbol.ID"])
    pc = t1.loc[t1["gene_biotype"] == "protein_coding", "Symbol.ID"]
    add("SCZ_prioritised_pc", pc)
    for name, col in SUBSETS.items():
        if col in t1.columns:
            add(name, t1.loc[t1[col] == 1, "Symbol.ID"])

    # The 685-gene pool the 120 were chosen from: a same-loci control.  If the
    # ABCD signal is enriched in the pool just as much as in the prioritised
    # subset, the enrichment is about being near a SCZ locus rather than about
    # the gene prioritisation adding anything.
    if "ST12 all criteria" in x.sheet_names:
        st12 = x.parse("ST12 all criteria")
        add("SCZ_locus_pool", st12["Symbol.ID"])
        # The pool CONTAINS the prioritised genes, so testing both and comparing
        # p-values compares overlapping sets.  The disjoint remainder is what
        # makes the comparison answer a question: if the remainder is enriched
        # and the prioritised subset is not, the signal lives at SCZ loci
        # generally rather than in the genes prioritisation nominated as causal.
        prio = set(t1["Symbol.ID"].dropna())
        add("SCZ_pool_not_prio",
            [s for s in st12["Symbol.ID"] if isinstance(s, str) and s not in prio])

    with open(out, "w") as fh:
        for gid, name in rows:
            fh.write(f"{gid}\t{name}\n")
    n_sets = len({n for _, n in rows})
    print(f"\nwrote {out}: {len(rows)} gene-set memberships across {n_sets} sets")
    print("use as: magma --gene-results <x>.genes.raw "
          f"--set-annot {out} col=1,2")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(*sys.argv[1:4]))
