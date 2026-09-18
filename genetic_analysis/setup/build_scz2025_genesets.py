#!/usr/bin/env python
"""MAGMA gene sets from the 2025 SCZ GWAS, with PGC3 analogues built the SAME way.

The legacy `SCZ_locus_pool` (462 genes) is Trubetskoy 2022's ST12 candidate
pool -- genes in GWS loci passing any of the paper's mapping criteria.  No such
table exists for the 2025 release, so two definitions that can be computed for
both GWAS from the files we have:

  <G>_locus_pool   genes whose annotated SNP window (35 kb up / 10 kb down,
                   union-panel annotation) contains at least one genome-wide
                   significant SNP (p < 5e-8) in the GWAS.  The natural
                   "genes under the GWS peaks" pool.
  <G>_genesig      genes Bonferroni-significant in the GWAS's own MAGMA gene
                   analysis (P < 0.05 / n genes tested).  A gene-level
                   definition that does not depend on a single SNP.

built for G = SCZ25 (the meta-analysed gene result and the meta .ma) and for
G = PGC3 (PGC3 european gene result; PGC3 primary .ma for the SNP pool, the
file step 7's SCZ_locus_pool test was read against).  The legacy sets
(SCZ_locus_pool, SCZ_prioritised, ...) are appended unchanged so every set is
tested in one MAGMA call and the conditional tests can name them.

Usage: build_scz2025_genesets.py <annot> <scz25_meta.ma> <SCZ25_META.genes.out>
                                 <pgc3_primary.ma> <PGC3_EUR.genes.out> <legacy_sets.txt> <out.txt>
"""
import sys

annot, ma25, go25, ma3, go3, legacy, out = sys.argv[1:8]

def gws(ma, thr=5e-8):
    s = set()
    with open(ma) as f:
        next(f)
        for line in f:
            v = line.split()
            if float(v[6]) < thr:
                s.add(v[0])
    return s

def genesig(go):
    rows = []
    with open(go) as f:
        hdr = f.readline().split(); ip = hdr.index("P"); ig = hdr.index("GENE")
        for line in f:
            v = line.split(); rows.append((v[ig], float(v[ip])))
    thr = 0.05 / len(rows)
    return {g for g, p in rows if p < thr}, len(rows), thr

def locus_pool(annot, snps):
    genes = set()
    with open(annot) as f:
        for line in f:
            if line.startswith("#"):
                continue
            v = line.rstrip("\n").split("\t")
            if any(s in snps for s in v[2:]):
                genes.add(v[0])
    return genes

g25, g3 = gws(ma25), gws(ma3)
lp25, lp3 = locus_pool(annot, g25), locus_pool(annot, g3)
gs25, n25, t25 = genesig(go25)
gs3, n3, t3 = genesig(go3)
print(f"SCZ25: {len(g25):,} GWS SNPs -> locus pool {len(lp25):,} genes; genesig {len(gs25):,} of {n25:,} (P<{t25:.2e})")
print(f"PGC3 : {len(g3):,} GWS SNPs -> locus pool {len(lp3):,} genes; genesig {len(gs3):,} of {n3:,} (P<{t3:.2e})")
print(f"overlap: locus pools {len(lp25 & lp3):,}; genesig {len(gs25 & gs3):,}; SCZ25 locus pool not in PGC3 pool {len(lp25 - lp3):,}")
with open(out, "w") as o:
    for name, genes in (("SCZ25_locus_pool", lp25), ("SCZ25_genesig", gs25),
                        ("PGC3_locus_pool", lp3), ("PGC3_genesig", gs3),
                        ("SCZ25_locus_pool_new", lp25 - lp3), ("SCZ25_genesig_new", gs25 - gs3)):
        for g in sorted(genes, key=int):
            o.write(f"{g}\t{name}\n")
    with open(legacy) as f:
        o.write(f.read())
print(f"wrote {out}")
