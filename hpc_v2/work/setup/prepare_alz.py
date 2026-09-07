"""Map PGC-ALZ2 chr:pos to rsID, so it can join the rsID-keyed pipeline.

PGCALZ2sumstatsExcluding23andMe.txt (Wightman et al. 2021, PGC-ALZ stage,
excluding 23andMe) carries chr / PosGRCh37 / testedAllele / otherAllele / z /
p / N and NO marker ID.  Everything downstream here -- MAGMA, LDSC and the PRS
-- joins on rsID (hpc/README_HPC.md 8.16.7), so the file is unusable until the
IDs are restored.

Mapping target is the PRS genotype fileset's .bim, deliberately: that is the
exact variant space the SCZ and MDD scores were built in, so the Alzheimer's
control is scored over the same SNPs rather than a different subset.  Both
sides are GRCh37 (the .bim after v1's rsID restore, and the file's own
PosGRCh37 column), so positions are directly comparable.

Alleles are checked, not assumed: a position match is kept only if the tested/
other pair equals the .bim pair, in either order or as its strand complement.
Ambiguous A/T and C/G SNPs are dropped -- with no allele frequency column in
this file there is no way to resolve their strand.
"""
import sys, gzip

bim_path, alz_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIG = {frozenset("AT"), frozenset("CG")}

bim = {}
with open(bim_path) as f:
    for line in f:
        c, rs, _, pos, a1, a2 = line.split()
        bim[(c, pos)] = (rs, a1.upper(), a2.upper())
print(f"bim: {len(bim):,} positions", flush=True)

kept = no_pos = allele_mismatch = ambiguous = 0
with open(alz_path) as fi, open(out_path, "w") as fo:
    fo.write("SNP\tA1\tA2\tZ\tP\tN\n")
    hdr = fi.readline().split()
    idx = {name.strip(): i for i, name in enumerate(hdr)}
    ic, ip = idx["chr"], idx["PosGRCh37"]
    ita, ioa, iz, ipv, iN = idx["testedAllele"], idx["otherAllele"], idx["z"], idx["p"], idx["N"]
    for line in fi:
        f = line.rstrip("\n").split("\t")
        hit = bim.get((f[ic].strip(), f[ip].strip()))
        if hit is None:
            no_pos += 1
            continue
        rs, b1, b2 = hit
        ta, oa = f[ita].strip().upper(), f[ioa].strip().upper()
        if frozenset((ta, oa)) in AMBIG:
            ambiguous += 1
            continue
        same = {ta, oa} == {b1, b2}
        flipped = {COMP.get(ta, "?"), COMP.get(oa, "?")} == {b1, b2}
        if not (same or flipped):
            allele_mismatch += 1
            continue
        fo.write(f"{rs}\t{ta}\t{oa}\t{f[iz]}\t{f[ipv]}\t{f[iN]}\n")
        kept += 1
print(f"kept {kept:,} | no position match {no_pos:,} | allele mismatch {allele_mismatch:,} "
      f"| ambiguous dropped {ambiguous:,}", flush=True)
