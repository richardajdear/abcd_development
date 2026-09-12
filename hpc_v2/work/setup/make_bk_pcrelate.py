"""Rewrite the Zaitlen second GRM so its non-zero pairs come from PC-Relate.

Reads the dense pooled GRM (GCTA binary, lower triangle including diagonal,
float32, row-major) and writes a copy in which an off-diagonal entry survives
ONLY IF PC-Relate calls that pair related.  Values and diagonal are taken
unchanged from the pooled GRM, so the component keeps Zaitlen's construction
(genotype-based sharing among close relatives); all that changes is WHICH pairs
are declared close, which is the thing the pooled GRM gets wrong.

Usage: make_bk_pcrelate.py SRC_PREFIX PAIRS_TSV DST_PREFIX
"""
import sys
import numpy as np

src, pairs_tsv, dst = sys.argv[1], sys.argv[2], sys.argv[3]

ids = [l.split() for l in open(src + ".grm.id")]
iid_index = {r[1]: i for i, r in enumerate(ids)}
n = len(ids)

from collections import defaultdict

cols_by_row = defaultdict(list)          # row -> columns kept, row > col
n_pairs = n_unmapped = 0
with open(pairs_tsv) as fh:
    fh.readline()                        # header
    for line in fh:
        a, b = line.split("\t")[:2]
        ia, ib = iid_index.get(a.strip()), iid_index.get(b.strip())
        if ia is None or ib is None:
            n_unmapped += 1
            continue
        cols_by_row[max(ia, ib)].append(min(ia, ib))
        n_pairs += 1
print(f"{n} GRM ids; {n_pairs} related pairs mapped into the GRM"
      + (f" ({n_unmapped} unmapped)" if n_unmapped else ""))

kept = dropped = 0
with open(src + ".grm.bin", "rb") as fi, open(dst + ".grm.bin", "wb") as fo:
    for i in range(n):
        row = np.fromfile(fi, dtype=np.float32, count=i + 1)
        if i:
            cols = cols_by_row.get(i)
            off = row[:i]
            if cols:
                idx = np.asarray(cols, dtype=np.int64)
                vals = off[idx].copy()
                dropped += int(np.count_nonzero(off)) - int(np.count_nonzero(vals))
                kept += int(np.count_nonzero(vals))
                off[:] = 0.0
                off[idx] = vals
            else:
                dropped += int(np.count_nonzero(off))
                off[:] = 0.0
        row.tofile(fo)
print(f"off-diagonal entries kept: {kept}   zeroed: {dropped}")
