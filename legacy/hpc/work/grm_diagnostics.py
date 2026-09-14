"""Describe a GCTA GRM: diagonal, off-diagonal, and who --grm-cutoff removes.

    python grm_diagnostics.py <grm-prefix> [--strata strata.tsv] [--cutoff 0.05]

WHY

`--grm-cutoff 0.05` is conventional for SNP-h2 and is what 01_grm and the
all-ancestry merge apply.  In a SINGLE-ancestry sample it does what it says:
drop one member of each pair related above 0.05, so shared environment cannot
load onto the additive term.

In a MULTI-ancestry sample it is not so simple, and that is the point of this
script.  Across ancestry groups the GRM off-diagonal is systematically negative
(individuals from different groups are less alike than the sample average);
within a group it is systematically positive.  A single global threshold
therefore does not treat the groups alike: it can prune preferentially inside
whichever groups are most homogeneous or most densely sampled, changing the
ancestry composition of the unrelated subset relative to the full one.

That is a real effect on what the pooled h2 describes, and it is invisible in
the h2 itself.  So it is measured here -- per-stratum retention through the
cutoff, and the within- vs between-stratum off-diagonal distributions -- rather
than left as a caveat nobody quantified.

Reads the lower triangle GCTA writes to .grm.bin as float32, n(n+1)/2 entries
in row-major order with the diagonal last in each row.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

TOKEN = re.compile(r"([A-Z0-9]{8})$")


def load_grm(prefix: str):
    ids = [ln.split()[1] for ln in Path(f"{prefix}.grm.id").read_text().splitlines() if ln.strip()]
    n = len(ids)
    v = np.fromfile(f"{prefix}.grm.bin", dtype=np.float32)
    exp = n * (n + 1) // 2
    if v.size != exp:
        raise SystemExit(f"{prefix}.grm.bin has {v.size} entries, expected {exp} for n={n}")
    return ids, n, v


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--strata", type=Path, default=None,
                    help="strata.tsv from assign_ancestry.py")
    ap.add_argument("--cutoff", type=float, default=0.05)
    ap.add_argument("--unrel", default=None,
                    help="prefix of the pruned GRM, to report who survived")
    a = ap.parse_args()

    ids, n, v = load_grm(a.prefix)
    idx = np.arange(n)
    diag_pos = idx * (idx + 1) // 2 + idx
    diag = v[diag_pos]
    off = np.delete(v, diag_pos)

    print(f"{a.prefix}: n={n}, {v.size} stored entries")
    print(f"  diagonal     mean {diag.mean():+.4f}  sd {diag.std():.4f}  "
          f"min {diag.min():+.4f}  max {diag.max():+.4f}")
    print(f"  off-diagonal mean {off.mean():+.4f}  sd {off.std():.4f}  "
          f"min {off.min():+.4f}  max {off.max():+.4f}")
    for t in (0.05, 0.1, 0.2, 0.4):
        print(f"  pairs > {t:<4}: {int((off > t).sum()):>12,}  "
              f"({100 * (off > t).mean():.4f}% of pairs)")

    # A diagonal far from 1 is the standard smell test for a GRM built on the
    # wrong allele frequencies or a badly-called variant set; say so explicitly
    # rather than leaving the reader to know that 1.0 is the expectation.
    if abs(diag.mean() - 1.0) > 0.05:
        print(f"  WARNING: mean diagonal {diag.mean():.4f} is far from 1.0")

    if a.strata and a.strata.exists():
        lab = {}
        for ln in a.strata.read_text().splitlines()[1:]:
            f = ln.split("\t")
            if len(f) >= 3:
                lab[f[0]] = f[2]
        names = sorted(set(lab.values()))
        code = np.array([names.index(lab[i]) if i in lab else -1 for i in ids])
        print(f"\n  strata: {', '.join(f'{s} n={(code == j).sum()}' for j, s in enumerate(names))}")

        # Off-diagonal within vs between strata.  Built by rows to avoid
        # materialising an n x n matrix for 11,670 subjects.
        wi, be = [], []
        pos = 0
        for i in range(n):
            row = v[pos:pos + i]            # entries (i,0..i-1)
            pos += i + 1
            if i == 0 or code[i] < 0:
                continue
            same = code[:i] == code[i]
            wi.append(row[same & (code[:i] >= 0)])
            be.append(row[~same & (code[:i] >= 0)])
        wi = np.concatenate(wi) if wi else np.array([])
        be = np.concatenate(be) if be else np.array([])
        print(f"  within-stratum  off-diagonal: mean {wi.mean():+.4f}  sd {wi.std():.4f}  n={wi.size:,}")
        print(f"  between-stratum off-diagonal: mean {be.mean():+.4f}  sd {be.std():.4f}  n={be.size:,}")
        print("  (a large gap is population structure in the GRM, which is what "
              "the in-sample PCs are covarying out)")

        if a.unrel:
            kept = {ln.split()[1] for ln in Path(f"{a.unrel}.grm.id").read_text().splitlines() if ln.strip()}
            print(f"\n  retention through --grm-cutoff {a.cutoff}, by stratum:")
            for j, s in enumerate(names):
                mem = [ids[i] for i in np.flatnonzero(code == j)]
                k = sum(1 for m in mem if m in kept)
                print(f"    {s:<12} {k:>6} of {len(mem):>6}  ({100 * k / max(len(mem), 1):.1f}%)")
            print("  (materially different retention means the unrelated subset "
                  "has a different ancestry composition from the full sample)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
