"""Assign ancestry strata to the all-ancestry sample from in-sample genetic PCs.

    python assign_ancestry.py <eigenvec> <out.tsv> [--eur-anchor FAM] [--k N]

WHY THIS EXISTS
---------------
README_HPC.md section 5.3.2 requires an ancestry-stratified REML alongside the
pooled one, because a single GRM across ancestries assumes a common allele
frequency and LD structure that does not hold.  Running that needs a stratum
label per subject, and **ABCD 7.0 does not ship one**.  Checked, not assumed:
the only genetics columns in the 7.0 tabulated release are
`ab_g_stc__gen_pc__01..32` (32 PCs) and `gn_y_genrel_*` (pihat/zygosity); there
is no ancestry-proportion or ancestry-group variable anywhere under
`derivatives/tabulated`.  Nor is a multi-ancestry 1000 Genomes panel available
on this account -- `hpc-work` holds `g1000_eur` and `g1000_eas` only -- so the
usual "project onto 1000G and assign by nearest reference centroid" is not
available without downloading a reference panel.

So the labels here are derived, and they are derived two ways so that neither
has to be trusted alone:

1. **EUR_anchor** -- membership of the `abcd_eur` fileset that produced every
   published result in this project.  This is somebody else's genetic-ancestry
   call, made with reference panels we do not have, and it is the *only* label
   here that is not our own inference.  It is what makes the stratified EUR h2
   directly comparable to the published EUR h2: same ancestry definition, and
   (because it runs on this GRM) the same variant set.

2. **kmeans clusters** on the top in-sample PCs, for the strata the anchor does
   not name.  Reported with each cluster's overlap against the anchor, so the
   reader can see how well the unsupervised split reproduces a label that was
   made independently.  A cluster that is 97 % anchor-EUR is believable; one
   that is 60 % is a warning that the split is not clean.

WHAT THIS IS NOT
----------------
It is not a continental-ancestry classifier.  ABCD contains a large admixed
group, and k-means assigns an admixed individual to whichever centroid happens
to be nearest -- there is no partial membership.  Strata built this way are
"regions of PC space", which is enough to ask whether h2 is stable across them
and is NOT enough to report an h2 "in African-ancestry participants".  Label
the results accordingly.

k-means is implemented here rather than imported: the env has numpy but adding
a scikit-learn dependency for 30 lines of Lloyd's algorithm is not worth it,
and a fixed seed with k-means++ init makes the assignment reproducible.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

TOKEN = re.compile(r"([A-Z0-9]{8})$")
SEED = 20260817


def token(s: str) -> str | None:
    m = TOKEN.search(s.strip())
    return m.group(1) if m else None


def read_eigenvec(path: Path, n_pcs: int):
    """GCTA's headerless `FID IID PC1..PCk` -> (iids, tokens, scaled PC matrix).

    Columns are scaled by sqrt(eigenvalue) when the matching `.eigenval` is
    present, and THAT MATTERS MORE THAN IT LOOKS.  GCTA writes unit-norm
    eigenVECTORS, so every column comes out with roughly the same spread --
    measured here, SD 0.0091, 0.0088, 0.0088, 0.0093 for PC1-4 -- even though
    their eigenvalues are 727, 175, 49 and 15, i.e. 6.04 %, 1.45 %, 0.41 % and
    0.13 % of variance.

    Clustering the raw columns therefore weights a PC explaining 0.13 % exactly
    as heavily as one explaining 6.04 %, and the noise PCs drive the split.
    Observed: on the unscaled matrix, k=4 put 8,385 subjects in one cluster at
    54.8 % anchor-EUR and named a separate 1,477-subject cluster "EURlike" at
    72 % -- the EUR mass split in two and neither piece clean.  Scaling to PC
    SCORES (eigenvector * sqrt(eigenvalue)) is the standard construction and is
    what makes the geometry mean what k-means assumes it means.
    """
    iids, toks, mat = [], [], []
    for line in path.read_text().splitlines():
        f = line.split()
        if len(f) < 2 + n_pcs:
            continue
        iids.append(f[1])
        toks.append(token(f[1]))
        mat.append([float(v) for v in f[2:2 + n_pcs]])
    if not mat:
        raise SystemExit(f"{path}: no rows with at least {n_pcs} PCs")
    X = np.asarray(mat, dtype=float)

    eigenval = path.with_suffix(".eigenval")
    if eigenval.exists():
        ev = np.array([float(v) for v in eigenval.read_text().split()][:n_pcs])
        if ev.size == n_pcs and (ev > 0).all():
            X = X * np.sqrt(ev)
            print(f"  scaled PC1..PC{n_pcs} by sqrt(eigenvalue) "
                  f"({', '.join(f'{v:.1f}' for v in ev)})")
        else:
            print(f"  WARNING: {eigenval.name} unusable; clustering UNSCALED "
                  "eigenvectors, which over-weights the low-variance PCs")
    else:
        print(f"  WARNING: no {eigenval.name}; clustering UNSCALED eigenvectors, "
              "which over-weights the low-variance PCs")
    return iids, toks, X


def kmeans(X: np.ndarray, k: int, seed: int = SEED, iters: int = 300):
    """Lloyd's algorithm with k-means++ init.  Deterministic given `seed`."""
    rng = np.random.default_rng(seed)
    n = X.shape[0]

    # k-means++ seeding: first centre uniform, each later centre sampled with
    # probability proportional to squared distance from the nearest chosen one.
    # Plain random init on genetic PCs routinely puts two centres inside the
    # same dense cluster and leaves a real group unsplit.
    centres = [X[rng.integers(n)]]
    for _ in range(1, k):
        d2 = np.min(((X[:, None, :] - np.asarray(centres)[None, :, :]) ** 2).sum(-1), axis=1)
        total = d2.sum()
        if total <= 0:
            centres.append(X[rng.integers(n)])
            continue
        centres.append(X[rng.choice(n, p=d2 / total)])
    C = np.asarray(centres)

    lab = np.zeros(n, dtype=int)
    for _ in range(iters):
        d2 = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1)
        new = d2.argmin(1)
        if np.array_equal(new, lab):
            break
        lab = new
        for j in range(k):
            m = lab == j
            if m.any():
                C[j] = X[m].mean(0)
    return lab, C


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("eigenvec", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--eur-anchor", type=Path, default=None,
                    help="a .fam whose subjects are the reference EUR set")
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--n-pcs", type=int, default=4,
                    help="PCs entering k-means (raw, so variance-weighted)")
    a = ap.parse_args(argv)

    iids, toks, X = read_eigenvec(a.eigenvec, a.n_pcs)
    print(f"{a.eigenvec}: {len(iids)} subjects, using PC1..PC{a.n_pcs}")

    anchor: set[str] = set()
    if a.eur_anchor:
        for line in a.eur_anchor.read_text().splitlines():
            f = line.split()
            if len(f) >= 2 and (t := token(f[1])):
                anchor.add(t)
        print(f"{a.eur_anchor}: {len(anchor)} EUR-anchor subjects")

    lab, C = kmeans(X, a.k)

    # Name the cluster holding the most anchor members EUR; the rest keep an
    # index.  Naming the others AFR/AMR/EAS would be an assertion this method
    # cannot support (see the module docstring).
    names = {}
    if anchor:
        # By anchor COUNT, not anchor SHARE.  Share picks whichever cluster is
        # purest, which on a bad split is a small offshoot: the unscaled run
        # named a 1,477-subject cluster at 72 % "EURlike" over an
        # 8,385-subject one at 54.8 %, even though the latter held four times
        # as many anchor subjects.  The EUR cluster is the one where the
        # anchor's 5,656 members actually are.
        count = {}
        for j in range(a.k):
            m = lab == j
            count[j] = sum(1 for i in np.flatnonzero(m) if toks[i] in anchor)
        eur_j = max(count, key=count.get)
        names[eur_j] = "EURlike"
        nxt = 1
        for j in sorted(set(range(a.k)) - {eur_j}, key=lambda j: -(lab == j).sum()):
            names[j] = f"cluster{nxt}"
            nxt += 1
    else:
        for j in range(a.k):
            names[j] = f"cluster{j + 1}"

    print(f"\n{'stratum':<12}{'n':>7}{'anchor-EUR':>12}{'% of stratum':>14}   PC centroid")
    for j in sorted(range(a.k), key=lambda j: -(lab == j).sum()):
        m = lab == j
        n_anchor = sum(1 for i in np.flatnonzero(m) if toks[i] in anchor)
        pct = 100 * n_anchor / max(m.sum(), 1)
        cen = "  ".join(f"{v:+.4f}" for v in C[j][:3])
        print(f"{names[j]:<12}{m.sum():>7}{n_anchor:>12}{pct:>13.1f}%   {cen}")

    if anchor:
        assigned = sum(1 for i, t in enumerate(toks) if t in anchor and names[lab[i]] == "EURlike")
        print(f"\nanchor EUR subjects landing in EURlike: {assigned} of "
              f"{sum(1 for t in toks if t in anchor)} present in the eigenvec")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w") as fh:
        fh.write("IID\ttoken\tstratum\tanchor_eur\t"
                 + "\t".join(f"PC{i+1}" for i in range(a.n_pcs)) + "\n")
        for i, iid in enumerate(iids):
            fh.write(f"{iid}\t{toks[i]}\t{names[lab[i]]}\t"
                     f"{int(toks[i] in anchor)}\t"
                     + "\t".join(f"{v:.6f}" for v in X[i]) + "\n")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
