"""Leave-fold-out test data from the reference panel itself (no ABCD data).

  python make_test_data.py --panel ../resources/MHC_...GRCh37.vcf.gz --out out/fold1_thin3 \
      --fold 1 --nfold 5 --thin 3 [--flip-frac 0.05] [--seed 1]

Splits the panel's 111 individuals into nfold folds.  The held-out fold becomes
a TARGET that looks like array data: C4 marker removed, every `thin`-th SNP
kept, genotypes unphased, a fraction of non-ambiguous SNPs strand-flipped, and
50 decoy SNPs at positions the panel does not have.  The other folds become
the REFERENCE.  Writes:
  <out>.ref.vcf.gz         reference panel without the held-out individuals
  <out>.target.vcf         target genotypes (VCF; run_smoketest.sh turns it into a PLINK fileset)
  <out>.truth.tsv          true C4 structures and GREx of the held-out individuals
Sample IDs become sub-T<NAxxxxx> so the 8-character token rule of the ABCD
pipeline (README_HPC.md rule 1) is exercised.
"""
from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location("c4grex", Path(__file__).resolve().parents[1] / "03_c4_grex.py")
c4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c4)  # noqa: E702

COMP = str.maketrans("ACGT", "TGCA")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fold", type=int, default=1)
    ap.add_argument("--nfold", type=int, default=5)
    ap.add_argument("--thin", type=int, default=3)
    ap.add_argument("--flip-frac", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bcftools", default="bcftools")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)

    header, rows = [], []
    with gzip.open(a.panel, "rt") as fh:
        for line in fh:
            (header if line.startswith("#") else rows).append(line.rstrip("\n"))
    samples = header[-1].split("\t")[9:]
    order = rng.permutation(len(samples))
    held = sorted(order[(a.fold - 1)::a.nfold])
    held_ids = [samples[i] for i in held]
    ref_ids = [s for s in samples if s not in held_ids]

    subprocess.run([a.bcftools, "view", "-Oz", "-o", f"{a.out}.ref.vcf.gz", "-s", ",".join(ref_ids), a.panel], check=True)
    subprocess.run([a.bcftools, "index", "-f", "-t", f"{a.out}.ref.vcf.gz"], check=True)

    # ---- truth -----------------------------------------------------------------
    c4row = next(r for r in rows if r.split("\t")[2] == "C4").split("\t")
    alleles = ["REF"] + [x.strip("<>") for x in c4row[4].split(",")]
    truth = []
    for i in held:
        h = [alleles[int(x)] for x in c4row[9 + i].split("|")]
        e = {s: c4.seg_counts(h[0])[s] + c4.seg_counts(h[1])[s] for s in c4.SEGMENTS}
        t = dict(sample=samples[i], IID=f"sub-T{samples[i]}", true_hap1=c4.structure(h[0]),
                 true_hap2=c4.structure(h[1]), **{f"true_{s}": v for s, v in e.items()})
        for g, w in c4.GREX_WEIGHTS.items():
            t[f"true_{g}"] = sum(wt * e[s] for s, wt in w.items())
        truth.append(t)
    pd.DataFrame(truth).to_csv(f"{a.out}.truth.tsv", sep="\t", index=False)

    # ---- target: thinned, unphased, some strand flips, decoys ------------------
    snps = [r for r in rows if r.split("\t")[2] != "C4"]
    keep = snps[::a.thin]
    ambiguous = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
    out = [l for l in header if l.startswith("##fileformat") or l.startswith("##contig")]
    out.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    out.append("\t".join(header[-1].split("\t")[:9] + [f"sub-T{s}" for s in held_ids]))
    nflip = 0
    body = []
    for r in keep:
        f = r.split("\t")
        if "," in f[4] or len(f[3]) != 1 or len(f[4]) != 1:
            continue
        ref, alt = f[3], f[4]
        if (ref, alt) not in ambiguous and rng.random() < a.flip_frac:
            ref, alt = ref.translate(COMP), alt.translate(COMP); nflip += 1
        gts = []
        for i in held:
            g = sorted(f[9 + i].split("|"))
            gts.append("/".join(g))
        body.append([f[0], f[1], f[2], ref, alt, ".", ".", ".", "GT"] + gts)
    pos = {int(b[1]) for b in body}
    for k in range(50):                                        # decoys: positions not in the panel
        p = int(rng.integers(24_900_000, 33_880_000))
        while p in pos:
            p += 1
        pos.add(p)
        g = ["/".join(sorted(rng.choice(["0", "1"], 2))) for _ in held]
        body.append(["6", str(p), f"decoy{k}", "A", "G", ".", ".", ".", "GT"] + g)
    body.sort(key=lambda b: int(b[1]))
    Path(f"{a.out}.target.vcf").write_text("\n".join(out + ["\t".join(b) for b in body]) + "\n")
    print(f"fold {a.fold}/{a.nfold}: {len(held_ids)} held out, {len(ref_ids)} reference; "
          f"target {len(body) - 50} panel SNPs (thin {a.thin}, {nflip} strand-flipped) + 50 decoys")


if __name__ == "__main__":
    main()
