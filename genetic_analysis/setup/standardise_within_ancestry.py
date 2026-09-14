"""Rewrite PLINK .profile scores as z-scores WITHIN ancestry cluster.

WHY.  The pooled-sample polygenic associations are inflated by ancestry: a
subject's PRS correlates with their ancestry, and cortical thickness/slope
varies with ancestry too, so part of the association is population structure
rather than biology.  The evidence is that restricting to Europeans roughly
HALVES every beta while the standard errors FALL -- a coefficient that is both
larger and less precisely estimated in the pooled arm is the signature of a
predictor collinear with the covariates it is adjusted for, not of a real
effect diluted by noise.

Centring and scaling the score within each ancestry cluster removes the
between-cluster component of the score's variance, so the association is
driven only by within-cluster variation.  This keeps the full n = 8,082 rather
than discarding the 3,966 non-European children, which is what the EUR-only
arm does.

It is not a complete fix -- residual within-cluster structure survives, and
k-means clusters are a coarse summary of a continuum -- but it is the cheapest
correction that uses the whole sample, and it is directly comparable to the
uncorrected numbers because nothing else changes.

Clusters come from v1's k-means on in-sample ancestry PCs (strata_k4.tsv, the
file legacy/hpc/config.sh names as the recommended STRATA_FILE): EURlike 7,448,
cluster1 2,309, cluster2 1,390, cluster3 523 of 11,670 genotyped.

Usage: standardise_within_ancestry.py <strata.tsv> <in_dir> <out_dir>
"""
import csv, glob, os, re, sys, math

strata_path, in_dir, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

def norm(x):
    return re.sub(r"_", "", re.sub(r"^sub-", "", str(x))).upper()

clus = {}
with open(strata_path) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        clus[norm(row["IID"])] = row["stratum"]
print(f"clusters: {len(clus):,} subjects, {len(set(clus.values()))} groups", flush=True)

for src in sorted(glob.glob(os.path.join(in_dir, "score_*.profile"))):
    rows = []
    with open(src) as f:
        header = f.readline().split()
        for line in f:
            v = line.split()
            if len(v) < 6:
                continue
            rows.append(v)
    idx = header.index("SCORESUM")
    # group scores by cluster, then z-score within each
    groups = {}
    for v in rows:
        g = clus.get(norm(v[1]), "UNASSIGNED")
        groups.setdefault(g, []).append(float(v[idx]))
    stats = {}
    for g, vals in groups.items():
        n = len(vals)
        m = sum(vals) / n
        sd = math.sqrt(sum((x - m) ** 2 for x in vals) / (n - 1)) if n > 1 else 0.0
        stats[g] = (m, sd)
    out = os.path.join(out_dir, os.path.basename(src))
    kept = dropped = 0
    with open(out, "w") as o:
        o.write("           FID            IID  PHENO    CNT   CNT2 SCORESUM\n")
        for v in rows:
            g = clus.get(norm(v[1]), "UNASSIGNED")
            m, sd = stats[g]
            if sd <= 0:
                dropped += 1
                continue
            z = (float(v[idx]) - m) / sd
            o.write(f"{v[0]:>14} {v[1]:>14} {v[2]:>6} {v[3]:>6} {v[4]:>6} {z:.6f}\n")
            kept += 1
    print(f"  {os.path.basename(src)}: {kept} subjects"
          + (f" ({dropped} dropped, zero-variance cluster)" if dropped else ""), flush=True)
