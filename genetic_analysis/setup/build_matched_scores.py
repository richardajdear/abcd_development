#!/usr/bin/env python
"""Ancestry-matched composite polygenic scores for the pooled arm.

WHY.  Until now the pooled (multi-ancestry, n = 8,596) arm used one weight set
for every child -- the multi-ancestry meta-analysis -- and standardised the
score within ancestry cluster (setup/standardise_within_ancestry.py) to strip
the between-cluster component.  With ancestry-specific discovery GWAS we can
also give each child the weights estimated in the discovery population closest
to their own, which is what PRS transferability theory says to do (LD and
allele frequencies match; effect-size attenuation across ancestry is avoided).

MAPPING (k-means cluster on in-sample PCs, strata_k4.tsv; legacy/hpc/README
validates them blind to self-report as 91 % White / 80 % Black / 95 % Hispanic
/ 45 % Asian):
    EURlike  -> EUR   weights
    cluster1 -> AFR   weights          (the African-American-like cluster)
    cluster3 -> EAS   weights          (the closest available; only ~45 % Asian,
                                        so read the per-stratum table, not this)
    cluster2 -> META  weights          (Hispanic/admixed-American-like; the 2025
                                        release has NO Latino/AMR cohort, so the
                                        trans-ancestry meta is the best available)
Each child's chosen score is then z-scored WITHIN their cluster, so the
composite has the same within-cluster-standardised form as the `_zanc`
profiles and can be read against them directly.

Built for every method that produced all four population weight sets:
    CT      -> CT/SCZ25_MATCHED/score_SCZ25MATCHED_<thr>.profile   (8 thresholds)
    PRSCSX  -> PRSCSX/SCZ25_MATCHED/score_SCZ25MATCHED_csx.profile (the per-
               population posteriors from the joint PRS-CSx model)

Usage: build_matched_scores.py <strata_k4.tsv> <scores_root>
"""
import csv, glob, math, os, re, sys

strata_path, root = sys.argv[1], sys.argv[2]
MAP = {"EURlike": "EUR", "cluster1": "AFR", "cluster3": "EAS", "cluster2": "META"}
POPS = ["EUR", "AFR", "EAS", "META"]

def norm(x):
    return re.sub(r"_", "", re.sub(r"^sub-", "", str(x))).upper()

clus = {}
with open(strata_path) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        clus[norm(row["IID"])] = row["stratum"]
print(f"clusters: {len(clus):,} subjects; mapping {MAP}", flush=True)

def read_profile(p):
    out = {}
    with open(p) as f:
        hdr = f.readline().split()
        idx = hdr.index("SCORESUM") if "SCORESUM" in hdr else hdr.index("SCORE")
        for line in f:
            v = line.split()
            if len(v) <= idx:
                continue
            out[v[1]] = (v[0], v[2], v[3], v[4], float(v[idx]))
    return out

def build(meth, token_suffix_re):
    dirs = {pop: os.path.join(root, meth, f"SCZ25_{pop}") for pop in POPS}
    if not all(os.path.isdir(d) for d in dirs.values()):
        print(f"{meth}: not all four population dirs present, skipping", flush=True)
        return
    # thresholds/tokens available in every population dir
    toks = None
    for pop, d in dirs.items():
        t = {re.match(rf"score_SCZ25{pop}_(.+)\.profile$", os.path.basename(p)).group(1)
             for p in glob.glob(os.path.join(d, "score_*.profile"))}
        toks = t if toks is None else toks & t
    toks = sorted(t for t in toks if re.fullmatch(token_suffix_re, t))
    if not toks:
        print(f"{meth}: no common score tokens", flush=True); return
    outd = os.path.join(root, meth, "SCZ25_MATCHED"); os.makedirs(outd, exist_ok=True)
    for tok in toks:
        prof = {pop: read_profile(os.path.join(dirs[pop], f"score_SCZ25{pop}_{tok}.profile")) for pop in POPS}
        ids = set.intersection(*(set(p) for p in prof.values()))
        chosen, groups = {}, {}
        for iid in ids:
            g = clus.get(norm(iid))
            if g not in MAP:
                continue
            pop = MAP[g]
            chosen[iid] = (g, pop, prof[pop][iid])
            groups.setdefault(g, []).append(prof[pop][iid][4])
        stats = {}
        for g, vals in groups.items():
            n = len(vals); m = sum(vals) / n
            sd = math.sqrt(sum((x - m) ** 2 for x in vals) / (n - 1)) if n > 1 else 0.0
            stats[g] = (m, sd)
        out = os.path.join(outd, f"score_SCZ25MATCHED_{tok}.profile")
        kept = 0
        with open(out, "w") as o:
            o.write("           FID            IID  PHENO    CNT   CNT2 SCORESUM\n")
            for iid in sorted(chosen):
                g, pop, (fid, ph, cnt, cnt2, s) = chosen[iid]
                m, sd = stats[g]
                if sd <= 0:
                    continue
                o.write(f"{fid:>14} {iid:>14} {ph:>6} {cnt:>6} {cnt2:>6} {(s - m) / sd:.6f}\n")
                kept += 1
        # n_snps is not one number for a composite, and prs_assoc.R would read
        # a <disorder>_<thr>.weights file's row count as n_snps -- so record the
        # per-population counts under a name it does not look for.
        with open(os.path.join(outd, f"nsnps_SCZ25MATCHED_{tok}.tsv"), "w") as w:
            for pop in POPS:
                wf = os.path.join(dirs[pop], f"SCZ25{pop}_{tok}.weights")
                n = sum(1 for _ in open(wf)) if os.path.exists(wf) else "NA"
                w.write(f"{pop}\t{n}\n")
        print(f"  {meth}/{tok}: {kept} subjects; per-cluster n "
              + ", ".join(f"{g}={len(v)}" for g, v in sorted(groups.items())), flush=True)
    with open(os.path.join(outd, "matched_mapping.tsv"), "w") as f:
        f.write("stratum\tweights\n" + "".join(f"{g}\t{p}\n" for g, p in MAP.items()))

build("CT", r"(5e-8|1e-5|0p001|0p01|0p05|0p1|0p5|1)")
build("PRSCSX", r"csx")
