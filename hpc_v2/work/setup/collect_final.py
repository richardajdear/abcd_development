#!/usr/bin/env python
"""Assemble the final PRS tables: four methods x two target arms.

ARM MATCHING is the organising principle.  A polygenic score is only as good as
the match between the ancestry of the discovery GWAS and the ancestry of the
target sample, so each arm is read off the GWAS built for it:

  pooled arm  n = 8,082 multi-ancestry target  <- SCZ primary / MDD div
  EUR arm     n = 4,116 European target        <- SCZ european / MDD eur

prs_assoc.R fits both strata for every score, so we select the (score, stratum)
pairs that are matched and ignore the rest.  ASD and ALZ have no stratified
release: one file serves both arms, and that is marked in the table.

For the pooled arm there are two versions of every number:
  raw    score used as PLINK wrote it
  zanc   score z-standardised within ancestry cluster before the regression,
         which removes the between-cluster component of the score variance
Differences between raw and zanc are ancestry stratification, not biology.

Output, all in results_v2/prs_final/:
  table_main.tsv       global_slope and baseline_thickness, every method x arm
  table_all.tsv        every phenotype, unfiltered
  table_family.tsv     Fulker between/within-family decomposition
"""
import csv, glob, os, sys
from collections import OrderedDict

FINAL = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/rajd2/rds/hpc-work/abcd_development/hpc_v2/work/results_v2/prs_final"

METHODS  = ["CT", "PRSCS", "SBayesR", "SBayesRC"]
# trait-arm -> (display trait, arm, discovery GWAS, stratum to read)
ARMS = OrderedDict([
    ("SCZ_pooled", ("SCZ", "pooled", "PGC3 primary (EUR+EAS+AFR+LAT)", "full")),
    ("SCZ_eur",    ("SCZ", "EUR",    "PGC3 european",                  "EUR")),
    ("MDD_pooled", ("MDD", "pooled", "PGC MDD2025 div (trans-anc)",    "full")),
    ("MDD_eur",    ("MDD", "EUR",    "PGC MDD2025 eur",                "EUR")),
    ("ASD",        ("ASD", "pooled", "SPARK+iPSYCH+PGC (no strata)",   "full")),
    ("ALZ",        ("ALZ", "pooled", "PGC-ALZ2 Wightman (no strata)",   "full")),
    # The APOE-excluded score is the one the specificity claim rests on: with
    # APOE in, ALZ looks as strong as SCZ, and that is one large-effect locus
    # rather than polygenic AD risk.
    ("ALZ_noAPOE", ("ALZnoAPOE", "pooled", "PGC-ALZ2 minus chr19:44.4-46.5Mb", "full")),
])
HEADLINE = ["global_slope", "baseline_thickness"]

def rd(path):
    if not os.path.exists(path): return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None

def best(rows, stratum):
    """One row per phenotype.  C+T scans 8 thresholds, so take the threshold
    with the smallest p and report the Bonferroni-adjusted p that prs_assoc.R
    already computed across those thresholds -- reporting the raw minimum would
    be selection on the outcome."""
    out = {}
    for r in rows:
        if r.get("stratum") != stratum: continue
        p = fnum(r.get("p"))
        if p is None: continue
        ph = r["phenotype"]
        if ph not in out or p < fnum(out[ph]["p"]):
            out[ph] = r
    return out

# (method, trait-arm, phenotype, threshold) -> n_snps, filled from the raw
# tables so the standardised rows can borrow it: same weights, same SNPs.
NSNP = {}

main, allrows, famrows = [], [], []
for meth in METHODS:
    for key, (trait, arm, gwas, stratum) in ARMS.items():
        for ver in ("raw", "zanc"):
            # the EUR arm has no ancestry-standardised version: one cluster
            if ver == "zanc" and arm == "EUR": continue
            suf = "" if ver == "raw" else "_zanc"
            rows = rd(f"{FINAL}/assoc_{meth}_{key}{suf}.tsv")
            if not rows: continue
            sel = best(rows, stratum)
            for ph, r in sorted(sel.items()):
                rec = OrderedDict([
                    ("method", meth), ("trait", trait), ("arm", arm),
                    ("score", ver), ("discovery_gwas", gwas),
                    ("phenotype", ph), ("threshold", r.get("threshold", "")),
                    ("n", r.get("n", "")), ("n_families", r.get("n_families", "")),
                    ("n_snps", r.get("n_snps", "")),
                    ("beta", r.get("beta", "")), ("se", r.get("se", "")),
                    ("p", r.get("p", "")), ("p_adj", r.get("p_adj", r.get("p", ""))),
                ])
                k = (meth, key, ph, rec["threshold"])
                if ver == "raw" and rec["n_snps"]:
                    NSNP[k] = rec["n_snps"]
                elif not rec["n_snps"]:
                    rec["n_snps"] = NSNP.get(k, "")
                allrows.append(rec)
                if ph in HEADLINE: main.append(rec)
        for ver in ("raw", "zanc"):
            if ver == "zanc" and arm == "EUR": continue
            suf = "" if ver == "raw" else "_zanc"
            for r in rd(f"{FINAL}/fam_{meth}_{key}{suf}.tsv"):
                rec = OrderedDict([("method", meth), ("trait", trait),
                                   ("arm", arm), ("score", ver)])
                rec.update(r)
                famrows.append(rec)

def write(path, rows):
    if not rows:
        print(f"  (nothing for {os.path.basename(path)})"); return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t",
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"  {os.path.basename(path)}: {len(rows)} rows")

order = {m: i for i, m in enumerate(METHODS)}
main.sort(key=lambda r: (r["phenotype"] != "global_slope", r["trait"],
                         r["arm"], order[r["method"]], r["score"]))
write(f"{FINAL}/table_main.tsv", main)
write(f"{FINAL}/table_all.tsv", allrows)
write(f"{FINAL}/table_family.tsv", famrows)

# ---- readable summary of the headline phenotype -----------------------------
print("\nglobal_slope, beta (SD per SD of score), p  [* p<0.05]")
hdr = f"{'trait':5s} {'arm':7s} {'score':5s} " + " ".join(f"{m:>18s}" for m in METHODS)
print(hdr); print("-" * len(hdr))
for key, (trait, arm, gwas, stratum) in ARMS.items():
    for ver in ("raw", "zanc"):
        if ver == "zanc" and arm == "EUR": continue
        cells = []
        for m in METHODS:
            hit = [r for r in main if r["method"] == m and r["trait"] == trait
                   and r["arm"] == arm and r["score"] == ver
                   and r["phenotype"] == "global_slope"]
            if not hit: cells.append(f"{'-':>18s}"); continue
            b, p = fnum(hit[0]["beta"]), fnum(hit[0]["p_adj"])
            star = "*" if p is not None and p < 0.05 else " "
            cells.append(f"{b:>+9.4f} {p:8.2e}{star}"[:18].rjust(18))
        print(f"{trait:5s} {arm:7s} {ver:5s} " + " ".join(cells))
