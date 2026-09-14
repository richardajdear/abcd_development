#!/usr/bin/env python
"""Assemble the final PRS tables: four methods x every trait-arm.

ANCESTRY MATCHING is the organising principle, and the table now reports it
rather than resolving it silently.

prs_assoc.R fits two target strata for every score: `full` (n = 8,082,
multi-ancestry) and `EUR` (n = 4,116).  A score is ancestry-MATCHED when the
discovery GWAS's ancestry matches the target stratum:

  multi-ancestry discovery (SCZ primary, MDD div)  -> matched to `full`
  European discovery (everything else)             -> matched to `EUR`

Every trait except SCZ_pooled and MDD_pooled has a European-only discovery GWAS,
including all four controls -- ASD, both ALZ releases, and EA.  An earlier
version of this script assigned those to the pooled arm because they have no
ancestry-stratified release, which is true but is not a reason to read them in
the mismatched stratum.  Both strata are now emitted for every cell with a
`matched` flag, so the reader chooses and nothing is hidden; the console summary
shows matched cells only.

Mismatch is not symmetric and the two directions fail differently:
  European GWAS -> pooled target: the score correlates with ancestry and so does
    the phenotype, so the estimate is confounded.  This is the bad one.  It is
    visible in the output as an SE that inflates with threshold density (ALZ_IGAP
    pooled runs 0.014 -> 0.046), the signature of collinearity with the PCs.
  multi-ancestry GWAS -> EUR target: loses power and some discovery-side LD
    match, but is not confounded.

Output, all in results_v2/prs_final/:
  table_main.tsv       headline phenotypes, every method x trait-arm x stratum
  table_all.tsv        every phenotype, unfiltered
  table_family.tsv     Fulker between/within-family decomposition
"""
import csv, os, sys
from collections import OrderedDict

FINAL = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/results/prs_final"

METHODS = ["CT", "PRSCS", "SBayesR", "SBayesRC"]

# trait-arm -> (display trait, discovery ancestry, discovery GWAS description)
ARMS = OrderedDict([
    ("SCZ_pooled", ("SCZ",        "multi", "PGC3 primary (EUR+EAS+AFR+LAT)")),
    ("SCZ_eur",    ("SCZ",        "EUR",   "PGC3 european")),
    ("MDD_pooled", ("MDD",        "multi", "PGC MDD2025 div (trans-ancestry)")),
    ("MDD_eur",    ("MDD",        "EUR",   "PGC MDD2025 eur")),
    ("ASD",        ("ASD",        "EUR",   "SPARK+iPSYCH+PGC (no strata)")),
    ("ALZ",        ("ALZ",        "EUR",   "PGC-ALZ2 Wightman (incl. UKB proxy)")),
    ("ALZ_noAPOE", ("ALZnoAPOE",  "EUR",   "Wightman minus APOE region")),
    ("ALZ_IGAP",   ("ALZigap",    "EUR",   "Kunkle 2019 IGAP, diagnosed only")),
    ("ALZ_IGAP_noAPOE", ("ALZigapNoAPOE", "EUR", "Kunkle minus APOE region")),
    ("EA",         ("EA",         "EUR",   "Okbay 2016 EduYears, N=405,072")),
])
HEADLINE = ["global_slope", "baseline_thickness"]
MATCHED = {"multi": "full", "EUR": "EUR"}

def rd(path):
    if not os.path.exists(path): return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None

def best_per_pheno(rows, stratum):
    """One row per phenotype: the threshold with the smallest p.  C+T scans 8
    thresholds, so the reported p_adj is the Bonferroni-adjusted value
    prs_assoc.R already computed across them -- taking the raw minimum would be
    selecting on the outcome."""
    out = {}
    for r in rows:
        if r.get("stratum") != stratum: continue
        p = fnum(r.get("p"))
        if p is None: continue
        ph = r["phenotype"]
        if ph not in out or p < fnum(out[ph]["p"]):
            out[ph] = r
    return out

NSNP = {}   # (method, arm, pheno, threshold) -> n_snps, to fill the zanc rows
main, allrows, famrows = [], [], []

for meth in METHODS:
    for key, (trait, gwas_anc, gwas) in ARMS.items():
        for ver in ("raw", "zanc"):
            suf = "" if ver == "raw" else "_zanc"
            rows = rd(f"{FINAL}/assoc_{meth}_{key}{suf}.tsv")
            if not rows: continue
            for stratum in ("full", "EUR"):
                # within-ancestry standardisation is meaningless inside the EUR
                # stratum: it is one cluster, so the z-score is the raw score
                if ver == "zanc" and stratum == "EUR": continue
                for ph, r in sorted(best_per_pheno(rows, stratum).items()):
                    rec = OrderedDict([
                        ("method", meth), ("trait", trait), ("trait_arm", key),
                        ("discovery_ancestry", gwas_anc), ("target_stratum", stratum),
                        ("matched", "yes" if MATCHED[gwas_anc] == stratum else "no"),
                        ("score", ver), ("discovery_gwas", gwas),
                        ("phenotype", ph), ("threshold", r.get("threshold", "")),
                        ("n", r.get("n", "")), ("n_families", r.get("n_families", "")),
                        ("n_snps", r.get("n_snps", "")),
                        ("beta", r.get("beta", "")), ("se", r.get("se", "")),
                        ("p", r.get("p", "")),
                        ("p_adj", r.get("p_adj", r.get("p", ""))),
                    ])
                    k = (meth, key, ph, rec["threshold"], stratum)
                    if ver == "raw" and rec["n_snps"]:
                        NSNP[k] = rec["n_snps"]
                    elif not rec["n_snps"]:
                        rec["n_snps"] = NSNP.get(k, "")
                    allrows.append(rec)
                    if ph in HEADLINE: main.append(rec)
            for r in rd(f"{FINAL}/fam_{meth}_{key}{suf}.tsv"):
                rec = OrderedDict([("method", meth), ("trait", trait),
                                   ("trait_arm", key), ("score", ver)])
                rec.update(r)
                famrows.append(rec)

def write(path, rows):
    if not rows:
        print(f"  (nothing for {os.path.basename(path)})"); return
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"  {os.path.basename(path)}: {len(rows)} rows")

order = {m: i for i, m in enumerate(METHODS)}
arm_order = {k: i for i, k in enumerate(ARMS)}
main.sort(key=lambda r: (r["phenotype"] != "global_slope", arm_order[r["trait_arm"]],
                         r["target_stratum"] != MATCHED[r["discovery_ancestry"]],
                         order[r["method"]], r["score"]))
write(f"{FINAL}/table_main.tsv", main)
write(f"{FINAL}/table_all.tsv", allrows)
write(f"{FINAL}/table_family.tsv", famrows)

print("\nglobal_slope, ANCESTRY-MATCHED cells only, beta (SD per SD) and")
print("threshold-adjusted p.  * = p_adj < 0.05.  'zanc' = score z-standardised")
print("within ancestry cluster (pooled arm only).\n")
hdr = f"{'trait_arm':17s} {'strat':5s} {'score':5s} " + " ".join(f"{m:>17s}" for m in METHODS)
print(hdr); print("-" * len(hdr))
for key, (trait, gwas_anc, gwas) in ARMS.items():
    st = MATCHED[gwas_anc]
    for ver in ("raw", "zanc"):
        if ver == "zanc" and st == "EUR": continue
        cells = []
        for m in METHODS:
            hit = [r for r in main
                   if r["method"] == m and r["trait_arm"] == key
                   and r["target_stratum"] == st and r["score"] == ver
                   and r["phenotype"] == "global_slope"]
            if not hit:
                cells.append(f"{'-':>17s}"); continue
            b, p = fnum(hit[0]["beta"]), fnum(hit[0]["p_adj"])
            star = "*" if p is not None and p < 0.05 else " "
            cells.append(f"{b:>+8.4f} {p:7.1e}{star}".rjust(17))
        print(f"{key:17s} {st:5s} {ver:5s} " + " ".join(cells))
print("\nMismatched cells are in the tables with matched=no; for a European-only")
print("discovery GWAS the pooled-target row is ancestry-confounded, not just")
print("noisier, so it is not a robustness check on the matched row.")
