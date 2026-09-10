"""Assemble every polygenic-score result into two tidy tables.

Writes, under results_v2/prs_tables_v2/:
  prs_all_methods.tsv        every (method, trait, phenotype, stratum) row
  prs_all_withinfamily.tsv   the between/within-family decomposition of each

Methods covered: C+T (best-of-8 thresholds and per-threshold), PRS-CS on the
1000G reference, PRS-CS on the UKB reference, PRS-CSx (EUR/EAS/sum), SBayesR,
SBayesRC, and the within-ancestry-standardised version of each.
"""
import csv, glob, os, re, sys

REPO = "/home/rajd2/rds/hpc-work/abcd_development"
V1 = f"{REPO}/hpc/work/results"
V2 = f"{REPO}/hpc_v2/work/results_v2"
OUT = f"{V2}/prs_tables_v2"
os.makedirs(OUT, exist_ok=True)

# (method label, ancestry-standardised?, path, how the disorder is named)
SRC = [
    ("C+T",            0, f"{V1}/prs_imp/prs_association.tsv"),
    ("C+T",            0, f"{V2}/control_ASD/prs_association_ASD.tsv"),
    ("C+T",            0, f"{V2}/control_ALZ/prs_association_ALZ.tsv"),
    ("PRS-CS/1000G",   0, f"{V2}/prscs/prs_association_*_prscs.tsv"),
    ("PRS-CS/UKB",     0, f"{V2}/prscs_ukbb/prs_association_*_prscs_ukbb.tsv"),
    ("PRS-CS/UKB",     0, f"{V2}/prscs_ukbb/prs_association_ALZnoAPOE_ukbb.tsv"),
    ("PRS-CSx",        0, f"{V2}/prscsx/prs_association_SCZcsx*.tsv"),
    ("SBayesR",        0, f"{V2}/sbayesr/prs_association_*_sbayesr.tsv"),
    ("SBayesRC",       0, f"{V2}/sbayesrc/prs_association_*_sbayesrc.tsv"),
    ("*_within-ancestry", 1, f"{V2}/within_ancestry/prs_association_*_wanc.tsv"),
]
FAM = [(m, z, p.replace("prs_association", "prs_withinfamily")) for m, z, p in SRC]

def rows_from(pattern, method, zflag):
    out = []
    for f in sorted(glob.glob(pattern)):
        base = os.path.basename(f)
        # the within-ancestry files encode the underlying method in the name
        meth = method
        if zflag:
            m = re.match(r"prs_(?:association|withinfamily)_([A-Za-z0-9]+)_", base)
            meth = (m.group(1) if m else "?") + " (within-ancestry z)"
        try:
            with open(f) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    r["_method"] = meth
                    r["_source"] = os.path.relpath(f, REPO)
                    out.append(r)
        except Exception as e:
            print(f"  skip {f}: {e}", file=sys.stderr)
    return out

def write(rows, path, cols):
    with open(path, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])
    print(f"{path}: {len(rows)} rows")

assoc, fam = [], []
for m, z, p in SRC:
    assoc += rows_from(p, m, z)
for m, z, p in FAM:
    fam += rows_from(p, m, z)

A = ["_method", "disorder", "threshold", "phenotype", "stratum", "n", "n_snps",
     "beta", "se", "p", "p_adj", "r2_partial", "_source"]
F = ["_method", "disorder", "threshold", "phenotype", "stratum", "n", "n_pairs",
     "beta_between", "se_between", "p_between", "beta_within", "se_within",
     "p_within", "beta_diff", "p_diff", "_source"]
write(assoc, f"{OUT}/prs_all_methods.tsv", A)
write(fam,   f"{OUT}/prs_all_withinfamily.tsv", F)
