"""Assemble the v1-vs-v2 comparison tables (GWAS, LDSC rg, MAGMA, PRS).

Every number is read from a summary file on disk; nothing is hard-coded except
the file locations.  Run after the EUR chain completes.
"""
import csv, os, sys, re

REPO = "/home/rajd2/rds/hpc-work/abcd_development"
V1 = f"{REPO}/hpc/work/results"
V2 = f"{REPO}/hpc_v2/work/results_v2"
PH = ["baseline_thickness", "global_slope", "slope_PC1", "slope_PC2", "slope_PC3"]

def read_tsv(p):
    if not os.path.exists(p): return []
    with open(p) as f: return list(csv.DictReader(f, delimiter="\t"))

def fmt(x, n=4):
    try: return f"{float(x):.{n}f}"
    except (TypeError, ValueError): return "NA"

def sci(x):
    try: return f"{float(x):.2e}"
    except (TypeError, ValueError): return "NA"

def section(t): print(f"\n{'='*100}\n{t}\n{'='*100}")

# ---- 1. GWAS summaries -------------------------------------------------------
section("1. GWAS: v1 (fastGWA, GCTA GRM) vs v2 (GENESIS, PC-AiR + PC-Relate)")
for arm, v1p, v2p in [("POOLED", f"{V1}/gwas_imp/gwas_summary.tsv", f"{V2}/assoc/gwas_summary.tsv"),
                      ("EUR",    f"{V1}/gwas_imp_eur/gwas_summary.tsv", f"{V2}/assoc_eur/gwas_summary.tsv")]:
    a = {r["phenotype"]: r for r in read_tsv(v1p)}
    b = {r["phenotype"]: r for r in read_tsv(v2p)}
    if not a or not b: print(f"\n[{arm}] missing ({v1p if not a else v2p})"); continue
    print(f"\n[{arm}]  {'phenotype':20} {'n_v1':>6} {'n_v2':>6} | {'lam_v1':>7} {'lam_v2':>7} | "
          f"{'hits_v1':>7} {'hits_v2':>7} | {'minp_v1':>9} {'minp_v2':>9}")
    for p in PH:
        if p not in a or p not in b: continue
        print(f"{'':6} {p:20} {a[p]['n_mean']:>6} {b[p]['n_mean']:>6} | "
              f"{fmt(a[p]['lambda_gc'],3):>7} {fmt(b[p]['lambda_gc'],3):>7} | "
              f"{a[p]['n_p5e8']:>7} {b[p]['n_p5e8']:>7} | {sci(a[p]['min_p']):>9} {sci(b[p]['min_p']):>9}")

# ---- 2. LDSC -----------------------------------------------------------------
section("2. LDSC genetic correlation with SCZ / MDD")
for arm, v1p, v2p in [("POOLED", f"{V1}/ldsc_imp/ldsc_rg_summary.tsv", None),
                      ("EUR",    f"{V1}/ldsc_imp_eur/ldsc_rg_summary.tsv", f"{V2}/ldsc_eur/ldsc_rg_summary.tsv")]:
    a = {(r["disorder"], r["phenotype"]): r for r in read_tsv(v1p)}
    b = {(r["disorder"], r["phenotype"]): r for r in read_tsv(v2p)} if v2p else {}
    if not a: continue
    print(f"\n[{arm}]  {'disorder':5} {'phenotype':20} | {'rg_v1':>9} {'se_v1':>7} {'h2z_v1':>7} | "
          f"{'rg_v2':>9} {'se_v2':>7} {'h2z_v2':>7}")
    for d in ("SCZ", "MDD"):
        for p in PH:
            r1, r2 = a.get((d, p)), b.get((d, p))
            if not r1 and not r2: continue
            g = lambda r, k: (r.get(k) if r else None)
            print(f"{'':6} {d:5} {p:20} | {fmt(g(r1,'rg')):>9} {fmt(g(r1,'se')):>7} {fmt(g(r1,'h2_z'),2):>7} | "
                  f"{fmt(g(r2,'rg')):>9} {fmt(g(r2,'se')):>7} {fmt(g(r2,'h2_z'),2):>7}")

# ---- 3. MAGMA: AHBA components ----------------------------------------------
section("3. MAGMA gene covariates: AHBA C1-C3 (EUR arm)")
def magma_map(p):
    out = {}
    for r in read_tsv(p):
        m = re.match(r"(.+)_ahba_components_posneg$", r["comparison"])
        if m: out[(m.group(1), r["variable"])] = r
    return out
a, b = magma_map(f"{V1}/magma_imp_eur/magma_summary.tsv"), magma_map(f"{V2}/magma_eur/magma_summary.tsv")
if a:
    print(f"\n{'phenotype':20} {'set':4} | {'beta_v1':>9} {'p_v1':>9} | {'beta_v2':>9} {'p_v2':>9}")
    for p in PH:
        for s in ("C1+", "C1-", "C2+", "C2-", "C3+", "C3-"):
            r1, r2 = a.get((p, s)), b.get((p, s))
            if not r1 and not r2: continue
            g = lambda r, k: (r.get(k) if r else None)
            print(f"{p:20} {s:4} | {fmt(g(r1,'beta')):>9} {sci(g(r1,'p')):>9} | "
                  f"{fmt(g(r2,'beta')):>9} {sci(g(r2,'p')):>9}")

# ---- 4. MAGMA: SCZ / MDD gene sets ------------------------------------------
section("4. MAGMA gene sets: SCZ / MDD (EUR arm)")
for tag, sets in [("sczprio", ["SCZ_locus_pool", "SCZ_prioritised", "SCZ_prio_finemap"]),
                  ("mddhc",   ["MDD_hc_finemap", "MDD_highconf", "MDD_pool"])]:
    a = {(r["trait"], r["set"]): r for r in read_tsv(f"{V1}/magma_prio_imp_eur/{tag}_gsa_summary.tsv")
         if r.get("model") == "marginal"}
    b = {(r["trait"], r["set"]): r for r in read_tsv(f"{V2}/magma_prio_eur/{tag}_gsa_summary.tsv")
         if r.get("model") == "marginal"}
    if not a: continue
    print(f"\n[{tag}]  {'phenotype':20} {'set':18} | {'beta_v1':>9} {'p_v1':>9} | {'beta_v2':>9} {'p_v2':>9}")
    for p in PH:
        for s in sets:
            r1, r2 = a.get((p, s)), b.get((p, s))
            if not r1 and not r2: continue
            g = lambda r, k: (r.get(k) if r else None)
            print(f"{'':8} {p:20} {s:18} | {fmt(g(r1,'beta')):>9} {sci(g(r1,'p')):>9} | "
                  f"{fmt(g(r2,'beta')):>9} {sci(g(r2,'p')):>9}")

# ---- 5. PRS ------------------------------------------------------------------
section("5. PRS: v1 population association vs v2 between/within-family (pooled sample)")
pop = read_tsv(f"{V1}/prs_imp/prs_association.tsv")
fam = read_tsv(f"{V2}/prs_family/prs_withinfamily.tsv")
peak = {}
for r in pop:
    if r["stratum"] != "full": continue
    k = (r["disorder"], r["phenotype"])
    if k not in peak or float(r["p"]) < float(peak[k]["p"]): peak[k] = r
fam_i = {(r["disorder"], r["threshold"], r["phenotype"]): r for r in fam if r["stratum"] == "full"}
print(f"\n{'disorder':5} {'phenotype':20} {'thr':6} | {'beta_pop_v1':>11} {'p_v1':>8} | "
      f"{'b_between':>10} {'p_btw':>8} | {'b_within':>9} {'p_wth':>7} | {'p_diff':>7} {'pairs':>5}")
for d in ("SCZ", "MDD"):
    for p in PH:
        r1 = peak.get((d, p))
        if not r1: continue
        r2 = fam_i.get((d, r1["threshold"], p))
        g = lambda r, k: (r.get(k) if r else None)
        print(f"{d:5} {p:20} {r1['threshold']:6} | {fmt(r1['beta']):>11} {sci(r1['p']):>8} | "
              f"{fmt(g(r2,'beta_between')):>10} {sci(g(r2,'p_between')):>8} | "
              f"{fmt(g(r2,'beta_within')):>9} {sci(g(r2,'p_within')):>7} | "
              f"{sci(g(r2,'p_diff')):>7} {g(r2,'n_pairs') or 'NA':>5}")
