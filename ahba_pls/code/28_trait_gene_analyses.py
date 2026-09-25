"""
28_trait_gene_analyses.py -- MAGMA gene analyses (SNP -> gene) for the extra
traits of the summary / mechanism slides, run LOCALLY on the 1000G EUR panel.

Why local: none of these traits had a gene analysis under
genetic_analysis/inputs/magma/, the public summary statistics are small enough
to download, and the gene analysis is single-threaded MAGMA (~20-40 min/trait).

Traits (European-ancestry GWAS, to match the g1000_eur LD panel):
  BIP     O'Connell et al. 2025, PGC bip2024_eur_no23andMe        figshare 27216117
  ADHD    Demontis et al. 2023, ADHD2022_iPSYCH_deCODE_PGC.meta   figshare 22564390
  ALZ     Bellenguez et al. 2022, GCST90027158 (already local)
  EA      Okbay et al. 2016 (EA2), EduYears_Main, excl. 23andMe   GWAS Catalog GCST003676
  INT     Savage, Jansen et al. 2018 intelligence meta-analysis   CTG (vu.data.surf.nl)
  HEIGHT  Yengo et al. 2022, GIANT EUR (HapMap3 SNPs)             non-brain negative control
EA3 (Lee 2018) is not used: its public file is a Dropbox link that did not resolve
here, and EA4 needs an SSGAC data-use agreement. EA2 is the EA GWAS the cluster
PRS analysis also used.

Sample size: MAGMA's --pval model needs N per SNP. Case-control traits use
Ncases + Ncontrols per SNP (BIP, ADHD, ALZ); EA2's file carries no N, so the
published N of the EduYears_Main meta-analysis (328,917) is used for every SNP.

Settings match the ones the cluster used for SCZ/MDD (step 10 of
genetic_analysis/README_HPC.md): NCBI37.3 genes, window 35 kb up / 10 kb down,
g1000_eur; the local annotation is ~/Git/AHBA/magma/snp_gene_annotation/
ncbi37.window35-10.genes.annot (the cluster's `union3` annotation differs only
in how multi-mapping SNPs are merged).

Inputs:  ~/Git/AHBA/magma/gwas/<file>          (downloaded; not in any repo)
Outputs: genetic_analysis/inputs/magma/<TRAIT>.genes.{raw,out}   (gitignored inputs)
         ahba_pls/results/trait_gene_analyses.tsv                 (provenance: N genes, SNPs used)
Usage:   python code/28_trait_gene_analyses.py [TRAIT ...]        (default: all, one process per trait)
"""
import gzip, subprocess, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import magma_utils  # noqa: E402

AHBA = Path.home() / "Git" / "AHBA" / "magma"
GWAS, REF = AHBA / "gwas", AHBA / "reference_data" / "g1000_eur"
ANNOT = AHBA / "snp_gene_annotation" / "ncbi37.window35-10.genes.annot"
OUT = REPO / "genetic_analysis" / "inputs" / "magma"
TMP = REPO / "ahba_pls" / "results" / "magma_runs" / "trait_prep"
RES = REPO / "ahba_pls" / "results"

# name -> (file, separator, SNP col, P col, N spec: column name, list of columns to sum, or a constant)
TRAITS = {
    "BIP":    ("bip2024_eur_no23andMe.gz", r"\s+", "SNP", "P", ["Nca", "Nco"]),
    "ADHD":   ("ADHD2022_iPSYCH_deCODE_PGC.meta.gz", r"\s+", "SNP", "P", ["Nca", "Nco"]),
    "ALZ":    ("GCST90027158_buildGRCh38.tsv", "\t", "variant_id", "p_value", ["n_cases", "n_controls"]),
    "EA":     ("Okbay_27225129-EduYears_Main.txt.gz", "\t", "MarkerName", "Pval", 328917),
    "INT":    ("SavageJansen_2018_intelligence_metaanalysis.txt", "\t", "SNP", "P", "N_analyzed"),
    "HEIGHT": ("GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz", "\t", "RSID", "P", "N"),
}


def prep(name: str) -> Path:
    f, sep, snp, p, n = TRAITS[name]
    cols = [snp, p] + (n if isinstance(n, list) else [n] if isinstance(n, str) else [])
    parts = []
    for ch in pd.read_csv(GWAS / f, sep=sep, usecols=cols, chunksize=2_000_000, engine="c"):
        N = ch[n].sum(axis=1) if isinstance(n, list) else ch[n] if isinstance(n, str) else n
        # MAGMA rejects denormal p (5e-324 prints as "not a number"); clip at 1e-300
        d = pd.DataFrame({"SNP": ch[snp], "P": ch[p].clip(lower=1e-300), "N": N})
        parts.append(d[d.SNP.astype(str).str.startswith("rs") & d.P.between(0, 1, inclusive="right")])
    d = pd.concat(parts).drop_duplicates("SNP")
    TMP.mkdir(parents=True, exist_ok=True)
    out = TMP / f"{name}.pval.tsv"
    d.to_csv(out, sep="\t", index=False)
    return out


def run(name: str) -> dict:
    pv = prep(name)
    OUT.mkdir(parents=True, exist_ok=True)
    cmd = [str(magma_utils.magma_bin()), "--bfile", str(REF), "--gene-annot", str(ANNOT),
           "--pval", str(pv), "use=SNP,P", "ncol=N", "--out", str(OUT / name)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    log = (OUT / f"{name}.log").read_text() if (OUT / f"{name}.log").exists() else r.stdout
    if r.returncode != 0 or not (OUT / f"{name}.genes.raw").exists():
        raise RuntimeError(f"{name}: MAGMA failed\n{log[-1500:]}")
    g = pd.read_csv(OUT / f"{name}.genes.out", sep=r"\s+")
    used = [l for l in log.splitlines() if "SNPs in file" in l or "remaining" in l.lower() or "synonymous" in l.lower()]
    return dict(trait=name, source=TRAITS[name][0], n_snps_input=sum(1 for _ in open(pv)) - 1,
                n_genes=len(g), n_genes_bonf=int((g.P < 0.05 / len(g)).sum()),
                median_nsnps=float(g.NSNPS.median()), log_note=" | ".join(used)[:300])


if __name__ == "__main__":
    # one trait per process: a process pool is not permitted in some sandboxes, so
    # with several traits this re-invokes itself once per trait and collects the rows
    names = sys.argv[1:] or list(TRAITS)
    if len(names) > 1:
        procs = [subprocess.Popen([sys.executable, __file__, n]) for n in names]
        bad = [n for n, pr in zip(names, procs) if pr.wait() != 0]
        R = pd.read_csv(RES / "trait_gene_analyses.tsv", sep="\t")
        print(R.drop(columns="log_note").to_string(index=False))
        sys.exit(f"failed: {bad}" if bad else 0)
    rows = [run(names[0])]
    R = pd.DataFrame(rows)
    f = RES / "trait_gene_analyses.tsv"
    import fcntl
    lock = open(RES / ".trait_gene_analyses.lock", "w"); fcntl.flock(lock, fcntl.LOCK_EX)
    if f.exists():
        R = pd.concat([pd.read_csv(f, sep="\t").query("trait not in @names"), R], ignore_index=True)
    R.to_csv(f, sep="\t", index=False)
    print(R.drop(columns="log_note").to_string(index=False))
