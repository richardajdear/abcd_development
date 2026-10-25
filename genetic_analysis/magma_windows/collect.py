"""MAGMA window robustness, step 3: reproduction and per-window/trait checks -> work/magma_windows/README.md."""
import subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

REPO = Path("/home/rajd2/rds/hpc-work/abcd_development")
W = REPO / "genetic_analysis/work/magma_windows"
GA = REPO / "genetic_analysis/work"
TAGS = ["w0_0", "w10_10", "w35_10", "w50_50", "w100_100"]
TRAITS = ["ABCD_thin_EUR", "ABCD_base_EUR", "SCZ25_EUR", "MDD_EUR"]
OLD = {"ABCD_thin_EUR": GA / "results_70tab_hcp/magma_panel/genes/global_slope_1lmm.genes.out",
       "ABCD_base_EUR": GA / "results_70tab_hcp/magma_panel/genes/baseline_thickness_1lmm.genes.out",
       "SCZ25_EUR": GA / "magma_scz2025/genes/SCZ25_EUR.genes.out",
       "MDD_EUR": GA / "magma_scz2025/genes/MDD_EUR.genes.out"}
loc = pd.read_csv("/rds/user/rajd2/hpc-work/magma/gene_locations/NCBI37.3.gene.loc", sep=r"\s+", header=None)
SYM = dict(zip(loc[0], loc[5]))

def read_out(p):
    return pd.read_csv(p, sep=r"\s+").set_index("GENE")

def lam(p):
    return np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)

def annot_stats(tag):
    f = W / f"ref/ncbi37.{tag}.union3.genes.annot"
    nsnp = []
    with open(f) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            nsnp.append(line.count("\t") - 1)
    nsnp = np.array(nsnp)
    # genes per SNP over all annotated SNPs
    awk = "!/^#/{for(i=3;i<=NF;i++) c[$i]++} END{for(k in c) print c[k]}"
    out = subprocess.run(["awk", awk, str(f)], capture_output=True, text=True, check=True).stdout
    gps = np.fromstring(out, sep="\n", dtype=int)
    return dict(window=tag, genes=len(nsnp), snps_annot=len(gps), med_snps_per_gene=np.median(nsnp),
                p90_snps_per_gene=np.percentile(nsnp, 90), med_genes_per_snp=np.median(gps),
                mean_genes_per_snp=round(gps.mean(), 3), frac_snps_multi_gene=round((gps > 1).mean(), 3))

def md(df):
    f = lambda v: f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    rows = [list(df.columns), ["---"] * df.shape[1]] + [[f(v) for v in r] for r in df.itertuples(index=False)]
    return "\n".join("| " + " | ".join(r) + " |" for r in rows)

lines = ["# MAGMA gene analysis under five SNP-to-gene windows", "",
         "Task: `transcriptional_maturation/gwas_cells/csd3/MAGMA_WINDOWS_INSTRUCTIONS.md`. "
         "Scripts: `genetic_analysis/magma_windows/` (common.sh, annotate.sbatch, genes.sbatch, collect.py).", ""]
ver = subprocess.run([str(REPO / "legacy/hpc/work/bin/magma"), "--version"], capture_output=True, text=True)
lines += [f"MAGMA: `{(ver.stdout + ver.stderr).strip().splitlines()[0] if (ver.stdout + ver.stderr).strip() else 'n/a'}`. "
          "LD ref g1000_eur, NCBI37.3 gene loc, union-3 snploc, default gene model (snp-wise=mean).", ""]

# reproduction
old_ann = GA / "magma_scz2025/ref/ncbi37.window35-10.union3.genes.annot"
same = subprocess.run(["cmp", "-s", str(old_ann), str(W / "ref/ncbi37.w35_10.union3.genes.annot")]).returncode == 0
rep = []
for t in TRAITS:
    a, b = read_out(OLD[t]), read_out(W / f"genes/w35_10/{t}.genes.out")
    g = a.index.intersection(b.index); dz = (a.loc[g, "ZSTAT"] - b.loc[g, "ZSTAT"]).abs()
    rep.append(dict(trait=t, n_old=len(a), n_new=len(b), n_common=len(g),
                    pearson=stats.pearsonr(a.loc[g, "ZSTAT"], b.loc[g, "ZSTAT"])[0],
                    spearman=stats.spearmanr(a.loc[g, "ZSTAT"], b.loc[g, "ZSTAT"])[0], max_abs_dZ=dz.max()))
lines += ["## Reproduction (w35_10 vs existing)", "",
          f"- Annotation identical to `magma_scz2025/ref/ncbi37.window35-10.union3.genes.annot`: **{'yes' if same else 'NO'}** (cmp).", "",
          md(pd.DataFrame(rep)), ""]

lines += ["## Per window (annotation)", "", md(pd.DataFrame([annot_stats(t) for t in TAGS])), ""]

ref = {t: read_out(W / f"genes/w35_10/{t}.genes.out") for t in TRAITS}
rows, tops = [], []
for tag in TAGS + ["w35_10_multi"]:
    for t in TRAITS:
        f = W / f"genes/{tag}/{t}.genes.out"
        if not f.exists():
            continue
        d = read_out(f); r = ref[t]; g = d.index.intersection(r.index)
        top_d = set(d.loc[g, "ZSTAT"].nlargest(1000).index); top_r = set(r.loc[g, "ZSTAT"].nlargest(1000).index)
        p = d["P_MULTI"] if "P_MULTI" in d else d["P"]
        z = d["ZSTAT"] if "ZSTAT" in d else stats.norm.isf(p)
        rows.append(dict(window=tag, trait=t, genes_tested=len(d), lambdaGC=lam(p), n_P_lt_2_5e6=int((p < 2.5e-6).sum()),
                         spearman_vs_w35_10=stats.spearmanr(z.loc[g], r.loc[g, "ZSTAT"])[0],
                         top1000_overlap_vs_w35_10=len(top_d & top_r)))
        top = z.nlargest(10)
        tops.append(dict(window=tag, trait=t, top10=", ".join(f"{SYM.get(i, i)} ({v:.2f})" for i, v in top.items())))
lines += ["## Per window x trait", "", "λGC = median χ²₁(P) / 0.4549. Overlap and Spearman over genes tested in both windows. "
          "For w35_10_multi, P = P_MULTI.", "", md(pd.DataFrame(rows)), "", "### Top 10 genes by ZSTAT", "", md(pd.DataFrame(tops)), ""]
lines += ["## Files", "", "`genes/<window>/<trait>.genes.{out,raw}` and `<trait>.log` (MAGMA names the log `<out>.log`). "
          "Annotations in `ref/` (not copied back; large).", "",
          "Rsync to local (run on the local machine):", "", "```bash",
          "rsync -av --include='*/' --include='*.genes.out' --include='*.genes.raw' --include='*.log' --include='README.md' "
          "--exclude='*' rajd2@login.hpc.cam.ac.uk:/home/rajd2/rds/hpc-work/abcd_development/genetic_analysis/work/magma_windows/ "
          "~/Git/transcriptional_maturation/gwas_cells/work/magma_windows/ --prune-empty-dirs", "```", ""]
(W / "README.md").write_text("\n".join(lines))
print("\n".join(lines))
