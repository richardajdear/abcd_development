"""
08_magma_gene_property.py -- H2, test B: MAGMA gene-property regression of
disorder gene-level Z (SCZ: PGC3 2022; MDD: Adams 2025) on continuous gene
weights, using the gene-gene correlation structure in the .genes.raw files.

Marginal models: one --gene-covar file per weight vector (its own universe =
genes with a weight AND a MAGMA gene result), covariates = thinning-oriented Z
and a top-decile indicator.  Default MAGMA gene-property test is two-sided.

Conditional models (ds25 universe, genes with ALL vectors): lead | C1,
lead | static control, C3 | lead, lead | C3, NSPN_PLS2 | lead -- to ask whether
the thinning signature carries disorder signal beyond the static axis, and
whether it accounts for the C3 result.

Sign convention: positive BETA = genes with higher thinning-oriented weight
(expressed more where adolescent thinning is faster) have larger disorder Z.

Output: results/magma_gene_property.tsv (marginal), results/magma_conditional.tsv,
        results/magma_runs/<...>.gsa.out and logs
"""
from __future__ import annotations
import subprocess, sys, shutil
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parent
RES, REF = ROOT / "results", ROOT / "data" / "reference"
GS = REF / "gene_sets"
RUN = RES / "magma_runs"; RUN.mkdir(exist_ok=True)
def _magma_bin(repo):
    # native arm64 build first (tools/bin/magma_src, compiled from the v1.10
    # source), then the x86_64 macOS build (needs Rosetta), then the Linux one
    for c in ("magma_src/magma", "magma_mac/magma", "magma"):
        if (repo / "tools" / "bin" / c).exists():
            return repo / "tools" / "bin" / c
    raise FileNotFoundError("no MAGMA binary under tools/bin; see tools/bin/README.md")
MAGMA = _magma_bin(REPO)
RAW = {"SCZ": REPO / "genetic_analysis/inputs/magma/SCZ.genes.raw", "MDD": REPO / "genetic_analysis/inputs/magma/MDD.genes.raw"}
for p in RAW.values():
    assert p.exists(), p

sym2ent = (pd.read_csv(GS / "magma_SCZ_genes.tsv", sep="\t").dropna(subset=["symbol"])
           .drop_duplicates("symbol").set_index("symbol")["GENE"].astype(int))

def load_Z(opt, ds, comp, sign=-1):
    t = pd.read_csv(RES / "pls_weights" / f"{opt}_{ds}.tsv", sep="\t", index_col=0)
    return sign * t[f"{comp}_Z"].dropna()

lead25 = load_Z("opt2_dCT_CT", "ds25", "PLS2")
opt5 = load_Z("opt5_slopePCs", "ds25", "PLS2", sign=1)
sh = opt5.index.intersection(lead25.index)
if stats.spearmanr(opt5.loc[sh], lead25.loc[sh]).statistic < 0:
    opt5 = -opt5
nspn = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")
c123 = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)
vectors = {
    "lead_opt2_PLS2_ds0": load_Z("opt2_dCT_CT", "ds0", "PLS2"),
    "lead_opt2_PLS2_ds25": lead25,
    "lead_opt2_PLS2_ds50": load_Z("opt2_dCT_CT", "ds50", "PLS2"),
    "opt1_dCT_PLS1_ds25": load_Z("opt1_dCT", "ds25", "PLS1"),
    "opt3_dCT_dT1T2_PLS2_ds25": load_Z("opt3_dCT_dT1T2", "ds25", "PLS2"),
    "opt4_full4_PLS2_ds25": load_Z("opt4_full4", "ds25", "PLS2"),
    "opt5_slopePCs_PLS2_ds25": opt5,
    "control_opt2_PLS1_static_ds25": load_Z("opt2_dCT_CT", "ds25", "PLS1"),
    "AHBA_C1": c123["C1"], "AHBA_C2": c123["C2"], "AHBA_C3": c123["C3"],
    "NSPN_PLS1_z": nspn["PLS1_z"].dropna(), "NSPN_PLS2_z": nspn["PLS2_z"].dropna(),
}

def covar_frame(cols: dict[str, pd.Series], how="outer") -> pd.DataFrame:
    df = pd.concat(cols, axis=1, join=how)
    df = df.loc[df.index.intersection(sym2ent.index)]
    df.index = sym2ent.loc[df.index].values
    df = df[~pd.Index(df.index).duplicated()]
    df.index.name = "GENE"
    return df

def run_magma(covar: pd.DataFrame, tag: str, extra: list[str]) -> pd.DataFrame:
    out = []
    for dis, raw in RAW.items():
        f = RUN / f"{tag}.covar.txt"
        covar.to_csv(f, sep="\t", na_rep="NA", float_format="%.6g")
        prefix = RUN / f"{tag}.{dis}"
        cmd = [str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(f), *extra, "--out", str(prefix)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 or not Path(f"{prefix}.gsa.out").exists():
            print(r.stdout[-2000:], r.stderr[-2000:], file=sys.stderr)
            raise RuntimeError(f"MAGMA failed: {tag} {dis}")
        t = pd.read_csv(f"{prefix}.gsa.out", sep=r"\s+", comment="#")
        if "FULL_NAME" in t.columns:   # MAGMA truncates long names into VARIABLE
            t["VARIABLE"] = t["FULL_NAME"].fillna(t["VARIABLE"])
        t.insert(0, "disorder", dis); t.insert(0, "run", tag)
        out.append(t)
    return pd.concat(out, ignore_index=True)

# ---- marginal ------------------------------------------------------------------
marg = []
for name, w in vectors.items():
    w = w.dropna()
    cov = covar_frame({name: w, f"{name}_top10": (w >= w.quantile(0.9)).astype(int)})
    marg.append(run_magma(cov, f"marginal_{name}", []))
    print(name, len(cov), file=sys.stderr)
M = pd.concat(marg, ignore_index=True)
M.to_csv(RES / "magma_gene_property.tsv", sep="\t", index=False, float_format="%.5g")

# ---- conditional (ds25 universe, genes present in all vectors) ------------------
base = {"lead": lead25, "static": vectors["control_opt2_PLS1_static_ds25"],
        "C1": c123["C1"], "C3": c123["C3"], "NSPN_PLS2": nspn["PLS2_z"].dropna()}
cov = covar_frame(base, how="inner")
conds = {"none": [], "cond_C1": ["--model", "condition=C1"], "cond_static": ["--model", "condition=static"],
         "cond_lead": ["--model", "condition=lead"], "cond_C3": ["--model", "condition=C3"]}
C = pd.concat([run_magma(cov, f"conditional_{k}", v) for k, v in conds.items()], ignore_index=True)
C.to_csv(RES / "magma_conditional.tsv", sep="\t", index=False, float_format="%.5g")

pd.set_option("display.width", 250)
Mm = M[~M.VARIABLE.str.endswith("_top10")]
print("MARGINAL (continuous Z):")
print(Mm.pivot(index="VARIABLE", columns="disorder", values=["NGENES", "BETA_STD", "P"]).round(4).to_string())
print("\nMARGINAL (top-decile indicator):")
print(M[M.VARIABLE.str.endswith("_top10")].pivot(index="VARIABLE", columns="disorder", values=["BETA", "P"]).round(4).to_string())
print(f"\nCONDITIONAL (n genes = {len(cov)}):")
Cc = C[(C.run == "conditional_none") | ((C.run.str.replace("conditional_cond_", "") != C.VARIABLE))]
print(Cc[["run", "disorder", "VARIABLE", "BETA_STD", "P"]].round(4).to_string(index=False))
