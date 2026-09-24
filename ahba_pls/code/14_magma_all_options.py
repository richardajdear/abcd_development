"""
14_magma_all_options.py -- one MAGMA gene-property table for every major gene
ranking in this project, on a SINGLE shared universe, plus the pairwise joint
models that say which ranking's signal survives when another is in the model.

Why re-run rather than reuse magma_gene_property.tsv and hcp_vs_dk_enrichment.tsv:
those were separate runs, each on its own universe (a vector's genes with a
MAGMA result), so their betas are not strictly comparable. Here every vector is
a column of one --gene-covar file, so all marginal and joint models share the
same genes.

Vectors (thinning orientation throughout: positive = expressed where adolescent
thinning is faster; C1-C3 and NSPN keep their published sign):
  ABCD_PLS2_HCP     option 2 PLS2 in HCP-MMP, 137 parcels          [12_hcp_pls]
  ABCD_PLS2_DK      option 2 PLS2 in DK, 33 regions, ds25          [04_fit_pls]
  ABCD_dCT_HCP      dCT alone in HCP-MMP
  ABCD_dCT_DK       dCT alone in DK (ds25)
  ABCD_dCT_dT1T2_DK option 3 PLS2 in DK (ds25)
  AHBA_C3, AHBA_C1  Dear et al. 2024
  NSPN_PLS2         Whitaker, Vertes et al. 2016 bootstrapped weights

Joint models: every ordered pair in PAIRS below, each fitted once with both
coefficients reported, so mutual attenuation is visible rather than inferred.

Outputs: results/magma_all_marginal.tsv, results/magma_all_joint.tsv
"""
from __future__ import annotations
import subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

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
RAW = {"SCZ": REPO / "genetic_analysis/inputs/magma/SCZ.genes.raw",
       "MDD": REPO / "genetic_analysis/inputs/magma/MDD.genes.raw"}

dk_w = pd.read_csv(RES / "pls_weights" / "opt2_dCT_CT_ds25.tsv", sep="\t", index_col=0)
dk_o1 = pd.read_csv(RES / "pls_weights" / "opt1_dCT_ds25.tsv", sep="\t", index_col=0)
dk_o3 = pd.read_csv(RES / "pls_weights" / "opt3_dCT_dT1T2_ds25.tsv", sep="\t", index_col=0)
hcp_w = pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)
c123 = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)
nspn = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")

# The HCP and AHBA_updated-DK vectors are fitted on DIFFERENT abagen builds and
# different DS gene sets, so a joint model of the two confounds parcellation
# with pipeline. ABCD_PLS2_DKmatched is the gene- and build-matched DK fit
# (dk_3d.csv, the same 7,973 genes as the HCP matrix; written by 12_hcp_pls.py)
# and is the vector the atlas pair below uses. ABCD_PLS2_DK stays in the
# marginal comparison because it is the main analysis's own signature.
dk_m = pd.read_csv(RES / "dk_matched_weights.tsv", sep="\t", index_col=0)
vectors = {
    "ABCD_PLS2_HCP": -hcp_w["hcp_opt2_dCT_CT_PLS2_Z"],
    "ABCD_PLS2_DKmatched": dk_m["DK_PLS2_matchedX"],
    "ABCD_PLS2_DK": -dk_w["PLS2_Z"],
    "ABCD_dCT_HCP": -hcp_w["hcp_opt1_dCT_PLS1_Z"],
    "ABCD_dCT_DK": -dk_o1["PLS1_Z"],
    "ABCD_dCT_dT1T2_DK": -dk_o3["PLS2_Z"],
    "AHBA_C3": c123["C3"],
    "AHBA_C1": c123["C1"],
    "NSPN_PLS2": nspn["PLS2_z"],
}
PAIRS = [("ABCD_PLS2_HCP", "ABCD_PLS2_DKmatched"),   # does the finer parcellation add?
         #  ^ gene- and build-matched, so this pair isolates parcellation
         ("ABCD_PLS2_HCP", "AHBA_C3"),          # does either ABCD vector add to C3?
         ("ABCD_PLS2_DK", "AHBA_C3"),
         ("ABCD_PLS2_HCP", "AHBA_C1"),          # is it just the static gradient?
         ("ABCD_PLS2_HCP", "NSPN_PLS2"),        # vs the original imaging-transcriptomic map
         ("ABCD_PLS2_DK", "NSPN_PLS2"),
         ("NSPN_PLS2", "AHBA_C3")]

sym2ent = (pd.read_csv(GS / "magma_SCZ_genes.tsv", sep="\t").dropna(subset=["symbol"])
           .drop_duplicates("symbol").set_index("symbol")["GENE"].astype(int))
cov = pd.DataFrame(vectors).dropna(how="any")
cov = cov.loc[cov.index.intersection(sym2ent.index)]
cov.index = sym2ent.loc[cov.index].values
cov = cov[~cov.index.duplicated()]
cov.index.name = "GENE"
covar = RUN / "all_options.covar"
cov.to_csv(covar, sep="\t", float_format="%.6g")
print(f"shared universe: {len(cov)} genes x {cov.shape[1]} vectors", file=sys.stderr)
assert len(cov) > 4000, len(cov)


def run(tag: str, cols: list[str], model: list[str]) -> pd.DataFrame:
    sub = RUN / f"all_{tag}.covar"
    cov[cols].to_csv(sub, sep="\t", float_format="%.6g")
    out = []
    for dis, raw in RAW.items():
        prefix = RUN / f"all_{tag}_{dis}"
        r = subprocess.run([str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(sub),
                            "--out", str(prefix), *model], capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"MAGMA failed ({tag} {dis}):\n{r.stdout[-1500:]}")
        t = pd.read_csv(f"{prefix}.gsa.out", sep=r"\s+", comment="#")
        if "FULL_NAME" in t.columns:
            t["VARIABLE"] = t["FULL_NAME"].fillna(t["VARIABLE"])
        t.insert(0, "disorder", dis); t.insert(0, "run", tag)
        out.append(t)
    return pd.concat(out, ignore_index=True)


# ---- marginal: one model per vector, alone --------------------------------
marg = pd.concat([run(f"marg_{v}", [v], []) for v in vectors], ignore_index=True)
marg["se_std"] = marg.SE * marg.BETA_STD / marg.BETA
marg.to_csv(RES / "magma_all_marginal.tsv", sep="\t", index=False, float_format="%.5g")

# ---- joint: both coefficients from one model ------------------------------
# With `--model condition=B` on a two-column covar file MAGMA fits ONE model
# containing both covariates and prints both coefficients (MODEL 1), so the
# pair's mutual attenuation is read off directly.
joint = []
for a, b in PAIRS:
    t = run(f"joint_{a}__{b}", [a, b], ["--model", f"condition={b}"])
    t = t[t.MODEL == 1].copy()
    t["pair"] = f"{a} + {b}"
    t["partner"] = np.where(t.VARIABLE == a, b, a)
    joint.append(t)
J = pd.concat(joint, ignore_index=True)
J["se_std"] = J.SE * J.BETA_STD / J.BETA
J.to_csv(RES / "magma_all_joint.tsv", sep="\t", index=False, float_format="%.5g")

pd.set_option("display.width", 220)
print(f"\nMARGINAL (shared universe, n = {int(marg.NGENES.iloc[0])} genes):")
print(marg.pivot(index="VARIABLE", columns="disorder", values=["BETA_STD", "P"])
      .loc[list(vectors)].round(4).to_string())
print("\nJOINT (both coefficients of each two-vector model):")
print(J[["pair", "disorder", "VARIABLE", "BETA_STD", "se_std", "P"]].round(4).to_string(index=False))
