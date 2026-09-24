"""
21_magma_disorder_panel.py -- MAGMA gene-property test of the thinning signature
(option 2, dCT + CT, PLS2) against EVERY available disorder gene analysis,
including the 2025 multi-ancestry SCZ GWAS.

Scripts 08 and 14 tested against the two gene analyses shipped in
genetic_analysis/inputs/magma/ -- SCZ.genes.raw = PGC3 *primary* and
MDD.genes.raw = MDD2025 *div*. Both are multi-ancestry GWAS whose gene analysis
was run against the 1000G EUR panel only, the LD mismatch step 10 of
genetic_analysis/README_HPC.md documents. Step 10/11 on CSD3 produced the
LD-matched replacements; this script uses whichever of them are present.

Disorder gene analyses (genetic_analysis/inputs/magma/<name>.genes.raw):
  SCZ25_META    2025 SCZ (AFR+EUR+EAS), per-ancestry gene analysis on the
                matched 1000G panel, then `magma --meta`         PRIMARY SCZ
  SCZ25_EUR     2025 SCZ, European GWAS x g1000_eur              LD-matched EUR
  PGC3_EUR      PGC3 european x g1000_eur                        old-GWAS comparator
  PGC3_primary  = SCZ.genes.raw (PGC3 primary x g1000_eur)       what 08/14 used
  MDD_EUR       MDD2025 eur x g1000_eur                          PRIMARY MDD
  MDD_div       = MDD.genes.raw (MDD2025 div x g1000_eur)        what 08/14 used
and, as a check of H3 (does the signature carry the genetics of thinning
itself?), the ABCD GENESIS GWAS of global slope, EUR arm x g1000_eur:
  global_slope_1lmm_{hcp,dk}       single-LMM whole-cortex slope (primary spec)
  global_slope_perregion_{hcp,dk}  mean of per-region slope BLUPs
Those GWAS have low SNP heritability (LDSC h2 z < 2 at n ~ 4,300 EUR), so a
null there is expected and uninformative; they are run because it costs nothing.
Missing files are reported and skipped; the table says which were run.

Gene-property analysis needs no LD reference of its own: the gene-gene
correlations it corrects for are carried inside each .genes.raw, computed at the
gene-analysis step from that analysis's panel. So the LD choice is made once per
disorder, upstream, and is recorded in the table's `ld` column.

Vectors (thinning orientation; positive = expressed where thinning is faster):
  ABCD_PLS2_DK         the lead signature, DK, AHBA_updated ds25
  ABCD_PLS2_HCP        the same design at 137 HCP-MMP parcels
  ABCD_PLS2_DKmatched  DK on the HCP gene basis (the parcellation-only contrast)
  ABCD_PLS1_HCP        PLS1 of the HCP fit, the static component, aligned with C1 (positive = thinner)
  AHBA_C3, NSPN_PLS2   the published benchmarks;  AHBA_C1  static-gradient control

Each vector is tested alone on its own universe (genes with both a weight and a
gene result), and the lead signature jointly with C3 (does it add to C3?) and
with C1 (is it the static gradient?) on the pair's shared universe. MAGMA's
default gene-property model is used, as in 08 and 14; for gene properties that
default is TWO-SIDED (manual v1.10, --model direction), so every p here is
two-sided and a negative coefficient can be significant.

Outputs
  results/magma_disorder_panel.tsv        marginal
  results/magma_disorder_panel_joint.tsv  joint models (both coefficients)
"""
from __future__ import annotations
import subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parent
RES, REF = ROOT / "results", ROOT / "data" / "reference"
RUN = RES / "magma_runs" / "disorder_panel"; RUN.mkdir(parents=True, exist_ok=True)
def _magma_bin(repo):
    # native arm64 build first (tools/bin/magma_src, compiled from the v1.10
    # source), then the x86_64 macOS build (needs Rosetta), then the Linux one
    for c in ("magma_src/magma", "magma_mac/magma", "magma"):
        if (repo / "tools" / "bin" / c).exists():
            return repo / "tools" / "bin" / c
    raise FileNotFoundError("no MAGMA binary under tools/bin; see tools/bin/README.md")
MAGMA = _magma_bin(REPO)
INP = REPO / "genetic_analysis" / "inputs" / "magma"

DISORDERS = {   # name -> (file stem, disorder, ld reference, role)
    "SCZ25_META":   ("SCZ25_META",   "SCZ", "g1000 EUR/AFR/EAS, per-ancestry + --meta", "primary"),
    "SCZ25_EUR":    ("SCZ25_EUR",    "SCZ", "g1000_eur (matched)",                      "EUR"),
    "PGC3_EUR":     ("PGC3_EUR",     "SCZ", "g1000_eur (matched)",                      "old GWAS, EUR"),
    "PGC3_primary": ("SCZ",          "SCZ", "g1000_eur (mismatched)",                   "old GWAS, as in 08/14"),
    "MDD_EUR":      ("MDD_EUR",      "MDD", "g1000_eur (matched)",                      "primary"),
    "MDD_div":      ("MDD",          "MDD", "g1000_eur (mismatched)",                   "as in 08/14"),
    "slope_1lmm_HCP":      ("global_slope_1lmm_hcp",      "ABCD slope", "g1000_eur (EUR arm)", "primary spec"),
    "slope_1lmm_DK":       ("global_slope_1lmm_dk",       "ABCD slope", "g1000_eur (EUR arm)", "DK"),
    "slope_perregion_HCP": ("global_slope_perregion_hcp", "ABCD slope", "g1000_eur (EUR arm)", "per-region BLUP"),
    "slope_perregion_DK":  ("global_slope_perregion_dk",  "ABCD slope", "g1000_eur (EUR arm)", "per-region BLUP, DK"),
}
present = {k: v for k, v in DISORDERS.items() if (INP / f"{v[0]}.genes.raw").exists()}
missing = sorted(set(DISORDERS) - set(present))
print("present:", list(present), "| missing:", missing, file=sys.stderr)
assert present, f"no gene analyses in {INP}"

# ---- vectors ----------------------------------------------------------------
dk_w = pd.read_csv(RES / "pls_weights" / "opt2_dCT_CT_ds25.tsv", sep="\t", index_col=0)
hcp_w = pd.read_csv(RES / "hcp_pls_weights.tsv", sep="\t", index_col=0)
dk_m = pd.read_csv(RES / "dk_matched_weights.tsv", sep="\t", index_col=0)
c123 = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)
nspn = pd.read_csv(REF / "nspn_pls_gene_weights.csv").set_index("gene")
vectors = {
    "ABCD_PLS2_DK": -dk_w["PLS2_Z"],
    "ABCD_PLS2_HCP": -hcp_w["hcp_opt2_dCT_CT_PLS2_Z"],
    # PLS1 of the same fit is the static (baseline-thickness) component: CT
    # salience +0.99; flipped to align with AHBA C1, so positive = expressed
    # where cortex is THINNER
    "ABCD_PLS1_HCP": -hcp_w["hcp_opt2_dCT_CT_PLS1_Z"],
    "ABCD_PLS2_DKmatched": dk_m["DK_PLS2_matchedX"],
    "AHBA_C3": c123["C3"], "NSPN_PLS2": nspn["PLS2_z"], "AHBA_C1": c123["C1"],
}
JOINT = [("ABCD_PLS2_DK", "AHBA_C3"), ("ABCD_PLS2_HCP", "AHBA_C3"), ("ABCD_PLS1_HCP", "AHBA_C1"),
         ("ABCD_PLS2_DK", "AHBA_C1"), ("ABCD_PLS2_HCP", "ABCD_PLS2_DKmatched")]

# symbol -> Entrez from the gene-location file the gene analyses were annotated
# with (NCBI37.3), so every disorder file maps through the same table. Symbols
# that occur twice are dropped rather than guessed.
loc = pd.read_csv(REPO / "tools" / "bin" / "NCBI37.3.gene.loc", sep="\t", header=None,
                  names=["GENE", "chr", "start", "end", "strand", "symbol"])
loc = loc[~loc.symbol.duplicated(keep=False)]
sym2ent = loc.set_index("symbol").GENE.astype(int)


def covar(cols: list[str], tag: str) -> Path:
    c = pd.DataFrame({k: vectors[k] for k in cols}).dropna(how="any")
    c = c.loc[c.index.intersection(sym2ent.index)]
    c.index = sym2ent.loc[c.index].values
    c = c[~c.index.duplicated()]; c.index.name = "GENE"
    p = RUN / f"{tag}.covar"
    c.to_csv(p, sep="\t", float_format="%.6g")
    return p


def magma(raw: Path, cov: Path, out: Path, extra: list[str]) -> pd.DataFrame:
    r = subprocess.run([str(MAGMA), "--gene-results", str(raw), "--gene-covar", str(cov),
                        "--out", str(out), *extra], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"MAGMA failed for {out.name}:\n{r.stdout[-1500:]}")
    t = pd.read_csv(f"{out}.gsa.out", sep=r"\s+", comment="#")
    if "FULL_NAME" in t.columns:
        t["VARIABLE"] = t["FULL_NAME"].fillna(t["VARIABLE"])
    return t


rows, jrows = [], []
for name, (stem, dis, ld, role) in present.items():
    raw = INP / f"{stem}.genes.raw"
    meta = dict(gene_analysis=name, disorder=dis, ld=ld, role=role)
    for v in vectors:
        t = magma(raw, covar([v], f"m_{v}"), RUN / f"{name}__{v}", [])
        rows.append({**meta, **t.iloc[0][["VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P"]].to_dict()})
    for a, b in JOINT:
        t = magma(raw, covar([a, b], f"j_{a}__{b}"), RUN / f"{name}__{a}__{b}",
                  ["--model", f"condition={b}"])
        t = t[t.MODEL == 1] if "MODEL" in t.columns else t
        for _, r_ in t.iterrows():
            jrows.append({**meta, "pair": f"{a} + {b}",
                          **r_[["VARIABLE", "NGENES", "BETA", "BETA_STD", "SE", "P"]].to_dict()})

M = pd.DataFrame(rows); J = pd.DataFrame(jrows)
for T in (M, J):
    T["se_std"] = T.SE * T.BETA_STD / T.BETA
M.to_csv(RES / "magma_disorder_panel.tsv", sep="\t", index=False, float_format="%.5g")
J.to_csv(RES / "magma_disorder_panel_joint.tsv", sep="\t", index=False, float_format="%.5g")

pd.set_option("display.width", 220)
print("\nMARGINAL  beta_std (p), n genes")
M["cell"] = M.apply(lambda r: f"{r.BETA_STD:+.3f} ({r.P:.2g})", axis=1)
print(M.pivot(index="VARIABLE", columns="gene_analysis", values="cell")
      .loc[list(vectors), list(present)].to_string())
print("\nn genes:", M.groupby("gene_analysis").NGENES.max().to_dict())
print("\nJOINT  (both coefficients of each pair)")
J["cell"] = J.apply(lambda r: f"{r.BETA_STD:+.3f} ({r.P:.2g})", axis=1)
print(J.pivot_table(index=["pair", "VARIABLE"], columns="gene_analysis", values="cell",
                    aggfunc="first")[list(present)].to_string())
if missing:
    print(f"\nNOT RUN (file absent from {INP.relative_to(REPO)}): {missing}")
