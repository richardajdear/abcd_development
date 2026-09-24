"""Figure 1 panel f input: PRS -> CT / ΔCT on the single-LMM HCP phenotypes.

Usage (repo root):  python genetic_analysis/fig1_prep_prs.py
Writes genetic_analysis/fig1_inputs/hcp70_prs_key_arms.tsv (coefficients only).

Reads every order-of-operations table under results_70tab_hcp/ (the per-region
vs single-LMM comparison each PRS step writes):
    prs_final_1lmm/table_order_of_operations_all.tsv      SCZ/MDD/ASD/ALZ/EA
    prs_scz2025/table_order_of_operations_scz2025.tsv     SCZ 2025 (primary cells)
    prs_*/table_order_of_operations_*.tsv                 any later trait, e.g.
                                                          bipolar / ADHD
A later table is picked up automatically if it has the same columns
(atlas phenotype trait_arm stratum method score beta_1lmm se_1lmm p_adj_1lmm n)
and its trait_arm starts with a name in NEW_TRAITS.  Cells kept: matched arms
only -- EUR stratum with the raw score, pooled ('full') stratum with the
within-ancestry z score ('zanc'), as for every other trait.
"""
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
R = REPO / "genetic_analysis/work/results_70tab_hcp"
OUT = REPO / "genetic_analysis/fig1_inputs/hcp70_prs_key_arms.tsv"
KEY = ["SCZ25_EUR", "SCZ25_META", "MDD_eur", "MDD_pooled", "ASD", "ALZ_noAPOE", "ALZ", "EA"]
NEW_TRAITS = ("BIP", "ADHD")          # trait_arm prefixes from the new GWAS downloads
METHODS = ["CT", "PRSCS", "SBayesR", "SBayesRC"]


def main():
    a = pd.read_csv(R / "prs_final_1lmm/table_order_of_operations_all.tsv", sep="\t")
    s = pd.read_csv(R / "prs_scz2025/table_order_of_operations_scz2025.tsv", sep="\t")
    s = s[s.cell.str.contains("primary") & s.method.isin(METHODS)
          & s.trait_arm.isin(["SCZ25_EUR", "SCZ25_META"])]
    parts = [a, s[a.columns]]
    known = {R / "prs_final_1lmm/table_order_of_operations_all.tsv",
             R / "prs_scz2025/table_order_of_operations_scz2025.tsv"}
    for f in sorted(R.glob("prs_*/table_order_of_operations_*.tsv")):
        if f in known:
            continue
        x = pd.read_csv(f, sep="\t")
        if not set(a.columns) <= set(x.columns):
            print(f"skip {f.relative_to(REPO)}: columns differ")
            continue
        x = x[x.trait_arm.astype(str).str.startswith(NEW_TRAITS) & x.method.isin(METHODS)]
        if "cell" in x.columns:
            x = x[x.cell.astype(str).str.contains("primary")]
        if len(x):
            print(f"add {f.relative_to(REPO)}: {sorted(x.trait_arm.unique())}")
            parts.append(x[a.columns])
    t = pd.concat(parts)
    keep = ((t.stratum == "EUR") & (t.score == "raw")) | ((t.stratum == "full") & (t.score == "zanc"))
    t = t[keep & (t.trait_arm.isin(KEY) | t.trait_arm.astype(str).str.startswith(NEW_TRAITS))]
    t.to_csv(OUT, sep="\t", index=False)
    print(len(t), "rows;", sorted(t.trait_arm.unique()))


if __name__ == "__main__":
    main()
