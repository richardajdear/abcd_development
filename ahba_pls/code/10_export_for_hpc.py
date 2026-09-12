"""
10_export_for_hpc.py -- export the lead signature in the formats the cluster
follow-up (FOLLOWUP_GENETICS.md; hypotheses H3 and H4) consumes.

  hpc/lead_pls2_dk_scores_68.csv        label (lh_/rh_), thinning_score  -- the
                                        33-region gene-side score map mirrored
                                        to both hemispheres; frontalpole = NA.
                                        Input to the H4 projection phenotype.
  hpc/lead_pls2_gene_covar_entrez.txt   MAGMA --gene-covar file (GENE = Entrez,
                                        NCBI37.3 ids as in hpc/): thinning_Z at
                                        ds0/ds25/ds50 + AHBA C3 for conditioning.
                                        Input to the H3 gene-property test.
  hpc/lead_pls2_gene_weights_symbol.tsv gene symbol, thinning_Z_ds25, U, rank.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES, REF, HPC = ROOT / "results", ROOT / "data" / "reference", ROOT / "hpc"
HPC.mkdir(exist_ok=True)

S = pd.read_csv(RES / "lead_signature_scores.csv", index_col=0)["thinning_scores_gene_side"]
labels = pd.read_csv(ROOT.parent / "data" / "region_labels.csv")
labels = labels[(labels.parcellation == "dsk") & (~labels.is_global)]
lab_col = "label"
assert len(labels) == 68, len(labels)
out = []
for l in labels[lab_col]:
    stem = l.split("_", 1)[1]
    out.append((l, S.get(f"lh_{stem}", float("nan"))))
dk = pd.DataFrame(out, columns=["label", "thinning_score"]).set_index("label")
assert dk.thinning_score.notna().sum() == 66, dk.thinning_score.notna().sum()
dk.to_csv(HPC / "lead_pls2_dk_scores_68.csv", float_format="%.6g")

W = pd.read_csv(RES / "lead_signature_weights.tsv", sep="\t", index_col=0)
W[["thinning_Z_ds25", "U_ds25_thinning", "rank_ds25", "decile_ds25"]].to_csv(
    HPC / "lead_pls2_gene_weights_symbol.tsv", sep="\t", float_format="%.6g")

sym2ent = (pd.read_csv(REF / "gene_sets" / "magma_SCZ_genes.tsv", sep="\t").dropna(subset=["symbol"])
           .drop_duplicates("symbol").set_index("symbol")["GENE"].astype(int))
c3 = pd.read_csv(REF / "ahba_c123_gene_weights.csv", index_col=0)["C3"]
cov = W[["thinning_Z_ds0", "thinning_Z_ds25", "thinning_Z_ds50"]].join(c3.rename("AHBA_C3"), how="left")
cov = cov.loc[cov.index.intersection(sym2ent.index)]
cov.index = sym2ent.loc[cov.index].values
cov = cov[~cov.index.duplicated()]
cov.index.name = "GENE"
cov.to_csv(HPC / "lead_pls2_gene_covar_entrez.txt", sep="\t", na_rep="NA", float_format="%.6g")
print(dk.shape, W.shape, cov.shape, cov.notna().sum().to_dict())
