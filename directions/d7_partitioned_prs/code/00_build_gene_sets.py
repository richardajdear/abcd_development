"""Step 0 (laptop): build the D7 gene sets.  Run from the repo root:

    python directions/d7_partitioned_prs/code/00_build_gene_sets.py

Writes to directions/d7_partitioned_prs/gene_sets/ (Entrez-keyed, committed;
gene lists are public and nothing here is per-subject):

    d7_gene_sets.tsv   set, entrez, symbol        one row per gene x set
    d7_sets_meta.tsv   set, tier, contrast, universe, n_genes, source
    d7_universes.tsv   universe, entrez           the AHBA and snRNA-seq universes

The GENOME universe is not built here. It is "every NCBI37.3 gene with at least
one assigned, weighted SNP", defined on CSD3 by 01_assign_snps.py, because it
depends on the genotype and weight files.

Sets, with the roles DIRECTIONS.md D7 pre-registered (README.md "Pre-registration"):
    primary    K1  brain-expressed genes vs the rest of the genome -- the POSITIVE
                   CONTROL: must enrich.  "Brain-expressed" = the 17.6k genes that
                   passed the expression filter of the developmental cortex
                   snRNA-seq analysis (data/velmeshev_PC1_gene_loadings.csv).
               K2  neuronal markers | glial markers (Seidlitz 2020; made disjoint)
               K3  SynGO 1.3 synaptic genes | oligodendrocyte markers (made disjoint)
               K4  adolescent-window maturation genes from D6 (slot, 0 genes)
    secondary  S0  SCZ locus pool (Trubetskoy 2022 ST12); also a mechanical check
                   that f is large for the SCZ score, not a biological test
               S1  snRNA-seq maturation PC1 top | bottom decile (Herring V3;
                   top = neuronal/synaptic pole, bottom = oligodendrocyte pole)
               S2  AHBA C3 top | bottom decile
               S3  AHBA C1 and C2 top | bottom deciles (controls for S2)
               S4  ABCD PLS2 (HCP-MMP, base AHBA) top | bottom decile
               S5  brain cell-type marker genes (union of Seidlitz classes), the
                   background that K2 is drawn from
S1 uses the snRNA-seq universe; S2-S4 use the AHBA universe (genes with AHBA
C1-C3 weights), which is brain-expressed and longer than average -- the reason
README_HPC 8.3 C3-D requires the null to be drawn from it.  All other sets use
the genome universe.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
OUT = HERE / "gene_sets"
SEIDLITZ = Path.home() / "Git" / "AHBA" / "data" / "seidlitz_cell_genes.csv"
SYNGO = OUT / "source" / "syngo1.3_genes.tsv"
SYM2ENT = REPO / "data" / "symbol2entrez.csv"
ENS2ENT = REPO / "data" / "ens2entrez_mygene.csv"
AHBA = REPO / "data" / "weights.csv"
PC1 = REPO / "data" / "velmeshev_PC1_gene_loadings.csv"
PLS2 = REPO / "ahba_pls" / "results" / "hcp_summary_base_gene_weights.tsv"
LOCUS = REPO / "ahba_pls" / "data" / "reference" / "gene_sets" / "SCZ_locus_pool.txt"

NEURONAL = {"Neuro-Ex", "Neuro-In", "Neuro"}
GLIAL = {"Astro", "Oligo", "OPC", "Micro"}
DECILE = 0.10


def seidlitz() -> pd.DataFrame:
    rows = []
    with open(SEIDLITZ, encoding="utf-8-sig") as f:
        r = csv.reader(f)
        next(r)
        for v in r:
            for g in v[4:]:
                if g.strip():
                    rows.append((v[3], g.strip()))
    return pd.DataFrame(rows, columns=["cls", "symbol"]).drop_duplicates()


def main() -> None:
    s2e = (pd.read_csv(SYM2ENT).dropna().drop_duplicates("symbol")
           .set_index("symbol").entrez.astype(int))
    e2e = (pd.read_csv(ENS2ENT).dropna().drop_duplicates("ensembl")
           .set_index("ensembl").entrez.astype(int))

    def from_symbols(syms) -> pd.Series:          # Series(symbol, index=entrez)
        e = s2e.reindex(pd.Index(sorted(set(syms)))).dropna().astype(int)
        return pd.Series(e.index.values, index=e.values)

    sets: dict[str, pd.Series] = {}
    meta: list[dict] = []

    def add(name, ent_sym, tier, contrast, universe, source):
        ent_sym = ent_sym[~ent_sym.index.duplicated()]
        sets[name] = ent_sym
        meta.append(dict(set=name, tier=tier, contrast=contrast, universe=universe,
                         n_genes=len(ent_sym), source=source))

    sd = seidlitz()
    neu = set(sd[sd.cls.isin(NEURONAL)].symbol)
    gli = set(sd[sd.cls.isin(GLIAL)].symbol)

    pc = pd.read_csv(PC1)
    pc["entrez"] = e2e.reindex(pc.ensembl).values
    pc = pc.dropna(subset=["entrez", "PC1_herringV3"]).drop_duplicates("entrez")
    pc["entrez"] = pc.entrez.astype(int)
    pc = pc.sort_values("PC1_herringV3")
    add("K1_brain_expressed", pd.Series(pc.feature_name.values, index=pc.entrez.values),
        "primary", "K1", "genome",
        "positive control (must enrich): genes passing the developmental-cortex snRNA-seq "
        "expression filter, velmeshev_PC1_gene_loadings.csv")
    add("K2_neuronal", from_symbols(neu - gli), "primary", "K2", "genome",
        f"Seidlitz neuronal classes, minus {len(neu & gli)} symbols also glial markers")
    add("K2_glial", from_symbols(gli - neu), "primary", "K2", "genome",
        "Seidlitz astro/oligo/OPC/micro, minus symbols also neuronal markers")

    sy = pd.read_csv(SYNGO, sep="\t").dropna(subset=["entrez_id"])
    sy_ent = pd.Series(sy.hgnc_symbol.values, index=sy.entrez_id.astype(int).values)
    oli_ent = from_symbols(set(sd[sd.cls == "Oligo"].symbol))
    shared = set(sy_ent.index) & set(oli_ent.index)
    add("K3_syngo", sy_ent[~sy_ent.index.isin(shared)], "primary", "K3", "genome",
        f"SynGO 1.3 (Koopmans 2019) annotated genes, minus {len(shared)} also oligodendrocyte markers")
    add("K3_oligodendrocyte", oli_ent[~oli_ent.index.isin(shared)], "primary", "K3", "genome",
        "Seidlitz Oligo class, minus genes also in SynGO")
    meta.append(dict(set="K4_adolescent_window", tier="primary", contrast="K4", universe="snrna",
                     n_genes=0, source="D6 (transcriptional_maturation): not built; append its rows "
                     "to d7_gene_sets.tsv, then re-run 03 and 05"))

    k = int(round(DECILE * len(pc)))
    add("S1_snrna_pc1_top", pd.Series(pc.feature_name.values[-k:], index=pc.entrez.values[-k:]),
        "secondary", "S1", "snrna", f"PC1_herringV3 top decile of {len(pc)} genes (neuronal/synaptic pole)")
    add("S1_snrna_pc1_bottom", pd.Series(pc.feature_name.values[:k], index=pc.entrez.values[:k]),
        "secondary", "S1", "snrna", f"PC1_herringV3 bottom decile of {len(pc)} genes (oligodendrocyte pole)")

    locus = [l.strip() for l in LOCUS.read_text().split() if l.strip()]
    add("S0_scz_locus_pool", from_symbols(locus), "secondary", "S0", "genome",
        "Trubetskoy 2022 ST12, ahba_pls/data/reference/gene_sets/SCZ_locus_pool.txt")
    add("S5_brain_celltype_markers", from_symbols(set(sd.symbol)), "secondary", "S5", "genome",
        "Seidlitz 2020 marker compilation, union of all classes")

    w = pd.read_csv(AHBA, index_col=0)
    w["entrez"] = s2e.reindex(w.index).values
    w = w.dropna(subset=["entrez"])
    w = w[~w.entrez.duplicated()]
    w["entrez"] = w.entrez.astype(int)
    w["PLS2"] = pd.read_csv(PLS2, sep="\t").set_index("gene")["PLS2"].reindex(w.index).values
    for comp, contrast in (("C3", "S2"), ("C1", "S3"), ("C2", "S3"), ("PLS2", "S4")):
        ww = w.dropna(subset=[comp]).sort_values(comp)
        k = int(round(DECILE * len(ww)))
        for pole, sl in (("top", slice(-k, None)), ("bottom", slice(None, k))):
            add(f"{contrast}_{comp.lower()}_{pole}",
                pd.Series(ww.index.values[sl], index=ww.entrez.values[sl]),
                "secondary", contrast, "ahba", f"{comp} {pole} decile of {len(ww)} AHBA genes")

    OUT.mkdir(parents=True, exist_ok=True)
    pd.concat([pd.DataFrame(dict(set=s, entrez=v.index.astype(int), symbol=v.values))
               for s, v in sets.items()], ignore_index=True).to_csv(
        OUT / "d7_gene_sets.tsv", sep="\t", index=False)
    pd.DataFrame(meta).to_csv(OUT / "d7_sets_meta.tsv", sep="\t", index=False)
    pd.concat([pd.DataFrame(dict(universe="ahba", entrez=w.entrez.values)),
               pd.DataFrame(dict(universe="snrna", entrez=pc.entrez.values))]).to_csv(
        OUT / "d7_universes.tsv", sep="\t", index=False)
    print(pd.DataFrame(meta)[["set", "tier", "universe", "n_genes"]].to_string(index=False))


if __name__ == "__main__":
    main()
