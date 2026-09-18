#!/usr/bin/env python
"""Collect step 10 (MAGMA on the 2025 SCZ GWAS) into tables.

Per atlas root, from results_70tab*/magma_scz2025/*.gsa.out:
  table_magma_scz2025.tsv   one row per test: kind (reverse gene-property =
                            disorder genes ~ phenotype gene Z; forward = phenotype
                            genes ~ disorder gene Z; gene-set marginal /
                            conditional; positive control), disorder result,
                            phenotype, variable, n_genes, beta, se, p
Disorder-side, parcellation-neutral, in work/magma_scz2025/ and copied into
both roots:
  table_gene_level_comparison.tsv   for each pair of disorder gene results:
                            genes in common, Spearman/Pearson of ZSTAT, and per
                            result the Bonferroni-significant gene count -- the
                            direct measurement of what the LD-panel mismatch
                            ("normal MAGMA" on the meta file) does to gene-level
                            inference compared with the per-ancestry --meta.
"""
from __future__ import annotations

import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"dk": REPO / "genetic_analysis/work/results_70tab",
         "hcp": REPO / "genetic_analysis/work/results_70tab_hcp"}
G = REPO / "genetic_analysis/work/magma_scz2025/genes"
DIS = ["SCZ25_EUR", "SCZ25_META", "SCZ25_METAnaive", "PGC3_EUR", "PGC3_primary"]
PHENOS = ["baseline_thickness", "global_slope", "slope_PC1", "slope_PC2", "slope_PC3"]


def read_gsa(p: Path) -> pd.DataFrame:
    rows, cols = [], None
    with open(p) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            if cols is None and v[0] == "VARIABLE":
                cols = v; continue
            if cols:
                rows.append(dict(zip(cols, v)))
    t = pd.DataFrame(rows)
    for c in ("NGENES", "BETA", "BETA_STD", "SE", "P"):
        if c in t:
            t[c] = pd.to_numeric(t[c], errors="coerce")
    return t


def collect(root: Path) -> pd.DataFrame:
    d = root / "magma_scz2025"
    out = []
    for p in sorted(d.glob("*.gsa.out")):
        stem = p.name[:-8]
        m = re.match(rf"({'|'.join(DIS)})_vs_({'|'.join(PHENOS)})$", stem)
        if m:
            kind, dis, ph = "reverse: disorder genes ~ phenotype gene Z", m.group(1), m.group(2)
        else:
            m = re.match(rf"({'|'.join(PHENOS)})_on_({'|'.join(DIS)})$", stem)
            if m:
                kind, dis, ph = "forward: phenotype genes ~ disorder gene Z", m.group(2), m.group(1)
            else:
                m = re.match(rf"({'|'.join(PHENOS)})_sets(_cond(\w+))?$", stem)
                if m:
                    ph = m.group(1); dis = "gene sets"
                    kind = "gene-set marginal" if not m.group(2) else f"gene-set conditional on {m.group(3)}"
                elif stem == "SCZ25_META_sets":
                    kind, dis, ph = "positive control: 2025 sets on SCZ25_META", "gene sets", "SCZ25_META"
                else:
                    continue
        t = read_gsa(p)
        if t.empty:
            continue
        for _, r in t.iterrows():
            out.append(dict(kind=kind, disorder_result=dis, phenotype=ph, variable=r.get("VARIABLE"),
                            type=r.get("TYPE", ""), n_genes=r.get("NGENES"), beta=r.get("BETA"),
                            beta_std=r.get("BETA_STD"), se=r.get("SE"), p=r.get("P"), file=p.name))
    t = pd.DataFrame(out)
    if not t.empty:
        t.to_csv(d / "table_magma_scz2025.tsv", sep="\t", index=False)
    return t


def gene_level() -> pd.DataFrame:
    z = {}
    for n in DIS:
        p = G / f"{n}.genes.out"
        if p.exists():
            t = pd.read_csv(p, sep=r"\s+")
            z[n] = t.set_index("GENE")[["ZSTAT", "P", "N"]]
    rows = []
    for n, t in z.items():
        thr = 0.05 / len(t)
        rows.append(dict(result=n, comparator="", n_genes=len(t), n_bonferroni_sig=int((t.P < thr).sum()),
                         median_N=float(t.N.median()), pearson_z=np.nan, spearman_z=np.nan, n_shared=np.nan,
                         sig_in_both=np.nan, sig_only_in_result=np.nan, sig_only_in_comparator=np.nan))
    for a, b in combinations(z, 2):
        j = z[a].join(z[b], lsuffix="_a", rsuffix="_b", how="inner")
        ta, tb = 0.05 / len(z[a]), 0.05 / len(z[b])
        sa, sb = j.P_a < ta, j.P_b < tb
        rows.append(dict(result=a, comparator=b, n_genes=np.nan, n_bonferroni_sig=np.nan, median_N=np.nan,
                         pearson_z=j.ZSTAT_a.corr(j.ZSTAT_b), spearman_z=j.ZSTAT_a.corr(j.ZSTAT_b, method="spearman"),
                         n_shared=len(j), sig_in_both=int((sa & sb).sum()), sig_only_in_result=int((sa & ~sb).sum()),
                         sig_only_in_comparator=int((~sa & sb).sum())))
    return pd.DataFrame(rows)


def main() -> int:
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
    gl = gene_level()
    for root in ROOTS.values():
        (root / "magma_scz2025").mkdir(parents=True, exist_ok=True)
        gl.to_csv(root / "magma_scz2025" / "table_gene_level_comparison.tsv", sep="\t", index=False)
    print("### gene-level: results and pairwise agreement")
    print(gl.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    for parc, root in ROOTS.items():
        t = collect(root)
        print(f"\n[{parc}] {len(t)} test rows")
        if t.empty:
            continue
        show = t[(t.phenotype.isin(["global_slope", "baseline_thickness"]))
                 & ((t.kind.str.startswith("reverse")) | (t.kind.str.startswith("forward"))
                    | (t.variable.isin(["SCZ25_locus_pool", "SCZ25_genesig", "PGC3_locus_pool", "PGC3_genesig",
                                        "SCZ_locus_pool", "SCZ_prioritised", "SCZ25_locus_pool_new", "SCZ25_genesig_new"])))]
        print(show[["kind", "disorder_result", "phenotype", "variable", "n_genes", "beta", "se", "p"]]
              .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
