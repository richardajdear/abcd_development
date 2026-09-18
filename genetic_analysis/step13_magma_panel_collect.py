#!/usr/bin/env python
"""Collect step 13: MAGMA panel (11 disorder results x 4 phenotype files x 2 atlases).
Writes results_70tab*/magma_panel/table_magma_panel.tsv with columns
kind, disorder_result, phenotype, construction, variable, n_genes, beta, se, p."""
from __future__ import annotations
import os, re, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from step10_scz2025_magma_collect import read_gsa  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
ROOTS = {"dk": REPO / "genetic_analysis/work/results_70tab", "hcp": REPO / "genetic_analysis/work/results_70tab_hcp"}
DIS = ["SCZ25_EUR", "SCZ25_META", "SCZ25_METAnaive", "PGC3_EUR", "PGC3_primary", "MDD_EUR", "MDD_div", "ASD", "ALZ", "ALZ_noAPOE", "EA"]
PH = ["global_slope", "baseline_thickness", "global_slope_1lmm", "baseline_thickness_1lmm"]
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400)
for parc, root in ROOTS.items():
    d = root / os.environ.get("MAGMA_PANEL_DIR", "magma_panel"); rows = []
    for p in sorted(d.glob("*.gsa.out")):
        stem = p.name[:-8]
        m = re.match(rf"({'|'.join(DIS)})_vs_({'|'.join(PH)})$", stem)
        if m: kind, dis, ph = "reverse", m.group(1), m.group(2)
        else:
            m = re.match(rf"({'|'.join(PH)})_on_({'|'.join(DIS)})$", stem)
            if m: kind, dis, ph = "forward", m.group(2), m.group(1)
            else:
                m = re.match(rf"({'|'.join(PH)})_sets_(scz|mdd)$", stem)
                if not m: continue
                kind, dis, ph = "gene-set", f"{m.group(2).upper()} sets", m.group(1)
        t = read_gsa(p)
        for _, r in t.iterrows():
            rows.append(dict(kind=kind, disorder_result=dis, phenotype=ph.replace("_1lmm", ""),
                             construction="1lmm" if ph.endswith("_1lmm") else "perregion", variable=r.get("VARIABLE"),
                             n_genes=r.get("NGENES"), beta=r.get("BETA"), se=r.get("SE"), p=r.get("P")))
    t = pd.DataFrame(rows); t.to_csv(d / f"table_{d.name}.tsv", sep="\t", index=False)
    print(f"\n[{parc}] {len(t)} rows; reverse gene-property p-values (disorder genes ~ phenotype gene Z):")
    r = t[t.kind == "reverse"].pivot_table(index="disorder_result", columns=["phenotype", "construction"], values="p")
    print(r.to_string(float_format=lambda x: f"{x:.3g}"))
    print(f"[{parc}] gene sets, p:")
    s = t[(t.kind == "gene-set") & (t.variable.isin(["SCZ_locus_pool", "SCZ25_locus_pool", "SCZ25_genesig", "MDD_pool", "MDD_highconf"]))]
    print(s.pivot_table(index="variable", columns=["phenotype", "construction"], values="p").to_string(float_format=lambda x: f"{x:.3g}"))
