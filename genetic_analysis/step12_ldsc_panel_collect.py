#!/usr/bin/env python
"""Parse the step-12 LDSC rg logs into one table.

Columns: disorder, atlas, phenotype, construction (perregion | 1lmm), rg, se, z,
p, h2_obs, h2_obs_se, h2_z, gcov_int, intercept_pheno, underpowered (h2 z < 4,
rule 13: rg not interpretable).  Written to work/ldsc_panel/table_ldsc_panel.tsv
and copied to results_70tab*/ldsc_1lmm/.
"""
from __future__ import annotations
import re, sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "genetic_analysis/work/ldsc_panel"
rows = []
for log in sorted(OUT.glob("rg_*.log")):
    D = log.name[3:-4]
    txt = log.read_text()
    m = re.search(r"Summary of Genetic Correlation Results\n(.*?)\n\n", txt, re.S)
    if not m:
        print(f"  {D}: no summary block"); continue
    lines = [l.split() for l in m.group(1).strip().splitlines()]
    hdr = lines[0]
    # the phenotype's own LDSC intercept comes from its "Heritability of phenotype 2/…" block
    ints = {}
    for blk in re.finditer(r"Heritability of phenotype (\d+)(?:/\d+)?\n-+\n(.*?)\n\n", txt, re.S):
        mi = re.search(r"Intercept: ([-\d.]+) \(([\d.]+)\)", blk.group(2))
        if mi: ints[int(blk.group(1))] = float(mi.group(1))
    for i, v in enumerate(lines[1:], start=2):
        r = dict(zip(hdr, v))
        ph = Path(r["p2"]).name.replace(".sumstats.gz", "")
        parc, pheno, cons = ph.split("__")
        def f(x):
            try: return float(x)
            except: return float("nan")
        h2, h2se = f(r["h2_obs"]), f(r["h2_obs_se"])
        z = h2 / h2se if h2se else float("nan")
        rows.append(dict(disorder=D, atlas=parc, phenotype=pheno, construction=cons,
                         rg=f(r["rg"]), se=f(r["se"]), z=f(r["z"]), p=f(r["p"]),
                         h2_obs=h2, h2_obs_se=h2se, h2_z=round(z, 2), gcov_int=f(r["gcov_int"]),
                         intercept_pheno=ints.get(i, float("nan")),
                         underpowered="yes" if z < 4 else "no"))
# a phenotype whose rg errored inside LDSC (jackknife h2 products negative --
# the slope's h2 is indistinguishable from zero) gets its h2 from a standalone
# --h2 run saved as h2__<atlas>__<pheno>__<construction>.log, and rg = NaN
for log in sorted(OUT.glob("h2__*.log")):
    _, parc, pheno, cons = log.stem.split("__")
    m = re.search(r"Total Observed scale h2: ([-\d.]+) \(([\d.]+)\)", log.read_text())
    mi = re.search(r"Intercept: ([-\d.]+)", log.read_text())
    if not m: continue
    h2, h2se = float(m.group(1)), float(m.group(2))
    for r in rows:   # rows LDSC emitted with NA h2 for this phenotype: fill h2, keep rg NaN
        if r["atlas"] == parc and r["phenotype"] == pheno and r["construction"] == cons and r["h2_obs"] != r["h2_obs"]:
            r.update(h2_obs=h2, h2_obs_se=h2se, h2_z=round(h2 / h2se, 2), underpowered="yes",
                     intercept_pheno=float(mi.group(1)) if mi else float("nan"))
    for D in sorted({r["disorder"] for r in rows}):
        if not any(r["atlas"] == parc and r["phenotype"] == pheno and r["construction"] == cons and r["disorder"] == D for r in rows):
            rows.append(dict(disorder=D, atlas=parc, phenotype=pheno, construction=cons, rg=float("nan"), se=float("nan"),
                             z=float("nan"), p=float("nan"), h2_obs=h2, h2_obs_se=h2se, h2_z=round(h2 / h2se, 2), gcov_int=float("nan"),
                             intercept_pheno=float(mi.group(1)) if mi else float("nan"), underpowered="yes"))
t = pd.DataFrame(rows)
if t.empty: sys.exit("no rows")
t.to_csv(OUT / "table_ldsc_panel.tsv", sep="\t", index=False)
for s in ("", "_hcp"):
    d = REPO / f"genetic_analysis/work/results_70tab{s}/ldsc_1lmm"; d.mkdir(parents=True, exist_ok=True)
    t.to_csv(d / "table_ldsc_panel.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(t.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
