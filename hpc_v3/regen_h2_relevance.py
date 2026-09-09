"""Rebuild `h2_vs_scz_relevance.csv` with the four Experiment A phenotypes.

`README.md` section 2 item 5 asks for four new rows on the h2-vs-relevance
table --- the anti-ranking of `SETUP_CONTEXT.md` section 4(c), whose whole
point is that this project's most heritable phenotypes are its least
disorder-relevant ones.  Experiment A is its first out-of-sample test: the
subset phenotypes were chosen for biological relevance, so if the pattern
holds they should sit low-h2 / high-association.

Every column is re-derived from a result file rather than copied, and the five
settled rows are rebuilt the same way as a check that the provenance below is
the real one (they were reverse-engineered from the committed CSV, not
documented anywhere):

| column | source |
|---|---|
| `h2_obs`, `h2_obs_se`, `h2_z`, `rg` | LDSC EUR arm, SCZ rows of `ldsc_rg_summary.tsv` |
| `p_between` | C+T SCZ p<0.5, EUR stratum, `beta_between`'s p in `prs_withinfamily*.tsv` |
| `p_v2` | MAGMA reverse direction, COVAR p in `SCZ_vs_<pheno>.gsa.out` (EUR arm) |
| `prs_evidence` | -log10(`p_between`) |
| `abs_rg` | abs(`rg`) |

    python hpc_v3/regen_h2_relevance.py [--check-only]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
V2 = REPO / "hpc_v2/work/results_v2"
SETTLED = ["baseline_thickness", "global_slope", "slope_PC3", "slope_PC2",
           "slope_PC1"]
NEW = ["slope_topDelta", "slope_topC3", "slope_projDelta", "slope_projC3"]


def ldsc_rows(*paths: Path) -> pd.DataFrame:
    """SCZ rows from one or more ldsc_rg_summary.tsv, later files winning."""
    frames = [pd.read_csv(p, sep="\t") for p in paths if p.exists()]
    if not frames:
        raise FileNotFoundError("no ldsc_rg_summary.tsv found")
    d = pd.concat(frames, ignore_index=True)
    d = d[d.disorder == "SCZ"].drop_duplicates("phenotype", keep="last")
    return d.set_index("phenotype")[["h2_obs", "h2_obs_se", "h2_z", "rg"]]


def prs_between(*paths: Path) -> pd.Series:
    frames = [pd.read_csv(p, sep="\t") for p in paths if p.exists()]
    d = pd.concat(frames, ignore_index=True)
    d = d[(d.disorder == "SCZ") & (d.threshold == "0p5") & (d.stratum == "EUR")]
    return d.drop_duplicates("phenotype", keep="last").set_index("phenotype").p_between


def magma_covar_p(pheno: str, *dirs: Path) -> float:
    """COVAR row's P from SCZ_vs_<pheno>.gsa.out, first directory that has it."""
    for d in dirs:
        f = d / f"SCZ_vs_{pheno}.gsa.out"
        if not f.exists():
            continue
        for line in f.read_text().splitlines():
            if line.startswith("#") or line.startswith("VARIABLE"):
                continue
            parts = re.split(r"\s+", line.strip())
            if len(parts) >= 6 and parts[1] == "COVAR":
                return float(parts[-1])
    return np.nan


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check-only", action="store_true",
                    help="rebuild and compare, but do not write")
    a = ap.parse_args()

    h = ldsc_rows(V2 / "ldsc_eur/ldsc_rg_summary.tsv",
                  V2 / "ldsc_eur_v3/ldsc_rg_summary.tsv")
    pb = prs_between(V2 / "prs_family/prs_withinfamily.tsv",
                     REPO / "hpc_v3/work/prs/prs_withinfamily_v3.tsv")

    rows = []
    for p in SETTLED + NEW:
        r = dict(phenotype=p)
        if p in h.index:
            r.update(h.loc[p].to_dict())
        r["p_between"] = float(pb.get(p, np.nan))
        r["p_v2"] = magma_covar_p(p, V2 / "magma_eur_v3", V2 / "magma_eur")
        r["prs_evidence"] = -np.log10(r["p_between"])
        r["abs_rg"] = abs(r.get("rg", np.nan))
        rows.append(r)
    out = pd.DataFrame(rows)[["phenotype", "h2_obs", "h2_obs_se", "h2_z", "rg",
                              "p_between", "p_v2", "prs_evidence", "abs_rg"]]

    old = pd.read_csv(REPO / "hpc_v3/h2_vs_scz_relevance.csv").set_index("phenotype")
    chk = out.set_index("phenotype").loc[SETTLED]
    print("reproduction of the five committed rows (rebuilt - committed):")
    for c in old.columns:
        d = (chk[c].astype(float) - old[c].astype(float)).abs().max()
        print(f"  {c:<14} max |diff| = {d:.3g}")
    print()
    print(out.to_string(index=False))

    if not a.check_only:
        out.to_csv(REPO / "hpc_v3/h2_vs_scz_relevance.csv", index=False)
        print(f"\nwrote hpc_v3/h2_vs_scz_relevance.csv ({len(out)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
