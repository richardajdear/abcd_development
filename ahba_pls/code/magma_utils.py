"""
magma_utils.py -- MAGMA gene-property helpers shared by the summary-slide
scripts: symbol -> Entrez mapping, --gene-covar file writing, and a runner
that returns the .gsa.out table. Mirrors 21_magma_disorder_panel.py (same
gene-location file, same duplicate-symbol rule, same default two-sided model).
"""
from __future__ import annotations
import subprocess, sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
INP = REPO / "genetic_analysis" / "inputs" / "magma"


def magma_bin() -> Path:
    for c in ("magma_src/magma", "magma_mac/magma", "magma"):
        if (REPO / "tools" / "bin" / c).exists():
            return REPO / "tools" / "bin" / c
    raise FileNotFoundError("no MAGMA binary under tools/bin; see tools/bin/README.md")


def sym2entrez() -> pd.Series:
    loc = pd.read_csv(REPO / "tools" / "bin" / "NCBI37.3.gene.loc", sep="\t", header=None,
                      names=["GENE", "chr", "start", "end", "strand", "symbol"])
    loc = loc[~loc.symbol.duplicated(keep=False)]          # ambiguous symbols dropped, not guessed
    return loc.set_index("symbol").GENE.astype(int)


def write_covar(vectors: dict[str, pd.Series], path: Path, s2e: pd.Series) -> Path:
    c = pd.DataFrame(vectors).dropna(how="any")
    c = c.loc[c.index.intersection(s2e.index)]
    c.index = s2e.loc[c.index].values
    c = c[~c.index.duplicated()]; c.index.name = "GENE"
    c.to_csv(path, sep="\t", float_format="%.6g")
    return path


def run(raw: Path, covar: Path, out: Path, extra: list[str] | None = None) -> pd.DataFrame:
    r = subprocess.run([str(magma_bin()), "--gene-results", str(raw), "--gene-covar", str(covar),
                        "--out", str(out), *(extra or [])], capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"MAGMA failed for {out.name}:\n{r.stdout[-1500:]}")
    t = pd.read_csv(f"{out}.gsa.out", sep=r"\s+", comment="#")
    if "FULL_NAME" in t.columns:
        t["VARIABLE"] = t["FULL_NAME"].fillna(t["VARIABLE"])
    t["se_std"] = t.SE * t.BETA_STD / t.BETA
    return t
