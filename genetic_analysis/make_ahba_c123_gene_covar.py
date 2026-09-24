#!/usr/bin/env python3
"""Write the AHBA C1-C3 gene-covariate file for MAGMA (step 15).

Why this exists rather than `python -m abcd.magma_export`:  that exporter's
preferred input is the Entrez-keyed loading table shipped by the *sibling*
AHBA repo (`ahba_weights.csv`, found via ABCD_AHBA_MAGMA_DIR or ~/Git/AHBA).
That file is not in this repo, so the exporter cannot be re-run on a checkout
that lacks the sibling -- which is every machine except the one step 7 was
built on.  Step 15 is a gene-property test whose whole point is comparability
across phenotype constructions, so its covariate file must be reproducible
from committed inputs alone.  This script builds it from `data/weights.csv`
(7,973 genes x C1..C3, symbol-keyed, the published Dear et al. 2024 loadings)
and two committed symbol->Entrez sources, and the result is committed.

The symbol->Entrez step is lossy in ways the Ensembl route is not (aliases,
withdrawn symbols, read-through genes sharing a prefix), so the mapping rate
is printed and a loss above 25% aborts.  `data/symbol2entrez.csv` alone maps
only 6,626/7,973: it was built for the snRNA-seq PC1 gene list, not for AHBA,
and is missing many current symbols outright (AAAS, AARS, ABR, ACTG1).  The
gap is closed with `ahba_pls/.../magma_entrez_mygene_raw.json`, which is a
mygene query keyed BY Entrez id -- inverting it is an exact id->symbol record,
not a string match against an alias table -- taking the rate to 7,289/7,973
(91.4%).  The 684 that remain are symbols retired since the AHBA probe set
was annotated.  The resulting gene set is NOT
guaranteed identical to step 7's covariate file; step 15 therefore re-runs the
per-region phenotypes through THIS file as well, so the single-LMM and
per-region columns are compared on one gene set rather than across two.

Columns written (MAGMA drops any covariate column containing NA, so every
column here is complete by construction -- `weights.csv` has no missing
values and unmapped genes are dropped as whole rows):

    C1 C2 C3           signed loadings, for the marginal single-component runs
    C1+ C1- ... C3-     positive/negative split, for the joint run that
                        replicates step 7's `ahba_components_posneg`

The split matters: a component's two poles are different biology (for C1,
sensory versus association cortex), and one signed covariate averages an
enrichment at one pole against the opposite at the other.

Usage:  python3 genetic_analysis/make_ahba_c123_gene_covar.py [-o OUT]
Stdlib only -- it runs on the CSD3 login node without an environment.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "genetic_analysis" / "hpc" / "ahba_c123_gene_covar_entrez.txt"
COMPONENTS = ("C1", "C2", "C3")


def load_symbol_to_entrez(path: Path) -> dict[str, str]:
    """symbol -> Entrez id, first mapping wins (the table is already unique)."""
    lut: dict[str, str] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            sym, ent = row["symbol"].strip(), row["entrez"].strip()
            if sym and ent and sym not in lut:
                lut[sym] = str(int(float(ent)))
    return lut


def add_entrez_keyed_symbols(lut: dict[str, str], path: Path) -> int:
    """Extend `lut` from a mygene query that was keyed BY Entrez id.

    Each record is `{"query": <entrez>, "_id": <entrez>, "symbol": <current>}`,
    so the symbol->id direction is read straight off an authoritative record
    rather than matched against an alias string.  Existing entries win, so
    `data/symbol2entrez.csv` stays the primary table.
    """
    added = 0
    for rec in json.loads(path.read_text()):
        sym, ent = rec.get("symbol"), rec.get("_id") or rec.get("query")
        if sym and ent and sym not in lut:
            lut[sym] = str(ent)
            added += 1
    return added


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default=str(DEFAULT_OUT))
    a = ap.parse_args(argv)

    weights_path = REPO / "data" / "weights.csv"
    lut = load_symbol_to_entrez(REPO / "data" / "symbol2entrez.csv")
    n_base = len(lut)
    n_added = add_entrez_keyed_symbols(
        lut, REPO / "ahba_pls" / "data" / "reference" / "gene_sets"
        / "magma_entrez_mygene_raw.json")
    print(f"symbol->Entrez: {n_base} from symbol2entrez.csv "
          f"+ {n_added} from the Entrez-keyed mygene records = {len(lut)}")

    rows: list[tuple[str, list[float]]] = []
    seen: set[str] = set()
    n_total = n_unmapped = n_dup = 0
    with weights_path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        # weights.csv has an unnamed index column; components follow in order.
        idx = [header.index(c) for c in COMPONENTS]
        for row in reader:
            n_total += 1
            entrez = lut.get(row[0].strip())
            if entrez is None:
                n_unmapped += 1
                continue
            if entrez in seen:
                # Two symbols collapsing onto one Entrez id: keep the first.
                n_dup += 1
                continue
            seen.add(entrez)
            rows.append((entrez, [float(row[i]) for i in idx]))

    kept = len(rows)
    dropped = n_total - kept
    print(f"weights: {n_total} genes x {len(COMPONENTS)} components")
    print(f"  Entrez mapping: kept {kept}/{n_total} "
          f"({n_unmapped} unmapped, {n_dup} duplicate Entrez)")
    if dropped > 0.25 * n_total:
        print(f"FATAL: {dropped}/{n_total} genes lost -- too many to proceed",
              file=sys.stderr)
        return 2

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = (list(COMPONENTS)
            + [f"{c}+" for c in COMPONENTS] + [f"{c}-" for c in COMPONENTS])
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["GENE", *cols])
        for entrez, vals in rows:
            pos = [v if v > 0 else 0.0 for v in vals]
            neg = [-v if v < 0 else 0.0 for v in vals]
            w.writerow([entrez, *(f"{v:.6g}" for v in vals + pos + neg)])
    print(f"wrote {out}  {kept} genes x {len(cols)} covariates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
