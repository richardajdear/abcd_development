"""Export gene-set / gene-covariate files for MAGMA.

Run as ``python -m abcd.magma_export --out-dir <dir>``.

MAGMA gene-set tests come in two flavours and the choice is not cosmetic:

* ``--set-annot``  hard sets: "these 1000 genes".  Throws away the magnitude
  of every loading and makes the result depend on an arbitrary cutoff.
* ``--gene-covar`` continuous covariates: every gene carries its loading.
  Uses all the information and has no cutoff to justify.

We write continuous covariates, and additionally split each signed component
into its positive and negative halves (``C1+``, ``C1-``, ...).  The split is
needed because a component's two tails are different biology -- for AHBA C1
the poles are sensory versus association cortex -- and a single signed
covariate would average a positive enrichment at one pole against a negative
one at the other, cancelling real signal.

Gene identifiers are converted to Entrez, because that is what the standard
MAGMA SNP-gene annotation files key on.  Genes that fail to map are dropped
and the count is reported: silent dropout of a third of the gene set would
change the answer without changing the output format.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import genemaps, paths


def _ens2entrez() -> pd.Series:
    """Ensembl -> Entrez, from the mapping table shipped with the repo."""
    p = paths.DATA_DIR / "ens2entrez_mygene.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found; it maps Ensembl ids to Entrez ids for MAGMA"
        )
    m = pd.read_csv(p).dropna()
    return (m.drop_duplicates(subset="ensembl")
              .set_index("ensembl")["entrez"].astype("Int64"))


def ahba_weights_ensembl() -> pd.DataFrame:
    """AHBA component loadings keyed by Ensembl id.

    Prefer this over symbol-keyed loadings: MAGMA's SNP-gene annotations are
    keyed on stable ids, and symbol-to-id conversion is lossy in both
    directions (aliases, withdrawn symbols, and the read-through genes that
    share a symbol prefix).  ``ABCD_AHBA_MAGMA_DIR`` or the sibling AHBA repo
    supplies the file.
    """
    import os
    env = os.environ.get("ABCD_AHBA_MAGMA_DIR")
    cands = [Path(env) / "ahba_weights_ensembl.txt"] if env else []
    cands += [Path.home() / "Git" / "AHBA" / "magma" / "ahba_weights_ensembl.txt",
              paths.DATA_DIR / "ahba_weights_ensembl.txt"]
    for p in cands:
        if p.exists():
            return pd.read_csv(p, sep=r"\s+", index_col=0)
    raise FileNotFoundError(
        "ahba_weights_ensembl.txt not found. Set ABCD_AHBA_MAGMA_DIR, or "
        f"place the file at {paths.DATA_DIR / 'ahba_weights_ensembl.txt'}. "
        f"(looked in: {[str(c) for c in cands]})"
    )


def _to_entrez(frame: pd.DataFrame, id_type: str) -> pd.DataFrame:
    """Attach Entrez ids, dropping unmapped genes and reporting the loss.

    Raises rather than proceeding if more than a quarter of genes fail to
    map: silent dropout of that size changes the gene-set test's answer
    without changing its output format.
    """
    if id_type == "ensembl":
        e = _ens2entrez()
        gid = frame.index.map(e)
    elif id_type == "entrez":
        gid = pd.Index(frame.index, dtype="Int64")
    else:
        raise ValueError(f"unsupported id_type {id_type!r}")
    out = frame.copy()
    out.insert(0, "geneID", gid)
    n_before = len(out)
    out = out.dropna(subset=["geneID"])
    out["geneID"] = out["geneID"].astype(int)
    out = out.drop_duplicates(subset="geneID")
    dropped = n_before - len(out)
    if dropped > 0.25 * n_before:
        raise ValueError(
            f"{dropped}/{n_before} genes failed Entrez mapping -- too many to "
            "proceed; check the mapping table matches the id type"
        )
    print(f"  Entrez mapping: kept {len(out)}/{n_before} genes ({dropped} dropped)")
    return out.reset_index(drop=True)


def _split_posneg(frame: pd.DataFrame) -> pd.DataFrame:
    """Split each signed column into a positive and a negative half."""
    pos = frame.clip(lower=0)
    neg = -frame.clip(upper=0)
    pos.columns = [f"{c}+" for c in frame.columns]
    neg.columns = [f"{c}-" for c in frame.columns]
    return pd.concat([pos, neg], axis=1)


def ahba_weights_entrez() -> pd.DataFrame:
    """AHBA component loadings already keyed by Entrez id.

    The AHBA repo ships these alongside the symbol-keyed version, so no
    identifier conversion is needed for the component covariates -- which is
    the preferred route, since every conversion step loses genes.
    """
    import os
    env = os.environ.get("ABCD_AHBA_MAGMA_DIR")
    cands = [Path(env) / "ahba_weights.csv"] if env else []
    cands += [Path.home() / "Git" / "AHBA" / "magma" / "ahba_weights.csv",
              paths.DATA_DIR / "ahba_weights_entrez.csv"]
    for p in cands:
        if p.exists():
            w = pd.read_csv(p, index_col=0)
            w.index.name = "entrez"
            return w
    raise FileNotFoundError(
        "Entrez-keyed AHBA weights not found. Set ABCD_AHBA_MAGMA_DIR to the "
        f"AHBA repo's magma/ directory. (looked in: {[str(c) for c in cands]})"
    )


def build_ahba() -> pd.DataFrame:
    """Component loadings as MAGMA gene covariates, positive/negative split."""
    try:
        w = ahba_weights_entrez()
        id_type = "entrez"
    except FileNotFoundError:
        print("  Entrez-keyed weights unavailable; converting from Ensembl")
        w = ahba_weights_ensembl()
        id_type = "ensembl"
    return _to_entrez(_split_posneg(w), id_type)


def build_pc1() -> pd.DataFrame:
    """snRNA-seq maturation PC1 loadings as MAGMA gene covariates.

    These loadings are symbol-keyed.  A symbol->Entrez table is required and
    is *not* interchangeable with the Ensembl table: converting symbols by
    string match to an Ensembl-keyed file would silently mis-assign aliases.
    If no symbol table is present this raises with instructions rather than
    guessing, because a mis-keyed gene set produces a plausible-looking but
    meaningless enrichment p-value.
    """
    frames = {}
    for which in ("PC1_herringV3", "PC1_U01V2"):
        try:
            frames[which] = genemaps.snrnaseq_pc1(which)
        except (FileNotFoundError, KeyError) as exc:
            print(f"  skipping {which}: {exc}")
    if not frames:
        raise FileNotFoundError("no snRNA-seq PC1 loadings available")
    loadings = pd.DataFrame(frames)

    sym_path = paths.DATA_DIR / "symbol2entrez.csv"
    if not sym_path.exists():
        raise FileNotFoundError(
            f"snRNA-seq PC1 loadings are symbol-keyed but {sym_path} is absent. "
            "Generate it with:\n"
            "    import mygene; mg = mygene.MyGeneInfo()\n"
            "    r = mg.querymany(symbols, scopes='symbol', fields='entrezgene',\n"
            "                     species='human', as_dataframe=True)\n"
            "and save columns [symbol, entrez]. Do NOT substitute the Ensembl "
            "table: symbol->Ensembl string matching mis-assigns aliases."
        )
    m = pd.read_csv(sym_path).dropna()
    sym_col = next(c for c in ("symbol", "query", "gene") if c in m)
    ent_col = next(c for c in ("entrez", "entrezgene") if c in m)
    lut = m.drop_duplicates(subset=sym_col).set_index(sym_col)[ent_col]
    loadings.index = loadings.index.map(lut)
    loadings = loadings[loadings.index.notna()]
    return _to_entrez(_split_posneg(loadings), "entrez")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=str(paths.OUT_DIR / "magma_genesets"))
    a = ap.parse_args(argv)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, fn in [("ahba_components_posneg", build_ahba),
                     ("snrnaseq_pc1", build_pc1)]:
        print(f"{name}:")
        try:
            frame = fn()
        except (FileNotFoundError, ValueError) as exc:
            print(f"  FAILED: {exc}")
            continue
        p = out / f"{name}.txt"
        frame.to_csv(p, sep="\t", index=False, float_format="%.6g")
        print(f"  wrote {p}  {frame.shape[0]} genes x {frame.shape[1]-1} covariates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
