"""Build a PRS-specific copy of the GCTA inputs whose FID is the family id.

Why a second copy rather than a fix in one place:

`tools/prs_assoc.R` fits `phenotype ~ scale(PRS) + ... + (1 | family_id)` and
takes family_id from FID, which is correct for `abcd.gcta_export`'s output
(FID = family_id).  But GCTA matches an individual on the FID+IID *pair*, and
the .fam on CSD3 spells FID = IID = sub-NDARINVxxxxxxxx, so align_ids.py has to
overwrite FID for 01-03 to see any subjects at all.  After that rewrite all
4,126 subjects look like singleton families, lmer errors on every model
("number of levels of each grouping factor must be < number of observations"),
prs_assoc.R skips each one via its try(), and the step exits with the
misleading "no models fitted -- check ID overlap between scores and
phenotypes".  The overlap was in fact complete.

The two consumers want genuinely different FID columns, so they get different
files.  This directory is for prs_assoc.R only; GCTA keeps reading pheno/.
Nothing about the IID column changes, so the merge against PLINK's .profile
output (which is by IID) is unaffected.

Sibling structure is not a rounding detail here: 657 of the 3,455 families in
the analysis sample contribute more than one child, and treating those as
independent would understate the standard error on every PRS coefficient.

    python make_prs_pheno.py <pheno-dir> <out-dir> [family-source]

family-source defaults to <pheno-dir>/family_map.tsv (written by align_ids.py).
It may also be a pre-alignment phenotypes_gcta.txt, whose FID column is the
family id -- regenerate one with `python -m abcd.gcta_export --out-dir ...` if
the pheno dir was aligned by a version of align_ids.py that predates the map.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

TOKEN = re.compile(r"([A-Z0-9]{8})$")

FILES = ("phenotypes_gcta.txt", "covar_quant.txt", "covar_categorical.txt")


def token(s: str) -> str | None:
    m = TOKEN.search(s.strip())
    return m.group(1) if m else None


def load_families(src: Path) -> dict[str, str]:
    """token -> family_id, from either family_map.tsv or a pre-align export."""
    lines = src.read_text().splitlines()
    if not lines:
        raise SystemExit(f"{src} is empty")
    sep = "\t" if "\t" in lines[0] else None
    out: dict[str, str] = {}
    for line in lines[1:]:
        f = line.split(sep)
        if len(f) < 2:
            continue
        # family_map.tsv is IID,family_id; a GCTA export is FID,IID.  In both
        # the family id and the subject id are the two columns, just swapped.
        (subj, fid) = (f[0], f[1]) if src.name == "family_map.tsv" else (f[1], f[0])
        t = token(subj)
        if t:
            out[t] = fid
    return out


def main(pheno_dir: str, out_dir: str, family_source: str | None) -> int:
    pheno_dir, out_dir = Path(pheno_dir), Path(out_dir)
    src = Path(family_source) if family_source else pheno_dir / "family_map.tsv"
    if not src.exists():
        raise SystemExit(
            f"no family source at {src}\n"
            "Re-run align_ids.py to write family_map.tsv, or pass a "
            "pre-alignment phenotypes_gcta.txt as the third argument."
        )
    fam = load_families(src)
    print(f"{src}: {len(fam)} subjects, {len(set(fam.values()))} families")

    out_dir.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        p = pheno_dir / name
        if not p.exists():
            print(f"  SKIP {name} (absent)")
            continue
        lines = p.read_text().splitlines()
        header, rows = lines[0], lines[1:]
        kept, unmapped = [], 0
        for row in rows:
            f = row.split(" ")
            t = token(f[1])
            if t is None or t not in fam:
                unmapped += 1
                continue
            f[0] = fam[t]          # FID <- family_id; IID left alone
            kept.append(" ".join(f))
        (out_dir / name).write_text("\n".join([header] + kept) + "\n")
        n_fam = len({r.split(" ")[0] for r in kept})
        note = f", {unmapped} unmapped (dropped)" if unmapped else ""
        print(f"  {name}: {len(kept)} rows in {n_fam} families{note}")

    manifest = pheno_dir / "phenotype_manifest.tsv"
    if manifest.exists():
        (out_dir / manifest.name).write_text(manifest.read_text())
    print(f"wrote {out_dir}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1], sys.argv[2],
                          sys.argv[3] if len(sys.argv) == 4 else None))
