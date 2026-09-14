"""Rewrite the GCTA input files' FID/IID to match the genotype .fam.

Site-specific glue, not a pipeline stage.  `abcd.gcta_export` writes

    FID = family_id (e.g. 10210)      IID = NDAR_INV005V6D2C

because its docstring records that "the genotype .fam follows the genetics
convention".  The .fam actually in use on CSD3 does not:

    FID = sub-NDARINV005V6D2C         IID = sub-NDARINV005V6D2C

Both columns therefore disagree, and GCTA identifies an individual by the
FID+IID *pair* -- so left alone this yields an analysis of zero subjects, with no
error.  00_check_inputs.sh would not have caught it either: it compares column 2
only, so an aligned IID and a mismatched FID reads as a clean intersection.

The join is on the 8-character NDAR token, which is invariant across all three
spellings ABCD uses (`NDAR_INVxxxxxxxx`, `sub-NDARINVxxxxxxxx`, and 7.0's bare
`sub-xxxxxxxx`), so this works whichever way either side is written.

    python align_ids.py <pheno-dir> <geno-prefix>

Rewrites phenotypes_gcta.txt, covar_quant.txt and covar_categorical.txt in
place, keeping only rows whose subject is in the .fam, and reports the counts.
Also writes family_map.tsv, because the rewrite overwrites the export's FID and
that column was the only record of which subjects are siblings.
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


def main(pheno_dir: str, geno_prefix: str) -> int:
    pheno_dir = Path(pheno_dir)
    fam = Path(f"{geno_prefix}.fam")

    # token -> (FID, IID) exactly as the .fam spells them
    ref: dict[str, tuple[str, str]] = {}
    for line in fam.read_text().splitlines():
        f = line.split()
        if len(f) < 2:
            continue
        t = token(f[1])
        if t:
            ref[t] = (f[0], f[1])
    print(f"{fam}: {len(ref)} genotyped subjects")

    # Overwriting FID destroys the only copy of the family structure: the export
    # writes FID = family_id, the .fam writes FID = IID, so after the rewrite
    # every subject looks like a singleton family.  GCTA does not care -- it
    # takes relatedness from the GRM -- but tools/prs_assoc.R reads family_id
    # straight off FID for its `(1 | family_id)` term, and with one observation
    # per level lmer errors on every model, so 06_prs died with "no models
    # fitted".  Record the mapping before it is lost; make_prs_pheno.py
    # reattaches it for the PRS step.
    family_map: dict[str, str] = {}

    for name in FILES:
        p = pheno_dir / name
        if not p.exists():
            print(f"  SKIP {name} (absent)")
            continue
        lines = p.read_text().splitlines()
        header, rows = lines[0], lines[1:]
        kept, dropped = [], 0
        for row in rows:
            f = row.split(" ")
            t = token(f[1])
            if t is None or t not in ref:
                dropped += 1
                continue
            family_map.setdefault(t, f[0])   # pre-rewrite FID *is* family_id
            f[0], f[1] = ref[t]
            kept.append(" ".join(f))
        p.write_text("\n".join([header] + kept) + "\n")
        print(f"  {name}: {len(kept)} rows kept, {dropped} not genotyped")

    fmap = pheno_dir / "family_map.tsv"
    with fmap.open("w") as fh:
        fh.write("IID\tfamily_id\n")
        for t, fid in sorted(family_map.items()):
            fh.write(f"{ref[t][1]}\t{fid}\n")
    print(f"  family_map.tsv: {len(family_map)} subjects in "
          f"{len(set(family_map.values()))} families")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1], sys.argv[2]))
