"""Replace the release ancestry PCs in covar_quant.txt with in-sample ones.

    python insample_pcs.py <pheno-dir> <eigenvec> [--n-pcs N] [--allow-missing]

WHY

`abcd.gcta_export` takes the PC block from the release's static genetics table
(`ab_g_stc__gen_pc__01..10`), which was computed on the FULL multi-ancestry ABCD
cohort.  Inside the EUR-only analysis sample those PCs are near-constant: PC1
has mean 0.00631 and SD 0.000523, so its mean is 12x its spread.  A PC that
separates ancestry groups looks exactly like that within one group, and it
cannot absorb within-EUR structure -- north-west vs south European gradients --
which is the structure that actually needs controlling in a GRM-based h2.

PCs computed *within* the analysis sample can.  They are the correct covariate
and they are free, so the pipeline uses them.

Do not expect the numbers to move.  Measured on baseline_thickness (README_HPC
section 15), the swap shifts h2 from 0.5754 to 0.5663 -- 0.06 of one standard
error -- because in-sample PC1 explains 0.13% of variance and there is no
residual structure inside EUR for a PC to correct.  This is a correctness fix,
not a results fix, and stating that plainly is the point: an unjustifiable
covariate that happens to be harmless is still worth removing, and the next
person should not have to re-derive that it was harmless.

BEHAVIOUR

The original file is preserved as covar_quant.txt.releasepcs on first run and
re-read as the source on every later run, so running this twice is the same as
running it once -- no compounding, and the release PCs stay recoverable.
Provenance (source path, md5, subject count) is written to covar_quant.pcsource
and printed by 00_check_inputs.sh, so which PCs a result used is a fact on disk
rather than an inference from when the file was written.

Order in the setup chain: gcta_export -> align_ids.py -> THIS -> make_prs_pheno.py.
After align_ids.py, because that one rewrites covar_quant.txt in place and the
.releasepcs backup taken here would otherwise hold the pre-alignment IDs.
Before make_prs_pheno.py, so that the PRS copy inherits the same covariates;
otherwise 02/03 and 06 silently adjust for different PCs.
If the PRS copy already exists (make_prs_pheno.py needs a family_map.tsv that
may be long gone), run this on that directory too instead of rebuilding it: the
join is on IID and only PC columns are rewritten, so its FID -- the family id,
which is the whole reason that copy exists -- is left alone.

A subject in covar_quant.txt with no row in the eigenvec is a hard error, not a
silent drop: dropping is how an analysis sample shrinks without anyone noticing.
Pass --allow-missing to drop them deliberately, and the count is reported.
"""
from __future__ import annotations

import hashlib
import re
import sys
from datetime import datetime
from pathlib import Path

TOKEN = re.compile(r"([A-Z0-9]{8})$")
BACKUP_SUFFIX = ".releasepcs"


def token(s: str) -> str | None:
    m = TOKEN.search(s.strip())
    return m.group(1) if m else None


def read_eigenvec(path: Path) -> tuple[dict[str, list[str]], int]:
    """IID -> [PC1..PCk], from GCTA's headerless FID IID PC1..PCk."""
    pcs: dict[str, list[str]] = {}
    n_pc = 0
    for line in path.read_text().splitlines():
        f = line.split()
        if len(f) < 3:
            continue
        n_pc = max(n_pc, len(f) - 2)
        pcs[f[1]] = f[2:]
    if not pcs:
        raise SystemExit(f"{path} holds no usable rows")
    return pcs, n_pc


def sd(raw: list[str]) -> float:
    """SD over the numeric entries; "NA" is a value this column really carries."""
    values = [float(v) for v in raw if v not in ("NA", "nan", "")]
    n = len(values)
    if n < 2:
        return float("nan")
    mean = sum(values) / n
    return (sum((v - mean) ** 2 for v in values) / (n - 1)) ** 0.5


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    flags = [a for a in argv if a.startswith("--")]
    if len(args) != 2 or any(
        f != "--allow-missing" and not f.startswith("--n-pcs=") for f in flags
    ):
        print(__doc__)
        return 2

    allow_missing = "--allow-missing" in flags
    # `--n-pcs=N` only, not `--n-pcs N`: the space form leaves the value looking
    # like a positional argument, and the arity check above would reject it with
    # a usage dump rather than the reason.
    want_pcs = next(
        (int(f.split("=", 1)[1]) for f in flags if f.startswith("--n-pcs=")), None
    )

    pheno_dir, eig_path = Path(args[0]), Path(args[1])
    target = pheno_dir / "covar_quant.txt"
    backup = pheno_dir / ("covar_quant.txt" + BACKUP_SUFFIX)
    if not eig_path.exists():
        raise SystemExit(f"no eigenvec at {eig_path}")
    if not target.exists() and not backup.exists():
        raise SystemExit(f"no covar_quant.txt in {pheno_dir}")

    # Re-read the pristine release file when one exists, so a second run cannot
    # merge in-sample PCs on top of in-sample PCs and call the result new.
    if backup.exists():
        source = backup
        print(f"source: {backup} (release PCs, preserved by an earlier run)")
    else:
        backup.write_text(target.read_text())
        source = backup
        print(f"source: {target} -> preserved as {backup.name}")

    lines = source.read_text().splitlines()
    header = lines[0].split(" ")
    rows = [ln.split(" ") for ln in lines[1:] if ln.strip()]

    pc_idx = [i for i, c in enumerate(header) if re.fullmatch(r"PC[0-9]+", c)]
    if not pc_idx:
        raise SystemExit(f"{source} has no PC columns to replace")
    # Match the existing block exactly unless overridden: downstream consumers
    # name these columns (GCTA --qcovar takes all of them, prs_assoc.R greps
    # ^PC[0-9]+$), so changing how many there are changes every model.
    n_pcs = want_pcs or len(pc_idx)
    if n_pcs != len(pc_idx):
        raise SystemExit(
            f"--n-pcs {n_pcs} but {source.name} carries {len(pc_idx)} PC columns; "
            "adding or removing PCs changes every downstream model, so do it in "
            "the export rather than here"
        )

    pcs, n_avail = read_eigenvec(eig_path)
    if n_avail < n_pcs:
        raise SystemExit(f"{eig_path} has {n_avail} PCs, need {n_pcs}")
    print(f"eigenvec: {len(pcs)} subjects, {n_avail} PCs, using the first {n_pcs}")

    iid_col = header.index("IID")
    hit = sum(1 for r in rows if r[iid_col] in pcs)
    by_token = False
    if hit < len(rows) // 2:
        # ABCD spells the same subject sub-NDARINVxxxxxxxx here and NDAR_INVxxxxxxxx
        # in the genetics tables; align_ids.py joins on the 8-character token for
        # exactly this reason.  Fall back rather than report an empty overlap.
        tokened = {token(k): v for k, v in pcs.items() if token(k)}
        hit_tok = sum(1 for r in rows if token(r[iid_col]) in tokened)
        if hit_tok > hit:
            pcs, by_token, hit = tokened, True, hit_tok
            print("  IID spellings differ; joined on the 8-character NDAR token")

    missing = [r[iid_col] for r in rows if (token(r[iid_col]) if by_token else r[iid_col]) not in pcs]
    if missing and not allow_missing:
        raise SystemExit(
            f"{len(missing)} of {len(rows)} subjects have no row in {eig_path.name} "
            f"(e.g. {', '.join(missing[:3])}).\n"
            "Those subjects would be dropped from every downstream model.  Use an "
            "eigenvec covering the analysis sample, or pass --allow-missing to drop "
            "them on purpose."
        )

    release_sd = [sd([r[i] for r in rows]) for i in pc_idx[:4]]
    # The release PC block is missing outright for a handful of subjects, and
    # GCTA and lmer both drop a row with any missing covariate.  Those subjects
    # come back when the PCs come from the in-sample eigenvec, so count them:
    # the analysis N changing is a consequence of this swap worth reporting.
    was_na = sum(1 for r in rows if any(r[i] == "NA" for i in pc_idx))

    kept = []
    for r in rows:
        key = token(r[iid_col]) if by_token else r[iid_col]
        if key not in pcs:
            continue
        for j, i in enumerate(pc_idx[:n_pcs]):
            r[i] = pcs[key][j]
        kept.append(r)

    insample_sd = [sd([r[i] for r in kept]) for i in pc_idx[:4]]
    still_na = sum(1 for r in kept if any(r[i] == "NA" for i in pc_idx))

    target.write_text("\n".join([" ".join(header)] + [" ".join(r) for r in kept]) + "\n")

    md5 = hashlib.md5(eig_path.read_bytes()).hexdigest()
    (pheno_dir / "covar_quant.pcsource").write_text(
        "\n".join(
            [
                f"source     {eig_path.resolve()}",
                f"md5        {md5}",
                f"n_pcs      {n_pcs}",
                f"subjects   {len(kept)}",
                f"joined_on  {'NDAR token' if by_token else 'IID'}",
                f"written    {datetime.now().isoformat(timespec='seconds')}",
                f"backup     {backup.name}",
            ]
        )
        + "\n"
    )

    dropped = len(rows) - len(kept)
    print(f"  {len(kept)} subjects written{f', {dropped} dropped' if dropped else ''}")
    if was_na:
        print(f"  {was_na} subjects had no release PCs at all; "
              f"{was_na - still_na} of them now have in-sample PCs and rejoin "
              "the analysis sample")
    print("  PC SDs   release : " + ", ".join(f"{v:.6f}" for v in release_sd))
    print("  PC SDs   in-sample: " + ", ".join(f"{v:.6f}" for v in insample_sd))
    print(f"wrote {target}")
    print(f"wrote {pheno_dir / 'covar_quant.pcsource'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
