"""Manifest of repo deliverables, for reconciling against the artifact tray.

    make audit-artifacts

Background: files edited on disk are not visible in the Claude Science artifact
tray until re-saved there, and committing does not update the tray.  This has
bitten the project three times -- most severely when the entire ``hpc/``
pipeline was correct on disk and in git but had never been saved as an artifact
at all, so the tray showed nothing.

**Why this script only does half the job.**  Comparing against the tray needs
``host.artifacts()``, and ``host`` is a *kernel global* injected into the Claude
Science python kernel -- not an importable module.  No subprocess can reach it,
so neither ``make audit-artifacts`` nor any script the user runs in a terminal
can query the tray.  This script therefore does the portable half: it enumerates
the deliverables and checksums them.  The comparison is done agent-side by
reading the manifest this writes.  See the ``audit-artifacts`` target in the
Makefile for the exact agent-side snippet.

The direction is deliberately one-way.  Repo deliverables must be current in the
tray; the converse does not hold -- the tray legitimately holds exploratory
figures, checkpoints, and one-off analyses that were never meant to land in the
repo, and flagging those as "missing from the repo" would be noise.

Artifact names are BASENAMES, because that is what ``save_artifacts`` keys on.
Where two deliverables share a basename (``README.md`` and ``hpc/README.md``)
the name is path-qualified, since otherwise the two files silently share one
artifact and an audit compares the wrong pair -- which is exactly how the
missing ``hpc/`` pipeline first presented as a merely "stale README".
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# What counts as a deliverable: something a collaborator or the user would read
# or receive.  Generated data under out/, the vendored release, and scratch/ are
# excluded by construction (not listed here) as well as by the gitignore check.
PATTERNS = [
    "README.md",
    "Makefile",
    "docs/*.md",
    "docs/*.csv",
    "docs/figures/*.png",
    "genetic_analysis/*.md",          # hpc/README.md -- the file whose absence from the tray
                         # prompted this script; omitting .md here would have
                         # made the audit blind to exactly that case.
    "genetic_analysis/*.sbatch",
    "genetic_analysis/*.sh",
    "genetic_analysis/*.example",
    "configs/*.yaml",
    "src/abcd/*.py",
    "tools/*.py",
    "tools/*.sh",
    "tools/*.R",
    "R/*.R",
    "tests/*.py",
]


def _gitignored(paths: list[Path]) -> set[Path]:
    """Which of ``paths`` git ignores.

    Batched through one ``check-ignore --stdin`` call rather than one per file.

    Two things this must get right, both learned the hard way:

    * ``check-ignore`` exits **1** when nothing matched (fine) but **128** on a
      real failure.  Under a sandboxed kernel it fails with 128 because
      ``~/.gitconfig`` is unreadable, and the previous version treated any
      non-zero exit as "nothing matched" -- so every gitignored file, including
      ``hpc/config.local.sh``, was silently offered for publishing.  A failure
      here must raise, because the safe-looking answer is the dangerous one.
    * Pointing the config vars at ``/dev/null`` stops git reading the
      unreadable user config in the first place, which is what makes the call
      work from a kernel at all.
    """
    if not paths:
        return set()
    rel = [str(p.relative_to(ROOT)) for p in paths]
    env = {**os.environ, "GIT_CONFIG_GLOBAL": "/dev/null",
           "GIT_CONFIG_SYSTEM": "/dev/null"}
    proc = subprocess.run(
        ["git", "-C", str(ROOT), "check-ignore", "--stdin"],
        input="\n".join(rel), capture_output=True, text=True, env=env,
    )
    if proc.returncode not in (0, 1):
        raise RuntimeError(
            f"git check-ignore failed (rc={proc.returncode}): "
            f"{proc.stderr.strip()[:200]}\n"
            "Refusing to continue: without it, gitignored files (local configs, "
            "credentials) would be reported as publishable deliverables."
        )
    return {ROOT / line for line in proc.stdout.splitlines() if line}


def collect() -> list[dict]:
    """Deliverables with their artifact name, size and checksum."""
    found: list[Path] = []
    for pat in PATTERNS:
        found += ROOT.glob(pat)
    found = sorted({p for p in found if p.is_file()})

    ignored = _gitignored(found)
    found = [p for p in found if p not in ignored]

    # Path-qualify only the basenames that actually collide, so the common case
    # keeps its natural name and diffs against existing artifacts stay clean.
    counts: dict[str, int] = {}
    for p in found:
        counts[p.name] = counts.get(p.name, 0) + 1

    rows = []
    for p in found:
        rel = p.relative_to(ROOT)
        name = p.name if counts[p.name] == 1 else str(rel).replace("/", "_")
        rows.append({
            "path": str(rel),
            "artifact_name": name,
            "size_bytes": p.stat().st_size,
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(ROOT / ".artifact_audit.json"),
                    help="where to write the manifest (default: repo root, gitignored)")
    a = ap.parse_args()

    rows = collect()
    Path(a.out).write_text(json.dumps(rows, indent=1) + "\n")

    collisions = [r for r in rows if "_" in r["artifact_name"]
                  and r["artifact_name"] != Path(r["path"]).name]
    print(f"{len(rows)} deliverables -> {a.out}")
    if collisions:
        print("\nbasename collisions, stored path-qualified:")
        for r in collisions:
            print(f"  {r['path']:44s} -> {r['artifact_name']}")
    print("\nThis manifest is only half the audit: comparing it against the tray "
          "needs host.artifacts(),\nwhich exists only inside a Claude Science "
          "kernel. Run `make help` for the agent-side step.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
