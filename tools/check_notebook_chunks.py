#!/usr/bin/env python3
"""Extract and execute every code chunk of a .qmd, in order.

Quarto is not always installed, and a notebook whose code fails will not render
no matter what else is right about it.  This runs the chunks as one script and
exits non-zero on the first failure, so it is usable in CI.

It does not verify prose, figures, or cross-references -- only that the code
executes.  That covers the common failure mode (a renamed column or a moved
file) and not the rare one (a broken @ref).

Usage:
    python tools/check_notebook_chunks.py notebooks/01_longitudinal_model.qmd
"""
from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

CHUNK = re.compile(r"^```\{(\w+)[^}]*\}\n(.*?)^```", re.S | re.M)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    qmd = Path(argv[1])
    if not qmd.exists():
        print(f"not found: {qmd}")
        return 2

    chunks = CHUNK.findall(qmd.read_text())
    if not chunks:
        print(f"no code chunks found in {qmd}")
        return 1
    langs = {lang for lang, _ in chunks}
    if len(langs) > 1:
        print(f"mixed-engine notebook ({langs}); this tool handles one engine")
        return 1
    lang = langs.pop()
    code = "\n".join(body for _, body in chunks)

    suffix = {"r": ".R", "python": ".py"}.get(lang)
    if suffix is None:
        print(f"unsupported engine: {lang}")
        return 1
    with tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False) as fh:
        fh.write(code)
        script = fh.name

    # Chunks use paths relative to the notebook, so run from its directory.
    cmd = ["Rscript", script] if lang == "r" else [sys.executable, script]
    print(f"{qmd.name}: {len(chunks)} {lang} chunks -> {' '.join(cmd[:1])}")
    proc = subprocess.run(cmd, cwd=qmd.parent, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stdout.write(proc.stdout[-4000:])
        sys.stderr.write(proc.stderr[-4000:])
        print(f"\nFAILED: {qmd} (exit {proc.returncode})")
        return proc.returncode
    print(f"OK: all {len(chunks)} chunks executed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
