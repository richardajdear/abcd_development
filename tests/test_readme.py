"""Guard the README's promises.

A README is the one document nobody re-reads critically once written, so its
claims rot silently. These tests check only claims that are mechanically
checkable: that every path it points at exists, and that the quickstart it
tells a new reader to run actually resolves a run directory. They deliberately
check no scientific numbers — those are guarded at source by
`test_docs_provenance.py`.

The `make help` test exists because the README's claim that it prints the
resolved run directory was false when written: GNU Make 3.81 (the macOS
default) does not propagate exported variables into `$(shell ...)`, so
`PYTHONPATH` was missing there and the error was swallowed by `2>/dev/null`,
leaving a blank line that read as "no run configured".
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"

# Roots an inline `code/path.py` span may be relative to.
_SEARCH = ["", "src/abcd", "configs", "tools", "tests", "R", "notebooks", "docs"]


def test_relative_links_resolve():
    """Every non-http markdown link points at a file that exists."""
    links = re.findall(r"\[[^\]]+\]\((?!https?:)([^)#]+)", README.read_text())
    missing = sorted({p for p in links if not (ROOT / p).exists()})
    assert not missing, f"README links to absent paths: {missing}"


def test_referenced_files_exist():
    """Inline code spans that look like repo files resolve somewhere sensible."""
    spans = set(re.findall(r"`([A-Za-z0-9_./]+\.(?:py|R|md|yaml|qmd))`", README.read_text()))
    missing = sorted(s for s in spans
                     if not any((ROOT / c / s).exists() for c in _SEARCH))
    assert not missing, f"README names files that do not exist: {missing}"


def test_no_stale_project_identity():
    """The README must not present itself as the snRNA-seq repo it was forked from."""
    head = "\n".join(README.read_text().splitlines()[:6])
    assert "# abcd_development" in head, "README title is not this project"
    assert "transcriptional_maturation" not in head, (
        "README opens as transcriptional_maturation; that project is an *input* "
        "to this one and belongs in the related-work section only"
    )


def test_every_config_the_readme_names_exists():
    named = set(re.findall(r"`(ct_[a-z0-9_]+)(?:\.yaml)?`", README.read_text()))
    # Brace expansions like ct_70_noglobal_mv{2,3,4} are prose, not literals.
    named = {n for n in named if "{" not in n}
    missing = sorted(n for n in named if not (ROOT / "configs" / f"{n}.yaml").exists())
    assert not missing, f"README names absent configs: {missing}"


def test_readme_test_count_is_current():
    """The README's advertised test count must match the suite as collected.

    A hardcoded count in prose is the archetype of a number that rots silently;
    this makes adding a test fail loudly until the claim is updated.

    The count comes from a subprocess collection of the whole tests/ directory,
    not from the running session, so the check gives the same answer whether the
    full suite or this one file was invoked.
    """
    claimed = set(re.findall(r"# (\d+) tests", README.read_text()))
    assert claimed, "README no longer states a test count"

    proc = subprocess.run(
        [os.environ.get("PYTEST_PY", "python"), "-m", "pytest",
         str(ROOT / "tests"), "-q", "--collect-only", "-p", "no:cacheprovider"],
        cwd=ROOT, capture_output=True, text=True, timeout=180,
        env={**os.environ, "PYTHONPATH": f"{ROOT / 'src'}:{os.environ.get('PYTHONPATH', '')}"},
    )
    m = re.search(r"(\d+) tests? collected", proc.stdout)
    assert m, f"could not read collected count:\n{proc.stdout[-2000:]}"
    actual = m.group(1)
    assert claimed == {actual}, (
        f"README claims {sorted(claimed)} tests; tests/ collects {actual}"
    )


@pytest.mark.skipif(shutil.which("make") is None, reason="make not available")
def test_make_help_resolves_a_run_dir():
    """`make help` must print a run directory, not a blank line.

    Requires no editable install: the Makefile is expected to put src/ on the
    path itself, which is what the README promises a new reader.
    """
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env["ABCD_CONFIG"] = "ct_70_noglobal_mv2_genetic"
    out = subprocess.run(["make", "help"], cwd=ROOT, env=env,
                         capture_output=True, text=True, timeout=120).stdout
    line = next((l for l in out.splitlines() if l.startswith("run dir")), None)
    assert line is not None, f"`make help` printed no run dir line:\n{out}"
    value = line.split("=", 1)[1].strip()
    assert value, (
        "`make help` printed an empty run dir. The README tells readers this "
        f"resolves the active config. Full output:\n{out}"
    )
    assert "thickness_dsk_70" in value, f"unexpected run dir: {value}"
