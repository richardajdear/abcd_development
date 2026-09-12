"""Guard the report's provenance invariants.

These do not check any *number*; they check that every number in the report can
be regenerated. The failure mode they exist to prevent is the one that actually
happened: figures and tables drawn in ad-hoc cells drifted out of sync with the
corrected numbers underneath them, and two cited figures did not exist at all.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
REPORT = DOCS / "REPORT_7.0.md"


def _cited(pattern: str) -> set[str]:
    return set(re.findall(pattern, REPORT.read_text()))


@pytest.mark.skipif(not REPORT.exists(), reason="report not present")
def test_every_cited_figure_exists():
    missing = sorted(f for f in _cited(r"([A-Za-z0-9_]+\.png)")
                     if not (DOCS / "figures" / f).exists())
    assert not missing, f"cited but absent from docs/figures/: {missing}"


@pytest.mark.skipif(not REPORT.exists(), reason="report not present")
def test_every_cited_table_exists():
    """Every table cited as a source must exist.

    The retirement paragraph at the end of section 14 deliberately names tables
    that no longer exist, so it is excluded -- naming a retired table is the
    opposite of the defect this guards against.
    """
    text = REPORT.read_text().split("have been retired as\nsuperseded")[0]
    cited = set(re.findall(r"`([A-Za-z0-9_]+\.csv)`", text))
    missing = sorted(f for f in cited if not (DOCS / f).exists())
    assert not missing, f"cited but absent from docs/: {missing}"


def test_every_figure_has_a_generator():
    """No figure may be an orphan -- an orphan cannot be re-derived."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(ROOT / "tools"))
    import regen_brain_maps
    import regen_report_figures

    generated = set(regen_report_figures.FIGURES) | set(regen_brain_maps.FIGURES)
    on_disk = {p.stem for p in (DOCS / "figures").glob("*.png")}
    assert not (on_disk - generated), \
        f"figures with no generator: {sorted(on_disk - generated)}"


def test_no_serial_level_scanner_claims_outside_the_withdrawal():
    """``mri_info_deviceserial`` is in neither release tree.

    Section 9 may *discuss* the withdrawal of serial-level results, but must not
    assert a serial-level count or correlation as a finding. The report
    previously did both, two paragraphs apart.
    """
    text = REPORT.read_text()
    section = text.split("## 9. Site and scanner effects")[1].split("\n## ")[0]
    offending = [ln for ln in section.splitlines()
                 if re.search(r"\d+ (?:distinct )?(?:scanner )?serials", ln)]
    assert not offending, f"serial-level claim stated as fact: {offending}"


HPC_README = ROOT / "hpc" / "README_HPC.md"


@pytest.mark.skipif(not HPC_README.exists(), reason="hpc/README_HPC.md not present")
def test_no_null_ahba_map_claim_in_hpc_readme():
    """The map-level AHBA result is positive; prose must not say otherwise.

    hpc/README_HPC.md asserted "no correlation between the group-mean developmental
    change map and AHBA C1-C3" long after that null had been superseded. The
    null was real but came from the GLOBAL-ADJUSTED age coefficient; the settled
    no-global specification gives 8 of 18 map x component pairs surviving the
    spin test.

    It sat undetected because no test guarded it: the provenance tests check
    that every cited figure and table EXISTS, not that prose agrees with what
    the tables say. This closes that specific gap the cheap way -- a direction
    check, not a numeric one -- because a stale sign is the error that misleads a
    reader, and the file is handover documentation for people who will not
    re-derive it.
    """
    import pandas as pd

    table = DOCS / "ahba_vs_maps_noglobal.csv"
    if not table.exists():
        pytest.skip("ahba_vs_maps_noglobal.csv not present")
    n_sig = int((pd.read_csv(table).p_spin < 0.05).sum())

    text = HPC_README.read_text()
    # Exclude the paragraph documenting the withdrawal, which quotes the old
    # wording deliberately.
    body = text.split("An earlier version of this file claimed")[0]
    # Collapse newlines before matching: markdown hard-wraps prose, so the
    # original claim put "no correlation between the" and "AHBA C1-C3" on
    # separate lines and a per-line regex missed it entirely.  Requiring the
    # subject within ~90 chars of the negation also avoids matching unrelated
    # sentences such as "No heritability estimate, no association statistic".
    flat = " ".join(body.split())
    offending = re.findall(
        r"no (?:correlation|association)(?:(?!\bno\b).){0,90}?(?:AHBA|\bC1\b|\bC2\b|\bC3\b)",
        flat, re.I)
    assert not offending or n_sig == 0, (
        f"hpc/README_HPC.md asserts a null map-level AHBA result, but "
        f"{n_sig} of 18 pairs survive the spin test: {offending}"
    )


@pytest.mark.skipif(not HPC_README.exists(), reason="hpc/README_HPC.md not present")
def test_hpc_readme_figure_references_exist():
    """Every ``figures/x.png`` hpc/README_HPC.md names must be in docs/figures/.

    It cited ``figures/abcd_spatial_nulls.png``, which no longer exists -- a
    dangling reference is the visible symptom of a claim that outlived its
    evidence, so it is worth failing on directly.
    """
    missing = [f for f in set(re.findall(r"figures/[\w.-]+\.png", HPC_README.read_text()))
               if not (DOCS / f).exists()]
    assert not missing, f"hpc/README_HPC.md cites missing figures: {sorted(missing)}"


def test_artifact_audit_covers_every_report_deliverable():
    """`make audit-artifacts` must not silently miss a deliverable.

    The audit is only as good as its glob list, and the first version of that
    list omitted ``hpc/*.md`` -- so it was blind to ``hpc/README.md``, the exact
    file whose absence from the artifact tray prompted writing it. This asserts
    the manifest covers everything the report cites plus the hpc/ handover docs,
    so adding a deliverable in a new location fails here rather than going
    unnoticed until someone reads a stale copy.
    """
    import sys
    sys.path.insert(0, str(ROOT / "tools"))
    import audit_artifacts

    covered = {r["path"] for r in audit_artifacts.collect()}

    required = {"README.md", "hpc/README_HPC.md", "docs/REPORT_7.0.md"}
    # Everything the report cites as a table or figure.
    required |= {f"docs/{t}" for t in _cited(r"`docs/([\w.-]+\.csv)`")}
    required |= {f"docs/figures/{f}" for f in
                 {Path(p).name for p in _cited(r"docs/figures/([\w.-]+\.png)")}}
    # Both pipeline entry points: these were the files never saved at all.
    required |= {"hpc/run_all.sh", "hpc/config.sh", "tools/local_test.sh"}

    present = {p for p in required if (ROOT / p).exists()}
    missing = sorted(present - covered)
    assert not missing, (
        f"audit_artifacts.PATTERNS does not cover: {missing} -- "
        "these deliverables would never be checked against the tray"
    )


def test_artifact_audit_qualifies_colliding_basenames():
    """Colliding basenames must be path-qualified, not silently shared.

    ``save_artifacts`` keys on basename, so ``README.md`` and ``hpc/README.md``
    would otherwise map to one artifact -- and an audit would compare the wrong
    pair, which is how the missing hpc/ pipeline first presented as a merely
    "stale README".
    """
    import sys
    sys.path.insert(0, str(ROOT / "tools"))
    import audit_artifacts

    rows = audit_artifacts.collect()
    names = [r["artifact_name"] for r in rows]
    dupes = sorted({n for n in names if names.count(n) > 1})
    assert not dupes, f"artifact names collide: {dupes}"
