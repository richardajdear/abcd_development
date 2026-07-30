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
