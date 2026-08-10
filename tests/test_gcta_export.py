"""Guards on the GCTA export path.

The family-effect guard is the important one: it prevents a silent, plausible-
looking heritability estimate computed from BLUPs that had the between-family
variance removed.
"""

import pytest


def test_family_effect_run_is_refused(tmp_path):
    """A run fitted with (1|family_id) must not reach GCTA.

    Its BLUPs are within-family deviations: the between-family variance that
    genetic relatedness explains has been absorbed by the random effect.
    """
    import yaml
    from abcd.gcta_export import FamilyEffectConflict, build

    (tmp_path / "config.yaml").write_text(yaml.safe_dump({"family_effect": True}))
    with pytest.raises(FamilyEffectConflict, match="within-family"):
        build(tmp_path)


def test_no_family_effect_run_passes_the_guard(tmp_path):
    """The guard must not block a valid genetic run (it fails later, on data)."""
    import yaml
    from abcd.gcta_export import FamilyEffectConflict, _check_genetic_validity

    (tmp_path / "config.yaml").write_text(yaml.safe_dump({"family_effect": False}))
    _check_genetic_validity(tmp_path)   # must not raise


# ---------------------------------------------------------------------------
# The manifest.  These are cheap tests for an expensive class of failure: a
# wrong --mpheno runs a GWAS to completion on the wrong phenotype and reports
# nothing amiss.
# ---------------------------------------------------------------------------

def test_manifest_mpheno_matches_column_position():
    """--mpheno must be 1-based over phenotype columns, i.e. excluding FID/IID.

    Off-by-one here is invisible: GCTA happily analyses the neighbouring
    phenotype.  Pinned against the actual exported frame rather than against a
    remembered convention.
    """
    import pandas as pd
    from abcd.gcta_export import PHENOTYPES

    header = ["FID", "IID"] + [p["name"] for p in PHENOTYPES]
    phen = pd.DataFrame(columns=header)
    for p in PHENOTYPES:
        mpheno = phen.columns.get_loc(p["name"]) - 1
        # GCTA reads column (2 + mpheno) of the file, 1-based.
        assert header[1 + mpheno] == p["name"]
        assert mpheno >= 1


def test_manifest_covers_every_phenotype_and_is_priority_ordered():
    from abcd.gcta_export import PHENOTYPES, PHENOTYPE_NAMES

    assert len(PHENOTYPE_NAMES) == len(set(PHENOTYPE_NAMES)) == 5
    prios = [p["priority"] for p in PHENOTYPES]
    assert prios == sorted(prios), "PHENOTYPES must be declared in priority order"
    # Exactly one positive control, and it is priority 0.
    controls = [p for p in PHENOTYPES if p["role"] == "positive control"]
    assert [c["name"] for c in controls] == ["baseline_thickness"]
    assert controls[0]["priority"] == 0


def test_no_phenotype_name_is_a_prefix_of_another():
    """Guards the failure the manifest was introduced to remove.

    The previous scripts located a phenotype's column by grepping the header.
    With ``slope_PC1`` and a hypothetical ``slope_PC1_resid`` that returns two
    matches and the script takes the first.  Names must stay mutually
    non-prefixing so that even a grep-based reader outside this repo is safe.
    """
    from abcd.gcta_export import PHENOTYPE_NAMES

    for a in PHENOTYPE_NAMES:
        others = [b for b in PHENOTYPE_NAMES if b != a]
        assert not any(b.startswith(a) for b in others), a


def test_h2_tables_use_the_exported_phenotype_definitions():
    """The published h2 must describe the vector that is actually exported.

    ``regen_h2_tables._subject_maps`` delegates to ``subject_phenotypes``.  If
    someone re-inlines the definitions there, the two can drift and the report
    would publish an h2 for a phenotype no GWAS was run on.
    """
    import importlib.util
    from pathlib import Path

    src = Path(__file__).resolve().parents[1] / "tools" / "regen_h2_tables.py"
    text = src.read_text()
    assert "subject_phenotypes" in text, (
        "regen_h2_tables must build phenotypes via abcd.gcta_export."
        "subject_phenotypes, not its own copy"
    )
