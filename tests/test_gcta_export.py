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
