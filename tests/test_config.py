"""Unit tests for RunConfig: round-trip fidelity and hash stability."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abcd.config import ConfigError, RunConfig  # noqa: E402


def test_dict_roundtrip_preserves_hash():
    cfg = RunConfig(metric="thickness", parcellation="dsk", min_visits=2)
    assert RunConfig.from_dict(cfg.to_dict()) == cfg
    assert RunConfig.from_dict(cfg.to_dict()).hash == cfg.hash


def test_yaml_roundtrip_preserves_hash(tmp_path):
    cfg = RunConfig(metric="t1t2_ratio", hemisphere="lh", age_basis="spline")
    p = cfg.to_yaml(tmp_path / "c.yaml")
    assert RunConfig.from_yaml(p) == cfg


def test_hash_is_literal_and_stable():
    """A fixed config must hash to a fixed value across sessions.

    If this fails after a deliberate schema change, update the literal - but
    every previously written out/<hash>/ directory now refers to a different
    config, so the change should be noted in the changelog.
    """
    cfg = RunConfig()
    assert cfg.hash == "5c104aff8f53", f"default-config hash drifted: {cfg.hash}"


def test_hash_ignores_note_but_tracks_values():
    base = RunConfig()
    assert base.replace(note="exploratory").hash == base.hash
    assert base.replace(min_visits=3).hash != base.hash
    assert base.replace(metric="volume").hash != base.hash


def test_hash_independent_of_field_insertion_order():
    d = RunConfig().to_dict()
    shuffled = dict(reversed(list(d.items())))
    assert RunConfig.from_dict(shuffled).hash == RunConfig.from_dict(d).hash


def test_frozen():
    cfg = RunConfig()
    with pytest.raises(Exception):
        cfg.metric = "volume"  # type: ignore[misc]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"metric": "cortical_thickness"},  # not in METRICS - the old name
        {"parcellation": "aparc"},         # old FreeSurfer-style name
        {"release": "5.0"},
        {"qc_policy": "euler"},
        {"min_visits": 0},
        {"age_basis": "spline", "spline_df": 1},
        {"visits": ()},
    ],
)
def test_invalid_values_rejected(kwargs):
    with pytest.raises(ConfigError):
        RunConfig(**kwargs)


def test_unknown_key_rejected():
    d = RunConfig().to_dict()
    d["lh_mean"] = True  # a 5.1-ism that is no longer a config field
    with pytest.raises(ConfigError):
        RunConfig.from_dict(d)


def test_hand_edited_yaml_detected(tmp_path):
    cfg = RunConfig()
    p = cfg.to_yaml(tmp_path / "c.yaml")
    text = p.read_text().replace("min_visits: 2", "min_visits: 3")
    p.write_text(text)
    with pytest.raises(ConfigError, match="recorded hash"):
        RunConfig.from_yaml(p)


def test_run_id_contains_hash():
    cfg = RunConfig()
    assert cfg.hash in cfg.run_id
