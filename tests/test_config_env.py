"""Guards on $ABCD_CONFIG resolution and multi-root release lookup.

Two classes of bug are pinned here, both of which have actually bitten:

1. **Silent wrong-run selection.** The whole point of ``$ABCD_CONFIG`` is that
   one export drives every step. If resolution ever falls back to a default
   config, or picks arbitrarily among candidates, the pipeline produces a
   real-looking run of the wrong phenotype. Every failure path must raise.

2. **Release lookup collapsing to one root.** 7.0 lives inside the repo
   (gitignored) and 5.1 outside it. A single-root resolver silently loses one
   of them -- which is how ``core/imaging/mri_y_smr_thk_dsk.csv`` (a 5.1 path)
   came to be looked up under the 7.0 root.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from abcd import paths
from abcd.config import (
    CONFIG_ENV_VAR,
    ConfigError,
    active_config,
    active_run_dir,
    resolve_config,
)


@pytest.fixture
def no_env(monkeypatch):
    monkeypatch.delenv(CONFIG_ENV_VAR, raising=False)
    monkeypatch.delenv("ABCD_ROOT", raising=False)


# ---------------------------------------------------------------- resolution
def test_env_var_resolves_bare_stem(monkeypatch, no_env):
    monkeypatch.setenv(CONFIG_ENV_VAR, "ct_70_genetic")
    assert resolve_config().name == "ct_70_genetic.yaml"


@pytest.mark.parametrize(
    "spec", ["ct_70_genetic", "ct_70_genetic.yaml", "configs/ct_70_genetic.yaml"]
)
def test_spellings_are_equivalent(spec, no_env):
    assert resolve_config(spec) == (paths.CONFIG_DIR / "ct_70_genetic.yaml").resolve()


def test_explicit_spec_beats_env(monkeypatch, no_env):
    """An argument must override the export, or --config is a lie."""
    monkeypatch.setenv(CONFIG_ENV_VAR, "ct_70_genetic")
    assert resolve_config("ct_70_baseline").stem == "ct_70_baseline"


def test_unset_raises_and_lists_options(no_env):
    with pytest.raises(ConfigError) as e:
        resolve_config()
    # the message must be actionable: name the variable and the choices
    assert CONFIG_ENV_VAR in str(e.value) and "ct_70_genetic" in str(e.value)


def test_unknown_config_raises_rather_than_defaulting(no_env):
    with pytest.raises(ConfigError):
        resolve_config("no_such_config")


def test_run_dir_matches_config_hash(monkeypatch, no_env):
    """active_run_dir must agree with the config's own run_id.

    These are two paths to the same identity; if they diverge, R (which shells
    out to abcd.run_dir) and Python fit different runs.
    """
    monkeypatch.setenv(CONFIG_ENV_VAR, "ct_70_genetic")
    assert active_run_dir(must_exist=False).name == active_config().run_id


def test_missing_run_dir_raises_with_next_step(monkeypatch, no_env, tmp_path):
    monkeypatch.setenv(CONFIG_ENV_VAR, "ct_70_genetic")
    monkeypatch.setattr(paths, "OUT_DIR", tmp_path / "nonexistent")
    with pytest.raises(ConfigError, match="assemble"):
        active_run_dir()


# ------------------------------------------------------------- release roots
def test_repo_root_searched_first(no_env):
    """A release vendored in the repo must resolve with no env var set."""
    assert paths.REPO_ROOT.resolve() in paths.abcd_roots()


def test_both_releases_resolve_from_different_roots(no_env):
    """The regression: 5.1 and 7.0 need not share a parent directory."""
    d51, d70 = paths.release_dir("5.1"), paths.release_dir("7.0")
    assert d51.exists() and d70.exists()
    assert "5.1" in d51.name and "7.0" in d70.name


def test_seven_zero_is_not_looked_up_under_a_five_one_layout(no_env):
    """Pin the exact bug: 7.0 tables are .tsv under y/, not core/imaging/*.csv."""
    d70 = paths.release_dir("7.0")
    assert not (d70 / "core" / "imaging").exists()
    assert paths.find_table(d70, "mr_y_smri__thk__dsk").suffix == ".tsv"


def test_bad_abcd_root_raises_rather_than_being_ignored(monkeypatch, no_env):
    """A typo'd ABCD_ROOT must fail loudly, not fall through to a fallback."""
    monkeypatch.setenv("ABCD_ROOT", "/nonexistent/abcd/tree")
    with pytest.raises(paths.DataRootError):
        paths.abcd_roots()


def test_explicit_abcd_root_is_searched(monkeypatch, no_env, tmp_path):
    monkeypatch.setenv("ABCD_ROOT", str(tmp_path))
    assert tmp_path.resolve() in paths.abcd_roots()


def test_loose_files_resolve_outside_the_preferred_root(no_env):
    """Regression: the legacy QC list lives beside 5.1, not in the repo.

    When the preferred root became the repo, a single-root lookup stopped
    finding this file -- and because the legacy policy *degrades* rather than
    failing, QC silently passed every scan (22,854 instead of 20,162). Any
    auxiliary input must be searched across all roots.
    """
    from abcd import qc

    assert paths.find_in_roots(qc.LEGACY_QC_FILE) is not None
    assert len(qc.load_legacy_list()) > 0


def test_find_in_roots_returns_none_when_absent(no_env):
    assert paths.find_in_roots("definitely_not_a_real_file.csv") is None
