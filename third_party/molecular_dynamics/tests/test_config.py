"""Tests for config loading and validation."""

import tempfile
from pathlib import Path

import pytest
import yaml

from lib.config import load_config, get_production_nsteps, get_save_nsteps


def test_load_default_config():
    cfg = load_config()
    assert cfg["md"]["force_field"] == "amber99sb-ildn"
    assert cfg["md"]["temperature"] == 310
    assert "H1" in cfg["analysis"]["cdr_regions"]


def test_load_missing_file():
    with pytest.raises(FileNotFoundError):
        load_config("/nonexistent/path.yaml")


def test_validation_missing_md_section():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        yaml.dump({"analysis": {"antibody_chains": ["A"], "antigen_chains": ["C"], "cdr_regions": {"H1": {}}}}, f)
        f.flush()
        with pytest.raises(ValueError, match="'md' section"):
            load_config(f.name)


def test_get_production_nsteps():
    cfg = load_config()
    nsteps = get_production_nsteps(cfg)
    # 50 ns / 0.002 ps = 25,000,000
    assert nsteps == 25_000_000


def test_get_save_nsteps():
    cfg = load_config()
    nsteps = get_save_nsteps(cfg)
    # 100 ps / 0.002 ps = 50,000
    assert nsteps == 50_000
