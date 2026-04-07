# third_party/molecular_dynamics/tests/test_gromacs_wrapper.py
"""Tests for GromacsWrapper — mock subprocess calls."""

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from lib.config import load_config
from lib.gromacs_wrapper import GromacsWrapper, GromacsError


@pytest.fixture
def cfg():
    return load_config()


@pytest.fixture
def wrapper(cfg, tmp_path):
    return GromacsWrapper(cfg, tmp_path)


@patch("lib.gromacs_wrapper.subprocess.run")
def test_run_success(mock_run, wrapper):
    mock_run.return_value = subprocess.CompletedProcess(
        args=["gmx", "help"], returncode=0, stdout="ok", stderr=""
    )
    result = wrapper._run(["help"])
    assert result.returncode == 0
    mock_run.assert_called_once()


@patch("lib.gromacs_wrapper.subprocess.run")
def test_run_failure_raises(mock_run, wrapper):
    mock_run.return_value = subprocess.CompletedProcess(
        args=["gmx", "pdb2gmx"], returncode=1, stdout="", stderr="Fatal error"
    )
    with pytest.raises(GromacsError, match="pdb2gmx failed"):
        wrapper._run(["pdb2gmx", "-f", "test.pdb"])


def test_render_mdp_minimization(wrapper):
    mdp_path = wrapper._render_minimization_mdp()
    assert mdp_path.exists()
    content = mdp_path.read_text()
    assert "integrator" in content
    assert "5000" in content  # max_steps from config


def test_render_mdp_production(wrapper):
    mdp_path = wrapper._render_mdp(
        "production.mdp.j2", "test_prod.mdp",
        dt=0.002, nsteps=25000000, save_nsteps=50000,
        temperature=310, cphmd_enabled=False,
    )
    content = mdp_path.read_text()
    assert "25000000" in content
    assert "lambda-dynamics" not in content  # CpHMD disabled


def test_render_mdp_production_with_cphmd(wrapper):
    mdp_path = wrapper._render_mdp(
        "production.mdp.j2", "test_cphmd.mdp",
        dt=0.002, nsteps=25000000, save_nsteps=50000,
        temperature=310,
        cphmd_enabled=True, cphmd_ph=6.0, cphmd_n_groups=1,
        cphmd_update_nst=100,
        cphmd_groups=[{
            "type": "histidine",
            "name": "HIS_42",
            "index_group": "HIS_42",
            "barrier": 5.0,
            "init_lambda": 0.5,
        }],
    )
    content = mdp_path.read_text()
    assert "lambda-dynamics" in content
    assert "6.0" in content
