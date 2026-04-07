# third_party/molecular_dynamics/tests/test_protonation.py
"""Tests for protonation manager."""

import tempfile
from pathlib import Path

import pytest

from lib.protonation import detect_his_residues, filter_titratable, build_cphmd_groups


# Minimal PDB with 2 His residues for testing
_MINI_PDB = """\
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   1.000   1.000  1.00  0.00           C
ATOM      3  N   HIS A  42       3.000   1.000   1.000  1.00  0.00           N
ATOM      4  CA  HIS A  42       4.000   1.000   1.000  1.00  0.00           C
ATOM      5  N   HIS B  50       5.000   1.000   1.000  1.00  0.00           N
ATOM      6  CA  HIS B  50       6.000   1.000   1.000  1.00  0.00           C
ATOM      7  N   GLU B  51       7.000   1.000   1.000  1.00  0.00           N
ATOM      8  CA  GLU B  51       8.000   1.000   1.000  1.00  0.00           C
END
"""


@pytest.fixture
def mini_pdb(tmp_path):
    p = tmp_path / "test.pdb"
    p.write_text(_MINI_PDB)
    return p


def test_detect_his_residues(mini_pdb):
    his_list = detect_his_residues(mini_pdb)
    assert len(his_list) == 2
    assert his_list[0] == {"chain": "A", "resid": 42, "resname": "HIS"}
    assert his_list[1] == {"chain": "B", "resid": 50, "resname": "HIS"}


def test_filter_titratable_empty_config(mini_pdb):
    his_list = detect_his_residues(mini_pdb)
    filtered = filter_titratable(his_list, [])
    assert len(filtered) == 2  # empty config = keep all


def test_filter_titratable_specific(mini_pdb):
    his_list = detect_his_residues(mini_pdb)
    filtered = filter_titratable(his_list, [{"chain": "A", "resid": 42}])
    assert len(filtered) == 1
    assert filtered[0]["resid"] == 42


def test_build_cphmd_groups():
    his_list = [
        {"chain": "A", "resid": 42, "resname": "HIS"},
        {"chain": "B", "resid": 50, "resname": "HIS"},
    ]
    groups = build_cphmd_groups(his_list, barrier_height=5.0)
    assert len(groups) == 2
    assert groups[0]["name"] == "HIS_A42"
    assert groups[0]["barrier"] == 5.0
    assert groups[1]["name"] == "HIS_B50"
