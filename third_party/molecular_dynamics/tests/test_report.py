# third_party/molecular_dynamics/tests/test_report.py
"""Tests for report generation."""

import json

import pytest

from lib.report import generate_report


@pytest.fixture
def mock_output_dir(tmp_path):
    for ph in ("7.4", "6.0"):
        d = tmp_path / "HE1H" / f"pH_{ph}" / "analysis"
        d.mkdir(parents=True)
        summary = {
            "rmsd_rmsd_mean": 2.0 if ph == "7.4" else 2.5,
            "hbond_n_hbond_mean": 5.0 if ph == "7.4" else 3.0,
        }
        (d / "summary.json").write_text(json.dumps(summary))

    comp = {
        "variant": "HE1H",
        "base_ph": 7.4,
        "target_ph": 6.0,
        "deltas": {"delta_rmsd": 0.5, "delta_n_hbond": -2.0},
    }
    (tmp_path / "HE1H" / "ph_comparison.json").write_text(json.dumps(comp))

    d2 = tmp_path / "WT" / "pH_7.4" / "analysis"
    d2.mkdir(parents=True)
    (d2 / "summary.json").write_text(json.dumps({"rmsd_rmsd_mean": 1.5}))
    return tmp_path


def test_generate_report(mock_output_dir):
    report_path = generate_report(mock_output_dir)
    assert report_path.exists()

    import pandas as pd
    df = pd.read_csv(report_path)
    assert len(df) == 3  # HE1H pH7.4, HE1H pH6.0, WT pH7.4
    assert "variant_id" in df.columns

    he1h_6 = df[(df["variant_id"] == "HE1H") & (df["ph"] == 6.0)]
    assert len(he1h_6) == 1
    assert he1h_6.iloc[0]["delta_rmsd"] == pytest.approx(0.5)
