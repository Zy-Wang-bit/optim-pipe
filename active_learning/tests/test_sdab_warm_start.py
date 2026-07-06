from pathlib import Path

import pandas as pd
import pytest

from active_learning.cli import main
from active_learning.tools import sequence_sha256


def _write_warm_start_inputs(tmp_path: Path):
    wt_sequence = "A" * 25 + "G" + "C"
    fasta = tmp_path / "wt.fasta"
    config = tmp_path / "config.yaml"
    legacy = tmp_path / "legacy.csv"
    direct = tmp_path / "direct.csv"
    out_csv = tmp_path / "training.csv"
    fasta.write_text(f">A\n{wt_sequence}\n", encoding="utf-8")
    config.write_text("system_id: sdab\nround_id: R5\nnaming:\n  fasta_record_id: A\n", encoding="utf-8")
    pd.DataFrame(
        [
            {
                "name": "WT",
                "mutations_unified": "",
                "log_pH74": 0.0,
                "log_pH6": 0.0,
                "source": "legacy",
                "round_id": "R4",
                "group": "WT",
                "collapsed": False,
                "leverage_flag": False,
            },
            {
                "name": "m1",
                "mutations_unified": "HG26H",
                "log_pH74": -1.0,
                "log_pH6": -2.0,
                "source": "legacy",
                "round_id": "R4",
                "group": "mutant",
                "collapsed": True,
                "leverage_flag": True,
            },
        ]
    ).to_csv(legacy, index=False)
    pd.DataFrame(
        [
            {
                "mutation_string": "HG26H",
                "A7": 1.2,
                "A6": -0.2,
                "confirmed": True,
                "round_id": "R5",
                "batch_id": "b1",
                "assay_id": "a1",
            }
        ]
    ).to_csv(direct, index=False)
    return fasta, config, legacy, direct, out_csv, wt_sequence


def test_sdab_warm_start_retains_legacy_rows_as_excluded(tmp_path: Path):
    fasta, config, legacy, _, out_csv, wt_sequence = _write_warm_start_inputs(tmp_path)

    main(["build-sdab-warm-start", "--config-yaml", str(config), "--wt-fasta", str(fasta), "--legacy-r4-training", str(legacy), "--out-csv", str(out_csv)])
    rows = pd.read_csv(out_csv)

    assert len(rows) == 2
    assert set(rows["label_source_type"]) == {"legacy_training_summary"}
    assert not rows["usable_for_training"].any()
    assert not rows["label_is_direct_measurement"].any()
    assert set(rows["experiment_status"]) == {"excluded"}
    assert set(rows["round_id"]) == {"R4"}
    assert sequence_sha256(wt_sequence) in set(rows["sequence_sha256"])
    assert {"system_id", "sequence", "sequence_sha256", "mutation_string", "round_id", "batch_id", "assay_id", "utility_recomputed"}.issubset(rows.columns)


def test_sdab_warm_start_infers_legacy_round_from_path(tmp_path: Path):
    fasta, config, legacy, _, out_csv, _ = _write_warm_start_inputs(tmp_path)
    legacy_dir = tmp_path / "experiments" / "sdab" / "R4" / "data"
    legacy_dir.mkdir(parents=True)
    legacy_r4 = legacy_dir / "training_data.csv"
    rows = pd.read_csv(legacy).drop(columns=["round_id"])
    rows.to_csv(legacy_r4, index=False)

    main(["build-sdab-warm-start", "--config-yaml", str(config), "--wt-fasta", str(fasta), "--legacy-r4-training", str(legacy_r4), "--out-csv", str(out_csv)])

    assert set(pd.read_csv(out_csv)["round_id"]) == {"R4"}


def test_sdab_warm_start_allows_confirmed_direct_wet_labels(tmp_path: Path):
    fasta, config, legacy, direct, out_csv, _ = _write_warm_start_inputs(tmp_path)

    main(
        [
            "build-sdab-warm-start",
            "--config-yaml",
            str(config),
            "--wt-fasta",
            str(fasta),
            "--legacy-r4-training",
            str(legacy),
            "--direct-wet-labels",
            str(direct),
            "--out-csv",
            str(out_csv),
        ]
    )
    direct_rows = pd.read_csv(out_csv).query("label_source_type == 'raw_elisa'")

    assert len(direct_rows) == 1
    assert direct_rows.iloc[0]["A7"] == 1.0
    assert direct_rows.iloc[0]["A6"] == 0.0
    assert direct_rows.iloc[0]["label_is_direct_measurement"]
    assert direct_rows.iloc[0]["usable_for_training"]
    assert direct_rows.iloc[0]["experiment_status"] == "measured"
    assert direct_rows.iloc[0]["round_id"] == "R5"


def test_sdab_warm_start_has_no_pseudo_label_promotion_flag(tmp_path: Path):
    fasta, config, legacy, _, out_csv, _ = _write_warm_start_inputs(tmp_path)

    with pytest.raises(SystemExit):
        main(
            [
                "build-sdab-warm-start",
                "--config-yaml",
                str(config),
                "--wt-fasta",
                str(fasta),
                "--legacy-r4-training",
                str(legacy),
                "--allow-pseudo-labels",
                "--out-csv",
                str(out_csv),
            ]
        )


def test_sdab_warm_start_refuses_to_overwrite_input(tmp_path: Path):
    fasta, config, legacy, _, _, _ = _write_warm_start_inputs(tmp_path)

    with pytest.raises(ValueError, match="Refusing to overwrite input file"):
        main(
            [
                "build-sdab-warm-start",
                "--config-yaml",
                str(config),
                "--wt-fasta",
                str(fasta),
                "--legacy-r4-training",
                str(legacy),
                "--out-csv",
                str(legacy),
            ]
        )
