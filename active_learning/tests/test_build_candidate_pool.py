from pathlib import Path

import pandas as pd
import pytest

from active_learning.cli import main


def test_build_candidate_pool_cli_writes_full_single_pool(tmp_path: Path):
    fasta = tmp_path / "tiny.fasta"
    config = tmp_path / "config.yaml"
    out_csv = tmp_path / "candidates.csv"
    fasta.write_text(">A\nACD\n", encoding="utf-8")
    config.write_text(
        """
naming:
  fasta_record_id: A
  mutation_chain_label: H
candidate_pool:
  excluded_positions: [1]
  excluded_mutations: [HC2A]
""".strip(),
        encoding="utf-8",
    )

    main(
        [
            "build-candidate-pool",
            "--config-yaml",
            str(config),
            "--wt-fasta",
            str(fasta),
            "--system-id",
            "sdab",
            "--round-id",
            "R5",
            "--candidate-mode",
            "full_single",
            "--out-csv",
            str(out_csv),
        ]
    )
    rows = pd.read_csv(out_csv)

    assert len(rows) == 37
    assert set(rows["system_id"]) == {"sdab"}
    assert set(rows["round_id"]) == {"R5"}
    assert "HC2A" not in set(rows["mutation_string"])
    assert rows.loc[rows["mutation_string"] == "HC2D", "sequence_id"].iloc[0] == "sdab:HC2D"
    assert not rows["mutation_string"].str.startswith("HA1").any()
    assert rows["available_pre_embedding"].all()
    assert rows["order"].tolist() == list(range(1, 38))


def test_build_candidate_pool_cli_rejects_unsupported_mode(tmp_path: Path):
    fasta = tmp_path / "tiny.fasta"
    config = tmp_path / "config.yaml"
    fasta.write_text(">A\nACD\n", encoding="utf-8")
    config.write_text("system_id: sdab\n", encoding="utf-8")

    with pytest.raises(ValueError, match="pairwise"):
        main(
            [
                "build-candidate-pool",
                "--config-yaml",
                str(config),
                "--wt-fasta",
                str(fasta),
                "--candidate-mode",
                "pairwise",
                "--out-csv",
                str(tmp_path / "out.csv"),
            ]
        )


def test_build_candidate_pool_refuses_to_overwrite_input(tmp_path: Path):
    fasta = tmp_path / "tiny.fasta"
    config = tmp_path / "config.yaml"
    fasta.write_text(">A\nACD\n", encoding="utf-8")
    config.write_text("system_id: sdab\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Refusing to overwrite input file"):
        main(
            [
                "build-candidate-pool",
                "--config-yaml",
                str(config),
                "--wt-fasta",
                str(fasta),
                "--out-csv",
                str(fasta),
            ]
        )
