from pathlib import Path

import json

import pandas as pd
import pytest

from active_learning.cli import main
from active_learning.tools import sequence_sha256, sha256_file


def _write_round_inputs(tmp_path: Path):
    fasta = tmp_path / "wt.fasta"
    config = tmp_path / "config.yaml"
    training = tmp_path / "training.csv"
    candidates = tmp_path / "candidates.csv"
    embeddings = tmp_path / "embeddings.csv"
    failed = tmp_path / "failed.csv"
    fasta.write_text(">A\nACD\n", encoding="utf-8")
    config.write_text(
        f"""
system_id: sdab
round_id: R5
wt_fasta_path: {fasta}
embedding:
  embedding_model_name: ESMC-600M
  embedding_model_version: test
  embedding_layer: final
  pooling: mean
  embedding_dim: 2
model:
  n_estimators: 9
  random_state: 1
acquisition:
  top_k: 4
""".strip(),
        encoding="utf-8",
    )
    training_rows = [
        {
            "sequence_id": "t1",
            "sequence": "ACD",
            "sequence_sha256": sequence_sha256("ACD"),
            "A7": 0.9,
            "A6": 0.1,
            "noise": 0.0,
            "utility_input": 999.0,
            "utility_recomputed": 999.0,
            "experiment_status": "measured",
            "label_source_type": "raw_elisa",
            "label_is_direct_measurement": True,
            "usable_for_training": True,
        },
        {
            "sequence_id": "t2",
            "sequence": "GCD",
            "sequence_sha256": sequence_sha256("GCD"),
            "A7": 0.6,
            "A6": 0.5,
            "noise": 0.0,
            "utility_input": 999.0,
            "utility_recomputed": 999.0,
            "experiment_status": "measured",
            "label_source_type": "raw_elisa",
            "label_is_direct_measurement": True,
            "usable_for_training": True,
        },
        {
            "sequence_id": "t3",
            "sequence": "AAD",
            "sequence_sha256": sequence_sha256("AAD"),
            "A7": 0.3,
            "A6": 0.2,
            "noise": 0.0,
            "utility_input": 999.0,
            "utility_recomputed": 999.0,
            "experiment_status": "measured",
            "label_source_type": "raw_elisa",
            "label_is_direct_measurement": True,
            "usable_for_training": True,
        },
        {
            "sequence_id": "pseudo",
            "sequence": "AAA",
            "sequence_sha256": sequence_sha256("AAA"),
            "A7": 1.0,
            "A6": 0.0,
            "noise": 0.0,
            "utility_input": 999.0,
            "utility_recomputed": 999.0,
            "experiment_status": "excluded",
            "label_source_type": "pseudo_label",
            "label_is_direct_measurement": False,
            "usable_for_training": False,
        },
    ]
    candidate_rows = [
        ("c1", "ACC", True, ""),
        ("c2", "ACA", True, ""),
        ("c3", "GCD", True, ""),
        ("c4", "CCD", False, "previously_selected"),
        ("c5", "DDD", True, ""),
    ]
    pd.DataFrame(training_rows).to_csv(training, index=False)
    pd.DataFrame(
        [
            {
                "system_id": "sdab",
                "round_id": "R5",
                "sequence_id": seq_id,
                "sequence": seq,
                "sequence_sha256": sequence_sha256(seq),
                "mutation_string": seq_id,
                "candidate_origin": "test",
                "available_pre_embedding": available,
                "unavailable_reason_pre_embedding": reason,
            }
            for seq_id, seq, available, reason in candidate_rows
        ]
    ).to_csv(candidates, index=False)
    embedding_rows = [
        ("t1", "ACD", 0.0, 0.0),
        ("t2", "GCD", 0.4, 0.4),
        ("t3", "AAD", 0.8, 0.2),
        ("c1", "ACC", 0.1, 0.1),
        ("c2", "ACA", 0.2, 0.3),
        ("c3", "GCD", 0.4, 0.4),
        ("c4", "CCD", 0.6, 0.6),
    ]
    pd.DataFrame(
        [
            {"sequence_id": seq_id, "sequence_sha256": sequence_sha256(seq), "emb_0": emb0, "emb_1": emb1}
            for seq_id, seq, emb0, emb1 in embedding_rows
        ]
    ).to_csv(embeddings, index=False)
    pd.DataFrame([{"sequence_id": "c2"}]).to_csv(failed, index=False)
    return config, training, candidates, embeddings, failed


def test_run_round_writes_reproducible_artifacts(tmp_path: Path):
    config, training, candidates, embeddings, failed = _write_round_inputs(tmp_path)
    out_dir = tmp_path / "out"

    main(
        [
            "run-round",
            "--config-yaml",
            str(config),
            "--training-csv",
            str(training),
            "--candidate-csv",
            str(candidates),
            "--embedding-csv",
            str(embeddings),
            "--failed-synthesis-csv",
            str(failed),
            "--round-id",
            "R5",
            "--top-k",
            "4",
            "--out-dir",
            str(out_dir),
        ]
    )

    expected = [
        "training_data_round_R5.csv",
        "candidate_pool_round_R5.csv",
        "model_config_round_R5.json",
        "round_provenance_round_R5.json",
        "prediction_table_round_R5.csv",
        "selected_candidates_round_R5.csv",
    ]
    assert all((out_dir / name).exists() for name in expected)

    training_out = pd.read_csv(out_dir / "training_data_round_R5.csv")
    candidates_out = pd.read_csv(out_dir / "candidate_pool_round_R5.csv").set_index("sequence_id")
    selected = pd.read_csv(out_dir / "selected_candidates_round_R5.csv")
    model_config = json.loads((out_dir / "model_config_round_R5.json").read_text(encoding="utf-8"))
    provenance = json.loads((out_dir / "round_provenance_round_R5.json").read_text(encoding="utf-8"))

    assert training_out["sequence_id"].tolist() == ["t1", "t2", "t3"]
    assert training_out.loc[0, "utility_recomputed"] == pytest.approx(0.81)
    assert selected["sequence_id"].tolist() == ["c1"]
    assert selected["sequence_id"].isin(["c2", "c3", "c4", "c5"]).sum() == 0
    assert len(selected) == 1
    assert candidates_out.loc["c2", "unavailable_reason_pre_embedding"] == "failed_synthesis"
    assert candidates_out.loc["c3", "unavailable_reason_pre_embedding"] == "tested"
    assert not candidates_out.loc["c4", "available"]
    assert candidates_out.loc["c5", "missing_embedding"]
    assert set(model_config) == {
        "model.type",
        "model.n_estimators",
        "model.criterion",
        "model.min_samples_leaf",
        "model.bootstrap",
        "model.random_state",
        "sklearn_version",
    }
    assert provenance["training_csv_sha256"]
    assert provenance["training_csv_sha256"] == sha256_file(training)
    assert provenance["objective_config_hash"]
    assert provenance["model_config_hash"]
    assert provenance["wt_sequence_sha256"] == sequence_sha256("ACD")
    assert provenance["embedding_model_name"] == "ESMC-600M"
    assert provenance["random_state"] == 1
    assert not any(column.startswith("emb_") for column in pd.read_csv(out_dir / "prediction_table_round_R5.csv").columns)


def test_run_round_rejects_zero_usable_direct_labels(tmp_path: Path):
    config, training, candidates, embeddings, _ = _write_round_inputs(tmp_path)
    rows = pd.read_csv(training)
    rows["usable_for_training"] = False
    rows.to_csv(training, index=False)

    with pytest.raises(ValueError, match="No usable direct wet-lab labels"):
        main(
            [
                "run-round",
                "--config-yaml",
                str(config),
                "--training-csv",
                str(training),
                "--candidate-csv",
                str(candidates),
                "--embedding-csv",
                str(embeddings),
                "--out-dir",
                str(tmp_path / "out"),
            ]
        )


def test_run_round_writes_prediction_schema_when_no_candidates_available(tmp_path: Path):
    config, training, candidates, embeddings, _ = _write_round_inputs(tmp_path)
    rows = pd.read_csv(candidates)
    rows["available_pre_embedding"] = False
    rows["unavailable_reason_pre_embedding"] = "invalid_sequence"
    rows.to_csv(candidates, index=False)
    out_dir = tmp_path / "out"

    main(
        [
            "run-round",
            "--config-yaml",
            str(config),
            "--training-csv",
            str(training),
            "--candidate-csv",
            str(candidates),
            "--embedding-csv",
            str(embeddings),
            "--out-dir",
            str(out_dir),
        ]
    )

    prediction_header = (out_dir / "prediction_table_round_R5.csv").read_text(encoding="utf-8").splitlines()[0]
    selected_header = (out_dir / "selected_candidates_round_R5.csv").read_text(encoding="utf-8").splitlines()[0]

    assert "predicted_utility" in prediction_header
    assert "selection_rank" in selected_header


def test_run_round_refuses_to_overwrite_input_training_csv(tmp_path: Path):
    config, training, candidates, embeddings, _ = _write_round_inputs(tmp_path)
    round_training = tmp_path / "training_data_round_R5.csv"
    pd.read_csv(training).to_csv(round_training, index=False)

    with pytest.raises(ValueError, match="Refusing to overwrite input training CSV"):
        main(
            [
                "run-round",
                "--config-yaml",
                str(config),
                "--training-csv",
                str(round_training),
                "--candidate-csv",
                str(candidates),
                "--embedding-csv",
                str(embeddings),
                "--out-dir",
                str(training.parent),
            ]
        )


def test_run_round_refuses_to_overwrite_optional_input_csv(tmp_path: Path):
    config, training, candidates, embeddings, _ = _write_round_inputs(tmp_path)
    selected_input = tmp_path / "selected_candidates_round_R5.csv"
    pd.DataFrame([{"sequence_id": "old"}]).to_csv(selected_input, index=False)

    with pytest.raises(ValueError, match="Refusing to overwrite input selected CSV"):
        main(
            [
                "run-round",
                "--config-yaml",
                str(config),
                "--training-csv",
                str(training),
                "--candidate-csv",
                str(candidates),
                "--embedding-csv",
                str(embeddings),
                "--previously-selected-csv",
                str(selected_input),
                "--out-dir",
                str(tmp_path),
            ]
        )


def test_run_round_refuses_to_overwrite_config_input(tmp_path: Path):
    config, training, candidates, embeddings, _ = _write_round_inputs(tmp_path)
    config_output_name = tmp_path / "model_config_round_R5.json"
    config_output_name.write_text(
        json.dumps(
            {
                "system_id": "sdab",
                "round_id": "R5",
                "wt_fasta_path": str(tmp_path / "wt.fasta"),
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Refusing to overwrite input model_config"):
        main(
            [
                "run-round",
                "--config-yaml",
                str(config_output_name),
                "--training-csv",
                str(training),
                "--candidate-csv",
                str(candidates),
                "--embedding-csv",
                str(embeddings),
                "--out-dir",
                str(tmp_path),
            ]
        )
