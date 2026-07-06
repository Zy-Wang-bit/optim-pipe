from pathlib import Path

import pytest

from active_learning.tools import (
    build_round_provenance,
    load_round_config,
    sha256_file,
    sha256_text,
)


def test_load_round_config_fills_defaults(tmp_path: Path):
    config_path = tmp_path / "round.yaml"
    config_path.write_text(
        "system_id: sdab\nround_id: R5\nwt_fasta_path: wt.fasta\n",
        encoding="utf-8",
    )

    config = load_round_config(config_path)

    assert config["objective"]["formula"] == "A7 * (1 - A6)"
    assert config["model"]["n_estimators"] == 100
    assert config["model"]["criterion"] == "friedman_mse"
    assert config["model"]["random_state"] == 42
    assert config["acquisition"]["policy"] == "topk"
    assert config["acquisition"]["top_k"] == 24


def test_sha256_file_is_stable(tmp_path: Path):
    path = tmp_path / "tiny.csv"
    path.write_text("sequence_id,A7\ns1,0.5\n", encoding="utf-8")

    digest = sha256_file(path)

    assert digest == sha256_file(path)
    assert digest == sha256_text(path.read_text(encoding="utf-8"))
    assert len(digest) == 64


def test_build_round_provenance_records_hashes(tmp_path: Path):
    config_path = tmp_path / "round.yaml"
    wt_path = tmp_path / "wt.fasta"
    training_path = tmp_path / "training.csv"
    candidate_path = tmp_path / "candidates.csv"
    embedding_path = tmp_path / "embeddings.csv"
    config_path.write_text(
        """
system_id: sdab
round_id: R5
wt_fasta_path: wt.fasta
embedding:
  embedding_model_name: ESMC-600M
  embedding_model_version: test
  embedding_layer: final
  pooling: mean
  embedding_dim: 2
""".strip(),
        encoding="utf-8",
    )
    wt_path.write_text(">A\nACD\n", encoding="utf-8")
    training_path.write_text("sequence_id,A7,A6\ns1,1,0\n", encoding="utf-8")
    candidate_path.write_text("sequence_id\ns2\n", encoding="utf-8")
    embedding_path.write_text("sequence_id,sequence_sha256,emb_0,emb_1\ns1,h,0,1\n", encoding="utf-8")

    config = load_round_config(config_path)
    provenance = build_round_provenance(
        config,
        {
            "wt_fasta": wt_path,
            "training_csv": training_path,
            "candidate_pool": candidate_path,
            "embedding_csv": embedding_path,
        },
    )

    assert provenance["system_id"] == "sdab"
    assert provenance["round_id"] == "R5"
    assert provenance["training_csv_sha256"] == sha256_file(training_path)
    assert provenance["candidate_pool_sha256"] == sha256_file(candidate_path)
    assert provenance["embedding_csv_sha256"] == sha256_file(embedding_path)
    assert len(provenance["objective_config_hash"]) == 64
    assert len(provenance["model_config_hash"]) == 64
    assert provenance["embedding_model_name"] == "ESMC-600M"


def test_build_round_provenance_hashes_configured_fasta_record(tmp_path: Path):
    wt_path = tmp_path / "wt.fasta"
    training_path = tmp_path / "training.csv"
    candidate_path = tmp_path / "candidates.csv"
    embedding_path = tmp_path / "embeddings.csv"
    wt_path.write_text(">A\nACD\n>B\nGGG\n", encoding="utf-8")
    training_path.write_text("sequence_id,A7,A6\ns1,1,0\n", encoding="utf-8")
    candidate_path.write_text("sequence_id\ns2\n", encoding="utf-8")
    embedding_path.write_text("sequence_id,sequence_sha256,emb_0\ns1,h,0\n", encoding="utf-8")
    config_path = tmp_path / "round.yaml"
    config_path.write_text("system_id: sdab\nround_id: R5\n", encoding="utf-8")
    config = load_round_config(config_path)
    config["wt_fasta_path"] = str(wt_path)
    config["naming"]["fasta_record_id"] = "A"

    provenance = build_round_provenance(
        config,
        {
            "wt_fasta": wt_path,
            "training_csv": training_path,
            "candidate_pool": candidate_path,
            "embedding_csv": embedding_path,
        },
    )

    assert provenance["wt_sequence_sha256"] == sha256_text("ACD")


def test_build_round_provenance_requires_existing_paths(tmp_path: Path):
    config_path = tmp_path / "round.yaml"
    config_path.write_text("system_id: sdab\nround_id: R5\n", encoding="utf-8")
    config = load_round_config(config_path)

    with pytest.raises(FileNotFoundError):
        build_round_provenance(config, {"training_csv": tmp_path / "missing.csv"})
