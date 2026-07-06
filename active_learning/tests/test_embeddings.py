from pathlib import Path

import pandas as pd
import pytest

from active_learning.tools import join_embeddings, load_embedding_table


def test_load_embedding_table_reads_csv(tmp_path: Path):
    path = tmp_path / "embeddings.csv"
    path.write_text("sequence_id,sequence_sha256,emb_1,emb_0\ns1,h1,2,1\n", encoding="utf-8")

    embeddings = load_embedding_table(path)

    assert embeddings.loc[0, "sequence_id"] == "s1"
    assert {"emb_0", "emb_1"}.issubset(embeddings.columns)


def test_join_embeddings_preserves_input_order_and_sorts_embedding_columns():
    rows = pd.DataFrame(
        [
            {"sequence_id": "s2", "sequence_sha256": "h2"},
            {"sequence_id": "s1", "sequence_sha256": "h1"},
        ]
    )
    embeddings = pd.DataFrame(
        [
            {"sequence_id": "s1", "sequence_sha256": "h1", "emb_1": 2.0, "emb_0": 1.0},
            {"sequence_id": "s2", "sequence_sha256": "h2", "emb_1": 4.0, "emb_0": 3.0},
        ]
    )

    joined, x = join_embeddings(rows, embeddings)

    assert joined["sequence_id"].tolist() == ["s2", "s1"]
    assert x.shape == (2, 2)
    assert x.tolist() == [[3.0, 4.0], [1.0, 2.0]]


def test_join_embeddings_rejects_hash_mismatch():
    rows = pd.DataFrame([{"sequence_id": "s1", "sequence_sha256": "expected"}])
    embeddings = pd.DataFrame([{"sequence_id": "s1", "sequence_sha256": "actual", "emb_0": 1.0}])

    with pytest.raises(ValueError, match="sequence_sha256 mismatch"):
        join_embeddings(rows, embeddings)


def test_join_embeddings_lists_missing_ids():
    rows = pd.DataFrame([{"sequence_id": "s1", "sequence_sha256": "h1"}])
    embeddings = pd.DataFrame(columns=["sequence_id", "sequence_sha256", "emb_0"])

    with pytest.raises(ValueError, match="Missing embeddings for sequence_id: s1"):
        join_embeddings(rows, embeddings)
