import pandas as pd
import pytest

from active_learning.tools import generate_single_substitution_pool, mark_available_candidates


def test_generate_single_substitution_pool_for_one_position():
    pool = generate_single_substitution_pool("ACD", allowed_positions=[2], mutation_chain_label="H")

    assert len(pool) == 19
    assert "HC2C" not in set(pool["mutation_string"])
    assert pool["sequence_sha256"].str.len().eq(64).all()
    assert pool["available_pre_embedding"].all()
    assert pool["order"].tolist() == list(range(1, 20))


def test_generate_single_substitution_pool_rejects_invalid_position():
    with pytest.raises(ValueError, match="4"):
        generate_single_substitution_pool("ACD", allowed_positions=[4])


def test_generate_single_substitution_pool_allows_empty_position_list():
    pool = generate_single_substitution_pool("ACD", allowed_positions=[])

    assert pool.empty


def test_mark_available_candidates_records_pre_embedding_reasons():
    candidates = pd.DataFrame(
        [
            {"sequence_id": "tested", "available_pre_embedding": True},
            {"sequence_id": "failed", "available_pre_embedding": True},
            {"sequence_id": "invalid", "available_pre_embedding": True},
            {"sequence_id": "selected", "available_pre_embedding": True},
            {"sequence_id": "open", "available_pre_embedding": True},
        ]
    )
    training = pd.DataFrame([{"sequence_id": "tested", "experiment_status": "measured"}])
    failed = pd.DataFrame([{"sequence_id": "failed"}])
    invalid = pd.DataFrame([{"sequence_id": "invalid"}])
    selected = pd.DataFrame([{"sequence_id": "selected"}])

    marked = mark_available_candidates(candidates, training, failed, invalid, selected).set_index("sequence_id")

    assert marked.loc["tested", "unavailable_reason_pre_embedding"] == "tested"
    assert marked.loc["failed", "unavailable_reason_pre_embedding"] == "failed_synthesis"
    assert marked.loc["invalid", "unavailable_reason_pre_embedding"] == "invalid_sequence"
    assert marked.loc["selected", "unavailable_reason_pre_embedding"] == "previously_selected"
    assert marked.loc["open", "unavailable_reason_pre_embedding"] == ""
    assert marked.loc["open", "available_pre_embedding"]
    assert not marked.loc["tested", "available_pre_embedding"]


def test_mark_available_candidates_treats_legacy_summary_as_tested():
    candidates = pd.DataFrame([{"sequence_id": "sdab:HC2A", "mutation_string": "HC2A"}])
    training = pd.DataFrame(
        [
            {
                "sequence_id": "sdab:HC2A",
                "mutation_string": "HC2A",
                "experiment_status": "excluded",
                "label_source_type": "legacy_training_summary",
            }
        ]
    )

    marked = mark_available_candidates(candidates, training)

    assert marked.loc[0, "unavailable_reason_pre_embedding"] == "tested"


def test_mark_available_candidates_preserves_string_false_availability():
    candidates = pd.DataFrame(
        [{"sequence_id": "bad", "available_pre_embedding": "False", "unavailable_reason_pre_embedding": "invalid_sequence"}]
    )

    marked = mark_available_candidates(candidates, pd.DataFrame())

    assert not marked.loc[0, "available_pre_embedding"]
    assert marked.loc[0, "unavailable_reason_pre_embedding"] == "invalid_sequence"
