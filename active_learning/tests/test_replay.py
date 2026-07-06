import pandas as pd

from active_learning.tools import sequence_sha256
from active_learning.replay import leave_round_out_enrichment


def _replay_training_rows():
    rows = []
    for sequence_id, sequence, round_id, utility, direct, source_type in [
        ("r1_high", "AAA", "R1", 0.9, True, "raw_elisa"),
        ("r1_high", "AAA", "R1", 0.7, True, "raw_elisa"),
        ("r1_low", "CCC", "R1", 0.1, True, "raw_elisa"),
        ("pseudo", "AAC", "R1", 1.0, False, "pseudo_label"),
        ("r2_high", "AAG", "R2", 0.85, True, "raw_elisa"),
        ("r2_mid", "ACC", "R2", 0.5, True, "raw_elisa"),
        ("r2_low", "CCA", "R2", 0.2, True, "raw_elisa"),
    ]:
        rows.append(
            {
                "sequence_id": sequence_id,
                "sequence_sha256": sequence_sha256(sequence),
                "round_id": round_id,
                "A7": utility,
                "A6": 0.0,
                "noise": 0.0,
                "utility_recomputed": utility,
                "experiment_status": "measured" if direct else "excluded",
                "label_source_type": source_type,
                "label_is_direct_measurement": direct,
                "usable_for_training": direct,
            }
        )
    return pd.DataFrame(rows)


def _replay_embeddings(training: pd.DataFrame):
    coords = {
        "r1_high": (1.0, 1.0),
        "r1_low": (0.0, 0.0),
        "pseudo": (1.0, 0.2),
        "r2_high": (0.9, 0.9),
        "r2_mid": (0.4, 0.4),
        "r2_low": (0.1, 0.0),
    }
    return pd.DataFrame(
        [
            {
                "sequence_id": row.sequence_id,
                "sequence_sha256": row.sequence_sha256,
                "emb_0": coords[row.sequence_id][0],
                "emb_1": coords[row.sequence_id][1],
            }
            for row in training.itertuples()
        ]
    ).drop_duplicates("sequence_id")


def test_leave_round_out_enrichment_recovers_high_utility_rows():
    training = _replay_training_rows()
    result = leave_round_out_enrichment(training, _replay_embeddings(training), heldout_round="R2", top_k=1)

    assert result["status"] == "ok"
    assert result["heldout_row_count"] == 3
    assert result["top_k_mean_utility"] >= result["random_baseline_mean_utility"]
    assert result["best_recovered_utility"] >= 0.5


def test_leave_round_out_enrichment_returns_insufficient_data_status():
    training = _replay_training_rows()
    result = leave_round_out_enrichment(training, _replay_embeddings(training), heldout_round="R3", top_k=1)

    assert result["status"] == "insufficient_data"


def test_leave_round_out_enrichment_excludes_pseudo_by_default():
    training = _replay_training_rows()
    result = leave_round_out_enrichment(training, _replay_embeddings(training), heldout_round="R2", top_k=1)

    assert not result["diagnostic_only"]
    assert result["training_row_count"] == 2


def test_leave_round_out_enrichment_can_include_pseudo_for_diagnostics_only():
    training = _replay_training_rows()
    result = leave_round_out_enrichment(
        training,
        _replay_embeddings(training),
        heldout_round="R2",
        top_k=1,
        include_pseudo_labels_for_diagnostic_only=True,
    )

    assert result["diagnostic_only"]
    assert result["training_row_count"] == 3


def test_leave_round_out_enrichment_uses_deterministic_top_k_random_baseline():
    training = _replay_training_rows()
    result = leave_round_out_enrichment(training, _replay_embeddings(training), heldout_round="R2", top_k=1)

    assert result["random_baseline_mean_utility"] in set(training.loc[training["round_id"].eq("R2"), "utility_recomputed"])
