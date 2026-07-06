import pandas as pd
import pytest

from active_learning.tools import (
    aggregate_training_replicates,
    clip01,
    compute_switch_utility,
    filter_training_rows,
    normalize_by_reference,
)


CORE_ROW = {
    "sequence_id": "s1",
    "sequence_sha256": "h1",
    "A7": 0.9,
    "A6": 0.1,
    "noise": 0.0,
    "utility_recomputed": 0.81,
    "experiment_status": "measured",
    "label_source_type": "raw_elisa",
    "label_is_direct_measurement": True,
    "usable_for_training": True,
}


def test_compute_switch_utility_rewards_ph_switch():
    assert compute_switch_utility(0.9, 0.1) > compute_switch_utility(0.9, 0.9)
    assert compute_switch_utility(0.2, 0.1) < compute_switch_utility(0.9, 0.1)
    assert compute_switch_utility(0.9, 0.1, noise=0.2, lambda_noise=0.5) == pytest.approx(0.71)


def test_clip_and_reference_normalization():
    assert clip01(-0.1) == 0.0
    assert clip01(1.2) == 1.0
    assert normalize_by_reference(5, blank=1, reference=9) == pytest.approx(0.5)
    assert normalize_by_reference(20, blank=1, reference=9) == 1.0


def test_aggregate_training_replicates_collapses_sequence_rows():
    rows = pd.DataFrame(
        [
            {**CORE_ROW, "A7": 0.8, "A6": 0.2, "noise": 0.05, "replicate_id": "r1"},
            {**CORE_ROW, "A7": 1.0, "A6": 0.4, "noise": 0.15, "replicate_id": "r2"},
            {**CORE_ROW, "sequence_id": "s2", "sequence_sha256": "h2", "A7": 0.4, "A6": 0.2},
        ]
    )

    aggregated = aggregate_training_replicates(rows)
    s1 = aggregated.set_index("sequence_id").loc["s1"]

    assert len(aggregated) == 2
    assert s1["A7_mean"] == pytest.approx(0.9)
    assert s1["A6_mean"] == pytest.approx(0.3)
    assert s1["A7"] == pytest.approx(0.9)
    assert s1["A6"] == pytest.approx(0.3)
    assert s1["n_replicates"] == 2
    assert s1["noise"] == pytest.approx(0.1)
    assert s1["utility_recomputed"] == pytest.approx(0.63)


def test_filter_training_rows_keeps_only_direct_measured_usable_rows():
    rows = pd.DataFrame(
        [
            CORE_ROW,
            {**CORE_ROW, "sequence_id": "legacy", "label_is_direct_measurement": False, "label_source_type": "legacy_training_summary"},
            {**CORE_ROW, "sequence_id": "pseudo", "label_is_direct_measurement": False, "label_source_type": "pseudo_label"},
            {**CORE_ROW, "sequence_id": "legacy_misflagged", "label_source_type": "legacy_training_summary"},
            {**CORE_ROW, "sequence_id": "pseudo_misflagged", "label_source_type": "pseudo_label"},
            {**CORE_ROW, "sequence_id": "unknown_misflagged", "label_source_type": "unknown"},
            {**CORE_ROW, "sequence_id": "unknown", "experiment_status": "unknown"},
            {**CORE_ROW, "sequence_id": "excluded", "experiment_status": "excluded"},
            {**CORE_ROW, "sequence_id": "not_usable", "usable_for_training": False},
            {**CORE_ROW, "sequence_id": "blank_source", "label_source_type": ""},
            {key: value for key, value in {**CORE_ROW, "sequence_id": "missing_source"}.items() if key != "label_source_type"},
        ]
    )

    filtered = filter_training_rows(rows)

    assert filtered["sequence_id"].tolist() == ["s1"]


def test_filter_training_rows_accepts_core_schema_without_source_type():
    row = {key: value for key, value in CORE_ROW.items() if key != "label_source_type"}

    filtered = filter_training_rows(pd.DataFrame([row]))

    assert filtered["sequence_id"].tolist() == ["s1"]


def test_filter_training_rows_requires_core_columns():
    with pytest.raises(ValueError, match="missing core training columns"):
        filter_training_rows(pd.DataFrame([{"sequence_id": "s1"}]))
