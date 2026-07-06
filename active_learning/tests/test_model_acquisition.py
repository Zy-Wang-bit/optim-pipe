import numpy as np
import pandas as pd
import pytest

from active_learning.model import fit_rf_regressor, predict_rf_with_uncertainty, score_candidates, select_top_k


def test_fit_rf_regressor_reads_config_and_predicts_uncertainty():
    x_train = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)
    y = np.array([0.0, 0.1, 0.8, 1.0], dtype=float)

    model = fit_rf_regressor(
        x_train,
        y,
        {
            "n_estimators": 7,
            "criterion": "friedman_mse",
            "random_state": 3,
            "min_samples_leaf": 1,
            "bootstrap": True,
        },
    )
    predictions = predict_rf_with_uncertainty(model, np.array([[1, 1], [0, 0]], dtype=float))

    assert model.n_estimators == 7
    assert model.criterion == "friedman_mse"
    assert list(predictions.columns) == ["predicted_utility", "prediction_uncertainty"]
    assert len(predictions) == 2


def test_score_candidates_topk_and_ucb():
    predictions = pd.DataFrame(
        {
            "sequence_id": ["a", "b"],
            "predicted_utility": [0.4, 0.5],
            "prediction_uncertainty": [0.3, 0.1],
        }
    )

    topk = score_candidates(predictions, policy="topk")
    ucb = score_candidates(predictions, policy="ucb", kappa=2.0)

    assert topk["acquisition_score"].tolist() == [0.4, 0.5]
    assert ucb["acquisition_score"].tolist() == [1.0, 0.7]


def test_select_top_k_marks_rank_and_selected_rows():
    scored = pd.DataFrame(
        {
            "sequence_id": ["low", "high", "mid"],
            "acquisition_score": [0.1, 0.9, 0.5],
        }
    )

    selected = select_top_k(scored, k=2)

    assert selected["sequence_id"].tolist() == ["high", "mid", "low"]
    assert selected["selection_rank"].tolist() == [1, 2, 3]
    assert selected["selected"].tolist() == [True, True, False]


def test_score_candidates_rejects_unknown_policy():
    with pytest.raises(ValueError, match="greedy"):
        score_candidates(pd.DataFrame({"predicted_utility": [1.0]}), policy="greedy")
