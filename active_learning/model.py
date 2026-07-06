from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestRegressor

from active_learning.tools import (
    aggregate_training_replicates,
    compute_switch_utility,
    filter_training_rows,
    join_embeddings,
    mark_available_candidates,
)


MODEL_CONFIG_KEYS = [
    "type",
    "n_estimators",
    "criterion",
    "min_samples_leaf",
    "bootstrap",
    "random_state",
]
PREDICTION_OUTPUT_COLUMNS = [
    "system_id",
    "sequence_id",
    "sequence",
    "sequence_sha256",
    "mutation_string",
    "candidate_origin",
    "predicted_utility",
    "prediction_uncertainty",
    "acquisition_score",
    "selection_rank",
    "selected",
]


def fit_rf_regressor(X: np.ndarray, y: np.ndarray, config: dict) -> RandomForestRegressor:
    model = RandomForestRegressor(
        n_estimators=int(config.get("n_estimators", 100)),
        criterion=config.get("criterion", "friedman_mse"),
        random_state=int(config.get("random_state", 42)),
        min_samples_leaf=int(config.get("min_samples_leaf", 1)),
        bootstrap=bool(config.get("bootstrap", True)),
    )
    return model.fit(X, y)


def predict_rf_with_uncertainty(model: RandomForestRegressor, X: np.ndarray) -> pd.DataFrame:
    tree_predictions = np.vstack([tree.predict(X) for tree in model.estimators_])
    return pd.DataFrame(
        {
            "predicted_utility": tree_predictions.mean(axis=0),
            "prediction_uncertainty": tree_predictions.std(axis=0),
        }
    )


def score_candidates(predictions: pd.DataFrame, policy: str = "topk", kappa: float = 0.0) -> pd.DataFrame:
    scored = predictions.copy()
    if policy == "topk":
        scored["acquisition_score"] = scored["predicted_utility"]
    elif policy == "ucb":
        uncertainty = scored.get("prediction_uncertainty", 0.0)
        scored["acquisition_score"] = scored["predicted_utility"] + float(kappa) * uncertainty
    else:
        raise ValueError(f"unsupported acquisition policy: {policy}")
    return scored


def select_top_k(scored: pd.DataFrame, k: int) -> pd.DataFrame:
    ranked = scored.sort_values("acquisition_score", ascending=False, kind="mergesort").reset_index(drop=True)
    ranked["selection_rank"] = range(1, len(ranked) + 1)
    ranked["selected"] = ranked["selection_rank"] <= int(k)
    return ranked


def recompute_utility(rows: pd.DataFrame, noise_lambda: float) -> pd.DataFrame:
    rows = rows.copy()
    if "noise" not in rows.columns:
        rows["noise"] = 0.0
    rows["utility_recomputed"] = [
        compute_switch_utility(a7, a6, noise, noise_lambda)
        for a7, a6, noise in zip(rows["A7"], rows["A6"], rows["noise"])
    ]
    return rows


def mark_final_candidate_availability(candidates: pd.DataFrame, embeddings: pd.DataFrame) -> pd.DataFrame:
    candidates = candidates.copy()
    embedding_ids = set(embeddings["sequence_id"].astype(str))
    sequence_ids = candidates["sequence_id"].astype(str)

    candidates["missing_embedding"] = ~sequence_ids.isin(embedding_ids)
    candidates["available"] = candidates["available_pre_embedding"] & ~candidates["missing_embedding"]
    candidates["unavailable_reason"] = candidates["unavailable_reason_pre_embedding"].fillna("").astype(str)
    missing_reason = candidates["unavailable_reason"].eq("") & candidates["missing_embedding"]
    candidates.loc[missing_reason, "unavailable_reason"] = "missing_embedding"
    return candidates


def model_config(config: dict) -> dict:
    model = config["model"]
    return {f"model.{key}": model.get(key) for key in MODEL_CONFIG_KEYS} | {"sklearn_version": sklearn.__version__}


def round_output_paths(out_dir: Path, round_id: str) -> dict[str, Path]:
    return {
        "training": out_dir / f"training_data_round_{round_id}.csv",
        "candidate": out_dir / f"candidate_pool_round_{round_id}.csv",
        "prediction": out_dir / f"prediction_table_round_{round_id}.csv",
        "selected": out_dir / f"selected_candidates_round_{round_id}.csv",
        "model_config": out_dir / f"model_config_round_{round_id}.json",
        "provenance": out_dir / f"round_provenance_round_{round_id}.json",
    }


def run_active_learning_round(
    config: dict,
    training: pd.DataFrame,
    candidates: pd.DataFrame,
    embeddings: pd.DataFrame,
    failed_synthesis: pd.DataFrame | None = None,
    invalid_candidates: pd.DataFrame | None = None,
    previously_selected: pd.DataFrame | None = None,
) -> dict[str, pd.DataFrame | dict]:
    noise_lambda = float(config["objective"].get("noise_lambda", 0.0))

    training = recompute_utility(training, noise_lambda)
    trainable = filter_training_rows(training)
    if trainable.empty:
        raise ValueError("No usable direct wet-lab labels")
    trainable = aggregate_training_replicates(trainable)
    trainable = recompute_utility(trainable, noise_lambda)
    _, x_train = join_embeddings(trainable, embeddings)

    candidates = mark_available_candidates(
        candidates,
        training,
        failed_synthesis=failed_synthesis,
        invalid_candidates=invalid_candidates,
        previously_selected=previously_selected,
    )
    candidates = mark_final_candidate_availability(candidates, embeddings)
    available_candidates = candidates.loc[candidates["available"]].copy()
    if available_candidates.empty:
        prediction_table = pd.DataFrame(columns=PREDICTION_OUTPUT_COLUMNS)
        selected = pd.DataFrame(columns=PREDICTION_OUTPUT_COLUMNS)
    else:
        joined_candidates, x_candidates = join_embeddings(available_candidates, embeddings)
        model = fit_rf_regressor(x_train, trainable["utility_recomputed"].to_numpy(dtype=float), config["model"])
        predictions = predict_rf_with_uncertainty(model, x_candidates)
        prediction_input = pd.concat([joined_candidates.reset_index(drop=True), predictions], axis=1)
        scored = score_candidates(
            prediction_input,
            policy=config["acquisition"].get("policy", "topk"),
            kappa=float(config["acquisition"].get("kappa", 0.0)),
        )
        prediction_table = select_top_k(scored, int(config["acquisition"].get("top_k", 24)))
        prediction_table = prediction_table[[column for column in PREDICTION_OUTPUT_COLUMNS if column in prediction_table.columns]]
        selected = prediction_table.loc[prediction_table["selected"]].copy()

    return {
        "training": trainable,
        "candidate": candidates,
        "prediction": prediction_table,
        "selected": selected,
        "model_config": model_config(config),
    }
