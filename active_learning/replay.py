from __future__ import annotations

import pandas as pd

from active_learning.model import fit_rf_regressor, predict_rf_with_uncertainty, score_candidates, select_top_k
from active_learning.tools import aggregate_training_replicates, filter_training_rows, join_embeddings


def leave_round_out_enrichment(
    training: pd.DataFrame,
    embeddings: pd.DataFrame,
    heldout_round: str,
    top_k: int,
    include_pseudo_labels_for_diagnostic_only: bool = False,
) -> dict:
    train_pool = training.loc[training["round_id"] != heldout_round].copy()
    heldout_pool = training.loc[training["round_id"] == heldout_round].copy()
    trainable = filter_training_rows(train_pool)
    heldout = filter_training_rows(heldout_pool)

    if include_pseudo_labels_for_diagnostic_only and "label_source_type" in train_pool.columns:
        pseudo = train_pool.loc[train_pool["label_source_type"].eq("pseudo_label")]
        trainable = pd.concat([trainable, pseudo], ignore_index=True)

    if trainable.empty or heldout.empty:
        return {
            "status": "insufficient_data",
            "diagnostic_only": bool(include_pseudo_labels_for_diagnostic_only),
            "training_row_count": int(len(trainable)),
            "heldout_row_count": int(len(heldout)),
        }

    trainable = aggregate_training_replicates(trainable)
    heldout = aggregate_training_replicates(heldout)
    _, x_train = join_embeddings(trainable, embeddings)
    heldout_joined, x_heldout = join_embeddings(heldout, embeddings)
    model = fit_rf_regressor(
        x_train,
        trainable["utility_recomputed"].to_numpy(dtype=float),
        {"n_estimators": 100, "criterion": "friedman_mse", "random_state": 42},
    )
    predictions = predict_rf_with_uncertainty(model, x_heldout)
    scored = score_candidates(pd.concat([heldout_joined.reset_index(drop=True), predictions], axis=1))
    selected = select_top_k(scored, min(int(top_k), len(scored))).query("selected")
    random_baseline = heldout.sample(n=min(int(top_k), len(heldout)), random_state=42)

    return {
        "status": "ok",
        "top_k_mean_utility": float(selected["utility_recomputed"].mean()),
        "random_baseline_mean_utility": float(random_baseline["utility_recomputed"].mean()),
        "best_recovered_utility": float(selected["utility_recomputed"].max()),
        "heldout_row_count": int(len(heldout)),
        "training_row_count": int(len(trainable)),
        "diagnostic_only": bool(include_pseudo_labels_for_diagnostic_only),
    }
