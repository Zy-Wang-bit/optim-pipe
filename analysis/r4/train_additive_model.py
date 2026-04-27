#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Train dual ElasticNet models (pH7.4, pH6.0) for sdab R4 His-mutation combos."""

import json
import os
import sys
import warnings
from itertools import combinations
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, r2_score

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_CSV = REPO_ROOT / "experiments/sdab/R4/data/training_data.csv"
OUT_DIR = REPO_ROOT / "experiments/sdab/R4/model"

CDR_POSITIONS = [
    ("G", 26), ("R", 27), ("T", 28), ("F", 29), ("S", 30), ("S", 31), ("Y", 32), ("A", 33),
    ("I", 51), ("N", 52), ("G", 53), ("N", 54), ("G", 55), ("G", 56), ("S", 57), ("T", 58),
    ("N", 97), ("A", 98), ("D", 99), ("Q", 100), ("R", 101), ("G", 102), ("W", 103),
    ("A", 104), ("V", 105), ("A", 106), ("V", 107), ("E", 108), ("G", 109), ("D", 110), ("Y", 111),
]
SINGLE_NAMES = [f"H{wt}{pos}H" for wt, pos in CDR_POSITIONS]
POS_TO_IDX = {name: i for i, name in enumerate(SINGLE_NAMES)}

ALPHAS = [0.001, 0.01, 0.1, 1.0]
L1_RATIOS = [0.3, 0.5, 0.7, 0.9]
RANDOM_STATE = 42
MAX_ITER = 10000
TOL = 1e-5


def mutations_to_indices(mutations_str):
    if not mutations_str or pd.isna(mutations_str):
        return []
    out = []
    for m in mutations_str.split(";"):
        if m not in POS_TO_IDX:
            raise ValueError(f"unknown mutation token {m!r}: not in CDR_POSITIONS whitelist")
        out.append(POS_TO_IDX[m])
    return sorted(out)


def build_pair_whitelist(index_lists, min_count=2):
    counts = {}
    for idxs in index_lists:
        for i, j in combinations(idxs, 2):
            counts[(i, j)] = counts.get((i, j), 0) + 1
    pairs = sorted(p for p, c in counts.items() if c >= min_count)
    return pairs


def build_feature_matrix(index_lists, pairs):
    n = len(index_lists)
    n_single = len(SINGLE_NAMES)
    pair_idx = {p: k for k, p in enumerate(pairs)}
    X = np.zeros((n, n_single + len(pairs)), dtype=np.float64)
    for r, idxs in enumerate(index_lists):
        for i in idxs:
            X[r, i] = 1.0
        for i, j in combinations(idxs, 2):
            if (i, j) in pair_idx:
                X[r, n_single + pair_idx[(i, j)]] = 1.0
    return X


def fit_elasticnet(X, y, alpha, l1_ratio):
    model = ElasticNet(
        alpha=alpha, l1_ratio=l1_ratio,
        max_iter=MAX_ITER, tol=TOL, random_state=RANDOM_STATE,
    )
    model.fit(X, y)
    return model


def loo_combo_predictions(X, y, combo_mask, alpha, l1_ratio):
    combo_indices = np.where(combo_mask)[0]
    preds = np.zeros(len(combo_indices))
    truths = np.zeros(len(combo_indices))
    for k, i in enumerate(combo_indices):
        train_mask = np.ones(len(y), dtype=bool)
        train_mask[i] = False
        model = fit_elasticnet(X[train_mask], y[train_mask], alpha, l1_ratio)
        preds[k] = model.predict(X[i:i + 1])[0]
        truths[k] = y[i]
    return truths, preds


def grid_search(X, y, combo_mask):
    best = None
    for alpha in ALPHAS:
        for l1 in L1_RATIOS:
            t, p = loo_combo_predictions(X, y, combo_mask, alpha, l1)
            rho, _ = spearmanr(t, p)
            if np.isnan(rho):
                rho = -1.0
            try:
                r2 = r2_score(t, p)
            except Exception:
                r2 = -np.inf
            score = (rho, r2)
            if best is None or score > best[0]:
                best = (score, alpha, l1, t, p)
    return best


def write_coefficients(path, single_names, pair_records, coef74, coef6):
    n_single = len(single_names)
    rows = []
    for i, name in enumerate(single_names):
        b74, b6 = coef74[i], coef6[i]
        rows.append({
            "feature": name, "type": "single",
            "beta_pH74": b74, "beta_pH6": b6, "delta_beta": b74 - b6,
        })
    for k, (i, j, name) in enumerate(pair_records):
        b74, b6 = coef74[n_single + k], coef6[n_single + k]
        rows.append({
            "feature": name, "type": "pair",
            "beta_pH74": b74, "beta_pH6": b6, "delta_beta": b74 - b6,
        })
    df = pd.DataFrame(rows)
    df["abs_delta"] = df["delta_beta"].abs()
    singles = df[df["type"] == "single"].sort_values("abs_delta", ascending=False)
    pairs = df[df["type"] == "pair"].sort_values("abs_delta", ascending=False)
    out = pd.concat([singles, pairs], ignore_index=True).drop(columns=["abs_delta"])
    out.to_csv(path, index=False)
    return out


def metrics(truth, pred):
    rho, _ = spearmanr(truth, pred)
    r2 = r2_score(truth, pred)
    try:
        pr, _ = pearsonr(truth, pred)
    except Exception:
        pr = np.nan
    mae = mean_absolute_error(truth, pred)
    return rho, r2, pr, mae


def write_report(path, best74, best6, m74, m6, coef_df, combos_df):
    top_delta = coef_df[coef_df["type"] == "single"].reindex(
        coef_df[coef_df["type"] == "single"]["delta_beta"].abs().sort_values(ascending=False).index
    ).head(10)
    pair_df = coef_df[coef_df["type"] == "pair"].copy()
    neg_pairs = pair_df[pair_df["beta_pH74"] < 0].copy()
    neg_pairs["abs74"] = neg_pairs["beta_pH74"].abs()
    top_neg_pairs = neg_pairs.sort_values("abs74", ascending=False).head(10).drop(columns=["abs74"])

    if "collapsed" in combos_df.columns:
        col_mask = combos_df["collapsed"].astype(bool)
        col_bias = (
            (combos_df.loc[col_mask, "y_pred_pH74"] - combos_df.loc[col_mask, "y_true_pH74"]).mean()
            if col_mask.any() else np.nan
        )
    else:
        col_bias = np.nan

    lines = []
    lines.append("# sdab R4 LOO-CV Training Report\n")
    lines.append("## Best Hyperparameters\n")
    lines.append(f"- pH7.4: alpha={best74[1]}, l1_ratio={best74[2]}")
    lines.append(f"- pH6.0: alpha={best6[1]}, l1_ratio={best6[2]}\n")
    lines.append("## LOO Metrics (combo hold-outs)\n")
    lines.append("| pH | Spearman rho | Pearson R | R2 | MAE |")
    lines.append("|----|--------------|-----------|-----|-----|")
    lines.append(f"| 7.4 | {m74[0]:.4f} | {m74[2]:.4f} | {m74[1]:.4f} | {m74[3]:.4f} |")
    lines.append(f"| 6.0 | {m6[0]:.4f} | {m6[2]:.4f} | {m6[1]:.4f} | {m6[3]:.4f} |\n")
    lines.append("## Top-10 singles by |delta_beta|\n")
    lines.append("| feature | beta_pH74 | beta_pH6 | delta_beta |")
    lines.append("|---------|-----------|----------|------------|")
    for _, r in top_delta.iterrows():
        lines.append(f"| {r['feature']} | {r['beta_pH74']:.4f} | {r['beta_pH6']:.4f} | {r['delta_beta']:.4f} |")
    lines.append("")
    lines.append("## Top-10 negative pairs by |beta_pH74|\n")
    lines.append("| feature | beta_pH74 | beta_pH6 | delta_beta |")
    lines.append("|---------|-----------|----------|------------|")
    for _, r in top_neg_pairs.iterrows():
        lines.append(f"| {r['feature']} | {r['beta_pH74']:.4f} | {r['beta_pH6']:.4f} | {r['delta_beta']:.4f} |")
    lines.append("")
    lines.append("## LOO Predictions (combos)\n")
    lines.append("| name | y_true_pH74 | y_pred_pH74 | y_true_pH6 | y_pred_pH6 | collapsed |")
    lines.append("|------|-------------|-------------|------------|------------|-----------|")
    for _, r in combos_df.iterrows():
        lines.append(
            f"| {r['name']} | {r['y_true_pH74']:.4f} | {r['y_pred_pH74']:.4f} | "
            f"{r['y_true_pH6']:.4f} | {r['y_pred_pH6']:.4f} | {bool(r.get('collapsed', False))} |"
        )
    lines.append("")
    lines.append("## Trend\n")
    if not np.isnan(col_bias):
        direction = "under-predicted" if col_bias < 0 else "over-predicted"
        lines.append(
            f"Collapsed samples are on average {direction} at pH7.4 by {abs(col_bias):.4f} log-units."
        )
    else:
        lines.append("No collapsed samples in training set (bias undefined).")
    path.write_text("\n".join(lines))


def main():
    csv_path = Path(os.environ.get("TRAIN_CSV", str(DEFAULT_CSV)))
    df = pd.read_csv(csv_path)
    df = df[df["source"].isin(["singles", "combos"])].reset_index(drop=True)

    idx_lists = [mutations_to_indices(s) for s in df["mutations_unified"].fillna("")]
    pairs = build_pair_whitelist(idx_lists, min_count=2)
    pair_records = [(i, j, f"{SINGLE_NAMES[i]}-{SINGLE_NAMES[j]}") for i, j in pairs]

    X = build_feature_matrix(idx_lists, pairs)
    y74 = df["log_pH74"].to_numpy()
    y6 = df["log_pH6"].to_numpy()
    combo_mask = (df["source"] == "combos").to_numpy()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        best74 = grid_search(X, y74, combo_mask)
        best6 = grid_search(X, y6, combo_mask)

        model74 = fit_elasticnet(X, y74, best74[1], best74[2])
        model6 = fit_elasticnet(X, y6, best6[1], best6[2])

    joblib.dump(model74, OUT_DIR / "model_pH74.joblib")
    joblib.dump(model6, OUT_DIR / "model_pH6.joblib")

    with open(OUT_DIR / "features.json", "w") as f:
        json.dump({
            "singles": SINGLE_NAMES,
            "pairs": [[i, j, name] for (i, j, name) in pair_records],
        }, f, indent=2)

    coef_df = write_coefficients(
        OUT_DIR / "coefficients.csv", SINGLE_NAMES, pair_records,
        model74.coef_, model6.coef_,
    )

    combo_names = df.loc[combo_mask, "name"].to_numpy()
    collapsed = df.loc[combo_mask, "collapsed"].to_numpy() if "collapsed" in df.columns else np.zeros(combo_mask.sum(), dtype=bool)
    combos_df = pd.DataFrame({
        "name": combo_names,
        "y_true_pH74": best74[3], "y_pred_pH74": best74[4],
        "y_true_pH6": best6[3], "y_pred_pH6": best6[4],
        "collapsed": collapsed,
    })

    m74 = metrics(best74[3], best74[4])
    m6 = metrics(best6[3], best6[4])

    write_report(OUT_DIR / "loo_cv_report.md", best74, best6, m74, m6, coef_df, combos_df)

    rho74, rho6 = m74[0], m6[0]
    pair_names = {name: (coef74, coef6) for name, coef74, coef6 in zip(
        [r["feature"] for _, r in coef_df[coef_df["type"] == "pair"].iterrows()],
        coef_df.loc[coef_df["type"] == "pair", "beta_pH74"].to_numpy(),
        coef_df.loc[coef_df["type"] == "pair", "beta_pH6"].to_numpy(),
    )}
    key_pairs = ["HV105H-HE108H", "HV105H-HY111H"]
    best_key = None
    for kp in key_pairs:
        if kp in pair_names:
            b = pair_names[kp][0]
            if best_key is None or b < best_key[1]:
                best_key = (kp, b)

    gate1 = rho74 >= 0.5
    gate2 = rho6 >= 0.5
    gate3 = best_key is not None and best_key[1] < -0.05

    print("=== train_additive_model.py acceptance gates ===")
    print(f"[{'PASS' if gate1 else 'FAIL'}] pH74 LOO Spearman rho >= 0.5: {rho74:.4f}")
    print(f"[{'PASS' if gate2 else 'FAIL'}] pH6  LOO Spearman rho >= 0.5: {rho6:.4f}")
    if best_key is None:
        print("[FAIL] at least one of {V105H-E108H, V105H-Y111H} has beta_pH74 < -0.05: neither pair in whitelist")
    else:
        print(f"[{'PASS' if gate3 else 'FAIL'}] at least one of {{V105H-E108H, V105H-Y111H}} has beta_pH74 < -0.05: {best_key[0]}={best_key[1]:.4f}")

    if not (gate1 and gate2 and gate3):
        sys.exit(1)


if __name__ == "__main__":
    main()
