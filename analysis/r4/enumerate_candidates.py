#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Enumerate and score sdab R4 candidates (2/3/4/5联 CDR His combos)."""

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

MODEL_DIR = REPO_ROOT / "experiments/sdab/R4/model"
DATA_DIR = REPO_ROOT / "experiments/sdab/R4/data"
OUT_DIR = REPO_ROOT / "experiments/sdab/R4/predictions"

TRAIN_CSV = DATA_DIR / "training_data.csv"
WT_JSON = DATA_DIR / "wt_baseline.json"
FEATURES_JSON = MODEL_DIR / "features.json"
MODEL_PH74 = MODEL_DIR / "model_pH74.joblib"
MODEL_PH6 = MODEL_DIR / "model_pH6.joblib"

OUT_ALL = OUT_DIR / "candidates_all.csv"
OUT_TOP50 = OUT_DIR / "candidates_top50.csv"
OUT_DIST = OUT_DIR / "score_distribution.txt"

# Scoring weights (tweakable; plan allows grid search if top-50 is uninspiring)
SCORE_ALPHA = 1.0   # pH6 penalty
SCORE_LAMBDA = 1.0  # pH74 collapse hinge
SCORE_MU = 0.5      # order-5 extrapolation penalty (raised to cap top-50 share ≤ 30%)
PH74_FLOOR_DROP = 0.8  # drop if log_pH74_pred < log_pH74_WT - 0.8

ORDERS = [2, 3, 4, 5]

# 15 covered positions — from hs32-91 combos
COVERED_NAMES = {
    "HG26H", "HR27H", "HT28H", "HS30H", "HS31H",
    "HG53H", "HN54H", "HG56H", "HS57H",
    "HQ100H", "HG102H", "HV105H", "HE108H", "HD110H", "HY111H",
}

# CDR region by position number (for diversity reporting)
def cdr_of(pos):
    if 26 <= pos <= 33:
        return "CDR1"
    if 51 <= pos <= 58:
        return "CDR2"
    if 97 <= pos <= 111:
        return "CDR3"
    return "?"


def load_artifacts():
    with open(FEATURES_JSON) as f:
        feats = json.load(f)
    singles = feats["singles"]                                # 31 names
    pairs = feats["pairs"]                                    # list of [i, j, name]
    pos_to_idx = {name: i for i, name in enumerate(singles)}
    # pair_pos: 2 arrays of length 40 for vectorized co-occurrence check
    pair_p = np.array([p[0] for p in pairs], dtype=np.int64)  # first index
    pair_q = np.array([p[1] for p in pairs], dtype=np.int64)  # second index
    pair_names = [p[2] for p in pairs]

    model_pH74 = joblib.load(MODEL_PH74)
    model_pH6 = joblib.load(MODEL_PH6)

    with open(WT_JSON) as f:
        wt = json.load(f)
    log_pH74_WT = wt["log_pH74_WT"]
    log_pH6_WT = wt["log_pH6_WT"]

    return singles, pos_to_idx, pair_p, pair_q, pair_names, model_pH74, model_pH6, log_pH74_WT, log_pH6_WT


def load_measured(pos_to_idx):
    """Return set of frozenset(indices) for all measured non-WT samples (to exclude)."""
    df = pd.read_csv(TRAIN_CSV)
    out = set()
    for s in df["mutations_unified"].fillna(""):
        if not s:
            continue
        idxs = tuple(sorted(pos_to_idx[m] for m in s.split(";")))
        out.add(idxs)
    return out


def enumerate_order(k, pos_to_idx, pair_p, pair_q, measured, model_pH74, model_pH6, log_pH74_WT):
    """Enumerate all C(31, k), filter, predict, score. Returns a DataFrame."""
    n_singles = len(pos_to_idx)
    v105 = pos_to_idx["HV105H"]
    d110 = pos_to_idx["HD110H"]
    e108 = pos_to_idx["HE108H"]
    y111 = pos_to_idx["HY111H"]

    combos = np.array(list(combinations(range(n_singles), k)), dtype=np.int64)
    N = len(combos)

    # Hard filter 1 + 2: banned co-occurrences (vectorized)
    has_v = (combos == v105).any(axis=1)
    has_d = (combos == d110).any(axis=1)
    has_e = (combos == e108).any(axis=1)
    has_y = (combos == y111).any(axis=1)
    banned = (has_v & has_d) | (has_v & (has_e | has_y))

    # Exclude already-measured
    measured_arr = np.zeros(N, dtype=bool)
    for i, row in enumerate(combos):
        if tuple(row) in measured:
            measured_arr[i] = True

    keep_mask = ~(banned | measured_arr)
    combos = combos[keep_mask]
    N = len(combos)

    if N == 0:
        return pd.DataFrame(), 0

    # Feature matrix: singles one-hot + pair one-hot
    n_pairs = len(pair_p)
    X = np.zeros((N, n_singles + n_pairs), dtype=np.float32)

    # Singles: for each row, X[row, combos[row, :]] = 1
    rows = np.repeat(np.arange(N), k)
    cols = combos.ravel()
    X[rows, cols] = 1.0

    # Pairs: for each whitelisted (p, q), combos containing both get their feature set
    for idx, (p, q) in enumerate(zip(pair_p, pair_q)):
        both = (combos == p).any(axis=1) & (combos == q).any(axis=1)
        X[both, n_singles + idx] = 1.0

    # Predict both pH
    y74 = model_pH74.predict(X)
    y6 = model_pH6.predict(X)

    # Hard filter 3: pH74 floor drop
    ph74_ok = y74 >= (log_pH74_WT - PH74_FLOOR_DROP)
    combos = combos[ph74_ok]
    y74 = y74[ph74_ok]
    y6 = y6[ph74_ok]
    X = X[ph74_ok]
    N = len(combos)

    if N == 0:
        return pd.DataFrame(), 0

    # Score
    is_order_5 = 1.0 if k == 5 else 0.0
    hinge = np.maximum(0.0, log_pH74_WT - y74)
    score = y74 - SCORE_ALPHA * y6 - SCORE_LAMBDA * hinge - SCORE_MU * is_order_5

    # Build output rows
    singles_names = list(pos_to_idx.keys())
    mutations_str = [";".join(singles_names[i] for i in row) for row in combos]

    # Confidence stratification
    covered_idx = {i for i, nm in enumerate(singles_names) if nm in COVERED_NAMES}
    pair_lookup = {(int(p), int(q)): idx for idx, (p, q) in enumerate(zip(pair_p, pair_q))}

    confidences = []
    for row in combos:
        idx_set = set(int(i) for i in row)
        all_in_covered = idx_set.issubset(covered_idx)
        if k <= 4 and all_in_covered:
            confidences.append("high")
        elif k == 5:
            # all pairs present in whitelist?
            pair_iter = combinations(sorted(row), 2)
            all_pairs_known = all((int(a), int(b)) in pair_lookup for a, b in pair_iter)
            if all_pairs_known:
                confidences.append("medium")
            else:
                confidences.append("low")
        else:
            # order 2-4 with at least one uncovered position — single data exists for all 31
            confidences.append("medium")

    # CDR presence per candidate + D110H flag
    cdr1_count, cdr2_count, cdr3_count, has_d110h = [], [], [], []
    # Map idx → CDR bucket once
    idx_to_cdr = {}
    for name, idx in pos_to_idx.items():
        pos = int("".join(c for c in name[2:-1] if c.isdigit()))
        idx_to_cdr[idx] = cdr_of(pos)

    for row in combos:
        c1 = c2 = c3 = 0
        for i in row:
            tag = idx_to_cdr[int(i)]
            if tag == "CDR1":
                c1 += 1
            elif tag == "CDR2":
                c2 += 1
            elif tag == "CDR3":
                c3 += 1
        cdr1_count.append(c1)
        cdr2_count.append(c2)
        cdr3_count.append(c3)
        has_d110h.append(d110 in row)

    df = pd.DataFrame({
        "mutations_unified": mutations_str,
        "order": k,
        "log_pH74_pred": y74,
        "log_pH6_pred": y6,
        "log_ratio_pred": y74 - y6,
        "score": score,
        "confidence": confidences,
        "positions_cdr1": cdr1_count,
        "positions_cdr2": cdr2_count,
        "positions_cdr3": cdr3_count,
        "contains_D110H": has_d110h,
    })
    return df, N


def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    (singles, pos_to_idx, pair_p, pair_q, pair_names,
     model_pH74, model_pH6, log_pH74_WT, log_pH6_WT) = load_artifacts()
    measured = load_measured(pos_to_idx)
    print(f"Loaded: {len(singles)} singles, {len(pair_p)} pair features, "
          f"{len(measured)} measured combos to exclude.")
    print(f"WT baseline: log_pH74={log_pH74_WT:.4f}, log_pH6={log_pH6_WT:.4f}")

    all_frames = []
    print(f"\n{'Order':>5} {'EnumTotal':>10} {'AfterFilter':>12} {'Elapsed':>8}")
    for k in ORDERS:
        ts = time.time()
        df_k, n_kept = enumerate_order(
            k, pos_to_idx, pair_p, pair_q, measured,
            model_pH74, model_pH6, log_pH74_WT,
        )
        elapsed = time.time() - ts
        total_k = 1
        for i in range(k):
            total_k = total_k * (len(singles) - i) // (i + 1)
        print(f"{k:>5} {total_k:>10} {n_kept:>12} {elapsed:>7.2f}s")
        if n_kept > 0:
            all_frames.append(df_k)

    df_all = pd.concat(all_frames, ignore_index=True)
    df_all = df_all.sort_values("score", ascending=False).reset_index(drop=True)

    df_all.to_csv(OUT_ALL, index=False)
    df_top = df_all.head(50)
    df_top.to_csv(OUT_TOP50, index=False)

    # Distribution report
    lines = [f"sdab R4 candidate distribution (total {len(df_all)})", ""]
    lines.append("By order:")
    lines.append(df_all.groupby("order").size().to_string())
    lines.append("")
    lines.append("By confidence:")
    lines.append(df_all.groupby("confidence").size().to_string())
    lines.append("")
    lines.append("By order × confidence:")
    lines.append(df_all.groupby(["order", "confidence"]).size().unstack(fill_value=0).to_string())
    lines.append("")
    lines.append("Score summary:")
    lines.append(df_all["score"].describe().to_string())
    lines.append("")
    lines.append("Top-50 composition:")
    lines.append("  order counts: " + df_top["order"].value_counts().sort_index().to_string().replace("\n", " | "))
    lines.append("  confidence:   " + df_top["confidence"].value_counts().to_string().replace("\n", " | "))
    lines.append(f"  D110H-containing: {int(df_top['contains_D110H'].sum())} / 50")
    OUT_DIST.write_text("\n".join(lines) + "\n")

    print(f"\nWrote {OUT_ALL} ({len(df_all)} candidates)")
    print(f"Wrote {OUT_TOP50} (top 50)")
    print(f"Wrote {OUT_DIST}")

    # Sanity checks
    print(f"\n=== enumerate_candidates.py sanity checks ===")
    checks = []
    v105 = pos_to_idx["HV105H"]
    d110 = pos_to_idx["HD110H"]
    e108 = pos_to_idx["HE108H"]
    y111 = pos_to_idx["HY111H"]
    singles_names = list(pos_to_idx.keys())

    def contains_pair(mut_str, p, q):
        idxs = {pos_to_idx[m] for m in mut_str.split(";")}
        return p in idxs and q in idxs

    bad_vd = df_top["mutations_unified"].apply(lambda s: contains_pair(s, v105, d110)).sum()
    bad_ve = df_top["mutations_unified"].apply(lambda s: contains_pair(s, v105, e108)).sum()
    bad_vy = df_top["mutations_unified"].apply(lambda s: contains_pair(s, v105, y111)).sum()
    checks.append((bad_vd == 0, f"top50 no V105H+D110H: {bad_vd}"))
    checks.append((bad_ve == 0, f"top50 no V105H+E108H: {bad_ve}"))
    checks.append((bad_vy == 0, f"top50 no V105H+Y111H: {bad_vy}"))

    non_d110h_in_top5 = (~df_top.head(5)["contains_D110H"]).sum()
    checks.append((non_d110h_in_top5 >= 1,
                   f"top5 has ≥1 non-D110H candidate (D110H decoupling): {non_d110h_in_top5} non-D110H in top5"))

    order5_share = (df_top["order"] == 5).sum() / 50.0
    checks.append((order5_share <= 0.30,
                   f"top50 order-5 share ≤ 30%: {order5_share:.2%}"))

    orders_present = set(df_all["order"].unique())
    checks.append((orders_present == {2, 3, 4, 5},
                   f"candidates_all covers orders {{2,3,4,5}}: got {sorted(orders_present)}"))
    confidence_present = set(df_all["confidence"].unique())
    checks.append((confidence_present == {"high", "medium", "low"},
                   f"candidates_all covers confidence {{high,medium,low}}: got {sorted(confidence_present)}"))

    wall = time.time() - t0
    checks.append((wall <= 120.0, f"total wall time ≤ 120s: {wall:.1f}s"))

    all_pass = True
    for ok, msg in checks:
        tag = "[PASS]" if ok else "[FAIL]"
        print(f"{tag} {msg}")
        if not ok:
            all_pass = False

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
