#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent verification of sdab R4 training data + additive model.

Does NOT import from build_training_data.py or train_additive_model.py.
Re-derives computations from raw CSVs and compares to on-disk artifacts.
"""

import json
import sys
import warnings
from itertools import combinations
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNet

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from analysis.naming.convert import parse_unified  # noqa: E402

RAW_ELISA = REPO_ROOT / "experiments/sdab/R2/data/sdab_elisa_raw.csv"
COMBOS_CSV = REPO_ROOT / "experiments/sdab/R2/data/hs32-92_8ng_ml_Elisa_results.csv"
TRAIN_CSV = REPO_ROOT / "experiments/sdab/R4/data/training_data.csv"
MODEL_DIR = REPO_ROOT / "experiments/sdab/R4/model"
FEATURES_JSON = MODEL_DIR / "features.json"
COEF_CSV = MODEL_DIR / "coefficients.csv"

OD_FLOOR = 0.05
W4, W20 = 0.569, 0.431
ABS_EPS = 1e-9


def interp_8ng(raw, sample_key):
    """Reproduce log-interp at 8 ng/mL for a single hs_key (e.g. 'hs1', 'wt')."""
    sub = raw[raw["sample"].str.split("-").str[-1] == sample_key]
    sub = sub[sub["concentration_ng_ml"].isin([4.0, 20.0])]
    out = {}
    for ph_col in ("od_ph74", "od_ph60"):
        v4 = sub.loc[sub["concentration_ng_ml"] == 4.0, ph_col].mean()
        v20 = sub.loc[sub["concentration_ng_ml"] == 20.0, ph_col].mean()
        log4 = np.log(max(v4, OD_FLOOR))
        log20 = np.log(max(v20, OD_FLOOR))
        out[ph_col] = W4 * log4 + W20 * log20
    return out["od_ph74"], out["od_ph60"]


def spot_check_row(name, got74, got6, exp74, exp6):
    d74, d6 = abs(got74 - exp74), abs(got6 - exp6)
    ok = d74 < ABS_EPS and d6 < ABS_EPS
    detail = ""
    if not ok:
        detail = (f" (pH74 got={got74:.9f} exp={exp74:.9f} d={d74:.2e}; "
                  f"pH6 got={got6:.9f} exp={exp6:.9f} d={d6:.2e})")
    return ok, detail


def t1_spot_checks(df):
    raw = pd.read_csv(RAW_ELISA)
    results = {}

    # WT
    g74, g6 = interp_8ng(raw, "wt")
    row = df[df["name"] == "WT"].iloc[0]
    results["WT"] = spot_check_row("WT", g74, g6, row["log_pH74"], row["log_pH6"])

    # hs1 single
    g74, g6 = interp_8ng(raw, "hs1")
    row = df[df["name"] == "hs1"].iloc[0]
    results["hs1"] = spot_check_row("hs1", g74, g6, row["log_pH74"], row["log_pH6"])

    # hs69 combo (fully collapsed)
    combos = pd.read_csv(COMBOS_CSV, encoding="utf-8-sig")
    hit = combos[combos["ID"].astype(str).str.contains("hs69")].iloc[0]
    od74 = float(str(hit["D-pH74_avg"]).strip())
    od6 = float(str(hit["D-pH60_Avg"]).strip())
    g74 = np.log(max(od74, OD_FLOOR))
    g6 = np.log(max(od6, OD_FLOOR))
    row = df[df["name"] == "hs69"].iloc[0]
    results["hs69"] = spot_check_row("hs69", g74, g6, row["log_pH74"], row["log_pH6"])
    return results


def t2_integrity(df):
    failures = []
    if len(df) != 93:
        failures.append(f"row count {len(df)} != 93")
    counts = df["source"].value_counts().to_dict()
    if counts.get("WT", 0) != 1:
        failures.append(f"source=WT count {counts.get('WT', 0)} != 1")
    if counts.get("singles", 0) != 31:
        failures.append(f"source=singles count {counts.get('singles', 0)} != 31")
    if counts.get("combos", 0) != 61:
        failures.append(f"source=combos count {counts.get('combos', 0)} != 61")
    if df[["log_pH74", "log_pH6", "log_ratio"]].isna().any().any():
        failures.append("NaN present in log columns")
    diff = (df["log_ratio"] - (df["log_pH74"] - df["log_pH6"])).abs().max()
    if diff >= 1e-9:
        failures.append(f"log_ratio identity violated (max diff {diff:.2e})")
    combos = df[df["source"] == "combos"]
    missing_d110h = [n for n, m in zip(combos["name"], combos["mutations_unified"])
                     if "HD110H" not in m.split(";")]
    if missing_d110h:
        failures.append(f"{len(missing_d110h)} combos missing HD110H: {missing_d110h[:3]}")
    for _, r in df.iterrows():
        m = r["mutations_unified"]
        if not isinstance(m, str) or not m:
            continue
        for tok in m.split(";"):
            try:
                parse_unified(tok)
            except Exception as e:
                failures.append(f"parse_unified fail on {tok!r} ({r['name']}): {e}")
                break
    return failures


def _mut_to_idx(m, pos_to_idx):
    if not isinstance(m, str) or not m:
        return []
    return sorted(pos_to_idx[t] for t in m.split(";") if t in pos_to_idx)


def _build_X(idx_lists, pairs, n_single):
    pair_idx = {tuple(p[:2]): k for k, p in enumerate(pairs)}
    X = np.zeros((len(idx_lists), n_single + len(pairs)), dtype=np.float64)
    for r, idxs in enumerate(idx_lists):
        for i in idxs:
            X[r, i] = 1.0
        for i, j in combinations(idxs, 2):
            if (i, j) in pair_idx:
                X[r, n_single + pair_idx[(i, j)]] = 1.0
    return X


def t3_in_sample(df):
    with open(FEATURES_JSON) as f:
        feat = json.load(f)
    singles = feat["singles"]
    pairs = feat["pairs"]
    pos_to_idx = {s: i for i, s in enumerate(singles)}

    sub = df[df["source"].isin(["singles", "combos"])].reset_index(drop=True)
    idx_lists = [_mut_to_idx(m, pos_to_idx) for m in sub["mutations_unified"].fillna("")]
    X = _build_X(idx_lists, pairs, len(singles))

    m74 = joblib.load(MODEL_DIR / "model_pH74.joblib")
    m6 = joblib.load(MODEL_DIR / "model_pH6.joblib")
    p74 = m74.predict(X)
    p6 = m6.predict(X)
    r74 = pearsonr(sub["log_pH74"].to_numpy(), p74)[0] ** 2
    r6 = pearsonr(sub["log_pH6"].to_numpy(), p6)[0] ** 2
    return r74, r6, m74, m6, X, sub


def t4_loo_one(df, X, sub, m74, target="hs60"):
    idx = sub.index[sub["name"] == target]
    if len(idx) == 0:
        return None
    i = int(idx[0])
    y = sub["log_pH74"].to_numpy()
    mask = np.ones(len(y), dtype=bool)
    mask[i] = False
    alpha = float(getattr(m74, "alpha", 0.001))
    l1 = float(getattr(m74, "l1_ratio", 0.9))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        model = ElasticNet(alpha=alpha, l1_ratio=l1, max_iter=10000,
                           tol=1e-5, random_state=42)
        model.fit(X[mask], y[mask])
    pred = float(model.predict(X[i:i + 1])[0])
    actual = float(y[i])
    return target, pred, actual, abs(pred - actual)


def t5_coef(coef_df):
    failures = []
    row_ve = coef_df[coef_df["feature"] == "HV105H-HE108H"]
    if row_ve.empty:
        failures.append("HV105H-HE108H pair missing from coefficients.csv")
        ve_ok = False
    else:
        ve_ok = float(row_ve.iloc[0]["beta_pH74"]) < 0
        if not ve_ok:
            failures.append(f"HV105H-HE108H beta_pH74={row_ve.iloc[0]['beta_pH74']:.4f} not < 0")

    row_d = coef_df[coef_df["feature"] == "HD110H"]
    if row_d.empty:
        failures.append("HD110H missing")
        d_ok = False
    else:
        b74 = float(row_d.iloc[0]["beta_pH74"])
        b6 = float(row_d.iloc[0]["beta_pH6"])
        d_ok = (b74 != 0.0) or (b6 != 0.0)
        if not d_ok:
            failures.append(f"HD110H both beta zero (b74={b74}, b6={b6})")

    singles = coef_df[coef_df["type"] == "single"].copy()
    singles["abs_d"] = singles["delta_beta"].abs()
    top_singles = singles.sort_values("abs_d", ascending=False).head(3)

    pairs = coef_df[coef_df["type"] == "pair"].copy()
    neg_pairs = pairs[pairs["beta_pH74"] < 0].sort_values("beta_pH74").head(3)

    return ve_ok, d_ok, failures, top_singles, neg_pairs


def main():
    df = pd.read_csv(TRAIN_CSV)
    coef_df = pd.read_csv(COEF_CSV)

    t1 = t1_spot_checks(df)
    t2 = t2_integrity(df)
    r74, r6, m74, m6, X, sub = t3_in_sample(df)
    t4 = t4_loo_one(df, X, sub, m74, target="hs60")
    ve_ok, d_ok, t5_fail, top_singles, neg_pairs = t5_coef(coef_df)

    all_pass = True
    print("=== verify_training_and_model.py summary ===")
    for key in ("WT", "hs1", "hs69"):
        ok, det = t1[key]
        tag = "PASS" if ok else "FAIL"
        print(f"[T1] {key} spot-check: {tag}{det}")
        if not ok:
            all_pass = False

    if t2:
        all_pass = False
        print(f"[T2] training_data.csv shape + integrity: FAIL ({'; '.join(t2)})")
    else:
        print("[T2] training_data.csv shape + integrity: PASS")

    gate_r = 0.7
    tag74 = "PASS" if r74 >= gate_r else "FAIL"
    tag6 = "PASS" if r6 >= gate_r else "FAIL"
    print(f"[T3] in-sample R^2 pH74: {r74:.4f}  [gate: >={gate_r}] {tag74}")
    print(f"[T3] in-sample R^2 pH6:  {r6:.4f}  [gate: >={gate_r}] {tag6}")
    if r74 < gate_r or r6 < gate_r:
        all_pass = False

    if t4 is None:
        print("[T4] LOO cross-check (hs60): FAIL (sample not found)")
        all_pass = False
    else:
        name, pred, actual, adiff = t4
        tagt4 = "PASS" if adiff <= 0.5 else "FAIL"
        print(f"[T4] LOO cross-check ({name}): predicted={pred:.4f}, "
              f"actual={actual:.4f}, abs_diff={adiff:.4f}  [gate: <=0.5] {tagt4}")
        if adiff > 0.5:
            all_pass = False

    print(f"[T5] V105H-E108H beta_pH74 < 0: {'PASS' if ve_ok else 'FAIL'}")
    print(f"[T5] D110H coefficients nonzero: {'PASS' if d_ok else 'FAIL'}")
    if not ve_ok or not d_ok:
        all_pass = False
    if t5_fail:
        for msg in t5_fail:
            print(f"       {msg}")

    print("[T5] top 3 |delta_beta| singles: " + ", ".join(
        f"{r['feature']}(Δβ={r['delta_beta']:+.3f})" for _, r in top_singles.iterrows()
    ))
    print("[T5] top 3 most-negative beta_pH74 pairs: " + ", ".join(
        f"{r['feature']}(β74={r['beta_pH74']:+.3f})" for _, r in neg_pairs.iterrows()
    ))

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
