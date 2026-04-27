#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build sdab R4 training data (91 ELISA + WT)."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from analysis.naming.convert import parse_unified, r1_to_unified  # noqa: E402

RAW_ELISA = REPO_ROOT / "experiments/sdab/R2/data/sdab_elisa_raw.csv"
VARIANTS = REPO_ROOT / "experiments/sdab/R2/data/sdab_variants.csv"
COMBOS = REPO_ROOT / "experiments/sdab/R2/data/hs32-92_8ng_ml_Elisa_results.csv"
OUT_DIR = REPO_ROOT / "experiments/sdab/R4/data"
OUT_CSV = OUT_DIR / "training_data.csv"
OUT_JSON = OUT_DIR / "wt_baseline.json"

OD_FLOOR = 0.05
# log-linear weights for 8 ng/mL between 4 and 20 ng/mL:
# log(8) lies at 0.569 of the way from log(4) to log(20).
W4, W20 = 0.569, 0.431

# Samples flagged for leverage inspection in LOO diagnostics.
# hs24 (HA104H): max OD across the dose curve is 0.24 — likely non-expressed or mis-plated.
LEVERAGE_HS = {"hs24"}


def build_singles_and_wt():
    raw = pd.read_csv(RAW_ELISA)
    raw["hs_key"] = raw["sample"].str.split("-").str[-1]  # SNb512-hs1 -> hs1
    raw = raw[raw["concentration_ng_ml"].isin([4.0, 20.0])].copy()
    mean = raw.groupby(["hs_key", "concentration_ng_ml"], as_index=False)[
        ["od_ph74", "od_ph60"]
    ].mean()

    rows = []
    for hs_key, grp in mean.groupby("hs_key"):
        vals = {"hs_key": hs_key}
        for ph_col in ("od_ph74", "od_ph60"):
            v4 = grp.loc[grp["concentration_ng_ml"] == 4.0, ph_col].iloc[0]
            v20 = grp.loc[grp["concentration_ng_ml"] == 20.0, ph_col].iloc[0]
            log4 = np.log(max(v4, OD_FLOOR))
            log20 = np.log(max(v20, OD_FLOOR))
            vals[ph_col] = W4 * log4 + W20 * log20
        rows.append(vals)
    interp_df = pd.DataFrame(rows)

    variants = pd.read_csv(VARIANTS)
    hs_to_mut = {
        row["hs_num"] if str(row["hs_num"]).startswith("hs") else f"hs{int(row['hs_num'])}":
        r1_to_unified(row["mutant_id"])
        for _, row in variants.iterrows()
    }

    out = []
    for _, r in interp_df.iterrows():
        hs_key = r["hs_key"]
        log74 = float(r["od_ph74"])
        log60 = float(r["od_ph60"])
        if hs_key == "wt":
            out.append({
                "name": "WT",
                "mutations_unified": "",
                "order": 0,
                "log_pH74": log74,
                "log_pH6": log60,
                "log_ratio": log74 - log60,
                "source": "WT",
                "group": "WT",
                "collapsed": False,
                "leverage_flag": False,
            })
        else:
            if hs_key not in hs_to_mut:
                raise KeyError(f"{hs_key!r} in {RAW_ELISA.name} has no mutation mapping in {VARIANTS.name}")
            mut = hs_to_mut[hs_key]
            out.append({
                "name": hs_key,
                "mutations_unified": mut,
                "order": 1,
                "log_pH74": log74,
                "log_pH6": log60,
                "log_ratio": log74 - log60,
                "source": "singles",
                "group": "interp",
                "collapsed": False,
                "leverage_flag": (hs_key in LEVERAGE_HS),
            })
    return out


def build_combos():
    df = pd.read_csv(COMBOS, encoding="utf-8-sig")
    out = []
    for _, r in df.iterrows():
        sample_id = str(r["ID"])
        hs_key = sample_id.split("-")[-1]  # 5125-Fcm3A-hs32 -> hs32
        raw_muts = str(r["Mutations"])
        tokens = [t.strip() for t in raw_muts.split(";") if t.strip()]
        mutations_unified = ";".join(f"H{t}" for t in tokens)
        order = len(tokens)
        od74 = float(str(r["D-pH74_avg"]).strip())
        od60 = float(str(r["D-pH60_Avg"]).strip())
        log74 = np.log(max(od74, OD_FLOOR))
        log60 = np.log(max(od60, OD_FLOOR))
        out.append({
            "name": hs_key,
            "mutations_unified": mutations_unified,
            "order": order,
            "log_pH74": log74,
            "log_pH6": log60,
            "log_ratio": log74 - log60,
            "source": "combos",
            "group": "D",
            "collapsed": log74 < -1.5,  # combos only; final value set in main()
            "leverage_flag": False,
        })
    return out


def sort_rows(rows):
    def key(r):
        if r["name"] == "WT":
            return (0, 0)
        return (1, int(r["name"][2:]))
    return sorted(rows, key=key)


def self_check(df, wt_baseline):
    checks = []

    row_ok = len(df) == 93
    checks.append((row_ok, f"row count == 93: {len(df)}"))

    ph74_raw = wt_baseline["pH74_raw_interp_8ng"]
    checks.append((3.45 <= ph74_raw <= 3.65,
                   f"WT pH74_raw within 3.45-3.65: {ph74_raw:.4f}"))

    ph6_raw = wt_baseline["pH6_raw_interp_8ng"]
    checks.append((3.28 <= ph6_raw <= 3.48,
                   f"WT pH6_raw within 3.28-3.48: {ph6_raw:.4f}"))

    collapsed_df = df.loc[df["collapsed"], ["name", "mutations_unified", "log_pH74"]]
    collapsed_names = set(collapsed_df["name"])
    must_include = {"hs69", "hs70", "hs89", "hs92"}
    print("[INFO] collapsed combos:")
    for _, r in collapsed_df.iterrows():
        print(f"       {r['name']:>6}  {r['mutations_unified']:<40}  log_pH74={r['log_pH74']:+.3f}")
    checks.append((must_include.issubset(collapsed_names),
                   f"collapsed set must include {sorted(must_include)}: got {sorted(collapsed_names)}"))

    combos = df[df["source"] == "combos"]
    all_d110h = all("HD110H" in m.split(";") for m in combos["mutations_unified"])
    checks.append((all_d110h, f"all {len(combos)} combos contain HD110H: {int(all_d110h)}"))

    parse_ok = True
    parse_err = None
    for m in df["mutations_unified"]:
        if not m:
            continue
        for tok in m.split(";"):
            try:
                parse_unified(tok)
            except Exception as e:
                parse_ok = False
                parse_err = f"{tok}: {e}"
                break
        if not parse_ok:
            break
    checks.append((parse_ok,
                   f"every mutations_unified parses with analysis.naming.convert.parse_unified: {parse_ok}"
                   + (f" (err={parse_err})" if parse_err else "")))

    no_nan = not df[["log_pH74", "log_pH6", "log_ratio"]].isna().any().any()
    checks.append((no_nan, f"no NaN in log_pH74, log_pH6, log_ratio: {no_nan}"))

    print("=== build_training_data.py self-check ===")
    all_pass = True
    for ok, msg in checks:
        tag = "[PASS]" if ok else "[FAIL]"
        print(f"{tag} {msg}")
        if not ok:
            all_pass = False
    return all_pass


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    singles_wt = build_singles_and_wt()
    combos = build_combos()
    all_rows = sort_rows(singles_wt + combos)
    df = pd.DataFrame(all_rows, columns=[
        "name", "mutations_unified", "order",
        "log_pH74", "log_pH6", "log_ratio",
        "source", "group", "collapsed", "leverage_flag",
    ])
    # collapsed = pair-induced pH7.4 collapse in combos (Lasso negative-pair diagnostic).
    # Low-signal singles are weak binders, not collapsed pairs — gate on source == "combos".
    df["collapsed"] = (df["source"] == "combos") & (df["log_pH74"] < -1.5)
    df.to_csv(OUT_CSV, index=False)

    wt_row = df[df["name"] == "WT"].iloc[0]
    wt_baseline = {
        "log_pH74_WT": float(wt_row["log_pH74"]),
        "log_pH6_WT": float(wt_row["log_pH6"]),
        "log_ratio_WT": float(wt_row["log_ratio"]),
        "pH74_raw_interp_8ng": float(np.exp(wt_row["log_pH74"])),
        "pH6_raw_interp_8ng": float(np.exp(wt_row["log_pH6"])),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(wt_baseline, f, indent=2)

    print(f"Wrote {OUT_CSV} ({len(df)} rows)")
    print(f"Wrote {OUT_JSON}")
    print(f"WT baseline: {wt_baseline}")

    ok = self_check(df, wt_baseline)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
