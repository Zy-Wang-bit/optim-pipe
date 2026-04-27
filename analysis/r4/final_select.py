#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 5 — Composite z-score ranking from MD deltas only, then diversity/QC
filtered selection of the top 12 candidates for wet-lab delivery.

Rationale for MD-only composite:
- The 16 candidates fed into Step 4 MD were already evaluated by FoldX + Rosetta
  (they sit at the top of the Step 3 ElasticNet model). Reusing delta_delta or
  dddG_elec in Step 5 double-counts the same signals.
- MD is the orthogonal physics layer: explicit solvent, flexible backbone,
  fixed-His protonation at two pH points. The 4 MD-side interface signals and
  per-CDR RMSF are the new information Step 5 should pivot on.

Signals and signs (positive composite = more pH6-selective dissociation):
    delta_rmsf_h3     (+, w=2.0) CDR3 flexibility increase at pH6 — dominant
                                  because CDR3 is the antigen-contact core.
    delta_n_hbond     (-, w=1.5) interface H-bonds lost at pH6.
    delta_buried_sasa (-, w=1.5) interface opened at pH6.
    delta_n_contacts  (-, w=1.5) interface contacts lost at pH6.
    delta_rmsf_h1     (+, w=0.5) CDR1 flexibility — supporting signal.
    delta_rmsf_h2     (+, w=0.5) CDR2 flexibility — supporting signal.

Not used (kept in output for context):
    model_score (Step 3, anti-correlated with delta_delta), delta_delta (FoldX,
    already used to reach these 16), dddG_elec/ph_score (Rosetta, same),
    delta_rmsd / delta_rmsf_global / delta_rmsd_h* (noisy or redundant).

QC:
    ddG_pH7_4 < 1.5 kcal/mol  — soft cap to avoid destabilizing the pH7.4 bound
                                 state that the assay measures.

Diversity (greedy fill):
    Track A >= 4, Track B >= 4
    D110H-containing >= 3, non-D110H >= 3
    Target 12 candidates.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SRC = REPO_ROOT / "experiments/sdab/R4/structures/all_metrics.csv"
OUT_DIR = REPO_ROOT / "experiments/sdab/R4/final"
OUT_CSV = OUT_DIR / "R4_delivery.csv"

TARGET_N = 10
QC_DDG_PH74_MAX = 1.5

COMPOSITE_WEIGHTS = {
    "delta_rmsf_h3":     ("+",  2.0),
    "delta_n_hbond":     ("-",  1.5),
    "delta_buried_sasa": ("-",  1.5),
    "delta_n_contacts":  ("-",  1.5),
    "delta_rmsf_h1":     ("+",  0.5),
    "delta_rmsf_h2":     ("+",  0.5),
}

CONSTRAINTS = {
    "track_A_min": 4,
    "track_B_min": 4,
    "d110h_yes_min": 3,
    "d110h_no_min": 3,
}

DISPLAY_COLS = [
    "rank", "name", "track", "mutations_unified", "contains_D110H", "confidence",
    "composite_z", "qc_pass", "qc_note",
    "ddG_pH7_4", "delta_delta", "dddG_elec", "ph_score",
    "delta_rmsf_h1", "delta_rmsf_h2", "delta_rmsf_h3",
    "delta_n_hbond", "delta_buried_sasa", "delta_n_contacts",
    "model_score", "log_pH74_pred", "log_pH6_pred",
]


def zscore(series: pd.Series) -> pd.Series:
    mu, sd = series.mean(), series.std(ddof=0)
    if sd == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - mu) / sd


def compute_composite(df: pd.DataFrame) -> pd.DataFrame:
    for col, (sign, _) in COMPOSITE_WEIGHTS.items():
        if col not in df.columns:
            sys.exit(f"[err] missing column: {col}")
        df[f"z_{col}"] = zscore(df[col])
    total = np.zeros(len(df))
    for col, (sign, w) in COMPOSITE_WEIGHTS.items():
        s = 1.0 if sign == "+" else -1.0
        total = total + s * w * df[f"z_{col}"].values
    df["composite_z"] = total
    return df


def apply_qc(df: pd.DataFrame) -> pd.DataFrame:
    notes = []
    passed = []
    for _, row in df.iterrows():
        note_parts = []
        ok = True
        if row["ddG_pH7_4"] >= QC_DDG_PH74_MAX:
            ok = False
            note_parts.append(f"ddG_pH7.4={row['ddG_pH7_4']:.2f}>{QC_DDG_PH74_MAX}")
        if str(row.get("foldx_status", "ok")) != "ok":
            ok = False
            note_parts.append("foldx_fail")
        if str(row.get("rosetta_status", "ok")) != "ok":
            note_parts.append("rosetta_fail")  # warn but don't fail; MD is primary
        passed.append(ok)
        notes.append(";".join(note_parts) if note_parts else "")
    df["qc_pass"] = passed
    df["qc_note"] = notes
    return df


def greedy_select(df: pd.DataFrame, target: int) -> pd.DataFrame:
    cand = df[df["qc_pass"]].copy().sort_values("composite_z", ascending=False)
    picked = []
    counts = {"A": 0, "B": 0, "d110h_yes": 0, "d110h_no": 0}
    min_A = CONSTRAINTS["track_A_min"]
    min_B = CONSTRAINTS["track_B_min"]
    min_yes = CONSTRAINTS["d110h_yes_min"]
    min_no = CONSTRAINTS["d110h_no_min"]

    for _, row in cand.iterrows():
        if len(picked) >= target:
            break
        picked.append(row)
        counts[row["track"]] += 1
        counts["d110h_yes" if row["contains_D110H"] else "d110h_no"] += 1

    picked_df = pd.DataFrame(picked)

    def shortfall(counts):
        return max(0, min_A - counts["A"]) + max(0, min_B - counts["B"]) \
             + max(0, min_yes - counts["d110h_yes"]) + max(0, min_no - counts["d110h_no"])

    # Force-fix quota shortfalls by swapping lowest-z picked items with best
    # remaining candidates that fill the missing bucket.
    if shortfall(counts) > 0:
        remaining = cand[~cand["name"].isin(picked_df["name"])].copy()
        while shortfall(counts) > 0 and len(remaining) > 0:
            need_track = None
            if counts["A"] < min_A:
                need_track = "A"
            elif counts["B"] < min_B:
                need_track = "B"
            need_d110h = None
            if counts["d110h_yes"] < min_yes:
                need_d110h = True
            elif counts["d110h_no"] < min_no:
                need_d110h = False

            def matches(r):
                return (need_track is None or r["track"] == need_track) \
                   and (need_d110h is None or bool(r["contains_D110H"]) == need_d110h)

            cand_fill = remaining[remaining.apply(matches, axis=1)]
            if len(cand_fill) == 0:
                break
            add_row = cand_fill.iloc[0]

            def drop_candidate(df_in, skip_name):
                mask = df_in["name"] != skip_name
                return df_in[mask]

            # Drop lowest-z picked that has the OPPOSITE surplus attribute.
            drop_candidates = picked_df.copy()
            if need_track:
                other = "B" if need_track == "A" else "A"
                drop_candidates = drop_candidates[drop_candidates["track"] == other]
            if need_d110h is not None:
                drop_candidates = drop_candidates[drop_candidates["contains_D110H"] != need_d110h]
            if len(drop_candidates) == 0:
                drop_candidates = picked_df
            drop_row = drop_candidates.sort_values("composite_z", ascending=True).iloc[0]

            picked_df = drop_candidate(picked_df, drop_row["name"])
            picked_df = pd.concat([picked_df, add_row.to_frame().T], ignore_index=True)
            remaining = drop_candidate(remaining, add_row["name"])

            counts[drop_row["track"]] -= 1
            counts["d110h_yes" if drop_row["contains_D110H"] else "d110h_no"] -= 1
            counts[add_row["track"]] += 1
            counts["d110h_yes" if add_row["contains_D110H"] else "d110h_no"] += 1

    picked_df = picked_df.sort_values("composite_z", ascending=False).reset_index(drop=True)
    picked_df["rank"] = np.arange(1, len(picked_df) + 1)
    return picked_df


def main():
    if not SRC.exists():
        sys.exit(f"[err] missing {SRC}")
    df = pd.read_csv(SRC)
    print(f"[load] {len(df)} rows from {SRC}")

    df = compute_composite(df)
    df = apply_qc(df)

    qc_fail = df[~df["qc_pass"]]
    if len(qc_fail) > 0:
        print(f"[qc] {len(qc_fail)} rows failed QC:")
        for _, r in qc_fail.iterrows():
            print(f"  - {r['name']}: {r['qc_note']}")

    picked = greedy_select(df, TARGET_N)
    print(f"\n[select] picked {len(picked)}/{TARGET_N}")
    print(f"  Track A: {(picked['track']=='A').sum()}  "
          f"Track B: {(picked['track']=='B').sum()}")
    print(f"  D110H+: {picked['contains_D110H'].sum()}  "
          f"D110H-: {(~picked['contains_D110H'].astype(bool)).sum()}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    display_cols = [c for c in DISPLAY_COLS if c in picked.columns]
    picked[display_cols].to_csv(OUT_CSV, index=False, float_format="%.4f")
    print(f"[OK] wrote {OUT_CSV}")

    print("\n--- Top 12 ---")
    show = picked[["rank", "name", "track", "mutations_unified", "contains_D110H",
                   "composite_z", "ddG_pH7_4", "delta_delta", "delta_rmsf_h3",
                   "delta_n_hbond", "delta_buried_sasa", "delta_n_contacts"]]
    print(show.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
