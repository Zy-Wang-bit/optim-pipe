#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Select 16 candidates (Track A=8, Track B=8) from candidates_combined_100.csv.

Output: experiments/sdab/R4/structures/mutlist.csv
Constraints:
- Track A = 8 by score (allow at most 2 `confidence=low`)
- Track B = 8 by score with quotas:
    - contains_D110H (global) >= 3
    - non-D110H (global)      >= 3
    - CDR1 / CDR2 / CDR3 (global) >= 2 each
If greedy picks fail a quota, swap the lowest-score selection in that track with a qualifying backup.
"""

import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

from analysis.naming.convert import unified_to_foldx

IN_CSV = REPO_ROOT / "experiments/sdab/R4/predictions/candidates_combined_100.csv"
OUT_DIR = REPO_ROOT / "experiments/sdab/R4/structures"
OUT_CSV = OUT_DIR / "mutlist.csv"

TRACK_A_K = 8
TRACK_B_K = 8
TRACK_A_MAX_LOW = 2
MIN_D110H = 3
MIN_NON_D110H = 3
MIN_PER_CDR = 2

POS_RE = re.compile(r"^H[A-Z](\d+)[A-Z]$")


def cdr_of(pos: int) -> str:
    if 26 <= pos <= 33:
        return "CDR1"
    if 51 <= pos <= 58:
        return "CDR2"
    if 97 <= pos <= 111:
        return "CDR3"
    raise ValueError(f"position {pos} outside sdab CDRs")


def cdrs_of_row(mutations_unified: str) -> set[str]:
    cdrs = set()
    for tok in mutations_unified.split(";"):
        m = POS_RE.match(tok.strip())
        if not m:
            raise ValueError(f"unexpected mutation token: {tok}")
        cdrs.add(cdr_of(int(m.group(1))))
    return cdrs


def greedy_track_a(df_a: pd.DataFrame) -> pd.DataFrame:
    picks = []
    n_low = 0
    for _, row in df_a.sort_values("score", ascending=False).iterrows():
        if len(picks) == TRACK_A_K:
            break
        if row["confidence"] == "low" and n_low >= TRACK_A_MAX_LOW:
            continue
        picks.append(row)
        if row["confidence"] == "low":
            n_low += 1
    if len(picks) < TRACK_A_K:
        raise RuntimeError(f"Track A yielded only {len(picks)} picks (need {TRACK_A_K})")
    return pd.DataFrame(picks)


def greedy_track_b(df_b: pd.DataFrame, existing: pd.DataFrame) -> pd.DataFrame:
    sorted_b = df_b.sort_values("score", ascending=False).reset_index(drop=True)
    picks = sorted_b.head(TRACK_B_K).copy()

    combined = pd.concat([existing, picks], ignore_index=True)
    n_d110h = int(combined["contains_D110H"].sum())
    n_non_d110h = int((~combined["contains_D110H"]).sum())
    cdr_counts = _cdr_counts(combined)

    backups = sorted_b.iloc[TRACK_B_K:].copy()

    def ensure(condition_fn, need_fn, description):
        """Swap lowest-score picks in Track B for backups that satisfy `need_fn`
        until the deficit is closed."""
        nonlocal picks
        while not condition_fn():
            candidates = backups[backups.apply(need_fn, axis=1)]
            if candidates.empty:
                raise RuntimeError(f"Track B quota unmet: {description}; backups exhausted")
            new_row = candidates.iloc[0]
            drop_idx = picks[picks.apply(lambda r: not need_fn(r), axis=1)]
            if drop_idx.empty:
                raise RuntimeError(f"Track B: nothing to drop for quota {description}")
            drop_idx = drop_idx.sort_values("score").head(1).index
            picks = pd.concat(
                [picks.drop(drop_idx), new_row.to_frame().T],
                ignore_index=True,
            )
            backups = backups[backups["mutations_unified"] != new_row["mutations_unified"]]

    def check_d110h():
        total_d = int(existing["contains_D110H"].sum() + picks["contains_D110H"].sum())
        return total_d >= MIN_D110H

    def check_non_d110h():
        total_nd = int((~existing["contains_D110H"]).sum() + (~picks["contains_D110H"]).sum())
        return total_nd >= MIN_NON_D110H

    def need_d110h(r):
        return bool(r["contains_D110H"])

    def need_non_d110h(r):
        return not bool(r["contains_D110H"])

    ensure(check_d110h, need_d110h, "contains_D110H >= 3")
    ensure(check_non_d110h, need_non_d110h, "non-D110H >= 3")

    for cdr in ("CDR1", "CDR2", "CDR3"):
        def check_cdr(cdr=cdr):
            combined2 = pd.concat([existing, picks], ignore_index=True)
            return _cdr_counts(combined2).get(cdr, 0) >= MIN_PER_CDR

        def need_cdr(r, cdr=cdr):
            return cdr in cdrs_of_row(r["mutations_unified"])

        ensure(check_cdr, need_cdr, f"{cdr} coverage >= {MIN_PER_CDR}")

    if len(picks) != TRACK_B_K:
        raise RuntimeError(f"Track B ended with {len(picks)} picks (need {TRACK_B_K})")
    return picks.reset_index(drop=True)


def _cdr_counts(df: pd.DataFrame) -> dict[str, int]:
    counts = {"CDR1": 0, "CDR2": 0, "CDR3": 0}
    for _, row in df.iterrows():
        for cdr in cdrs_of_row(row["mutations_unified"]):
            counts[cdr] += 1
    return counts


def to_foldx(mutations_unified: str) -> str:
    return ",".join(unified_to_foldx(tok.strip()) for tok in mutations_unified.split(";"))


def main():
    df = pd.read_csv(IN_CSV)
    if list(df["track"].unique()) != ["A", "B"] and set(df["track"].unique()) != {"A", "B"}:
        raise RuntimeError(f"unexpected track values: {df['track'].unique()}")
    df["contains_D110H"] = df["contains_D110H"].astype(bool)

    df_a = df[df["track"] == "A"].copy()
    df_b = df[df["track"] == "B"].copy()

    picks_a = greedy_track_a(df_a)
    picks_b = greedy_track_b(df_b, picks_a)

    picks = pd.concat([picks_a, picks_b], ignore_index=True)
    picks = picks.sort_values(["track", "score"], ascending=[True, False]).reset_index(drop=True)

    picks["name"] = [f"r4_{i+1:02d}" for i in range(len(picks))]
    picks["mutations_foldx"] = picks["mutations_unified"].apply(to_foldx)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = picks[[
        "name", "track", "mutations_unified", "mutations_foldx", "order",
        "contains_D110H", "confidence", "score", "log_pH74_pred", "log_pH6_pred",
    ]]
    out.to_csv(OUT_CSV, index=False)

    # Summary
    print(f"[OK] wrote {OUT_CSV} ({len(out)} rows)")
    print(f"  Track A: {(out['track']=='A').sum()}, Track B: {(out['track']=='B').sum()}")
    print(f"  contains_D110H: {int(out['contains_D110H'].sum())}")
    print(f"  non-D110H: {int((~out['contains_D110H']).sum())}")
    cdr_counts = _cdr_counts(out)
    print(f"  CDR coverage: {cdr_counts}")
    print(f"  confidence: {out['confidence'].value_counts().to_dict()}")
    print(f"  order distribution: {out['order'].value_counts().sort_index().to_dict()}")
    print()
    print(out[["name", "track", "mutations_unified", "mutations_foldx", "score"]].to_string(index=False))


if __name__ == "__main__":
    main()
