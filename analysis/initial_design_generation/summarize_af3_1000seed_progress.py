#!/usr/bin/env python3
"""Summarize AF3 1000-seed validation progress by entry.

The user-facing unit for this stage is one AF3 entry with 1000 seeds,
not the internal 50-seed execution batch.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import pandas as pd


DEFAULT_ROOT = Path("results/initial_design_generation/af3_1000seed_validation")
SEEDS_PER_BATCH = 50
TARGET_SEEDS_PER_ENTRY = 1000
SAMPLES_PER_SEED = 5
SEED_DIR_RE = re.compile(r"seed-(\d+)_sample-\d+")


def load_expected_entries(root: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []

    panel = root / "validation_panel_64.csv"
    if panel.exists():
        df = pd.read_csv(panel)
        if "variant_id" not in df.columns:
            raise ValueError(f"{panel} is missing variant_id")
        if "target" not in df.columns:
            df["target"] = "unknown"
        df = df.copy()
        df["entry_type"] = "mutant"
        frames.append(df[["variant_id", "target", "entry_type"]])

    parents = root / "validation_panel_parent_baselines.csv"
    if parents.exists():
        df = pd.read_csv(parents)
        if "variant_id" not in df.columns:
            raise ValueError(f"{parents} is missing variant_id")
        if "target" not in df.columns:
            df["target"] = "parent_baseline"
        df = df.copy()
        df["entry_type"] = "parent_baseline"
        frames.append(df[["variant_id", "target", "entry_type"]])

    if not frames:
        return pd.DataFrame(columns=["variant_id", "target", "entry_type"])

    expected = pd.concat(frames, ignore_index=True)
    expected["variant_id"] = expected["variant_id"].astype(str)
    expected["target"] = expected["target"].astype(str)
    return expected.drop_duplicates(subset=["variant_id"], keep="first")


def load_success_seed_counts(root: Path) -> pd.DataFrame:
    successes: set[tuple[str, str]] = set()

    for status_file in sorted(root.glob("af3_inference*.csv")):
        if not status_file.exists() or status_file.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(status_file)
        except Exception:
            continue
        required = {"variant_id", "shard_id", "status"}
        if not required.issubset(df.columns):
            continue
        for _, row in df[df["status"].eq("success")].iterrows():
            successes.add((str(row["variant_id"]), str(row["shard_id"])))

    counts: dict[str, int] = {}
    for variant_id, _batch_id in successes:
        counts[variant_id] = counts.get(variant_id, 0) + SEEDS_PER_BATCH

    file_counts = load_file_verified_seed_counts(root)
    for variant_id, file_verified_seeds in file_counts.items():
        counts[variant_id] = max(counts.get(variant_id, 0), file_verified_seeds)

    return pd.DataFrame(
        [{"variant_id": variant_id, "completed_seeds": seeds} for variant_id, seeds in counts.items()]
    )


def count_unique_seeds_with_summary(output_dir: Path) -> int:
    """Count unique random seeds with at least one summary file.

    Some resumed AF3 jobs are run directly outside the queue controller, so
    status CSVs can undercount successful shards. File-level verification is
    the source of truth for completed seed outputs. Use os.walk(followlinks=True)
    because the remaining sdAb outputs are symlinked into /home/ziyang.
    """
    if not output_dir.exists():
        return 0
    seeds: set[str] = set()
    for dirpath, _dirnames, filenames in os.walk(output_dir, followlinks=True):
        if "summary_confidences.json" not in filenames:
            continue
        match = SEED_DIR_RE.search(dirpath)
        if match:
            seeds.add(match.group(1))
    return len(seeds)


def load_file_verified_seed_counts(root: Path) -> dict[str, int]:
    shard_path = root / "af3_shard_schedule.csv"
    if not shard_path.exists() or shard_path.stat().st_size == 0:
        return {}

    try:
        shard = pd.read_csv(shard_path, usecols=["variant_id", "output_dir"], low_memory=False)
    except Exception:
        return {}

    counts: dict[str, int] = {}
    for _, row in shard.iterrows():
        variant_id = str(row["variant_id"])
        output_dir = root.parent.parent.parent / str(row["output_dir"])
        unique_seed_count = count_unique_seeds_with_summary(output_dir)
        if unique_seed_count >= SEEDS_PER_BATCH:
            counts[variant_id] = counts.get(variant_id, 0) + SEEDS_PER_BATCH
    return counts


def build_summary(root: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    expected = load_expected_entries(root)
    seed_counts = load_success_seed_counts(root)

    if expected.empty:
        progress = seed_counts.copy()
        progress["target"] = "unknown"
        progress["entry_type"] = "unknown"
    else:
        progress = expected.merge(seed_counts, on="variant_id", how="left")

    if "completed_seeds" not in progress.columns:
        progress["completed_seeds"] = 0
    progress["completed_seeds"] = progress["completed_seeds"].fillna(0).astype(int)
    progress["target_seeds"] = TARGET_SEEDS_PER_ENTRY
    progress["completion_fraction"] = progress["completed_seeds"] / TARGET_SEEDS_PER_ENTRY
    progress["status"] = "not_started"
    progress.loc[progress["completed_seeds"].between(1, TARGET_SEEDS_PER_ENTRY - 1), "status"] = "partial"
    progress.loc[progress["completed_seeds"].ge(TARGET_SEEDS_PER_ENTRY), "status"] = "complete"

    totals = {
        "expected_entries": int(len(progress)),
        "complete_entries": int(progress["status"].eq("complete").sum()),
        "partial_entries": int(progress["status"].eq("partial").sum()),
        "not_started_entries": int(progress["status"].eq("not_started").sum()),
        "completed_seeds": int(progress["completed_seeds"].sum()),
        "target_seeds": int(len(progress) * TARGET_SEEDS_PER_ENTRY),
    }
    return progress.sort_values(["status", "target", "variant_id"]), totals


def print_summary(progress: pd.DataFrame, totals: dict[str, int]) -> None:
    print(
        "AF3 1000-seed progress: "
        f"{totals['complete_entries']} / {totals['expected_entries']} entries complete; "
        f"{totals['completed_seeds']} / {totals['target_seeds']} seeds complete"
    )
    print(
        f"partial entries: {totals['partial_entries']}; "
        f"not started: {totals['not_started_entries']}"
    )

    if not progress.empty and {"target", "entry_type"}.issubset(progress.columns):
        print("\nBy target:")
        by_target = (
            progress.groupby(["target", "entry_type"], dropna=False)
            .agg(
                entries=("variant_id", "count"),
                complete=("status", lambda x: int((x == "complete").sum())),
                partial=("status", lambda x: int((x == "partial").sum())),
                not_started=("status", lambda x: int((x == "not_started").sum())),
                completed_seeds=("completed_seeds", "sum"),
            )
            .reset_index()
        )
        for _, row in by_target.iterrows():
            print(
                f"- {row['target']} / {row['entry_type']}: "
                f"{row['complete']}/{row['entries']} complete, "
                f"{row['partial']} partial, {row['not_started']} not started, "
                f"{row['completed_seeds']} seeds"
            )

    partial = progress[progress["status"].eq("partial")].sort_values(
        ["completed_seeds", "variant_id"], ascending=[False, True]
    )
    if not partial.empty:
        print("\nPartial entries:")
        for _, row in partial.iterrows():
            print(f"- {row['variant_id']}: {row['completed_seeds']} / {TARGET_SEEDS_PER_ENTRY} seeds")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-csv", type=Path, default=None)
    args = parser.parse_args()

    progress, totals = build_summary(args.root)
    print_summary(progress, totals)

    if args.output_csv is not None:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        progress.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
