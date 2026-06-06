#!/usr/bin/env python
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import pandas as pd

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis.window_selection.common import (
    load_capacity_config,
    make_arg_parser,
    read_table,
    tables_dir,
    write_csv,
)


DISTANCE_COLUMNS = [
    "distance_to_glycan_envelope",
    "n146_glycan_min_distance",
    "min_distance_to_n146_glycan",
    "distance_to_n146_glycan",
    "n146_min_distance",
    "distance_to_n146",
]

PARENT_OVERLAP_COLUMNS = [
    "parent_overlap_status",
    "n146_parent_overlap_status",
    "glycan_parent_overlap_status",
]


def split_classes(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return [item for item in str(value).split("|") if item]


def pipe_join(values: list[str]) -> str:
    return "|".join(values)


def as_float(value: Any) -> float | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(out):
        return None
    return out


def first_present(row: pd.Series, columns: list[str]) -> Any:
    for column in columns:
        if column in row.index and not pd.isna(row[column]) and str(row[column]).strip():
            return row[column]
    return None


def normalize_feature_columns(features: pd.DataFrame) -> pd.DataFrame:
    features = features.copy()
    rename: dict[str, str] = {}
    if "pos" not in features.columns:
        for candidate in ["local_pos", "antibody_pos", "residue_pos", "position"]:
            if candidate in features.columns:
                rename[candidate] = "pos"
                break
    if "chain" not in features.columns:
        for candidate in ["antibody_chain", "chain_id"]:
            if candidate in features.columns:
                rename[candidate] = "chain"
                break
    if rename:
        features = features.rename(columns=rename)
    if "pos" in features.columns:
        features["pos"] = pd.to_numeric(features["pos"], errors="coerce").astype("Int64")
    return features


def load_structure_features(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    features = pd.read_csv(path, dtype=str, keep_default_na=False)
    if features.empty:
        return None
    return normalize_feature_columns(features)


def first_nonempty(series: pd.Series) -> Any:
    for value in series:
        if value is not None and not pd.isna(value) and str(value).strip():
            return value
    return ""


def numeric_min(series: pd.Series) -> float:
    return pd.to_numeric(series, errors="coerce").min()


def collapse_duplicate_features(features: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    if not features.duplicated(keys).any():
        return features

    aggregations = {}
    for column in features.columns:
        if column in keys:
            continue
        aggregations[column] = numeric_min if column in DISTANCE_COLUMNS else first_nonempty
    return features.groupby(keys, as_index=False).agg(aggregations)


def merge_features(base: pd.DataFrame, features: pd.DataFrame | None) -> pd.DataFrame:
    if features is None:
        return base

    merge_keys = ["chain", "pos"]
    if "target" in features.columns and "target" in base.columns:
        merge_keys = ["target", *merge_keys]
    if "window_id" in features.columns:
        merge_keys = ["window_id", *[key for key in merge_keys if key != "target"]]

    available_keys = [key for key in merge_keys if key in base.columns and key in features.columns]
    if len(available_keys) < 2:
        return base
    if "pos" in available_keys:
        base = base.copy()
        features = features.copy()
        base["pos"] = pd.to_numeric(base["pos"], errors="coerce").astype("Int64")
        features["pos"] = pd.to_numeric(features["pos"], errors="coerce").astype("Int64")
    features = collapse_duplicate_features(features, available_keys)
    return base.merge(features, on=available_keys, how="left", suffixes=("", "_structure"))


def compute_row_guardrail(
    row: pd.Series,
    *,
    mode: str,
    radius: float,
    class_order: list[str],
    bulky_or_positive_classes: set[str],
    no_structure_features: bool,
) -> dict[str, Any]:
    allowed_classes = split_classes(row.get("allowed_mutation_classes"))
    forbidden_classes = split_classes(row.get("forbidden_mutation_classes"))
    allowed_set = set(allowed_classes)
    forbidden_set = set(forbidden_classes)

    if no_structure_features:
        return {
            "window_id": row["window_id"],
            "chain": row["chain"],
            "pos": int(row["pos"]),
            "aa": row["aa"],
            "glycan_guardrail_mode": "heuristic",
            "distance_to_glycan_envelope": "",
            "parent_overlap_status": "not_evaluated_no_structure_features",
            "mutation_expansion_risk": "low",
            "allowed_mutation_classes": pipe_join(allowed_classes),
            "forbidden_mutation_classes": pipe_join(forbidden_classes),
            "control_needed": False,
            "guardrail_reason": f"heuristic_dry_run_no_structure_features; configured_mode={mode}",
        }

    distance = as_float(first_present(row, DISTANCE_COLUMNS))
    parent_overlap_status = first_present(row, PARENT_OVERLAP_COLUMNS)
    close_to_envelope = distance is not None and distance <= radius

    if parent_overlap_status is None:
        if distance is None:
            parent_overlap_status = "not_evaluated_no_distance"
        elif close_to_envelope:
            parent_overlap_status = "near_n146_envelope_proxy"
        else:
            parent_overlap_status = "outside_n146_envelope_proxy"

    control_needed = False
    reason = "outside_configured_glycan_radius"
    mutation_expansion_risk = "low"

    if row.get("mask_status") == "hard_protected":
        reason = "hard_protected_position_no_mutation_expansion"
    elif distance is None:
        reason = "no_glycan_distance_available"
    elif close_to_envelope and allowed_set.intersection(bulky_or_positive_classes):
        mutation_expansion_risk = "high"
        forbidden_set.update(bulky_or_positive_classes)
        allowed_set.difference_update(bulky_or_positive_classes)
        control_needed = True
        reason = "within_configured_glycan_radius_forbid_bulky_or_positive_classes"
    elif close_to_envelope:
        mutation_expansion_risk = "medium"
        control_needed = True
        reason = "within_configured_glycan_radius_no_extra_expansion_classes"

    allowed_out = [cls for cls in class_order if cls in allowed_set]
    forbidden_out = [cls for cls in class_order if cls in forbidden_set]

    return {
        "window_id": row["window_id"],
        "chain": row["chain"],
        "pos": int(row["pos"]),
        "aa": row["aa"],
        "glycan_guardrail_mode": mode,
        "distance_to_glycan_envelope": "" if distance is None else distance,
        "parent_overlap_status": parent_overlap_status,
        "mutation_expansion_risk": mutation_expansion_risk,
        "allowed_mutation_classes": pipe_join(allowed_out),
        "forbidden_mutation_classes": pipe_join(forbidden_out),
        "control_needed": control_needed,
        "guardrail_reason": reason,
    }


def build_guardrails() -> pd.DataFrame:
    capacity_config = load_capacity_config()
    mode = str(capacity_config.get("glycan_envelope_mode", "heuristic")).strip() or "heuristic"
    guardrail_config = capacity_config.get("glycan_guardrails", {})
    radius = float(guardrail_config.get("heuristic_radius_angstrom", 8.0))
    bulky_or_positive_classes = set(guardrail_config.get("bulky_or_positive_classes", []))
    class_order = list(capacity_config["mutation_classes"].keys())

    candidates = read_table("candidate_windows.csv")
    mask = read_table("window_mutation_mask.csv")
    base = mask.merge(candidates[["window_id", "target"]], on="window_id", how="left")

    features = load_structure_features(tables_dir() / "residue_structure_features.csv")
    no_structure_features = features is None
    merged = merge_features(base, features)

    rows = [
        compute_row_guardrail(
            row,
            mode=mode,
            radius=radius,
            class_order=class_order,
            bulky_or_positive_classes=bulky_or_positive_classes,
            no_structure_features=no_structure_features,
        )
        for _, row in merged.iterrows()
    ]
    return pd.DataFrame(rows)


def main() -> None:
    make_arg_parser("Build N146 glycan guardrail table for candidate windows.").parse_args()
    guardrails = build_guardrails()
    write_csv(guardrails, "glycan_guardrail_table.csv")
    print(f"Wrote {len(guardrails)} glycan guardrail rows.")


if __name__ == "__main__":
    main()
