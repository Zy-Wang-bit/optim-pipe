#!/usr/bin/env python
from __future__ import annotations

import re
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis.window_selection.common import (
    count_variant_space,
    load_capacity_config,
    load_decision_rules,
    load_inputs,
    make_arg_parser,
    read_table,
    site_option_count,
    split_list,
    write_csv,
)


WINDOW_LENGTH = 40


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def pipe_join(values: list[str]) -> str:
    return "|".join(values)


def parse_position_tokens(text: Any) -> list[tuple[str, int, str]]:
    if text is None or (isinstance(text, float) and pd.isna(text)):
        return []
    tokens: list[tuple[str, int, str]] = []
    for match in re.finditer(r"([A-Za-z])?(\d+)([A-Za-z])?", str(text)):
        chain = match.group(1) or ""
        pos = int(match.group(2))
        suffix = match.group(3) or ""
        tokens.append((chain, pos, suffix))
    return tokens


def constraint_positions(row: pd.Series) -> list[int]:
    positions: list[int] = []
    position = row.get("position")
    if position is not None and not pd.isna(position) and str(position).strip():
        positions.append(int(float(position)))

    tokens = parse_position_tokens(row.get("position_pair"))
    if (
        row.get("constraint_type") == "soft_risk"
        and len(tokens) == 2
        and not tokens[0][2]
        and not tokens[1][2]
        and 0 < abs(tokens[1][1] - tokens[0][1]) <= 5
    ):
        start, end = sorted([tokens[0][1], tokens[1][1]])
        positions.extend(range(start, end + 1))
    else:
        positions.extend(pos for _, pos, _ in tokens)

    return sorted(set(positions))


def constraint_reason(row: pd.Series) -> str:
    pieces = [str(row.get("constraint_id", "")).strip()]
    rationale = str(row.get("rationale", "")).strip()
    if rationale and rationale.lower() != "nan":
        pieces.append(rationale)
    return ": ".join(piece for piece in pieces if piece)


def build_constraint_index(constraints: pd.DataFrame) -> dict[tuple[str, str, int], dict[str, list[str]]]:
    index: dict[tuple[str, str, int], dict[str, list[str]]] = defaultdict(lambda: {"hard": [], "soft": []})
    for _, row in constraints.iterrows():
        if not as_bool(row.get("applies_to_AeS_window_selection", True)):
            continue
        if not as_bool(row.get("applies_to_library_design", True)):
            continue

        allowed_usage = str(row.get("allowed_usage", "")).strip()
        if allowed_usage == "context_only":
            continue

        target = str(row.get("target", "")).strip()
        chain = str(row.get("chain", "")).strip()
        ctype = str(row.get("constraint_type", "")).strip()
        positions = constraint_positions(row)
        reason = constraint_reason(row)

        for pos in positions:
            key = (target, chain, pos)
            if ctype == "hard_protect" and allowed_usage == "hard_filter":
                index[key]["hard"].append(reason)
            elif ctype in {"soft_risk", "negative_pair"} and allowed_usage in {"risk_label", "mechanism_label"}:
                index[key]["soft"].append(reason)
    return index


def allowed_classes_for_status(mask_status: str, capacity_config: dict[str, Any]) -> list[str]:
    defaults = capacity_config["default_allowed_classes"]
    if mask_status == "hard_protected":
        return []
    if mask_status == "soft_risk":
        return list(defaults.get("soft_risk", []))
    return list(defaults.get("designable", []))


def mutation_options(classes: list[str], wt_aa: str, capacity_config: dict[str, Any]) -> set[str]:
    residues: set[str] = set()
    for cls in classes:
        residues.update(capacity_config["mutation_classes"].get(cls, []))
    residues.discard(wt_aa)
    return residues


def parse_forbidden_pair(pair: str) -> tuple[tuple[str, int, str], tuple[str, int, str]] | None:
    tokens = parse_position_tokens(pair)
    mutation_tokens = [(chain, pos, suffix) for chain, pos, suffix in tokens if chain and suffix]
    if len(mutation_tokens) != 2:
        return None
    return mutation_tokens[0], mutation_tokens[1]


def forbidden_pair_variant_count(
    target: str,
    site_options: dict[tuple[str, int], set[str]],
    capacity_config: dict[str, Any],
    max_order: int,
) -> tuple[int, list[str]]:
    if max_order < 2:
        return 0, []

    count = 0
    notes: list[str] = []
    option_counts = {site: len(options) for site, options in site_options.items() if options}

    for constraint in capacity_config.get("forbidden_pair_constraints", []):
        if str(constraint.get("target", "")).strip() != target:
            continue
        parsed = parse_forbidden_pair(str(constraint.get("pair", "")))
        if parsed is None:
            continue
        (chain_a, pos_a, mut_a), (chain_b, pos_b, mut_b) = parsed
        site_a = (chain_a, pos_a)
        site_b = (chain_b, pos_b)
        if mut_a not in site_options.get(site_a, set()) or mut_b not in site_options.get(site_b, set()):
            continue

        other_counts = [
            site_count
            for site, site_count in option_counts.items()
            if site not in {site_a, site_b}
        ]
        pair_count = 0
        for extra_order in range(0, min(max_order - 2, len(other_counts)) + 1):
            for combo in combinations(other_counts, extra_order):
                product = 1
                for site_count in combo:
                    product *= site_count
                pair_count += product

        count += pair_count
        notes.append(f"subtracted_forbidden_pair={constraint.get('pair')}")

    return count, notes


def make_window_id(target: str, segment: str, start: int, end: int) -> str:
    return f"{target}_{segment}_{start:03d}_{end:03d}"


def build_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inputs = load_inputs()
    capacity_config = load_capacity_config()
    decision_rules = load_decision_rules()
    ref = read_table("reference_sequence_map.csv")
    constraints = read_table("prior_constraints_table.csv")

    constraint_index = build_constraint_index(constraints)
    all_classes = list(capacity_config["mutation_classes"].keys())
    max_order = int(capacity_config["max_mutation_order"])
    control_budget = int(capacity_config.get("control_variant_budget", 0))
    min_capacity = int(
        capacity_config.get(
            "min_usable_noncontrol_variant_count",
            decision_rules.get("thresholds", {}).get("min_usable_noncontrol_variant_count", 10000),
        )
    )

    candidate_rows: list[dict[str, Any]] = []
    mask_rows: list[dict[str, Any]] = []
    capacity_rows: list[dict[str, Any]] = []

    for target_cfg in inputs["targets"]:
        target = target_cfg["name"]
        for segment_cfg in target_cfg["antibody_segments"]:
            segment = segment_cfg["segment"]
            chain = segment_cfg["chain"]
            segment_ref = (
                ref[(ref["target"] == target) & (ref["chain"] == chain)]
                .copy()
            )
            segment_ref["local_pos"] = pd.to_numeric(segment_ref["local_pos"], errors="raise")
            segment_ref = segment_ref.sort_values("local_pos").reset_index(drop=True)
            if len(segment_ref) < WINDOW_LENGTH:
                continue

            max_pos = int(segment_ref["local_pos"].max())
            for start_idx in range(0, len(segment_ref) - WINDOW_LENGTH + 1):
                window_ref = segment_ref.iloc[start_idx : start_idx + WINDOW_LENGTH].copy()
                positions = [int(x) for x in window_ref["local_pos"].tolist()]
                start = positions[0]
                end = positions[-1]
                window_id = make_window_id(target, segment, start, end)

                hard_exclusion_reasons: list[str] = []
                if positions != list(range(start, end + 1)):
                    hard_exclusion_reasons.append("non_contiguous_reference_positions")
                if end - start + 1 != WINDOW_LENGTH:
                    hard_exclusion_reasons.append("window_length_not_40")

                near_construct_edge = start <= 3 or end >= max_pos - 2
                soft_flags: list[str] = []
                if near_construct_edge:
                    soft_flags.append("near_construct_edge")

                option_counts: list[int] = []
                site_options: dict[tuple[str, int], set[str]] = {}
                designable_position_count = 0
                hard_protected_count = 0
                soft_risk_count = 0

                for _, residue in window_ref.iterrows():
                    pos = int(residue["local_pos"])
                    aa = str(residue["aa"])
                    region = str(residue.get("region", ""))
                    base_eligible = as_bool(residue.get("is_design_eligible_base", True))
                    hard_reasons = [] if base_eligible else ["reference_position_not_design_eligible_base"]
                    soft_reasons: list[str] = []

                    indexed = constraint_index.get((target, chain, pos), {"hard": [], "soft": []})
                    hard_reasons.extend(indexed["hard"])
                    soft_reasons.extend(indexed["soft"])

                    if hard_reasons:
                        mask_status = "hard_protected"
                        hard_protected_count += 1
                    elif soft_reasons:
                        mask_status = "soft_risk"
                        soft_risk_count += 1
                    else:
                        mask_status = "designable"

                    allowed_classes = allowed_classes_for_status(mask_status, capacity_config)
                    forbidden_classes = [cls for cls in all_classes if cls not in set(allowed_classes)]
                    option_count = site_option_count(pipe_join(allowed_classes), aa, capacity_config)
                    options = mutation_options(allowed_classes, aa, capacity_config)
                    option_counts.append(option_count)
                    site_options[(chain, pos)] = options
                    if option_count > 0:
                        designable_position_count += 1

                    mask_rows.append(
                        {
                            "window_id": window_id,
                            "chain": chain,
                            "pos": pos,
                            "aa": aa,
                            "region": region,
                            "mask_status": mask_status,
                            "allowed_mutation_classes": pipe_join(allowed_classes),
                            "forbidden_mutation_classes": pipe_join(forbidden_classes),
                            "hard_protect_reason": "; ".join(hard_reasons),
                            "soft_risk_reason": "; ".join(soft_reasons),
                        }
                    )

                if hard_protected_count:
                    soft_flags.append(f"contains_hard_protected_positions={hard_protected_count}")
                if soft_risk_count:
                    soft_flags.append(f"contains_soft_risk_positions={soft_risk_count}")

                raw_variant_count = count_variant_space(option_counts, max_order)
                forbidden_count, forbidden_notes = forbidden_pair_variant_count(
                    target=target,
                    site_options=site_options,
                    capacity_config=capacity_config,
                    max_order=max_order,
                )
                raw_variant_count = max(raw_variant_count - forbidden_count, 0)
                usable_count = max(raw_variant_count - control_budget, 0)
                if designable_position_count == 0:
                    capacity_status = "no_designable_positions"
                elif usable_count >= min_capacity:
                    capacity_status = "meets_minimum"
                else:
                    capacity_status = "below_minimum"

                candidate_rows.append(
                    {
                        "target": target,
                        "window_id": window_id,
                        "segment": segment,
                        "start": start,
                        "end": end,
                        "length": WINDOW_LENGTH,
                        "crosses_linker": False,
                        "near_construct_edge": near_construct_edge,
                        "hard_excluded": bool(hard_exclusion_reasons),
                        "hard_exclusion_reasons": "; ".join(hard_exclusion_reasons),
                        "soft_flags": "; ".join(soft_flags),
                    }
                )
                capacity_rows.append(
                    {
                        "window_id": window_id,
                        "designable_position_count": designable_position_count,
                        "raw_feasible_variant_count": raw_variant_count,
                        "control_variant_budget": control_budget,
                        "usable_noncontrol_variant_count": usable_count,
                        "capacity_status": capacity_status,
                        "capacity_notes": "; ".join(forbidden_notes),
                    }
                )

    return pd.DataFrame(candidate_rows), pd.DataFrame(mask_rows), pd.DataFrame(capacity_rows)


def main() -> None:
    make_arg_parser("Enumerate legal 40 aa antibody windows and design capacity tables.").parse_args()
    candidates, mask, capacity = build_outputs()
    write_csv(candidates, "candidate_windows.csv")
    write_csv(mask, "window_mutation_mask.csv")
    write_csv(capacity, "window_design_capacity_table.csv")
    print(
        "Wrote "
        f"{len(candidates)} candidate windows, "
        f"{len(mask)} mask rows, "
        f"{len(capacity)} capacity rows."
    )


if __name__ == "__main__":
    main()
