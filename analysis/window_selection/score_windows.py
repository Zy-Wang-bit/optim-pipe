#!/usr/bin/env python
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.window_selection.common import (
    load_capacity_config,
    load_decision_rules,
    load_inputs,
    load_source_policy,
    parse_constraint_positions,
    repo_path,
)


OUTPUT_COLUMNS = [
    "target",
    "window_id",
    "segment",
    "start",
    "end",
    "grade",
    "recommendation_status",
    "blocker_type",
    "blocker_reason",
    "decision_path",
    "evidence_route",
    "pH_mechanism_hypothesis_strength",
    "triggered_rules",
    "failed_rules",
    "soft_downgrade_rules",
    "hard_filter_reasons",
    "wet_prior_summary",
    "structure_summary",
    "glycan_risk_summary",
    "raw_feasible_variant_count",
    "usable_noncontrol_variant_count",
]


@dataclass
class TableBundle:
    candidate_windows: pd.DataFrame
    mask: pd.DataFrame
    capacity: pd.DataFrame
    prior_constraints: pd.DataFrame
    af3_manifest: pd.DataFrame
    residue_features: pd.DataFrame
    glycan_guardrails: pd.DataFrame
    wet_tables: dict[str, pd.DataFrame]
    missing_required_tables: list[str]
    table_errors: list[str]


def _tables_dir(output_root: str | Path) -> Path:
    return repo_path(output_root) / "tables"


def _read_optional_table(tables_dir: Path, filename: str) -> tuple[pd.DataFrame, str | None]:
    path = tables_dir / filename
    if not path.exists():
        return pd.DataFrame(), f"missing:{filename}"
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False), None
    except Exception as exc:  # noqa: BLE001 - validation report needs the parse reason.
        return pd.DataFrame(), f"unparseable:{filename}:{exc}"


def _truthy(value: Any) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _split_pipe(value: Any) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [part.strip() for part in re.split(r"[|;,]", str(value)) if part.strip()]


def _join(items: list[str]) -> str:
    return "|".join(dict.fromkeys(str(item) for item in items if str(item)))


def _numeric(value: Any, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _target_names() -> list[str]:
    return [target["name"] for target in load_inputs().get("targets", [])]


def load_tables(output_root: str | Path) -> TableBundle:
    td = _tables_dir(output_root)
    required = {
        "candidate_windows": "candidate_windows.csv",
        "mask": "window_mutation_mask.csv",
        "capacity": "window_design_capacity_table.csv",
        "prior_constraints": "prior_constraints_table.csv",
        "af3_manifest": "af3_model_manifest.csv",
        "residue_features": "residue_structure_features.csv",
        "glycan_guardrails": "glycan_guardrail_table.csv",
    }
    loaded: dict[str, pd.DataFrame] = {}
    missing_required: list[str] = []
    errors: list[str] = []
    for name, filename in required.items():
        df, error = _read_optional_table(td, filename)
        loaded[name] = df
        if error:
            if name in {"candidate_windows", "mask", "capacity"}:
                missing_required.append(filename)
            errors.append(error)

    wet_tables: dict[str, pd.DataFrame] = {}
    for target, filename in {"1E62": "wet_observation_table_1e62.csv", "sdAb": "wet_observation_table_sdab.csv"}.items():
        df, error = _read_optional_table(td, filename)
        wet_tables[target] = df
        if error:
            errors.append(error)

    return TableBundle(
        candidate_windows=loaded["candidate_windows"],
        mask=loaded["mask"],
        capacity=loaded["capacity"],
        prior_constraints=loaded["prior_constraints"],
        af3_manifest=loaded["af3_manifest"],
        residue_features=loaded["residue_features"],
        glycan_guardrails=loaded["glycan_guardrails"],
        wet_tables=wet_tables,
        missing_required_tables=missing_required,
        table_errors=errors,
    )


def wet_input_errors() -> dict[str, list[str]]:
    policy = load_source_policy()
    out: dict[str, list[str]] = {}
    for target, files in policy.get("scoring_allowed_files", {}).items():
        problems = []
        for file_name in files:
            path = repo_path(file_name)
            if not path.exists():
                problems.append(f"missing:{file_name}")
                continue
            try:
                pd.read_csv(path, nrows=1)
            except Exception as exc:  # noqa: BLE001
                problems.append(f"unparseable:{file_name}:{exc}")
        out[target] = problems
    return out


def source_policy_violations(bundle: TableBundle) -> list[str]:
    policy = load_source_policy()
    excluded = [str(pat).lower() for pat in policy.get("excluded_patterns", [])]
    violations: list[str] = []

    for target, files in policy.get("scoring_allowed_files", {}).items():
        for file_name in files:
            lower = str(file_name).lower()
            if any(pattern in lower for pattern in excluded):
                violations.append(f"{target}:scoring_allowed_file_matches_excluded_pattern:{file_name}")

    for target, wet in bundle.wet_tables.items():
        if wet.empty or "source_file" not in wet.columns:
            continue
        for source_file in wet["source_file"].dropna().astype(str).unique():
            lower = source_file.lower()
            if any(pattern in lower for pattern in excluded):
                violations.append(f"{target}:wet_source_matches_excluded_pattern:{source_file}")

    constraints = bundle.prior_constraints
    if not constraints.empty:
        hard = constraints[constraints.get("allowed_usage", "").astype(str).eq("hard_filter")]
        if not hard.empty:
            for _, row in hard.iterrows():
                antigen = str(row.get("source_antigen", "")).strip()
                applies = _truthy(row.get("applies_to_AeS_window_selection"))
                if antigen and antigen != "AeS" and not applies:
                    violations.append(
                        f"{row.get('target', '')}:non_AeS_hard_filter_without_AeS_apply:{row.get('constraint_id', '')}"
                    )

    return violations


def qualified_af3_counts(af3_manifest: pd.DataFrame) -> dict[str, int]:
    if af3_manifest.empty or "target" not in af3_manifest.columns:
        return {target: 0 for target in _target_names()}

    df = af3_manifest.copy()
    quality = df.get("quality_status", pd.Series([""] * len(df))).astype(str).str.lower()
    blocker = df.get("blocker_flag", pd.Series([False] * len(df))).map(_truthy)
    passed = df[quality.eq("pass") & ~blocker]
    return {target: int((passed["target"] == target).sum()) for target in _target_names()}


def target_capacity_blockers(candidates: pd.DataFrame, capacity: pd.DataFrame, min_capacity: int) -> dict[str, str]:
    blockers: dict[str, str] = {}
    if candidates.empty or capacity.empty or "window_id" not in candidates.columns or "window_id" not in capacity.columns:
        return blockers

    merged = candidates.merge(capacity, on="window_id", how="left", suffixes=("", "_capacity"))
    for target, group in merged.groupby("target", dropna=False):
        legal = group[
            (group.get("length", 0).map(_numeric).eq(40))
            & ~group.get("hard_excluded", False).map(_truthy)
            & ~group.get("crosses_linker", False).map(_truthy)
        ]
        if legal.empty:
            blockers[str(target)] = "no_legal_windows_available"
            continue
        usable = legal.get("usable_noncontrol_variant_count", pd.Series([0] * len(legal))).map(_numeric)
        if not usable.ge(min_capacity).any():
            blockers[str(target)] = f"no_legal_window_has_usable_noncontrol_variant_count_ge_{min_capacity}"
    return blockers


def collect_target_blockers(bundle: TableBundle) -> dict[str, tuple[str, str]]:
    rules = load_decision_rules()
    capacity_cfg = load_capacity_config()
    thresholds = rules.get("thresholds", {})
    min_structures = int(thresholds.get("min_qualified_structures_for_main_scoring", 15))
    min_capacity = int(capacity_cfg.get("min_usable_noncontrol_variant_count", thresholds.get("min_usable_noncontrol_variant_count", 10000)))

    wet_errors = wet_input_errors()
    af3_counts = qualified_af3_counts(bundle.af3_manifest)
    capacity_blockers = target_capacity_blockers(bundle.candidate_windows, bundle.capacity, min_capacity)
    policy_violations = source_policy_violations(bundle)

    blockers: dict[str, tuple[str, str]] = {}
    for target in _target_names():
        reasons: list[str] = []
        if bundle.missing_required_tables:
            blockers[target] = ("schema_error", f"required scoring table missing or unparseable: {_join(bundle.missing_required_tables)}")
            continue
        target_policy = [v for v in policy_violations if v.startswith(f"{target}:")]
        if target_policy:
            blockers[target] = ("source_policy_violation", _join(target_policy))
            continue
        if wet_errors.get(target):
            blockers[target] = ("missing_wet_data", _join(wet_errors[target]))
            continue
        if af3_counts.get(target, 0) < min_structures:
            reasons.append(f"qualified_AF3_models={af3_counts.get(target, 0)} < {min_structures}")
            blockers[target] = ("insufficient_AF3_models", _join(reasons))
            continue
        if target in capacity_blockers:
            blockers[target] = ("capacity_failure", capacity_blockers[target])
    return blockers


def _positions_for_window(window_id: str, mask: pd.DataFrame) -> pd.DataFrame:
    if mask.empty or "window_id" not in mask.columns:
        return pd.DataFrame()
    return mask[mask["window_id"].astype(str).eq(str(window_id))].copy()


def _capacity_for_window(window_id: str, capacity: pd.DataFrame) -> dict[str, Any]:
    if capacity.empty or "window_id" not in capacity.columns:
        return {}
    hit = capacity[capacity["window_id"].astype(str).eq(str(window_id))]
    if hit.empty:
        return {}
    return hit.iloc[0].to_dict()


def _prior_constraints_for_window(row: pd.Series, constraints: pd.DataFrame, mask_rows: pd.DataFrame) -> tuple[list[str], list[str]]:
    if constraints.empty:
        return [], []
    window_positions = set(mask_rows.get("pos", pd.Series(dtype=int)).map(lambda value: int(_numeric(value, -1))).tolist())
    window_chains = set(mask_rows.get("chain", pd.Series(dtype=str)).dropna().astype(str).tolist())
    hard: list[str] = []
    soft: list[str] = []
    for _, constraint in constraints.iterrows():
        if str(constraint.get("target", "")) != str(row.get("target", "")):
            continue
        allowed_usage = str(constraint.get("allowed_usage", ""))
        if allowed_usage == "context_only":
            continue
        chain = str(constraint.get("chain", "")).strip()
        if chain and window_chains and chain not in window_chains:
            continue
        positions = set(parse_constraint_positions(constraint))
        if positions and window_positions and not positions.intersection(window_positions):
            continue
        cid = str(constraint.get("constraint_id", constraint.get("constraint_type", "")))
        if allowed_usage == "hard_filter" and _truthy(constraint.get("applies_to_AeS_window_selection")):
            hard.append(cid)
        elif allowed_usage in {"risk_label", "mechanism_label"}:
            soft.append(cid)
    return hard, soft


def _wet_prior_summary(row: pd.Series, wet_table: pd.DataFrame, mask_rows: pd.DataFrame) -> str:
    if wet_table.empty:
        return "unknown_or_neutral:no_wet_observation_table"
    if "positions" not in wet_table.columns:
        return "unknown_or_neutral:wet_table_has_no_positions"
    positions = set(mask_rows.get("pos", pd.Series(dtype=int)).map(lambda value: int(_numeric(value, -1))).tolist())
    if not positions:
        return "unknown_or_neutral:no_window_positions"
    overlaps = 0
    for _, wet in wet_table.iterrows():
        wet_positions = {int(_numeric(pos, -1)) for pos in _split_pipe(wet.get("positions"))}
        if positions.intersection(wet_positions):
            overlaps += 1
    if overlaps == 0:
        return "unknown_or_neutral:no_direct_window_prior"
    return f"unknown_or_neutral:{overlaps}_direct_or_neighbor_wet_observations_not_directionally_scored"


def _structure_summary(row: pd.Series, features: pd.DataFrame, mask_rows: pd.DataFrame) -> tuple[bool, str]:
    if features.empty:
        return False, "insufficient_evidence:no_residue_structure_features"

    window_id = str(row.get("window_id", ""))
    if "window_id" in features.columns:
        subset = features[features["window_id"].astype(str).eq(window_id)].copy()
    else:
        positions = set(mask_rows.get("pos", pd.Series(dtype=int)).map(lambda value: int(_numeric(value, -1))).tolist())
        chains = set(mask_rows.get("chain", pd.Series(dtype=str)).dropna().astype(str).tolist())
        subset = features.copy()
        if "target" in subset.columns:
            subset = subset[subset["target"].astype(str).eq(str(row.get("target", "")))]
        if "pos" in subset.columns:
            subset = subset[subset["pos"].map(lambda value: int(_numeric(value, -1))).isin(positions)]
        if "chain" in subset.columns and chains:
            subset = subset[subset["chain"].astype(str).isin(chains)]

    if subset.empty:
        return False, "insufficient_evidence:no_features_for_window"

    rules = load_decision_rules()
    thresholds = rules.get("thresholds", {})
    min_seed = int(thresholds.get("min_seed_support_for_stable_contact", 3))
    min_cluster = int(thresholds.get("min_cluster_support_for_stable_contact", 1))

    seed_count = 0
    for col in ["dominant_cluster_seed_count", "seed_support_count", "contact_seed_count"]:
        if col in subset.columns:
            seed_count = max(seed_count, int(subset[col].map(_numeric).max()))
    cluster_support = 0
    for col in ["cluster_support_count", "dominant_cluster_support_count"]:
        if col in subset.columns:
            cluster_support = max(cluster_support, int(subset[col].map(_numeric).max()))
    if cluster_support == 0 and "cluster_id" in subset.columns:
        cluster_support = int(subset["cluster_id"].dropna().nunique())

    ag_loop_fraction = 0.0
    for col in ["dominant_cluster_ag_loop_contact_fraction", "ag_loop_contact_fraction", "antigenic_loop_contact_fraction"]:
        if col in subset.columns:
            ag_loop_fraction = max(ag_loop_fraction, float(subset[col].map(_numeric).max()))
    a_det_fraction = 0.0
    if "a_determinant_contact_fraction" in subset.columns:
        a_det_fraction = float(subset["a_determinant_contact_fraction"].map(_numeric).max())

    anomalous = False
    for col in ["dominant_cluster_is_anomalous", "dominant_cluster_tm_interface_flag", "tm_interface_flag"]:
        if col in subset.columns and subset[col].map(_truthy).any():
            anomalous = True

    quality_ok = True
    for col in ["dominant_cluster_quality_status", "model_quality_status", "quality_status"]:
        if col in subset.columns:
            values = set(subset[col].dropna().astype(str).str.lower())
            if values and not values.intersection({"pass", "ok", "interpretable"}):
                quality_ok = False

    stable = seed_count >= min_seed and cluster_support >= min_cluster and ag_loop_fraction > 0 and not anomalous and quality_ok
    status = "stable_contact" if stable else "insufficient_or_unstable_contact"
    summary = (
        f"{status}:seed_support={seed_count};cluster_support={cluster_support};"
        f"ag_loop_fraction={ag_loop_fraction:.3g};a_determinant_fraction={a_det_fraction:.3g};"
        f"quality_ok={quality_ok};anomalous_or_tm={anomalous}"
    )
    return stable, summary


def _glycan_summary(window_id: str, glycan: pd.DataFrame) -> tuple[str, list[str]]:
    if glycan.empty or "window_id" not in glycan.columns:
        return "unknown:no_glycan_guardrail_table", []
    subset = glycan[glycan["window_id"].astype(str).eq(str(window_id))]
    if subset.empty:
        return "low:no_window_glycan_guardrails", []

    mode = _join(subset.get("glycan_guardrail_mode", pd.Series(dtype=str)).dropna().astype(str).tolist()) or "unknown_mode"
    risk = subset.get("mutation_expansion_risk", pd.Series(dtype=str)).astype(str).str.lower()
    reasons = subset.get("guardrail_reason", pd.Series(dtype=str)).astype(str).str.lower()
    high = risk.str.contains("high|forbidden").any()
    forbidden = reasons.str.contains("forbid|forbidden|extra_expansion").any()
    control = subset.get("control_needed", pd.Series([False] * len(subset))).map(_truthy).any()
    if high or forbidden:
        level = "high"
        rules = ["glycan_high_or_forbidden_class"]
    elif control:
        level = "medium"
        rules = ["glycan_control_needed"]
    else:
        level = "low"
        rules = []
    return f"{level}:mode={mode};positions={len(subset)}", rules


def _mask_summaries(mask_rows: pd.DataFrame) -> tuple[int, int, float, float, list[str]]:
    if mask_rows.empty:
        return 0, 0, 1.0, 1.0, ["missing_mask_rows"]
    status = mask_rows.get("mask_status", pd.Series([""] * len(mask_rows))).astype(str).str.lower()
    hard_count = int(status.str.contains("hard").sum())
    soft_count = int(status.str.contains("soft|risk").sum())
    designable = int(status.str.contains("designable|soft").sum())
    his_candidates = int(mask_rows.get("allowed_mutation_classes", pd.Series([""] * len(mask_rows))).astype(str).str.contains("histidine").sum())
    denominator = max(len(mask_rows), 1)
    notes: list[str] = []
    if len(mask_rows) != 40:
        notes.append(f"mask_row_count={len(mask_rows)}")
    return designable, his_candidates, hard_count / denominator, soft_count / denominator, notes


def _score_window(row: pd.Series, bundle: TableBundle, target_blocker: tuple[str, str] | None) -> dict[str, Any]:
    target = str(row.get("target", ""))
    window_id = str(row.get("window_id", ""))
    mask_rows = _positions_for_window(window_id, bundle.mask)
    capacity = _capacity_for_window(window_id, bundle.capacity)
    capacity_cfg = load_capacity_config()
    rules = load_decision_rules()
    thresholds = rules.get("thresholds", {})
    min_capacity = int(capacity_cfg.get("min_usable_noncontrol_variant_count", thresholds.get("min_usable_noncontrol_variant_count", 10000)))
    min_designable = int(thresholds.get("min_designable_positions_for_main_window", 8))
    min_his = int(thresholds.get("min_his_candidates_for_A_or_B", 2))
    max_hard_fraction = float(thresholds.get("max_hard_protected_fraction_for_A_or_B", 0.35))
    max_high_risk_fraction = float(thresholds.get("max_high_risk_fraction_for_A_or_B", 0.35))

    raw_capacity = int(_numeric(capacity.get("raw_feasible_variant_count"), 0))
    usable_capacity = int(_numeric(capacity.get("usable_noncontrol_variant_count"), 0))
    designable, his_candidates, hard_fraction, soft_fraction, mask_notes = _mask_summaries(mask_rows)
    hard_constraints, soft_constraints = _prior_constraints_for_window(row, bundle.prior_constraints, mask_rows)
    wet_summary = _wet_prior_summary(row, bundle.wet_tables.get(target, pd.DataFrame()), mask_rows)
    stable_contact, structure_summary = _structure_summary(row, bundle.residue_features, mask_rows)
    glycan_summary, glycan_rules = _glycan_summary(window_id, bundle.glycan_guardrails)

    triggered: list[str] = []
    failed: list[str] = []
    soft: list[str] = []
    hard: list[str] = []

    if target_blocker:
        blocker_type, blocker_reason = target_blocker
        return {
            "target": target,
            "window_id": window_id,
            "segment": row.get("segment", ""),
            "start": row.get("start", ""),
            "end": row.get("end", ""),
            "grade": "not_scored",
            "recommendation_status": "blocked",
            "blocker_type": blocker_type,
            "blocker_reason": blocker_reason,
            "decision_path": f"target_blocker:{blocker_type}",
            "evidence_route": "insufficient_evidence",
            "pH_mechanism_hypothesis_strength": "not_scored_due_to_target_blocker",
            "triggered_rules": "",
            "failed_rules": blocker_type,
            "soft_downgrade_rules": "",
            "hard_filter_reasons": "",
            "wet_prior_summary": wet_summary,
            "structure_summary": structure_summary,
            "glycan_risk_summary": glycan_summary,
            "raw_feasible_variant_count": raw_capacity,
            "usable_noncontrol_variant_count": usable_capacity,
        }

    illegal_reasons: list[str] = []
    if int(_numeric(row.get("length"), 0)) != 40:
        illegal_reasons.append("window_length_not_40")
    if _truthy(row.get("hard_excluded")):
        illegal_reasons.append(str(row.get("hard_exclusion_reasons", "hard_excluded")))
    if _truthy(row.get("crosses_linker")):
        illegal_reasons.append("crosses_linker")
    if hard_constraints:
        illegal_reasons.extend(f"hard_constraint:{item}" for item in hard_constraints)
    if usable_capacity < min_capacity:
        illegal_reasons.append(f"window_capacity_below_{min_capacity}")
    if mask_notes:
        soft.extend(mask_notes)

    if hard_fraction > max_hard_fraction:
        illegal_reasons.append(f"hard_protected_fraction>{max_hard_fraction}")
    if soft_fraction > max_high_risk_fraction:
        soft.append(f"soft_or_high_risk_fraction>{max_high_risk_fraction}")
    if designable < min_designable:
        soft.append(f"designable_positions<{min_designable}")
    if soft_constraints:
        soft.extend(f"prior_risk:{item}" for item in soft_constraints)
    soft.extend(glycan_rules)

    if illegal_reasons:
        grade = "D"
        status = "exploratory_only"
        evidence_route = "insufficient_evidence"
        pH_strength = "insufficient"
        hard = illegal_reasons
        failed.extend(illegal_reasons)
        decision = "hard_filter_failed_or_unacceptable_risk"
    else:
        if stable_contact and his_candidates >= min_his:
            evidence_route = "antigen_contact_plus_CDR_His"
            pH_strength = "strong_hypothesis"
            triggered.append("direct_CDR_His_route")
        elif stable_contact and designable >= min_designable:
            evidence_route = "framework_CDR_control"
            pH_strength = "moderate_hypothesis"
            triggered.append("structure_control_route")
        elif "direct_or_neighbor_wet_observations" in wet_summary and designable >= min_designable:
            evidence_route = "wet_module_extension"
            pH_strength = "moderate_hypothesis"
            triggered.append("wet_module_extension_route")
        else:
            evidence_route = "insufficient_evidence"
            pH_strength = "insufficient"
            failed.append("no_stable_structure_or_wet_module_route")

        high_glycan = glycan_summary.startswith("high:")
        manageable_risk = not high_glycan and hard_fraction <= max_hard_fraction and soft_fraction <= max_high_risk_fraction
        if triggered and manageable_risk and designable >= min_designable:
            grade = "A" if pH_strength == "strong_hypothesis" else "B"
            status = "recommended"
            decision = _join(triggered)
        elif triggered and not manageable_risk:
            grade = "C"
            status = "exploratory_only"
            decision = "support_present_but_risk_high"
            soft.append("risk_prevents_A_or_B")
        else:
            grade = "C"
            status = "exploratory_only"
            decision = "insufficient_evidence_or_high_risk"

    return {
        "target": target,
        "window_id": window_id,
        "segment": row.get("segment", ""),
        "start": row.get("start", ""),
        "end": row.get("end", ""),
        "grade": grade,
        "recommendation_status": status,
        "blocker_type": "",
        "blocker_reason": "",
        "decision_path": decision,
        "evidence_route": evidence_route,
        "pH_mechanism_hypothesis_strength": pH_strength,
        "triggered_rules": _join(triggered),
        "failed_rules": _join(failed),
        "soft_downgrade_rules": _join(soft),
        "hard_filter_reasons": _join(hard),
        "wet_prior_summary": wet_summary,
        "structure_summary": structure_summary,
        "glycan_risk_summary": glycan_summary,
        "raw_feasible_variant_count": raw_capacity,
        "usable_noncontrol_variant_count": usable_capacity,
    }


def target_level_rows(blockers: dict[str, tuple[str, str]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for target, (blocker_type, blocker_reason) in blockers.items():
        rows.append(
            {
                "target": target,
                "window_id": "TARGET_LEVEL",
                "segment": "",
                "start": "",
                "end": "",
                "grade": "not_scored",
                "recommendation_status": "blocked",
                "blocker_type": blocker_type,
                "blocker_reason": blocker_reason,
                "decision_path": f"target_blocker:{blocker_type}",
                "evidence_route": "insufficient_evidence",
                "pH_mechanism_hypothesis_strength": "not_scored_due_to_target_blocker",
                "triggered_rules": "",
                "failed_rules": blocker_type,
                "soft_downgrade_rules": "",
                "hard_filter_reasons": "",
                "wet_prior_summary": "",
                "structure_summary": "",
                "glycan_risk_summary": "",
                "raw_feasible_variant_count": 0,
                "usable_noncontrol_variant_count": 0,
            }
        )
    return rows


def score_windows(output_root: str | Path) -> pd.DataFrame:
    bundle = load_tables(output_root)
    blockers = collect_target_blockers(bundle)

    rows: list[dict[str, Any]] = []
    if bundle.candidate_windows.empty:
        rows.extend(target_level_rows(blockers or {target: ("schema_error", "candidate_windows.csv missing or empty") for target in _target_names()}))
    else:
        for _, candidate in bundle.candidate_windows.iterrows():
            target = str(candidate.get("target", ""))
            rows.append(_score_window(candidate, bundle, blockers.get(target)))

    scores = pd.DataFrame(rows)
    for column in OUTPUT_COLUMNS:
        if column not in scores.columns:
            scores[column] = ""
    return scores[OUTPUT_COLUMNS]


def main() -> None:
    parser = argparse.ArgumentParser(description="Score candidate 40 aa windows for pH-sensitive antibody design.")
    parser.add_argument("--output-root", default=load_inputs()["output_root"], help="Output root containing tables/.")
    args = parser.parse_args()

    scores = score_windows(args.output_root)
    out_path = _tables_dir(args.output_root) / "window_scores.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(out_path, index=False)
    print(f"Wrote {len(scores)} rows to {out_path}")


if __name__ == "__main__":
    main()
