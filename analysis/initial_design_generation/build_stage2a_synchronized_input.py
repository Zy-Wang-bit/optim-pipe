#!/usr/bin/env python3
"""Build synchronized Stage-2A input lists and audit reports.

This script intentionally does not run Stage-2A compute. It constructs audited
input lists from existing Stage-1/Stage-1.5 and sdAb recovery outputs.
"""

from __future__ import annotations

import csv
import hashlib
import math
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(".")
OUT_DIR = ROOT / "results/initial_design_generation/stage2a_synchronized_input"
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"

SDAB_BANK = ROOT / "results/initial_design_generation/sdab_recovery_loop/round_02_passlike_supplement/sdab_candidate_bank.csv"
STAGE1_RISK = ROOT / "results/initial_design_generation/stage1_5_stage2a/stage1_refined_structure_risk.csv"
STAGE2A_OLD = ROOT / "results/initial_design_generation/stage1_5_stage2a/stage2a_candidate_list.csv"
ONEE62_CONTROL_DRAFT = ROOT / "results/initial_design_generation/sdab_recovery_loop/1e62_sync_waiting_bank/1e62_controls_anchors_draft.csv"
VALIDATOR_CONFIG = TASK_DIR / "stage2a_input_validator_config.yaml"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"
PROGRESS = TASK_DIR / "progress.md"


STANDARD_COLUMNS = [
    "variant_id",
    "target",
    "canonical_sequence_hash",
    "window_sequence_hash",
    "canonical_unique_key",
    "mutation_list",
    "source_bank",
    "source_round",
    "boundary_support_tags",
    "boundary_support_count",
    "boundary_support_level",
    "compute_eligible",
    "audit_only",
    "control_anchor",
    "final_slot_type",
    "input_stratum",
    "selection_rank",
    "selection_reason",
    "sequence_hash",
    "canonical_sequence_hash_full",
    "canonical_recovery_sequence_hash",
    "canonical_mutated_window_sequence",
    "normalized_mutation_list",
    "primary_generation_route",
    "all_source_routes",
    "his_seed_set",
    "near_duplicate_cluster_id",
    "mutation_count",
    "His_count",
    "refined_structure_risk_class",
    "bank_eligibility",
    "supported_boundary_status",
    "stage2a_list_action",
    "stage2_action",
    "risk_reason_codes",
    "selection_class_from_tier1_or_stage1",
    "control_type",
    "control_type_if_applicable",
    "rosetta_delta",
    "local_validity_score",
    "new_clash_count_total",
    "new_clash_within_mutation_6A_fraction",
    "new_bad_contact_count",
    "lost_parent_contact_count",
    "his_min_antigen_distance",
    "his_sasa_min",
    "his_sasa_median",
    "cdr_ca_rmsd",
    "window_ca_rmsd",
    "mutation_shell_ca_rmsd",
    "neutral_retention_t2_score",
    "neutral_retention_score",
    "mpnn_score_status",
    "mpnn_parent_delta",
    "mpnn_total_score_per_residue",
    "mpnn_score_raw",
    "recovery_round",
    "bank_source",
]


@dataclass
class AuditResult:
    target: str
    status: str
    reasons: list[str]
    metrics: dict[str, str | int | float]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def to_float(value: object, default: float | None = None) -> float | None:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def to_int(value: object, default: int = 0) -> int:
    try:
        if value is None or value == "":
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def first_nonempty(row: dict[str, str], names: Iterable[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value not in ("", None):
            return str(value)
    return ""


def canonical_hash(row: dict[str, str]) -> str:
    return first_nonempty(
        row,
        [
            "canonical_full_sequence_hash",
            "full_antibody_sequence_hash",
            "canonical_recovery_sequence_hash",
            "canonical_sequence_hash_full",
            "sequence_hash",
            "variant_id",
        ],
    )


def window_hash(row: dict[str, str]) -> str:
    value = first_nonempty(
        row,
        [
            "canonical_sequence_hash_short",
            "window_sequence_hash",
            "sequence_hash",
            "canonical_mutated_window_sequence",
        ],
    )
    if value == row.get("canonical_mutated_window_sequence") and value:
        return hashlib.sha256(value.encode()).hexdigest()[:12]
    return value


def target_label(target: str) -> str:
    if target in {"Ab_sdAb", "sdAb"}:
        return "sdAb"
    if target in {"Ab_1E62", "1E62"}:
        return "1E62"
    return target


def seed_contains(row: dict[str, str], token: str) -> bool:
    seed = row.get("his_seed_set", "") or row.get("mutation_list", "")
    muts = row.get("mutation_list", "")
    return token in seed or token in muts


def support_tags(row: dict[str, str]) -> list[str]:
    tags: set[str] = set()
    status = row.get("supported_boundary_status", "")
    route = ";".join(
        [
            row.get("primary_generation_route", ""),
            row.get("all_source_routes", ""),
            row.get("selection_class_from_tier1_or_stage1", ""),
            row.get("p0_evidence_class", ""),
        ]
    )
    risk = row.get("risk_reason_codes", "")

    if "low_clash_boundary" in status:
        tags.add("low_new_clash")
    if "neutral_retention" in status:
        tags.add("neutral_retention_proxy_supported")
    if "mpnn_provenance" in status:
        tags.add("mpnn_score_supported")
        tags.add("rescue_provenance_supported")
    if "tier1_support" in status:
        tags.add("route_diversity_supported")
    if "wetlab_seed" in status or "wetlab" in route:
        tags.add("wetlab_seed_supported")

    local_validity = to_float(row.get("local_validity_score"))
    if local_validity is not None and local_validity >= 0.65:
        tags.add("acceptable_local_validity")

    new_clash = to_float(row.get("new_clash_count_total"))
    old_clash = to_float(row.get("local_clash_count"))
    clash = new_clash if new_clash is not None else old_clash
    shell_fraction = to_float(row.get("new_clash_within_mutation_6A_fraction"), 0.0)
    if clash is not None and (clash <= 1 or (clash <= 5 and shell_fraction is not None and shell_fraction >= 0.8)):
        tags.add("low_new_clash")

    new_bad = to_float(row.get("new_bad_contact_count"))
    if new_bad is not None and new_bad <= 0:
        tags.add("no_new_bad_interface_contact")

    lost_parent = to_float(row.get("lost_parent_contact_count"))
    if lost_parent is not None and lost_parent <= 1:
        tags.add("limited_parent_contact_loss")

    his_dist = to_float(row.get("his_min_antigen_distance"))
    his_sasa = to_float(row.get("his_sasa_min"))
    if (his_dist is not None and his_dist <= 12.0) or (his_sasa is not None and 10.0 <= his_sasa <= 180.0):
        tags.add("his_near_interface_or_interpretable")

    neutral_t2 = to_float(row.get("neutral_retention_t2_score"))
    neutral_v0 = to_float(row.get("neutral_retention_score"))
    if (neutral_t2 is not None and neutral_t2 >= 0.70) or (neutral_v0 is not None and neutral_v0 >= 0.75):
        tags.add("neutral_retention_proxy_supported")

    mpnn_parent_delta = to_float(row.get("mpnn_parent_delta"))
    if row.get("mpnn_score_status") == "scored_by_mpnn" or (mpnn_parent_delta is not None and mpnn_parent_delta <= 0):
        tags.add("mpnn_score_supported")

    if row.get("rescue_count", "") not in ("", "0", "0.0") or "rescue" in route.lower():
        tags.add("rescue_provenance_supported")

    if "alternate_support" in risk:
        tags.add("route_diversity_supported")

    seed = row.get("his_seed_set", "")
    if target_label(row.get("target", "")) == "sdAb":
        if seed in {"AY111H", "AD110H", "AQ100H;AY111H", "AQ100H;AD110H"}:
            tags.add("wetlab_seed_supported")
    elif target_label(row.get("target", "")) == "1E62":
        if any(pos in seed for pos in ["LK24H", "LY31H", "LD34H", "LQ35H", "LY38H"]):
            tags.add("wetlab_seed_supported")

    return sorted(tags)


def refined_class(row: dict[str, str]) -> str:
    return row.get("refined_structure_risk_class") or row.get("tier2_class") or row.get("t2_class_current", "")


def boundary_level(row: dict[str, str], tags: list[str]) -> str:
    klass = refined_class(row)
    if klass == "T2_pass_like_structure":
        return "pass_like"
    if klass == "T2_manual_review" or row.get("primary_generation_route") == "control_anchor":
        return "control_anchor"
    if klass == "T2_severe_structure_risk":
        return "severe"
    if "unsupported" in row.get("supported_boundary_status", ""):
        return "unsupported_boundary"
    if "boundary" in klass:
        if len(tags) >= 2:
            return "high_support_boundary"
        if len(tags) == 1:
            return "medium_support_boundary"
        return "weak_boundary"
    return "unclassified"


def normalize_row(row: dict[str, str], source_bank: str, source_round: str = "") -> dict[str, object]:
    tags = support_tags(row)
    level = boundary_level(row, tags)
    canon = canonical_hash(row)
    win = window_hash(row)
    out: dict[str, object] = dict(row)
    out.update(
        {
            "target": row.get("target", ""),
            "variant_id": row.get("variant_id", ""),
            "canonical_sequence_hash": canon,
            "window_sequence_hash": win,
            "canonical_unique_key": f"{row.get('target', '')}|{canon}",
            "source_bank": source_bank,
            "source_round": source_round or row.get("recovery_round", "") or row.get("bank_source", ""),
            "boundary_support_tags": ";".join(tags),
            "boundary_support_count": len(tags),
            "boundary_support_level": level,
            "compute_eligible": "yes",
            "audit_only": "no",
            "control_anchor": "no",
            "final_slot_type": "design",
            "input_stratum": "",
            "selection_rank": "",
            "selection_reason": "",
        }
    )
    return out


def selection_score(row: dict[str, object]) -> tuple[float, float, float, float]:
    level = str(row.get("boundary_support_level", ""))
    class_score = {
        "pass_like": 100.0,
        "high_support_boundary": 70.0,
        "medium_support_boundary": 40.0,
        "control_anchor": 30.0,
    }.get(level, 0.0)
    stage2_rank = to_float(row.get("stage2a_rank_score"), 0.0) or 0.0
    local_validity = to_float(row.get("local_validity_score"), 0.0) or 0.0
    neutral = to_float(row.get("neutral_retention_t2_score"), to_float(row.get("neutral_retention_score"), 0.0)) or 0.0
    support_count = to_float(row.get("boundary_support_count"), 0.0) or 0.0
    return (class_score, support_count, stage2_rank + neutral, local_validity)


def deduplicate(rows: Iterable[dict[str, object]]) -> list[dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for row in rows:
        key = str(row.get("canonical_unique_key", ""))
        if not key:
            key = f"{row.get('target')}|{row.get('variant_id')}"
        if key not in best or selection_score(row) > selection_score(best[key]):
            best[key] = row
    return list(best.values())


def cap_value(total: int, fraction: float) -> int:
    return max(1, math.floor(total * fraction))


def can_add(row: dict[str, object], selected: list[dict[str, object]], caps: dict[str, int | float]) -> bool:
    target_total = int(caps.get("target_total", 999999))
    projected_total = max(target_total, len(selected) + 1)
    seed_cap = int(caps.get("seed_cap", cap_value(projected_total, 0.22)))
    cluster_cap = int(caps.get("cluster_cap", cap_value(projected_total, 0.10)))
    ag102_cap = int(caps.get("ag102_cap", 999999))
    supplement_cap = int(caps.get("supplement_cap", 999999))
    four_mut_cap = int(caps.get("four_mut_cap", 999999))

    seed = str(row.get("his_seed_set", ""))
    cluster = str(row.get("near_duplicate_cluster_id", ""))
    if seed and sum(1 for item in selected if item.get("his_seed_set") == seed) >= seed_cap:
        return False
    if cluster and sum(1 for item in selected if item.get("near_duplicate_cluster_id") == cluster) >= cluster_cap:
        return False
    if seed_contains(row, "AG102H") and sum(1 for item in selected if seed_contains(item, "AG102H")) >= ag102_cap:
        return False
    if str(row.get("bank_source", "")) == "round_02_passlike_supplement_stage1_5":
        if sum(1 for item in selected if str(item.get("bank_source", "")) == "round_02_passlike_supplement_stage1_5") >= supplement_cap:
            return False
    if to_int(row.get("mutation_count")) >= 4:
        if sum(1 for item in selected if to_int(item.get("mutation_count")) >= 4) >= four_mut_cap:
            return False
    return True


def add_rows(
    selected: list[dict[str, object]],
    pool: Iterable[dict[str, object]],
    count: int,
    caps: dict[str, int | float],
    final_slot_type: str,
    input_stratum: str,
    reason: str,
) -> None:
    selected_keys = {str(row.get("canonical_unique_key", "")) for row in selected}
    for row in sorted(pool, key=selection_score, reverse=True):
        if len(selected) >= count:
            break
        key = str(row.get("canonical_unique_key", ""))
        if key in selected_keys:
            continue
        if not can_add(row, selected, caps):
            continue
        row = dict(row)
        row["final_slot_type"] = final_slot_type
        row["input_stratum"] = input_stratum
        row["selection_reason"] = reason
        selected.append(row)
        selected_keys.add(key)


def mark_controls(selected: list[dict[str, object]], desired_count: int) -> None:
    existing = [row for row in selected if row.get("final_slot_type") in {"control", "audit"}]
    needed = max(0, desired_count - len(existing))
    if needed == 0:
        return
    candidates = sorted(
        [row for row in selected if row.get("final_slot_type") == "design"],
        key=lambda r: (
            str(r.get("primary_generation_route", "")) != "control_anchor",
            str(r.get("input_stratum", "")) != "diversity_representative",
            str(r.get("his_seed_set", "")),
            str(r.get("near_duplicate_cluster_id", "")),
        ),
    )
    for row in candidates[:needed]:
        row["final_slot_type"] = "control"
        row["control_anchor"] = "yes"
        if "selection_reason" in row and row["selection_reason"]:
            row["selection_reason"] = f"{row['selection_reason']};control_anchor_quota"
        else:
            row["selection_reason"] = "control_anchor_quota"


def mark_diversity(selected: list[dict[str, object]], desired_count: int) -> None:
    existing = [row for row in selected if row.get("final_slot_type") == "diversity"]
    needed = max(0, desired_count - len(existing))
    if needed == 0:
        return
    buckets: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in selected:
        if row.get("final_slot_type") != "design":
            continue
        key = (
            str(row.get("his_seed_set", "")),
            str(row.get("near_duplicate_cluster_id", "")),
            str(row.get("primary_generation_route", "")),
        )
        if key not in buckets or selection_score(row) > selection_score(buckets[key]):
            buckets[key] = row
    candidates = sorted(buckets.values(), key=selection_score, reverse=True)
    for row in candidates[:needed]:
        row["final_slot_type"] = "diversity"
        row["input_stratum"] = "diversity_representative"
        if row.get("selection_reason"):
            row["selection_reason"] = f"{row['selection_reason']};diversity_representative"
        else:
            row["selection_reason"] = "diversity_representative"


def finalize_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    for idx, row in enumerate(rows, 1):
        row["selection_rank"] = idx
        if row.get("final_slot_type") in {"control", "audit"}:
            row["control_anchor"] = "yes"
        if row.get("final_slot_type") == "audit":
            row["audit_only"] = "yes"
            row["compute_eligible"] = "no"
        if seed_contains(row, "AV105H") and row.get("final_slot_type") not in {"control", "audit"}:
            row["compute_eligible"] = "no"
            row["audit_only"] = "yes"
    return rows


def build_sdab_input() -> list[dict[str, object]]:
    bank = [normalize_row(row, "sdab_final_recovery_bank") for row in read_csv(SDAB_BANK)]
    controls = [
        normalize_row(row, "stage1_control_anchor", "stage1_5_control")
        for row in read_csv(STAGE1_RISK)
        if row.get("target") == "Ab_sdAb" and row.get("primary_generation_route") == "control_anchor"
    ]
    for row in controls:
        row["final_slot_type"] = "control"
        row["control_anchor"] = "yes"
        if seed_contains(row, "AV105H"):
            row["final_slot_type"] = "audit"
            row["audit_only"] = "yes"
            row["compute_eligible"] = "no"
        row["input_stratum"] = "control_anchor"
        row["selection_reason"] = "explicit_control_anchor"

    bank = deduplicate(bank)
    target_total = 820
    caps = {
        "target_total": target_total,
        "seed_cap": cap_value(target_total, 0.22),
        "cluster_cap": cap_value(target_total, 0.10),
        "ag102_cap": cap_value(target_total, 0.03),
        "supplement_cap": cap_value(target_total, 0.40),
        "four_mut_cap": cap_value(target_total, 0.05),
    }

    selected: list[dict[str, object]] = []
    add_rows(selected, controls, len(controls), caps, "control", "control_anchor", "explicit_control_anchor")

    pass_like = [
        row
        for row in bank
        if row.get("boundary_support_level") == "pass_like"
        and not (seed_contains(row, "AV105H") or "V105H;D110H" in str(row.get("mutation_list", "")))
    ]
    add_rows(selected, pass_like, len(selected) + 390, caps, "design", "pass_like", "pass_like_priority")

    nonboundary_count = lambda: sum(
        1
        for row in selected
        if row.get("boundary_support_level") not in {"high_support_boundary", "medium_support_boundary", "weak_boundary", "unsupported_boundary", "severe"}
    )
    boundary_max = min(math.floor(target_total * 0.55), math.floor(nonboundary_count() * 0.55 / 0.45))
    boundary_selected = lambda: sum(1 for row in selected if row.get("boundary_support_level") in {"high_support_boundary", "medium_support_boundary"})
    high_boundary = [
        row
        for row in bank
        if row.get("boundary_support_level") == "high_support_boundary"
        and not seed_contains(row, "AV105H")
        and "V105H;D110H" not in str(row.get("mutation_list", ""))
    ]
    while len(selected) < target_total and boundary_selected() < boundary_max:
        before = len(selected)
        add_rows(
            selected,
            high_boundary,
            len(selected) + 1,
            caps,
            "design",
            "high_support_boundary",
            "high_support_boundary_fill",
        )
        if len(selected) == before:
            break

    boundary_max = min(math.floor(target_total * 0.55), math.floor(nonboundary_count() * 0.55 / 0.45))
    diversity_target = min(120, max(80, target_total - len(selected)))
    medium_boundary = [
        row
        for row in bank
        if row.get("boundary_support_level") == "medium_support_boundary" and not seed_contains(row, "AV105H")
    ]
    if len(selected) < target_total and boundary_selected() < boundary_max:
        add_rows(
            selected,
            medium_boundary,
            min(target_total, len(selected) + diversity_target),
            caps,
            "diversity",
            "diversity_representative",
            "medium_boundary_diversity",
        )

    mark_diversity(selected, 80)
    mark_controls(selected, 50)
    return finalize_rows(selected)


def build_1e62_input() -> list[dict[str, object]]:
    stage1 = [
        normalize_row(row, "stage1_5_refined_structure_risk", "stage1_5")
        for row in read_csv(STAGE1_RISK)
        if row.get("target") == "Ab_1E62"
    ]
    controls = [
        normalize_row(row, "1e62_control_anchor_draft", "stage1_5_control")
        for row in read_csv(ONEE62_CONTROL_DRAFT)
    ]
    for row in controls:
        row["final_slot_type"] = "control"
        row["control_anchor"] = "yes"
        row["input_stratum"] = "control_anchor"
        row["selection_reason"] = "explicit_control_anchor"

    candidates = deduplicate(stage1)
    target_total = 550
    caps = {
        "target_total": target_total,
        "seed_cap": cap_value(target_total, 0.22),
        "cluster_cap": cap_value(target_total, 0.12),
        "four_mut_cap": cap_value(target_total, 0.10),
    }

    selected: list[dict[str, object]] = []
    add_rows(selected, controls, len(controls), caps, "control", "control_anchor", "explicit_control_anchor")

    pass_like = [row for row in candidates if row.get("boundary_support_level") == "pass_like"]
    add_rows(selected, pass_like, len(selected) + 250, caps, "design", "pass_like", "pass_like_priority")

    high_boundary = [
        row
        for row in candidates
        if row.get("boundary_support_level") == "high_support_boundary"
        and row.get("stage2a_list_action") in {"stage2a_include_low_quota", "stage2a_include_priority", "stage2a_hold"}
    ]
    add_rows(selected, high_boundary, target_total, caps, "design", "high_support_boundary", "high_support_boundary_fill")

    if len(selected) < target_total:
        medium_boundary = [
            row
            for row in candidates
            if row.get("boundary_support_level") == "medium_support_boundary"
            and row.get("stage2a_list_action") in {"stage2a_include_low_quota", "stage2a_include_priority", "stage2a_hold"}
        ]
        add_rows(selected, medium_boundary, target_total, caps, "diversity", "diversity_representative", "medium_boundary_diversity")

    mark_diversity(selected, 80)
    mark_controls(selected, 50)
    return finalize_rows(selected)


def metric_summary(rows: list[dict[str, object]], target: str) -> dict[str, str | int | float]:
    total = len(rows)
    unique = len({row.get("canonical_unique_key") for row in rows})
    duplicate_count = total - unique
    levels = Counter(str(row.get("boundary_support_level", "")) for row in rows)
    seeds = Counter(str(row.get("his_seed_set", "")) for row in rows if row.get("his_seed_set", "") != "")
    clusters = Counter(str(row.get("near_duplicate_cluster_id", "")) for row in rows if row.get("near_duplicate_cluster_id", "") != "")
    top_seed_count = seeds.most_common(1)[0][1] if seeds else 0
    top_seed = seeds.most_common(1)[0][0] if seeds else ""
    top_cluster_count = clusters.most_common(1)[0][1] if clusters else 0
    top_cluster = clusters.most_common(1)[0][0] if clusters else ""
    four_mut_count = sum(1 for row in rows if to_int(row.get("mutation_count")) >= 4)
    ag102_count = sum(1 for row in rows if seed_contains(row, "AG102H"))
    av105_main = sum(1 for row in rows if seed_contains(row, "AV105H") and row.get("final_slot_type") not in {"control", "audit"})
    supplement = sum(1 for row in rows if row.get("bank_source") == "round_02_passlike_supplement_stage1_5")
    controls = sum(1 for row in rows if row.get("control_anchor") == "yes" or row.get("final_slot_type") in {"control", "audit"})
    diversity = sum(1 for row in rows if row.get("final_slot_type") == "diversity" or row.get("input_stratum") == "diversity_representative")
    compute_eligible = sum(1 for row in rows if row.get("compute_eligible") == "yes")
    audit_only = sum(1 for row in rows if row.get("audit_only") == "yes")
    supported_boundary = levels.get("high_support_boundary", 0) + levels.get("medium_support_boundary", 0)
    pass_like = levels.get("pass_like", 0)
    metrics: dict[str, str | int | float] = {
        "target": target,
        "total_input": total,
        "canonical_unique_count": unique,
        "exact_duplicate_count": duplicate_count,
        "pass_like_count": pass_like,
        "pass_like_fraction": round(pass_like / total, 4) if total else 0,
        "supported_boundary_count": supported_boundary,
        "supported_boundary_fraction": round(supported_boundary / total, 4) if total else 0,
        "high_support_boundary_count": levels.get("high_support_boundary", 0),
        "medium_support_boundary_count": levels.get("medium_support_boundary", 0),
        "weak_boundary_count": levels.get("weak_boundary", 0),
        "unsupported_boundary_count": levels.get("unsupported_boundary", 0),
        "severe_count": levels.get("severe", 0),
        "control_anchor_count": controls,
        "diversity_representative_count": diversity,
        "top_seed": top_seed,
        "top_seed_count": top_seed_count,
        "top_seed_fraction": round(top_seed_count / total, 4) if total else 0,
        "top_near_duplicate_cluster": top_cluster,
        "top_near_duplicate_cluster_count": top_cluster_count,
        "top_near_duplicate_cluster_fraction": round(top_cluster_count / total, 4) if total else 0,
        "four_mut_count": four_mut_count,
        "four_mut_fraction": round(four_mut_count / total, 4) if total else 0,
        "AG102H_count": ag102_count,
        "AG102H_fraction": round(ag102_count / total, 4) if total else 0,
        "AV105H_main_count": av105_main,
        "supplement_derived_count": supplement,
        "supplement_derived_fraction": round(supplement / total, 4) if total else 0,
        "compute_eligible_count": compute_eligible,
        "audit_only_count": audit_only,
    }
    return metrics


def audit_target(target: str, rows: list[dict[str, object]]) -> AuditResult:
    m = metric_summary(rows, target)
    reasons: list[str] = []
    fail = False
    patch = False

    if m["exact_duplicate_count"] != 0:
        fail = True
        reasons.append("exact duplicate count is not zero")
    if m["severe_count"] != 0:
        fail = True
        reasons.append("severe-risk rows present")
    if m["unsupported_boundary_count"] != 0 or m["weak_boundary_count"] != 0:
        fail = True
        reasons.append("weak or unsupported boundary rows present")

    if target == "sdAb":
        if m["total_input"] < 600:
            fail = True
            reasons.append("sdAb total input < 600")
        elif m["total_input"] < 800:
            patch = True
            reasons.append("sdAb total input is pilot range 600-799")
        if m["pass_like_count"] < 300 or m["pass_like_fraction"] < 0.35:
            patch = True
            reasons.append("sdAb pass-like count/fraction below full PASS threshold")
        if m["supported_boundary_fraction"] > 0.55:
            patch = True
            reasons.append("sdAb supported-boundary fraction > 55%")
        if not 50 <= m["control_anchor_count"] <= 80:
            patch = True
            reasons.append("sdAb controls/anchors outside 50-80")
        if not 80 <= m["diversity_representative_count"] <= 120:
            patch = True
            reasons.append("sdAb diversity representatives outside 80-120")
        if m["top_seed_fraction"] > 0.25:
            fail = True
            reasons.append("sdAb top seed fraction > 25%")
        elif m["top_seed_fraction"] > 0.22:
            patch = True
            reasons.append("sdAb top seed fraction in PATCH range >22%")
        if m["top_near_duplicate_cluster_fraction"] > 0.12:
            fail = True
            reasons.append("sdAb top near-duplicate cluster fraction > 12%")
        elif m["top_near_duplicate_cluster_fraction"] > 0.10:
            patch = True
            reasons.append("sdAb top near-duplicate cluster fraction in PATCH range >10%")
        if m["AG102H_fraction"] > 0.05:
            fail = True
            reasons.append("sdAb AG102H fraction > 5%")
        elif m["AG102H_fraction"] > 0.03:
            patch = True
            reasons.append("sdAb AG102H fraction in PATCH range >3%")
        if m["AV105H_main_count"] > 0:
            fail = True
            reasons.append("sdAb AV105H appears as main candidate")
        if m["supplement_derived_fraction"] > 0.50:
            fail = True
            reasons.append("sdAb supplement-derived fraction > 50%")
        elif m["supplement_derived_fraction"] > 0.40:
            patch = True
            reasons.append("sdAb supplement-derived fraction in PATCH range >40%")
    elif target == "1E62":
        if m["total_input"] < 280:
            fail = True
            reasons.append("1E62 input count < 280")
        elif m["total_input"] < 500:
            patch = True
            reasons.append("1E62 input count 280-499")
        if m["top_seed_fraction"] > 0.25:
            fail = True
            reasons.append("1E62 top seed fraction > 25%")
        elif m["top_seed_fraction"] > 0.22:
            patch = True
            reasons.append("1E62 top seed fraction in PATCH range >22%")
        if m["top_near_duplicate_cluster_fraction"] > 0.12:
            fail = True
            reasons.append("1E62 top near-duplicate cluster fraction > 12%")
        if not 50 <= m["control_anchor_count"] <= 80:
            patch = True
            reasons.append("1E62 controls/anchors outside 50-80")

    status = "FAIL" if fail else "PATCH" if patch else "PASS"
    if not reasons:
        reasons.append("all configured checks passed")
    return AuditResult(target=target, status=status, reasons=reasons, metrics=m)


def summary_rows(rows: list[dict[str, object]], field: str) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row.get("target", "")), str(row.get(field, "")))].append(row)
    out = []
    for (target, value), items in sorted(grouped.items()):
        out.append(
            {
                "target": target,
                field: value,
                "count": len(items),
                "fraction_within_target": round(len(items) / max(1, sum(1 for row in rows if row.get("target") == target)), 4),
                "pass_like_count": sum(1 for row in items if row.get("boundary_support_level") == "pass_like"),
                "supported_boundary_count": sum(
                    1 for row in items if row.get("boundary_support_level") in {"high_support_boundary", "medium_support_boundary"}
                ),
                "control_anchor_count": sum(1 for row in items if row.get("control_anchor") == "yes"),
            }
        )
    return out


def write_audit_checks(audits: list[AuditResult], combined_status: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for audit in audits:
        for reason in audit.reasons:
            rows.append({"scope": audit.target, "status": audit.status, "check_or_reason": reason})
    rows.append({"scope": "overall", "status": combined_status, "check_or_reason": "combined target audit decision"})
    return rows


def write_compute_config(path: Path, input_path: Path, audit_report_path: Path) -> None:
    content = f"""stage: Stage-2A
status: draft_only_until_manual_unlock
input_list:
  path: {input_path.as_posix()}
  sha256: {sha256_file(input_path)}
audit_report:
  path: {audit_report_path.as_posix()}
  sha256: {sha256_file(audit_report_path)}
manual_unlock:
  required: true
  unlocked_by: null
  unlocked_at: null
targets:
  1E62:
    input_count: from_audit
    compute_scope: synchronized_stage2a
  sdAb:
    input_count: from_audit
    compute_scope: synchronized_stage2a
tools:
  pyrosetta_local_repack_minimize: planned_after_unlock
  pka_his_environment: planned_after_unlock
  simplefold_or_af3_heavy: hold
  md: hold
post_compute_gate:
  decide_targeted_tier2_heavy_planning: manual_review_required
"""
    path.write_text(content)


def md_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def build_report(audits: list[AuditResult], combined_status: str, outputs: dict[str, Path]) -> str:
    metrics_rows = [audit.metrics for audit in audits]
    check_rows = write_audit_checks(audits, combined_status)
    out_lines = [
        "# 当前阶段报告：Synchronized Stage-2A Input Construction + Audit",
        "",
        "日期：2026-05-29",
        "",
        "## 1. 阶段定位",
        "",
        "本阶段按照已确认方案执行 `Synchronized Stage-2A Input Construction + Audit`。本阶段只构建和审计 Stage-2A input list，不运行 Stage-2A compute，也不解锁 Tier2-heavy、AF3/SimpleFold heavy review、MD 或 final 10K。",
        "",
        "本轮执行目标：",
        "",
        "```text",
        "1. 冻结 sdAb final recovery bank；",
        "2. 审计并整理 1E62 waiting bank；",
        "3. 构建 sdAb 分层 input pool；",
        "4. 构建 synchronized Stage-2A input list；",
        "5. 用 validator 输出 PASS / PATCH / FAIL；",
        "6. 只在 PASS 且人工 unlock 后才允许后续 compute。",
        "```",
        "",
        f"本轮总体结论：`{combined_status}`。",
        "",
        "## 2. 输入来源",
        "",
        "| target | input source | role |",
        "| --- | --- | --- |",
        "| sdAb | final recovery bank after round 02 pass-like supplement | recovered candidate bank |",
        "| sdAb | Stage-1 control anchors | controls / audit anchors |",
        "| 1E62 | Stage-1.5 refined structure-risk table | waiting bank expansion source |",
        "| 1E62 | previous Stage-2A list and control-anchor draft | existing candidate/control reference |",
        "",
        "唯一候选口径固定为 `target + canonical_sequence_hash`。所有输出清单均保留 `variant_id`、`target`、`canonical_sequence_hash`、`window_sequence_hash`、`mutation_list`、`source_bank`、`source_round`、`boundary_support_tags`、`boundary_support_level`、`compute_eligible`、`audit_only`、`control_anchor` 和 `final_slot_type`。",
        "",
        "## 3. Output Package",
        "",
        "| output | path |",
        "| --- | --- |",
    ]
    for name, path in outputs.items():
        out_lines.append(f"| {name} | `{path.as_posix()}` |")
    out_lines.extend(
        [
            "",
            "## 4. Input List Metrics",
            "",
            md_table(
                metrics_rows,
                [
                    "target",
                    "total_input",
                    "canonical_unique_count",
                    "exact_duplicate_count",
                    "pass_like_count",
                    "pass_like_fraction",
                    "supported_boundary_count",
                    "supported_boundary_fraction",
                    "control_anchor_count",
                    "diversity_representative_count",
                    "top_seed",
                    "top_seed_fraction",
                    "top_near_duplicate_cluster",
                    "top_near_duplicate_cluster_fraction",
                    "four_mut_fraction",
                    "compute_eligible_count",
                    "audit_only_count",
                ],
            ),
            "",
            "## 5. sdAb 结果解读",
            "",
            "sdAb input list 从 1,380 条 final recovery bank 中分层抽取，并额外纳入 Stage-1 control/audit anchors。选择过程遵守以下限制：severe、weak boundary、unsupported boundary 不进入主 input；AV105H 只作为显式 audit/control；AG102H 严格限额；supplement-derived 候选不得主导 input list。",
            "",
        ]
    )
    sdab = next(a for a in audits if a.target == "sdAb")
    out_lines.append(f"sdAb audit status: `{sdab.status}`。")
    out_lines.append("")
    out_lines.append("sdAb 主要判定原因：")
    out_lines.append("")
    for reason in sdab.reasons:
        out_lines.append(f"- {reason}")
    out_lines.extend(
        [
            "",
            "## 6. 1E62 结果解读",
            "",
            "1E62 input list 没有按分数直接抽取，而是从已完成 Stage-1/Stage-1.5 结果中按 pass-like、supported-boundary、seed、near-duplicate cluster 和 controls/anchors 进行约束选择。当前 1E62 可以从已有计算结果中构建超过 500 条的同步 input list，不需要启动新的 heavy compute 来补数。",
            "",
        ]
    )
    one = next(a for a in audits if a.target == "1E62")
    out_lines.append(f"1E62 audit status: `{one.status}`。")
    out_lines.append("")
    out_lines.append("1E62 主要判定原因：")
    out_lines.append("")
    for reason in one.reasons:
        out_lines.append(f"- {reason}")
    out_lines.extend(
        [
            "",
            "## 7. Audit Checks",
            "",
            md_table(check_rows, ["scope", "status", "check_or_reason"]),
            "",
            "## 8. Decision Boundary",
            "",
            "本轮没有启动 Stage-2A compute。`stage2a_compute_config_draft.yaml` 只记录 future compute 的草案，并写入 input list 与 audit report 的 SHA256。真正 compute 仍需要人工 unlock。",
            "",
        ]
    )
    if combined_status == "PASS":
        out_lines.extend(
            [
                "当前 input package 达到 validator PASS，可进入人工 compute unlock 讨论。",
                "",
                "下一步建议：",
                "",
                "```text",
                "1. 人工审阅 stage2a_input_audit_report.md；",
                "2. 若认可 input list 和 compute draft，手动解锁 Stage-2A compute；",
                "3. compute 完成后再决定是否进入 targeted Tier2-heavy planning；",
                "4. Tier2-heavy / AF3 / MD / final 10K 仍保持锁定。",
                "```",
            ]
        )
    elif combined_status == "PATCH":
        out_lines.extend(
            [
                "当前 input package 方向正确，但仍需修补 input list。不得启动 Stage-2A compute。",
                "",
                "下一步建议：只修补 audit reason 中列出的局部问题，然后重新运行 input audit。",
            ]
        )
    else:
        out_lines.extend(
            [
                "当前 input package 未达到基本条件，不得启动 Stage-2A compute。",
                "",
                "下一步建议：回到 bank cleanup 或 targeted generation。",
            ]
        )
    out_lines.extend(
        [
            "",
            "## 9. Explicit Non-Goals",
            "",
            "本阶段仍不执行：",
            "",
            "```text",
            "Stage-2A compute",
            "Tier2-heavy",
            "AF3 / SimpleFold heavy review",
            "MD",
            "final 10K selection",
            "```",
            "",
        ]
    )
    return "\n".join(out_lines)


def append_progress(combined_status: str, audits: list[AuditResult]) -> None:
    lines = [
        "",
        "## 2026-05-29 - Synchronized Stage-2A input construction + audit execution",
        "",
        f"- Built synchronized Stage-2A input package; overall audit status = `{combined_status}`.",
    ]
    for audit in audits:
        m = audit.metrics
        lines.append(
            f"- {audit.target}: status `{audit.status}`, total={m['total_input']}, unique={m['canonical_unique_count']}, "
            f"pass-like={m['pass_like_count']} ({m['pass_like_fraction']}), "
            f"supported-boundary={m['supported_boundary_count']} ({m['supported_boundary_fraction']}), "
            f"top seed={m['top_seed']} ({m['top_seed_fraction']}), "
            f"top cluster={m['top_near_duplicate_cluster']} ({m['top_near_duplicate_cluster_fraction']})."
        )
    lines.extend(
        [
            "- Outputs written to `results/initial_design_generation/stage2a_synchronized_input/`.",
            "- `current_stage_report.md` overwritten with the new audit report.",
            "- Stage-2A compute remains locked pending manual unlock.",
            "",
        ]
    )
    with PROGRESS.open("a") as handle:
        handle.write("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if VALIDATOR_CONFIG.exists():
        shutil.copyfile(VALIDATOR_CONFIG, OUT_DIR / "stage2a_input_validator_config.yaml")

    sdab_rows = build_sdab_input()
    one_rows = build_1e62_input()
    combined = sdab_rows + one_rows

    output_fields = STANDARD_COLUMNS
    write_csv(OUT_DIR / "stage2a_input_sdAb.csv", sdab_rows, output_fields)
    write_csv(OUT_DIR / "stage2a_input_1E62.csv", one_rows, output_fields)
    write_csv(OUT_DIR / "stage2a_synchronized_input_list.csv", combined, output_fields)

    audits = [audit_target("sdAb", sdab_rows), audit_target("1E62", one_rows)]
    if any(audit.status == "FAIL" for audit in audits):
        combined_status = "FAIL"
    elif any(audit.status == "PATCH" for audit in audits):
        combined_status = "PATCH"
    else:
        combined_status = "PASS"

    metrics_rows = [audit.metrics for audit in audits]
    write_csv(OUT_DIR / "stage2a_input_audit_metrics.csv", metrics_rows, list(metrics_rows[0].keys()))
    check_rows = write_audit_checks(audits, combined_status)
    write_csv(OUT_DIR / "stage2a_input_audit_checks.csv", check_rows, ["scope", "status", "check_or_reason"])
    write_csv(OUT_DIR / "stage2a_input_seed_summary.csv", summary_rows(combined, "his_seed_set"), ["target", "his_seed_set", "count", "fraction_within_target", "pass_like_count", "supported_boundary_count", "control_anchor_count"])
    write_csv(OUT_DIR / "stage2a_input_cluster_summary.csv", summary_rows(combined, "near_duplicate_cluster_id"), ["target", "near_duplicate_cluster_id", "count", "fraction_within_target", "pass_like_count", "supported_boundary_count", "control_anchor_count"])
    write_csv(OUT_DIR / "stage2a_input_boundary_summary.csv", summary_rows(combined, "boundary_support_level"), ["target", "boundary_support_level", "count", "fraction_within_target", "pass_like_count", "supported_boundary_count", "control_anchor_count"])
    write_csv(OUT_DIR / "stage2a_input_control_anchor_summary.csv", summary_rows(combined, "final_slot_type"), ["target", "final_slot_type", "count", "fraction_within_target", "pass_like_count", "supported_boundary_count", "control_anchor_count"])

    outputs = {
        "stage2a_input_sdAb.csv": OUT_DIR / "stage2a_input_sdAb.csv",
        "stage2a_input_1E62.csv": OUT_DIR / "stage2a_input_1E62.csv",
        "stage2a_synchronized_input_list.csv": OUT_DIR / "stage2a_synchronized_input_list.csv",
        "stage2a_input_audit_checks.csv": OUT_DIR / "stage2a_input_audit_checks.csv",
        "stage2a_input_audit_metrics.csv": OUT_DIR / "stage2a_input_audit_metrics.csv",
        "stage2a_input_seed_summary.csv": OUT_DIR / "stage2a_input_seed_summary.csv",
        "stage2a_input_cluster_summary.csv": OUT_DIR / "stage2a_input_cluster_summary.csv",
        "stage2a_input_boundary_summary.csv": OUT_DIR / "stage2a_input_boundary_summary.csv",
        "stage2a_input_control_anchor_summary.csv": OUT_DIR / "stage2a_input_control_anchor_summary.csv",
        "stage2a_input_validator_config.yaml": OUT_DIR / "stage2a_input_validator_config.yaml",
        "stage2a_compute_config_draft.yaml": OUT_DIR / "stage2a_compute_config_draft.yaml",
        "stage2a_input_audit_report.md": OUT_DIR / "stage2a_input_audit_report.md",
    }
    report = build_report(audits, combined_status, outputs)
    audit_report_path = OUT_DIR / "stage2a_input_audit_report.md"
    audit_report_path.write_text(report)
    write_compute_config(OUT_DIR / "stage2a_compute_config_draft.yaml", OUT_DIR / "stage2a_synchronized_input_list.csv", audit_report_path)
    report = build_report(audits, combined_status, outputs)
    audit_report_path.write_text(report)
    CURRENT_STAGE_REPORT.write_text(report)
    append_progress(combined_status, audits)

    print(f"overall_status={combined_status}")
    for audit in audits:
        print(f"{audit.target}: {audit.status} {audit.metrics}")


if __name__ == "__main__":
    main()
