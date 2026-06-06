#!/usr/bin/env python3
"""Build final candidate-pool planning drafts from expanded Tier2 v2 pools.

This is a planning/audit stage. It does not create a synthesis-ready library
and does not run AF3, SimpleFold, Rosetta/FoldX, MD, glycan modeling, or any
new structural compute.
"""

import csv
import hashlib
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


ROOT = Path(__file__).resolve().parents[2]
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"
OUT_DIR = ROOT / "results/initial_design_generation/final_candidate_pool_planning"
INPUT_POOL = ROOT / "results/initial_design_generation/expanded_tier2_candidate_pool_v2/expanded_tier2_candidate_pool_v2.csv"
INPUT_AUDIT = ROOT / "results/initial_design_generation/expanded_tier2_candidate_pool_v2/expanded_pool_v2_audit_report.md"
CONTROL_SOURCE = ROOT / "results/initial_design_generation/tier2_route_limited_stage_b/tier2_heavy_route_controls.csv"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"


DEFAULT_SCENARIOS = {
    "10k": {"total": 10000, "quotas": {"Ab_1E62": 6000, "Ab_sdAb": 4000}},
    "15k": {"total": 15000, "quotas": {"Ab_1E62": 8500, "Ab_sdAb": 6500}},
}

ALT_SCENARIOS = {
    "10k_balanced": {"total": 10000, "quotas": {"Ab_1E62": 5000, "Ab_sdAb": 5000}},
    "10k_sdab_heavy": {"total": 10000, "quotas": {"Ab_1E62": 4500, "Ab_sdAb": 5500}},
    "15k_balanced": {"total": 15000, "quotas": {"Ab_1E62": 7500, "Ab_sdAb": 7500}},
    "15k_sdab_heavy": {"total": 15000, "quotas": {"Ab_1E62": 6500, "Ab_sdAb": 8500}},
}


CAPS = {
    "Ab_1E62": {
        "top_seed_soft": 0.20,
        "top_seed_hard": 0.25,
        "top_cluster_soft": 0.03,
        "top_cluster_hard": 0.05,
        "four_mut_soft": 0.50,
        "four_mut_hard": 0.55,
        "local_expansion_soft": 0.45,
        "local_expansion_hard": 0.55,
        "seed_only_soft": 0.10,
        "seed_only_hard": 0.12,
        "neutral_boundary_soft": 0.25,
        "neutral_boundary_hard": 0.30,
        "ay111h_soft": None,
        "ay111h_hard": None,
    },
    "Ab_sdAb": {
        "top_seed_soft": 0.25,
        "top_seed_hard": 0.30,
        "top_cluster_soft": 0.03,
        "top_cluster_hard": 0.05,
        "four_mut_soft": 0.08,
        "four_mut_hard": 0.10,
        "local_expansion_soft": 0.40,
        "local_expansion_hard": 0.45,
        "seed_only_soft": 0.15,
        "seed_only_hard": 0.18,
        "neutral_boundary_soft": 0.10,
        "neutral_boundary_hard": 0.15,
        "ay111h_soft": 0.30,
        "ay111h_hard": 0.35,
    },
}


TIER_SCORE = {
    "template_seed": 5.0,
    "A_strong": 4.4,
    "local_expansion_A": 4.0,
    "B_medium": 3.2,
    "local_expansion_B": 2.8,
    "C_seed_only": 1.0,
}

QUALITY_SCORE = {
    "bank_reviewed": 3.0,
    "high_support": 2.3,
    "supported": 1.8,
    "mpnn_favorable_pH_weak": 0.9,
    "neutral_boundary_or_high_risk": -0.9,
}

STATUS_SCORE = {
    "active_primary_template": 2.3,
    "active_primary_backfill_template": 1.5,
    "limited_boundary_template": -0.8,
    "secondary_template_complex_weak": 0.9,
    "secondary_template_complex_unchecked": 0.7,
    "low_priority_secondary_template": -0.4,
    "boundary_representative_only": -1.0,
}

SOURCE_SCORE = {
    "tier2_candidate_bank": 1.4,
    "stage2a_candidate_list": 1.0,
    "tier2_candidate_snapshot": 0.8,
    "sdab_recovery_passlike_supplement": 0.7,
    "tier1_ranked_production_pool": 0.35,
    "constrained_local_expansion_v2": -0.15,
    "constrained_low_frequency_seed_local_expansion_v2": -0.3,
}

EXTRA_FIELDS = [
    "selection_scenario",
    "selection_rank",
    "selection_phase",
    "final_selection_role",
    "final_planning_score",
    "score_template_component",
    "score_quality_component",
    "score_status_component",
    "score_source_component",
    "score_mechanism_component",
    "score_risk_penalty",
    "score_diversity_penalty",
    "is_local_expansion",
    "is_seed_only",
    "is_four_mut",
    "is_neutral_boundary_or_high_risk",
    "contains_AY111H",
    "contains_AG102H_or_AV105H",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: Optional[List[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: List[str] = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: Any, default: float = 0.0) -> float:
    if value is None:
        return default
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "<na>"}:
        return default
    try:
        out = float(text)
    except ValueError:
        return default
    if math.isnan(out) or math.isinf(out):
        return default
    return out


def as_int(value: Any, default: int = 0) -> int:
    return int(round(as_float(value, float(default))))


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def norm_key(value: Any, fallback: str) -> str:
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "<na>"}:
        return fallback
    return text


def mutation_tokens(row: Dict[str, str]) -> List[str]:
    text = row.get("normalized_mutation_list") or row.get("mutation_list") or ""
    return [tok.strip() for tok in str(text).split(";") if tok.strip()]


def mutation_introduces_cys(row: Dict[str, str]) -> bool:
    return any(tok.endswith("C") for tok in mutation_tokens(row))


def is_seed_only(row: Dict[str, str]) -> bool:
    return row.get("template_match_tier") == "C_seed_only"


def is_local_expansion(row: Dict[str, str]) -> bool:
    return row.get("template_match_tier", "").startswith("local_expansion") or row.get("source_pool", "").startswith("constrained")


def is_four_mut(row: Dict[str, str]) -> bool:
    return as_int(row.get("mutation_count")) >= 4


def is_neutral_boundary(row: Dict[str, str]) -> bool:
    return row.get("candidate_quality_tier") == "neutral_boundary_or_high_risk"


def contains_ay111h(row: Dict[str, str]) -> bool:
    return "AY111H" in (row.get("his_seed_set", "") + ";" + row.get("mutation_list", ""))


def contains_ag102h_or_av105h(row: Dict[str, str]) -> bool:
    text = row.get("his_seed_set", "") + ";" + row.get("mutation_list", "")
    return "AG102H" in text or "AV105H" in text


def hard_reject_reason(row: Dict[str, str]) -> Optional[str]:
    if not row.get("canonical_unique_key") and not row.get("canonical_sequence_hash_full"):
        return "missing_canonical_unique_key"
    if not row.get("sequence"):
        return "missing_sequence"
    if mutation_introduces_cys(row):
        return "new_cys_mutation"
    target = row.get("target")
    if target == "Ab_sdAb" and contains_ag102h_or_av105h(row):
        return "sdab_AG102H_or_AV105H_main_block"
    for field in ("hard_filter_status", "forbidden_pair_status", "buildability_light_status"):
        value = str(row.get(field, "")).lower()
        if "fail" in value or "violation" in value:
            return f"{field}_not_pass"
    return None


def cap_value(quota: int, target: str, cap_name: str, hardness: str) -> Optional[int]:
    fraction = CAPS[target].get(f"{cap_name}_{hardness}")
    if fraction is None:
        return None
    return max(1, int(math.floor(quota * float(fraction))))


def planning_score(row: Dict[str, str]) -> Dict[str, float]:
    tier = row.get("template_match_tier", "")
    quality = row.get("candidate_quality_tier", "")
    status = row.get("assigned_template_status", "")
    source = row.get("source_pool", "")
    template_component = TIER_SCORE.get(tier, 0.0)
    quality_component = QUALITY_SCORE.get(quality, 0.0)
    status_component = STATUS_SCORE.get(status, 0.0)
    source_component = SOURCE_SCORE.get(source, 0.0)

    mechanism = 0.0
    mechanism += 0.8 * as_float(row.get("hit_likelihood_score_v0"))
    mechanism += 0.6 * as_float(row.get("acidic_release_support_score"))
    mechanism += 0.6 * as_float(row.get("neutral_retention_score"))
    mechanism += 0.3 * max(0.0, 1.0 - as_float(row.get("mpnn_score_percentile_within_target"), 0.5))

    risk = 0.0
    risk += 0.7 * as_float(row.get("global_weakening_risk_score"))
    risk += 0.5 * as_float(row.get("display_or_expression_risk_score"))
    risk += 0.4 * as_float(row.get("glycan_or_epitope_risk_score"))
    risk += 0.12 * max(0, as_int(row.get("mutation_count")) - 2)
    if is_four_mut(row):
        risk += 0.35
    if is_neutral_boundary(row):
        risk += 0.9

    diversity_penalty = 0.0
    if is_seed_only(row):
        diversity_penalty += 1.0
    if is_local_expansion(row):
        diversity_penalty += 0.35

    total = template_component + quality_component + status_component + source_component + mechanism - risk - diversity_penalty
    total += as_float(row.get("expanded_pool_score")) * 0.02
    total += as_float(row.get("tier1_rank_score")) * 0.001
    return {
        "final_planning_score": total,
        "score_template_component": template_component,
        "score_quality_component": quality_component,
        "score_status_component": status_component,
        "score_source_component": source_component,
        "score_mechanism_component": mechanism,
        "score_risk_penalty": risk,
        "score_diversity_penalty": diversity_penalty,
    }


def add_features(rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in rows:
        enriched: Dict[str, Any] = dict(row)
        score_parts = planning_score(row)
        enriched.update({k: round(v, 6) for k, v in score_parts.items()})
        enriched["is_local_expansion"] = is_local_expansion(row)
        enriched["is_seed_only"] = is_seed_only(row)
        enriched["is_four_mut"] = is_four_mut(row)
        enriched["is_neutral_boundary_or_high_risk"] = is_neutral_boundary(row)
        enriched["contains_AY111H"] = contains_ay111h(row)
        enriched["contains_AG102H_or_AV105H"] = contains_ag102h_or_av105h(row)
        out.append(enriched)
    return out


def check_cap(
    row: Dict[str, Any],
    quota: int,
    target: str,
    selected_counts: Dict[str, Counter],
    hardness: str,
) -> Optional[str]:
    seed = norm_key(row.get("his_seed_set"), "missing_seed")
    cluster = norm_key(row.get("near_duplicate_cluster_id"), "missing_cluster")
    seed_cap = cap_value(quota, target, "top_seed", hardness)
    cluster_cap = cap_value(quota, target, "top_cluster", hardness)
    if seed_cap is not None and selected_counts["seed"][seed] >= seed_cap:
        return "seed_cap"
    if cluster_cap is not None and selected_counts["cluster"][cluster] >= cluster_cap:
        return "cluster_cap"
    if row["is_four_mut"]:
        cap = cap_value(quota, target, "four_mut", hardness)
        if cap is not None and selected_counts["four_mut"]["yes"] >= cap:
            return "four_mut_cap"
    if row["is_local_expansion"]:
        cap = cap_value(quota, target, "local_expansion", hardness)
        if cap is not None and selected_counts["local_expansion"]["yes"] >= cap:
            return "local_expansion_cap"
    if row["is_seed_only"]:
        cap = cap_value(quota, target, "seed_only", hardness)
        if cap is not None and selected_counts["seed_only"]["yes"] >= cap:
            return "seed_only_cap"
    if row["is_neutral_boundary_or_high_risk"]:
        cap = cap_value(quota, target, "neutral_boundary", hardness)
        if cap is not None and selected_counts["neutral_boundary"]["yes"] >= cap:
            return "neutral_boundary_or_high_risk_cap"
    if target == "Ab_sdAb" and row["contains_AY111H"]:
        cap = cap_value(quota, target, "ay111h", hardness)
        if cap is not None and selected_counts["ay111h"]["yes"] >= cap:
            return "AY111H_containing_cap"
    return None


def increment_counts(row: Dict[str, Any], counts: Dict[str, Counter]) -> None:
    seed = norm_key(row.get("his_seed_set"), "missing_seed")
    cluster = norm_key(row.get("near_duplicate_cluster_id"), "missing_cluster")
    counts["seed"][seed] += 1
    counts["cluster"][cluster] += 1
    if row["is_four_mut"]:
        counts["four_mut"]["yes"] += 1
    if row["is_local_expansion"]:
        counts["local_expansion"]["yes"] += 1
    if row["is_seed_only"]:
        counts["seed_only"]["yes"] += 1
    if row["is_neutral_boundary_or_high_risk"]:
        counts["neutral_boundary"]["yes"] += 1
    if row["contains_AY111H"]:
        counts["ay111h"]["yes"] += 1


def select_target(
    rows: List[Dict[str, Any]],
    target: str,
    quota: int,
    scenario: str,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    target_rows = [r for r in rows if r.get("target") == target]
    target_rows.sort(
        key=lambda r: (
            as_float(r["final_planning_score"]),
            as_float(r.get("expanded_pool_score")),
            -as_int(r.get("mutation_count")),
            r.get("variant_id", ""),
        ),
        reverse=True,
    )

    selected: List[Dict[str, Any]] = []
    selected_ids: Set[str] = set()
    selected_keys: Set[str] = set()
    counts: Dict[str, Counter] = defaultdict(Counter)
    hard_rejects: Dict[str, str] = {}

    for row in target_rows:
        reason = hard_reject_reason(row)
        if reason:
            hard_rejects[row["variant_id"]] = reason

    def try_phase(hardness: str) -> None:
        nonlocal selected
        for row in target_rows:
            if len(selected) >= quota:
                return
            vid = row["variant_id"]
            if vid in selected_ids or vid in hard_rejects:
                continue
            canonical = row.get("canonical_unique_key") or f"{row.get('target')}|{row.get('canonical_sequence_hash_full')}"
            if canonical in selected_keys:
                continue
            cap_reason = check_cap(row, quota, target, counts, hardness)
            if cap_reason:
                continue
            out = dict(row)
            out["selection_scenario"] = scenario
            out["selection_phase"] = hardness
            out["final_selection_role"] = "main_candidate"
            selected.append(out)
            selected_ids.add(vid)
            selected_keys.add(canonical)
            increment_counts(row, counts)

    try_phase("soft")
    if len(selected) < quota:
        try_phase("hard")

    for rank, row in enumerate(selected, start=1):
        row["selection_rank"] = rank

    excluded: List[Dict[str, Any]] = []
    for row in target_rows:
        if row["variant_id"] in selected_ids:
            continue
        reason = hard_rejects.get(row["variant_id"])
        if not reason:
            reason = check_cap(row, quota, target, counts, "hard") or "not_selected_by_rank_within_quota"
        excluded.append(
            {
                "selection_scenario": scenario,
                "target": target,
                "variant_id": row.get("variant_id"),
                "mutation_list": row.get("mutation_list"),
                "his_seed_set": row.get("his_seed_set"),
                "near_duplicate_cluster_id": row.get("near_duplicate_cluster_id"),
                "template_match_tier": row.get("template_match_tier"),
                "candidate_quality_tier": row.get("candidate_quality_tier"),
                "source_pool": row.get("source_pool"),
                "final_planning_score": row.get("final_planning_score"),
                "excluded_reason": reason,
            }
        )

    audit = audit_target(selected, target, quota, scenario)
    if len(selected) < quota:
        audit["verdict"] = "PATCH" if audit["verdict"] == "PASS" else audit["verdict"]
        audit["reasons_list"].append("selected_count_below_quota")
    if any(row.get("selection_phase") == "hard" for row in selected):
        audit["verdict"] = "PATCH" if audit["verdict"] == "PASS" else audit["verdict"]
        audit["reasons_list"].append("hard_cap_phase_used")
    audit["reasons"] = ";".join(audit["reasons_list"]) if audit["reasons_list"] else "all_checks_pass"
    return selected, excluded, audit


def frac(count: int, denom: int) -> float:
    if not denom:
        return 0.0
    return round(count / denom, 6)


def audit_target(selected: List[Dict[str, Any]], target: str, quota: int, scenario: str) -> Dict[str, Any]:
    count = len(selected)
    seed_counts = Counter(norm_key(r.get("his_seed_set"), "missing_seed") for r in selected)
    cluster_counts = Counter(norm_key(r.get("near_duplicate_cluster_id"), "missing_cluster") for r in selected)
    canonical = [r.get("canonical_unique_key") or f"{r.get('target')}|{r.get('canonical_sequence_hash_full')}" for r in selected]
    duplicate_count = len(canonical) - len(set(canonical))
    top_seed, top_seed_count = seed_counts.most_common(1)[0] if seed_counts else ("", 0)
    top_cluster, top_cluster_count = cluster_counts.most_common(1)[0] if cluster_counts else ("", 0)

    four_mut_count = sum(bool(r["is_four_mut"]) for r in selected)
    local_count = sum(bool(r["is_local_expansion"]) for r in selected)
    seed_only_count = sum(bool(r["is_seed_only"]) for r in selected)
    neutral_count = sum(bool(r["is_neutral_boundary_or_high_risk"]) for r in selected)
    ay111h_count = sum(bool(r["contains_AY111H"]) for r in selected)
    ag_av_count = sum(bool(r["contains_AG102H_or_AV105H"]) for r in selected)
    sdab_complex_diag_count = sum(as_bool(r.get("complex_diagnostic_required")) for r in selected)
    hard_rows = sum(1 for r in selected if hard_reject_reason(r))
    glycan_low_claim_count = 0

    metrics = {
        "selection_scenario": scenario,
        "target": target,
        "selected_count": count,
        "target_quota": quota,
        "canonical_duplicate_count": duplicate_count,
        "top_seed": top_seed,
        "top_seed_count": top_seed_count,
        "top_seed_fraction": frac(top_seed_count, count),
        "top_seed_soft_cap": cap_value(quota, target, "top_seed", "soft"),
        "top_seed_hard_cap": cap_value(quota, target, "top_seed", "hard"),
        "top_cluster": top_cluster,
        "top_cluster_count": top_cluster_count,
        "top_cluster_fraction": frac(top_cluster_count, count),
        "top_cluster_soft_cap": cap_value(quota, target, "top_cluster", "soft"),
        "top_cluster_hard_cap": cap_value(quota, target, "top_cluster", "hard"),
        "four_mut_count": four_mut_count,
        "four_mut_fraction": frac(four_mut_count, count),
        "four_mut_soft_cap": cap_value(quota, target, "four_mut", "soft"),
        "four_mut_hard_cap": cap_value(quota, target, "four_mut", "hard"),
        "local_expansion_count": local_count,
        "local_expansion_fraction": frac(local_count, count),
        "local_expansion_soft_cap": cap_value(quota, target, "local_expansion", "soft"),
        "local_expansion_hard_cap": cap_value(quota, target, "local_expansion", "hard"),
        "seed_only_count": seed_only_count,
        "seed_only_fraction": frac(seed_only_count, count),
        "seed_only_soft_cap": cap_value(quota, target, "seed_only", "soft"),
        "seed_only_hard_cap": cap_value(quota, target, "seed_only", "hard"),
        "neutral_boundary_or_high_risk_count": neutral_count,
        "neutral_boundary_or_high_risk_fraction": frac(neutral_count, count),
        "neutral_boundary_soft_cap": cap_value(quota, target, "neutral_boundary", "soft"),
        "neutral_boundary_hard_cap": cap_value(quota, target, "neutral_boundary", "hard"),
        "AY111H_containing_count": ay111h_count if target == "Ab_sdAb" else "",
        "AY111H_containing_fraction": frac(ay111h_count, count) if target == "Ab_sdAb" else "",
        "AY111H_soft_cap": cap_value(quota, target, "ay111h", "soft"),
        "AY111H_hard_cap": cap_value(quota, target, "ay111h", "hard"),
        "AG102H_AV105H_main_count": ag_av_count if target == "Ab_sdAb" else "",
        "sdAb_complex_diagnostic_required_count": sdab_complex_diag_count if target == "Ab_sdAb" else "",
        "hard_reject_selected_count": hard_rows,
        "glycan_low_risk_claim_count": glycan_low_claim_count,
    }
    verdict = "PASS"
    reasons: List[str] = []
    hard_checks = [
        ("canonical_duplicate_count", duplicate_count, 0),
        ("hard_reject_selected_count", hard_rows, 0),
        ("glycan_low_risk_claim_count", glycan_low_claim_count, 0),
    ]
    if target == "Ab_sdAb":
        hard_checks.append(("AG102H_AV105H_main_count", ag_av_count, 0))
        if sdab_complex_diag_count != count:
            verdict = "FAIL"
            reasons.append("sdAb_complex_diagnostic_required_not_retained")
    for name, value, allowed in hard_checks:
        if value > allowed:
            verdict = "FAIL"
            reasons.append(f"{name}_above_{allowed}")

    cap_checks = [
        ("top_seed", top_seed_count, cap_value(quota, target, "top_seed", "soft"), cap_value(quota, target, "top_seed", "hard")),
        ("top_cluster", top_cluster_count, cap_value(quota, target, "top_cluster", "soft"), cap_value(quota, target, "top_cluster", "hard")),
        ("four_mut", four_mut_count, cap_value(quota, target, "four_mut", "soft"), cap_value(quota, target, "four_mut", "hard")),
        ("local_expansion", local_count, cap_value(quota, target, "local_expansion", "soft"), cap_value(quota, target, "local_expansion", "hard")),
        ("seed_only", seed_only_count, cap_value(quota, target, "seed_only", "soft"), cap_value(quota, target, "seed_only", "hard")),
        ("neutral_boundary", neutral_count, cap_value(quota, target, "neutral_boundary", "soft"), cap_value(quota, target, "neutral_boundary", "hard")),
    ]
    if target == "Ab_sdAb":
        cap_checks.append(("AY111H_containing", ay111h_count, cap_value(quota, target, "ay111h", "soft"), cap_value(quota, target, "ay111h", "hard")))
    for name, value, soft, hard in cap_checks:
        if hard is not None and value > hard:
            verdict = "FAIL"
            reasons.append(f"{name}_hard_cap_exceeded")
        elif soft is not None and value > soft and verdict == "PASS":
            verdict = "PATCH"
            reasons.append(f"{name}_soft_cap_exceeded")

    if count != quota and verdict == "PASS":
        verdict = "PATCH"
        reasons.append("selected_count_not_equal_quota")
    metrics["verdict"] = verdict
    metrics["reasons"] = ";".join(reasons) if reasons else "all_checks_pass"
    metrics["reasons_list"] = reasons
    return metrics


def select_scenario(rows: List[Dict[str, Any]], name: str, scenario: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    selected_all: List[Dict[str, Any]] = []
    excluded_all: List[Dict[str, Any]] = []
    audits: List[Dict[str, Any]] = []
    for target, quota in scenario["quotas"].items():
        selected, excluded, audit = select_target(rows, target, int(quota), name)
        selected_all.extend(selected)
        excluded_all.extend(excluded)
        audits.append(audit)
    selected_all.sort(key=lambda r: (r["target"], as_int(r["selection_rank"])))
    return selected_all, excluded_all, audits


def counter_rows(rows: List[Dict[str, Any]], group_fields: List[str], scenario: Optional[str] = None) -> List[Dict[str, Any]]:
    counter: Counter[Tuple[str, ...]] = Counter()
    for row in rows:
        key = tuple(str(row.get(field, "")) for field in group_fields)
        counter[key] += 1
    out = []
    for key, count in sorted(counter.items(), key=lambda item: (item[0], item[1])):
        rec = {field: value for field, value in zip(group_fields, key)}
        if scenario is not None:
            rec = {"selection_scenario": scenario, **rec}
        rec["count"] = count
        out.append(rec)
    return out


def md_table(rows: List[Dict[str, Any]], fields: List[str]) -> str:
    if not rows:
        return "\n".join(["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"])
    lines = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(f, "")) for f in fields) + " |")
    return "\n".join(lines)


def write_manifest(input_rows: int) -> None:
    text = f"""stage: final_candidate_pool_planning
status: draft_planning_only
inputs:
  expanded_pool_v2:
    path: {INPUT_POOL.relative_to(ROOT)}
    rows: {input_rows}
    sha256: {sha256_file(INPUT_POOL)}
  expanded_pool_v2_audit:
    path: {INPUT_AUDIT.relative_to(ROOT)}
    sha256: {sha256_file(INPUT_AUDIT)}
identity:
  canonical_unique_key: target + canonical_sequence_hash_full
locks:
  final_synthesis_selection: locked
  wet_lab_library_construction: locked
  new_structure_compute: locked
  md: locked
  sdAb_structure_confirmed_upgrade: locked
  glycan_low_risk_claim: locked
"""
    (OUT_DIR / "final_candidate_input_manifest.yaml").write_text(text, encoding="utf-8")


def write_config() -> None:
    lines = [
        "stage: final_candidate_pool_planning",
        "status: draft_only_until_manual_review",
        "default_scenarios:",
    ]
    for name, scenario in DEFAULT_SCENARIOS.items():
        lines.append(f"  {name}:")
        lines.append(f"    total: {scenario['total']}")
        lines.append("    quotas:")
        for target, quota in scenario["quotas"].items():
            lines.append(f"      {target}: {quota}")
    lines.append("caps:")
    for target, caps in CAPS.items():
        lines.append(f"  {target}:")
        for key, value in caps.items():
            lines.append(f"    {key}: {value}")
    (OUT_DIR / "final_candidate_selection_config.yaml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_control_panel() -> None:
    if not CONTROL_SOURCE.exists():
        write_csv(OUT_DIR / "final_candidate_pool_control_anchor_panel.csv", [])
        return
    rows = read_csv(CONTROL_SOURCE)
    for row in rows:
        row["control_panel_source"] = str(CONTROL_SOURCE.relative_to(ROOT))
        row["final_selection_role"] = "control_anchor_panel"
    write_csv(OUT_DIR / "final_candidate_pool_control_anchor_panel.csv", rows)


def write_audit_report(
    default_audits: List[Dict[str, Any]],
    alt_audits: List[Dict[str, Any]],
    selected_by_scenario: Dict[str, List[Dict[str, Any]]],
) -> str:
    all_default_verdicts = [row["verdict"] for row in default_audits]
    if any(v == "FAIL" for v in all_default_verdicts):
        overall = "FINAL_CANDIDATE_POOL_DRAFTS_BUILT__AUDIT_FAIL__FINAL_SELECTION_LOCKED"
    elif any(v == "PATCH" for v in all_default_verdicts):
        overall = "FINAL_CANDIDATE_POOL_DRAFTS_BUILT__PATCH_REVIEW_REQUIRED__FINAL_SELECTION_LOCKED"
    else:
        overall = "FINAL_CANDIDATE_POOL_DRAFTS_BUILT__AUDIT_PASS__FINAL_SELECTION_LOCKED"

    audit_fields = [
        "selection_scenario",
        "target",
        "selected_count",
        "target_quota",
        "top_seed",
        "top_seed_count",
        "top_seed_fraction",
        "top_cluster",
        "top_cluster_count",
        "top_cluster_fraction",
        "four_mut_count",
        "four_mut_fraction",
        "local_expansion_count",
        "local_expansion_fraction",
        "seed_only_count",
        "seed_only_fraction",
        "neutral_boundary_or_high_risk_count",
        "neutral_boundary_or_high_risk_fraction",
        "AY111H_containing_count",
        "AY111H_containing_fraction",
        "AG102H_AV105H_main_count",
        "sdAb_complex_diagnostic_required_count",
        "canonical_duplicate_count",
        "verdict",
        "reasons",
    ]
    source_rows = [
        {
            "input": "expanded_tier2_candidate_pool_v2.csv",
            "rows": 40000,
            "sha256": sha256_file(INPUT_POOL),
        },
        {
            "input": "expanded_pool_v2_audit_report.md",
            "rows": "",
            "sha256": sha256_file(INPUT_AUDIT),
        },
    ]
    target_counts = []
    for name, selected in selected_by_scenario.items():
        for target, count in Counter(r["target"] for r in selected).items():
            target_counts.append({"selection_scenario": name, "target": target, "count": count})

    text = f"""# Final Candidate Pool Planning Audit

## Executive Summary

Verdict: `{overall}`.

This stage converts the validated 20K-per-target expanded Tier2 mother pools into 10K and 15K final candidate-pool drafts. It did not run AF3, SimpleFold, PyRosetta, FoldX, MD, explicit glycan modeling, or wet-lab synthesis selection.

The outputs are draft candidate pools for manual review. Final synthesis-ready selection remains locked.

## Source Inputs

{md_table(source_rows, ["input", "rows", "sha256"])}

## Default Draft Audit

{md_table(default_audits, audit_fields)}

## Alternative Allocation Feasibility

{md_table(alt_audits, audit_fields)}

## Draft Target Counts

{md_table(target_counts, ["selection_scenario", "target", "count"])}

## Interpretation

- `Ab_1E62` remains the primary branch.
- `Ab_sdAb` remains a secondary branch; selected sdAb rows must retain `complex_diagnostic_required`.
- A/B/local-expansion and C seed-only rows remain explicitly labeled.
- Glycan status remains unchecked; this report must not be cited as low glycan risk.
- The 10K and 15K tables are draft candidate pools, not synthesis-ready libraries.

## Next Gate

Allowed now:

- Review the 10K and 15K draft pools.
- Compare default, balanced, and sdAb-heavy target allocation feasibility.
- Decide whether to PATCH selection caps, request targeted refill, or manually unlock final library planning.

Still locked:

- Final synthesis-ready 10K / 15K library selection.
- Wet-lab order list generation.
- New AF3 / SimpleFold / PyRosetta / FoldX compute.
- MD or constant-pH MD.
- Treating sdAb as structure-confirmed.
- Reporting glycan low risk.

## Output Files

- `results/initial_design_generation/final_candidate_pool_planning/final_candidate_pool_10k_draft.csv`
- `results/initial_design_generation/final_candidate_pool_planning/final_candidate_pool_15k_draft.csv`
- `results/initial_design_generation/final_candidate_pool_planning/final_candidate_pool_10k_audit.md`
- `results/initial_design_generation/final_candidate_pool_planning/final_candidate_pool_15k_audit.md`
- `results/initial_design_generation/final_candidate_pool_planning/final_candidate_pool_planning_audit_report.md`
"""
    (OUT_DIR / "final_candidate_pool_planning_audit_report.md").write_text(text, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(text, encoding="utf-8")
    return overall


def write_scenario_audit(name: str, audits: List[Dict[str, Any]]) -> None:
    fields = [
        "selection_scenario",
        "target",
        "selected_count",
        "target_quota",
        "top_seed",
        "top_seed_count",
        "top_seed_fraction",
        "top_cluster",
        "top_cluster_count",
        "top_cluster_fraction",
        "four_mut_count",
        "four_mut_fraction",
        "local_expansion_count",
        "local_expansion_fraction",
        "seed_only_count",
        "seed_only_fraction",
        "verdict",
        "reasons",
    ]
    text = f"# Final Candidate Pool {name} Audit\n\n{md_table(audits, fields)}\n"
    (OUT_DIR / f"final_candidate_pool_{name}_audit.md").write_text(text, encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = add_features(read_csv(INPUT_POOL))
    original_fields = list(rows[0].keys()) if rows else []
    draft_fields = original_fields + [field for field in EXTRA_FIELDS if field not in original_fields]
    write_manifest(len(rows))
    write_config()
    copy_control_panel()

    selected_by_scenario: Dict[str, List[Dict[str, Any]]] = {}
    default_audits: List[Dict[str, Any]] = []
    excluded_all: List[Dict[str, Any]] = []
    for name, scenario in DEFAULT_SCENARIOS.items():
        selected, excluded, audits = select_scenario(rows, name, scenario)
        selected_by_scenario[name] = selected
        default_audits.extend(audits)
        excluded_all.extend(excluded)
        write_csv(OUT_DIR / f"final_candidate_pool_{name}_draft.csv", selected, draft_fields)
        write_scenario_audit(name, audits)

    alt_audits: List[Dict[str, Any]] = []
    for name, scenario in ALT_SCENARIOS.items():
        selected, excluded, audits = select_scenario(rows, name, scenario)
        alt_audits.extend(audits)
        excluded_all.extend(excluded[:1000])

    write_csv(OUT_DIR / "final_candidate_pool_excluded_candidates.csv", excluded_all)
    write_csv(OUT_DIR / "final_candidate_pool_allocation_sensitivity.csv", alt_audits)
    write_csv(OUT_DIR / "final_candidate_pool_audit_checks.csv", default_audits)

    combined_selected = []
    for selected in selected_by_scenario.values():
        combined_selected.extend(selected)

    write_csv(OUT_DIR / "final_candidate_pool_by_target.csv", counter_rows(combined_selected, ["selection_scenario", "target"]))
    write_csv(OUT_DIR / "final_candidate_pool_by_seed.csv", counter_rows(combined_selected, ["selection_scenario", "target", "his_seed_set"]))
    write_csv(OUT_DIR / "final_candidate_pool_by_cluster.csv", counter_rows(combined_selected, ["selection_scenario", "target", "near_duplicate_cluster_id"]))
    write_csv(OUT_DIR / "final_candidate_pool_by_match_tier.csv", counter_rows(combined_selected, ["selection_scenario", "target", "template_match_tier", "expanded_pool_role"]))
    write_csv(OUT_DIR / "final_candidate_pool_by_source.csv", counter_rows(combined_selected, ["selection_scenario", "target", "source_pool", "template_match_tier"]))
    write_csv(OUT_DIR / "final_candidate_pool_by_mutation_count.csv", counter_rows(combined_selected, ["selection_scenario", "target", "mutation_count"]))

    overall = write_audit_report(default_audits, alt_audits, selected_by_scenario)
    next_gate = f"""# Final Candidate Pool Planning Next Gate

Overall verdict: `{overall}`.

Review the draft pools before any final synthesis selection. New structural compute, MD, wet-lab order generation, sdAb structure-confirmed upgrade, and glycan low-risk claims remain locked.
"""
    (OUT_DIR / "final_candidate_pool_next_gate_report.md").write_text(next_gate, encoding="utf-8")
    print(overall)
    print(OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
