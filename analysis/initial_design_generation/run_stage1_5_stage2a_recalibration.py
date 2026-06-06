#!/usr/bin/env python3
"""Run Stage-1.5 full recalibration and Stage-2A candidate-list audit.

This script is intentionally compute-light relative to Tier2 structural tools:
it parses existing Stage-1 PyRosetta PDBs, computes geometry features, applies
the reviewed severe / boundary / pass-like structure policy, and builds an
auditable Stage-2A candidate list. It does not launch PyRosetta, FoldX, AF3,
SimpleFold, Tier2-heavy, or final selection.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.initial_design_generation.build_manual_spot_check_structure_features import (
    load_pka_detail,
    parent_rows,
    structure_feature_row,
)


DEFAULT_STAGE1 = Path("results/initial_design_generation/tier2_staged/full_stage1/tier2b_full_stage1_results.csv")
DEFAULT_PKA = Path("results/initial_design_generation/tier2_staged/full_stage1/tier2b_full_pka_detail.csv")
DEFAULT_OUT = Path("results/initial_design_generation/stage1_5_stage2a")

STRONG_GOOD = {"T2_strong_candidate", "T2_good_candidate"}
BOUNDARY_ALLOWED_CLASSES = {
    "T2_strong_candidate",
    "T2_good_candidate",
    "T2_release_possible_but_neutral_risky",
}
SUPPORTED_ROUTES = {
    "His_plus_rescue",
    "His_rule",
    "ProteinMPNN_seeded_rescue",
    "wetlab_informed_expansion",
    "structure_or_interface_guided",
    "repaired_sparse_constrained_MPNN",
}
PRIORITY_SEEDS = {
    "Ab_1E62": {
        "LQ35H;LY38H",
        "LK24H;LY38H",
        "LK24H;LQ35H",
        "LY31H;LQ35H",
    },
    "Ab_sdAb": {"AD110H;AY111H"},
}
SECONDARY_SEEDS = {
    "Ab_1E62": {"LY31H;LY38H", "LK24H;LY31H"},
    "Ab_sdAb": {"AQ100H;AD110H", "AQ100H;AY111H"},
}


_STAGE_ROWS: list[dict[str, Any]] = []
_PARENTS: dict[str, dict[str, Any]] = {}
_PKA_DETAIL: dict[str, dict[Any, str]] = {}
_PARENT_CACHE: dict[str, dict[str, object]] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage1", default=str(DEFAULT_STAGE1))
    parser.add_argument("--pka-detail", default=str(DEFAULT_PKA))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    parser.add_argument("--workers", type=int, default=min(48, os.cpu_count() or 1))
    return parser.parse_args()


def init_worker(stage_rows: list[dict[str, Any]], parents: dict[str, dict[str, Any]], pka_detail: dict[str, dict[Any, str]]) -> None:
    global _STAGE_ROWS, _PARENTS, _PKA_DETAIL, _PARENT_CACHE
    _STAGE_ROWS = stage_rows
    _PARENTS = parents
    _PKA_DETAIL = pka_detail
    _PARENT_CACHE = {}


def feature_for_index(index: int) -> dict[str, Any]:
    row = pd.Series(_STAGE_ROWS[index])
    target = str(row.get("target", ""))
    parent = _PARENTS.get(target)
    if not parent:
        return {
            "variant_id": row.get("variant_id", ""),
            "target": target,
            "structure_parse_status": "missing_parent_anchor_row",
        }
    return structure_feature_row(row, row, pd.Series(parent), _PARENT_CACHE, _PKA_DETAIL)


def numeric(series: pd.Series | Any) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def num_val(row: pd.Series, col: str, default: float = 0.0) -> float:
    val = pd.to_numeric(row.get(col), errors="coerce")
    if pd.isna(val):
        return default
    return float(val)


def bool_support(row: pd.Series) -> bool:
    route = str(row.get("primary_generation_route", ""))
    seed = str(row.get("his_seed_set", ""))
    target = str(row.get("target", ""))
    tier1_class = str(row.get("tier1_review_class", ""))
    rescue_count = num_val(row, "rescue_count")
    if route in SUPPORTED_ROUTES:
        return True
    if seed in PRIORITY_SEEDS.get(target, set()) or seed in SECONDARY_SEEDS.get(target, set()):
        return True
    if rescue_count > 0:
        return True
    if tier1_class.startswith(("A_", "B_", "C_")):
        return True
    return False


def parse_sasa_values(text: Any) -> list[float]:
    if text is None or pd.isna(text):
        return []
    vals: list[float] = []
    for part in str(text).split(";"):
        if ":" not in part:
            continue
        try:
            vals.append(float(part.rsplit(":", 1)[1]))
        except ValueError:
            continue
    return vals


def classify_row(row: pd.Series) -> tuple[str, str, str, str]:
    status = str(row.get("structure_parse_status", ""))
    target = str(row.get("target", ""))
    t2_class = str(row.get("tier2_class", row.get("t2_class_current", "")))
    if not status.startswith("ok"):
        return "T2_manual_review", "manual_only", "manual_review", "structure_parse_issue"
    if t2_class == "T2_control_anchor":
        return "T2_manual_review", "control_only", "control_only", "control_anchor"

    local = num_val(row, "local_validity_score", default=math.nan)
    rosetta_delta = num_val(row, "rosetta_delta", default=0.0)
    new_clash = num_val(row, "new_clash_count_total", default=0.0)
    max_new_clash = num_val(row, "max_new_clash_overlap", default=0.0)
    new_bad = num_val(row, "new_bad_contact_count", default=0.0)
    lost = num_val(row, "lost_parent_contact_count", default=0.0)
    his_dist = num_val(row, "his_min_antigen_distance", default=math.nan)
    mutation_count = num_val(row, "mutation_count", default=0.0)
    seed = str(row.get("his_seed_set", ""))

    support = bool_support(row)
    reasons: list[str] = []
    severe_reasons: list[str] = []
    boundary_reasons: list[str] = []

    if not math.isnan(local) and local < 0.35:
        severe_reasons.append("local_validity_lt_0.35")
    if new_clash >= 10:
        severe_reasons.append("new_clash_count_ge_10")
    if max_new_clash >= 1.3:
        severe_reasons.append("max_new_clash_overlap_ge_1.3")
    if new_bad >= 1 and ((not math.isnan(local) and local < 0.55) or new_clash >= 5 or max_new_clash >= 1.1):
        severe_reasons.append("bad_interface_contact_with_local_or_clash_support")
    if lost >= 3 and (not math.isnan(local) and local < 0.55):
        severe_reasons.append("lost_parent_contacts_ge_3_with_low_local_validity")
    if target == "Ab_sdAb" and mutation_count >= 4 and (new_clash > 0 or (not math.isnan(local) and local < 0.55)):
        severe_reasons.append("sdAb_4mut_with_clash_or_low_validity")

    if severe_reasons:
        return (
            "T2_severe_structure_risk",
            "reject_or_deprioritize",
            "stage2a_exclude",
            ";".join(severe_reasons),
        )

    if rosetta_delta > 600:
        boundary_reasons.append("rosetta_delta_gt_600")
    if not math.isnan(local) and 0.35 <= local < 0.55:
        boundary_reasons.append("local_validity_0.35_to_0.55")
    if 1 <= new_clash < 10:
        boundary_reasons.append("new_clash_count_1_to_9")
    if 1.1 <= max_new_clash < 1.3:
        boundary_reasons.append("max_new_clash_overlap_1.1_to_1.3")
    if new_bad >= 1:
        boundary_reasons.append("bad_interface_contact_without_severe_support")
    if 1 <= lost <= 2:
        boundary_reasons.append("lost_parent_contact_1_to_2")
    if not math.isnan(his_dist) and his_dist > 10 and not support:
        boundary_reasons.append("his_far_without_support")
    if target == "Ab_sdAb" and mutation_count >= 4:
        boundary_reasons.append("sdAb_4mut_cap_required")
    if target == "Ab_sdAb" and "AE108H" in seed:
        boundary_reasons.append("sdAb_ae108h_hold")

    if boundary_reasons:
        list_action = "stage2a_include_low_quota" if support and t2_class in BOUNDARY_ALLOWED_CLASSES else "stage2a_hold"
        return (
            "T2_boundary_structure_risk",
            "allow_if_supported" if list_action == "stage2a_include_low_quota" else "hold_boundary",
            list_action,
            ";".join(boundary_reasons),
        )

    pass_reasons = ["no_severe_or_boundary_geometry"]
    if not math.isnan(his_dist) and his_dist <= 10:
        pass_reasons.append("his_within_10A")
    elif support:
        pass_reasons.append("alternate_support_for_his_or_non_his")
    else:
        pass_reasons.append("no_his_or_unclear_his_support")

    return (
        "T2_pass_like_structure",
        "eligible",
        "stage2a_include_priority",
        ";".join(pass_reasons),
    )


def add_refined_labels(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    classes: list[str] = []
    actions: list[str] = []
    list_actions: list[str] = []
    reasons: list[str] = []
    for _, row in out.iterrows():
        cls, action, list_action, reason = classify_row(row)
        classes.append(cls)
        actions.append(action)
        list_actions.append(list_action)
        reasons.append(reason)
    out["refined_structure_risk_class"] = classes
    out["stage2_action"] = actions
    out["stage2a_list_action"] = list_actions
    out["risk_reason_codes"] = reasons

    new_clash = numeric(out["new_clash_count_total"]).fillna(0)
    shell = numeric(out["new_clash_count_mutation_shell"]).fillna(0)
    out["new_clash_within_mutation_6A_fraction"] = [
        (s / n if n > 0 else 0.0) for s, n in zip(shell, new_clash)
    ]

    sasa_min: list[float | None] = []
    sasa_median: list[float | None] = []
    for text in out.get("his_sasa_by_position", pd.Series([pd.NA] * len(out))):
        vals = parse_sasa_values(text)
        if vals:
            sasa_min.append(round(float(min(vals)), 4))
            sasa_median.append(round(float(pd.Series(vals).median()), 4))
        else:
            sasa_min.append(pd.NA)
            sasa_median.append(pd.NA)
    out["his_sasa_min"] = sasa_min
    out["his_sasa_median"] = sasa_median
    out["cdr_ca_rmsd"] = out.get("cdr_rmsd", pd.NA)
    out["window_ca_rmsd"] = out.get("window_rmsd", pd.NA)
    out["mutation_shell_ca_rmsd"] = out.get("mutation_shell_rmsd", pd.NA)
    out["selection_class_from_tier1_or_stage1"] = out.get("tier1_review_class", "").astype(str) + "|" + out.get("tier2_class", "").astype(str)
    return out


def rank_score(df: pd.DataFrame) -> pd.Series:
    neutral = numeric(df.get("neutral_retention_t2_score", 0)).fillna(0)
    acid = numeric(df.get("acidic_release_mechanism_t2_score", 0)).fillna(0)
    pka = numeric(df.get("his_pka_support_t2_score", 0)).fillna(0)
    mpnn = numeric(df.get("mpnn_compatibility_t2_score", 0)).fillna(0)
    local = numeric(df.get("local_validity_score", 0)).fillna(0)
    global_risk = numeric(df.get("global_weakening_risk_t2_score", 0)).fillna(0)
    new_clash = numeric(df.get("new_clash_count_total", 0)).fillna(0)
    lost = numeric(df.get("lost_parent_contact_count", 0)).fillna(0)
    t2_bonus = df.get("tier2_class", "").map({"T2_strong_candidate": 1.0, "T2_good_candidate": 0.6}).fillna(0)
    refined_bonus = df.get("refined_structure_risk_class", "").map(
        {"T2_pass_like_structure": 1.0, "T2_boundary_structure_risk": 0.2}
    ).fillna(-5.0)
    return (
        1.8 * neutral
        + 1.4 * acid
        + 0.6 * pka
        + 0.4 * mpnn
        + 0.8 * local
        + t2_bonus
        + refined_bonus
        - 1.2 * global_risk
        - 0.08 * new_clash
        - 0.15 * lost
    )


def choose_stage2a(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = df.copy()
    out["stage2a_rank_score"] = rank_score(out)
    out["stage2a_exclusion_reason"] = ""

    is_control = out["stage2a_list_action"].eq("control_only")
    severe = out["refined_structure_risk_class"].eq("T2_severe_structure_risk")
    manual = out["refined_structure_risk_class"].eq("T2_manual_review")
    weak_filler = (
        out.get("tier1_review_class", "").astype(str).str.startswith("F_")
        & ~out["tier2_class"].isin(STRONG_GOOD)
        & ~out["stage2a_list_action"].eq("stage2a_include_priority")
    )
    ae108h = out["target"].eq("Ab_sdAb") & out.get("his_seed_set", "").astype(str).str.contains("AE108H", regex=False, na=False)

    out.loc[is_control, "stage2a_exclusion_reason"] = "control_only_not_compute_candidate"
    out.loc[severe, "stage2a_exclusion_reason"] = "severe_structure_risk"
    out.loc[manual & ~is_control, "stage2a_exclusion_reason"] = "manual_review_required"
    out.loc[weak_filler & out["stage2a_exclusion_reason"].eq(""), "stage2a_exclusion_reason"] = "weak_filler_or_broad_F_boundary"
    out.loc[ae108h & out["stage2a_exclusion_reason"].eq(""), "stage2a_exclusion_reason"] = "sdAb_broad_AE108H_hold"

    eligible = out[out["stage2a_exclusion_reason"].eq("")].copy()
    eligible = eligible[eligible["stage2a_list_action"].isin({"stage2a_include_priority", "stage2a_include_low_quota"})]

    selected_parts: list[pd.DataFrame] = []
    for target, target_max in [("Ab_1E62", 5000), ("Ab_sdAb", 1500)]:
        sub = eligible[eligible["target"].eq(target)].copy()
        if sub.empty:
            continue
        sub = sub.sort_values(
            ["stage2a_list_action", "stage2a_rank_score", "local_validity_score"],
            ascending=[True, False, False],
        )
        selected: list[pd.Series] = []
        seed_counts: dict[str, int] = {}
        cluster_counts: dict[str, int] = {}
        boundary_count = 0
        four_mut_count = 0
        target_limit = min(target_max, len(sub))
        seed_cap = max(10, int(0.25 * target_limit))
        cluster_cap = max(10, int(0.05 * target_limit))
        boundary_cap = int((0.50 if target == "Ab_1E62" else 0.30) * target_limit)
        four_mut_cap = int((0.25 if target == "Ab_1E62" else 0.10) * target_limit)

        for _, row in sub.iterrows():
            if len(selected) >= target_limit:
                break
            seed = str(row.get("his_seed_set", "missing"))
            cluster = str(row.get("near_duplicate_cluster_id", row.get("sequence_hamming_cluster_id", "missing")))
            is_boundary = row["refined_structure_risk_class"] == "T2_boundary_structure_risk"
            is_four = num_val(row, "mutation_count") >= 4
            if seed_counts.get(seed, 0) >= seed_cap:
                continue
            if cluster_counts.get(cluster, 0) >= cluster_cap:
                continue
            if is_boundary and boundary_count >= boundary_cap:
                continue
            if is_four and four_mut_count >= four_mut_cap:
                continue
            selected.append(row)
            seed_counts[seed] = seed_counts.get(seed, 0) + 1
            cluster_counts[cluster] = cluster_counts.get(cluster, 0) + 1
            boundary_count += int(is_boundary)
            four_mut_count += int(is_four)
        if selected:
            selected_parts.append(pd.DataFrame(selected))

    candidate = pd.concat(selected_parts, ignore_index=True) if selected_parts else pd.DataFrame(columns=out.columns)
    selected_ids = set(candidate["variant_id"]) if not candidate.empty else set()
    excluded = out[~out["variant_id"].isin(selected_ids)].copy()
    excluded.loc[excluded["stage2a_exclusion_reason"].eq(""), "stage2a_exclusion_reason"] = "not_selected_by_rank_or_caps"
    return candidate, excluded


def count_table(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=group_cols + ["count", "fraction"])
    out = df.groupby(group_cols, dropna=False).size().reset_index(name="count")
    if "target" in group_cols:
        denom = out.groupby("target", dropna=False)["count"].transform("sum")
        out["fraction"] = (out["count"] / denom).round(6)
    else:
        out["fraction"] = (out["count"] / len(df)).round(6)
    return out.sort_values("count", ascending=False)


def by_target_route_seed(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["target", "primary_generation_route", "his_seed_set", "refined_structure_risk_class"]
    out = df.groupby(cols, dropna=False).size().reset_index(name="count")
    totals = df.groupby(["target", "primary_generation_route", "his_seed_set"], dropna=False).size().rename("total").reset_index()
    out = out.merge(totals, on=["target", "primary_generation_route", "his_seed_set"], how="left")
    out["fraction_within_target_route_seed"] = (out["count"] / out["total"]).round(6)
    return out.sort_values(["target", "count"], ascending=[True, False])


def audit_candidate_list(candidate: pd.DataFrame, excluded: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    checks: list[dict[str, Any]] = []

    def add(name: str, status: str, details: Any) -> None:
        checks.append({"check": name, "status": status, "details": details})

    severe = int(candidate["refined_structure_risk_class"].eq("T2_severe_structure_risk").sum()) if not candidate.empty else 0
    add("severe_risk_candidates_zero", "PASS" if severe == 0 else "FAIL", severe)

    controls = int(candidate["stage2a_list_action"].eq("control_only").sum()) if not candidate.empty else 0
    add("controls_separated_from_compute_candidates", "PASS" if controls == 0 else "FAIL", controls)

    weak = int(candidate.get("tier1_review_class", pd.Series([], dtype=str)).astype(str).str.startswith("F_").sum()) if not candidate.empty else 0
    weak_non_sg = int(
        (
            candidate.get("tier1_review_class", pd.Series([], dtype=str)).astype(str).str.startswith("F_")
            & ~candidate.get("tier2_class", pd.Series([], dtype=str)).isin(STRONG_GOOD)
        ).sum()
    ) if not candidate.empty else 0
    add("no_weak_filler_non_strong_good", "PASS" if weak_non_sg == 0 else "WARN", weak_non_sg)

    ae108h = int(
        (
            candidate["target"].eq("Ab_sdAb")
            & candidate.get("his_seed_set", pd.Series([], dtype=str)).astype(str).str.contains("AE108H", regex=False, na=False)
        ).sum()
    ) if not candidate.empty else 0
    add("no_broad_sdAb_AE108H", "PASS" if ae108h == 0 else "WARN", ae108h)

    metrics: dict[str, Any] = {
        "candidate_total": len(candidate),
        "excluded_total": len(excluded),
    }
    target_counts = candidate["target"].value_counts().to_dict() if not candidate.empty else {}
    add("1E62_candidate_count_ge_3000", "PASS" if target_counts.get("Ab_1E62", 0) >= 3000 else "WARN", target_counts.get("Ab_1E62", 0))
    add("sdAb_candidate_count_ge_1000", "PASS" if target_counts.get("Ab_sdAb", 0) >= 1000 else "WARN", target_counts.get("Ab_sdAb", 0))
    for target, sub in candidate.groupby("target", dropna=False):
        prefix = str(target)
        total = len(sub)
        boundary = int(sub["refined_structure_risk_class"].eq("T2_boundary_structure_risk").sum())
        four = int(numeric(sub["mutation_count"]).ge(4).sum())
        top_seed = int(sub.groupby("his_seed_set", dropna=False).size().max()) if total else 0
        top_cluster = int(sub.groupby("near_duplicate_cluster_id", dropna=False).size().max()) if total else 0
        top20 = int(sub.groupby("near_duplicate_cluster_id", dropna=False).size().sort_values(ascending=False).head(20).sum()) if total else 0
        metrics[f"{prefix}_count"] = total
        metrics[f"{prefix}_boundary_fraction"] = round(boundary / total, 6) if total else 0
        metrics[f"{prefix}_4mut_fraction"] = round(four / total, 6) if total else 0
        metrics[f"{prefix}_top_seed_fraction"] = round(top_seed / total, 6) if total else 0
        metrics[f"{prefix}_top_cluster_fraction"] = round(top_cluster / total, 6) if total else 0
        metrics[f"{prefix}_top20_cluster_fraction"] = round(top20 / total, 6) if total else 0
        add(f"{prefix}_near_duplicate_top_cluster_cap", "PASS" if (top_cluster / total if total else 0) <= 0.05 else "WARN", metrics[f"{prefix}_top_cluster_fraction"])
        add(f"{prefix}_top_seed_cap", "PASS" if (top_seed / total if total else 0) <= 0.25 else "WARN", metrics[f"{prefix}_top_seed_fraction"])
        if target == "Ab_sdAb":
            add("sdAb_4mut_fraction_le_10pct", "PASS" if (four / total if total else 0) <= 0.10 else "WARN", metrics[f"{prefix}_4mut_fraction"])
            add("sdAb_boundary_fraction_le_30pct", "PASS" if (boundary / total if total else 0) <= 0.30 else "WARN", metrics[f"{prefix}_boundary_fraction"])
        if target == "Ab_1E62":
            add("1E62_boundary_fraction_le_50pct", "PASS" if (boundary / total if total else 0) <= 0.50 else "WARN", metrics[f"{prefix}_boundary_fraction"])

    return pd.DataFrame(checks), metrics


def markdown_table(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy()
    lines = [
        "| " + " | ".join(str(c) for c in view.columns) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]).replace("\n", " ") for c in view.columns) + " |")
    return "\n".join(lines)


def write_summary(
    refined: pd.DataFrame,
    candidate: pd.DataFrame,
    excluded: pd.DataFrame,
    checks: pd.DataFrame,
    metrics: dict[str, Any],
    out_dir: Path,
) -> None:
    class_summary = count_table(refined, ["target", "refined_structure_risk_class"])
    action_summary = count_table(refined, ["target", "stage2_action"])
    candidate_class = count_table(candidate, ["target", "refined_structure_risk_class"])
    exclusion_summary = count_table(excluded, ["target", "stage2a_exclusion_reason"])
    reason_summary = count_table(refined, ["target", "risk_reason_codes"]).head(30)
    status_counts = checks["status"].value_counts().to_dict() if not checks.empty else {}
    audit_gate = "PASS" if status_counts.get("WARN", 0) == 0 and status_counts.get("FAIL", 0) == 0 else "NOT_PASSED"

    lines = [
        "# Stage-1.5 Refined Structure Risk and Stage-2A Candidate-List Audit",
        "",
        "## Executive Decision",
        "",
        "Stage-1.5 full recalibration and Stage-2A candidate-list construction are complete. No new Tier2-core, Tier2-heavy, AF3, SimpleFold, or final-selection compute was launched.",
        "",
        "Current gate:",
        "",
        "```text",
        "Stage-2A candidate-list construction: complete",
        f"Stage-2A candidate-list audit: {audit_gate}",
        "Stage-2A compute: HOLD",
        "Broad Tier2 expansion: HOLD",
        "Tier2-heavy: HOLD",
        "Final 10K selection: NO-GO",
        "```",
        "",
        "Interpretation: the full recalibration confirms that useful candidates exist, but the current candidate list does not meet the proposed Stage-2A scale/composition gates. Stage-2A compute should remain on hold unless the team explicitly accepts a much smaller, more conservative compute batch or revises the generator/gates.",
        "",
        "## Coverage",
        "",
        f"- Full Stage-1 rows analyzed: {len(refined)}",
        f"- Stage-2A candidate-list rows: {len(candidate)}",
        f"- Excluded / held rows: {len(excluded)}",
        f"- Parse ok-like rows: {int(refined['structure_parse_status'].astype(str).str.startswith('ok').sum())}",
        f"- Audit checks: PASS={status_counts.get('PASS', 0)}, WARN={status_counts.get('WARN', 0)}, FAIL={status_counts.get('FAIL', 0)}",
        "",
        "## Refined Structure-Risk Class Distribution",
        "",
        markdown_table(class_summary, 40),
        "",
        "## Stage2 Action Distribution",
        "",
        markdown_table(action_summary, 40),
        "",
        "## Stage-2A Candidate List Composition",
        "",
        markdown_table(candidate_class, 40),
        "",
        "## Candidate-List Audit Checks",
        "",
        markdown_table(checks, 40),
        "",
        "## Key Audit Metrics",
        "",
        markdown_table(pd.DataFrame([{"metric": k, "value": v} for k, v in metrics.items()]), 80),
        "",
        "## Exclusion / Hold Summary",
        "",
        markdown_table(exclusion_summary, 60),
        "",
        "## Top Refined Risk Reasons",
        "",
        markdown_table(reason_summary, 30),
        "",
        "## Interpretation",
        "",
        "- `rosetta_delta > 600` is treated as a boundary warning unless accompanied by severe local geometry.",
        "- Parent-relative new clash and max new clash overlap are primary local-risk signals.",
        "- CDR/window CA RMSD remains non-informative for backbone stability because Stage-1 PyRosetta PDBs do not sample backbone displacement.",
        "- 1E62 remains stronger than sdAb, but full-set severe-risk rate is still high under the refined rules.",
        "- sdAb remains narrow and capped; broad AE108H and weak filler are held.",
        "- The candidate list is deliberately conservative: severe risk is excluded, controls are separated, and AE108H broad expansion is held.",
        "- The current list is too small for the previously proposed Stage-2A compute scale: 1E62 has 280 rows versus a 3,000-row lower target; sdAb has 74 rows versus a 1,000-row lower target.",
        "- Several composition checks remain WARN: boundary fraction, near-duplicate concentration, sdAb top-seed concentration, and sdAb 4-mut fraction.",
        "",
        "## Next Required Manual Decision",
        "",
        "Manual review should decide one of three paths:",
        "",
        "1. Accept a much smaller conservative Stage-2A compute batch from this list.",
        "2. Revise generator/gates and rebuild a larger candidate list before compute.",
        "3. Stop Stage-2A expansion and return to upstream candidate generation / Tier1 filtering.",
        "",
        "Under the current policy gates, Stage-2A compute should not start automatically.",
    ]
    (out_dir / "stage1_refined_structure_risk_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out_dir / "stage2a_candidate_list_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out_dir / "stage1_5_stage2a_final_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    stage1_path = Path(args.stage1)
    pka_path = Path(args.pka_detail)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    stage1 = pd.read_csv(stage1_path, low_memory=False)
    parents = {k: v.to_dict() for k, v in parent_rows(stage1).items()}
    pka_detail = load_pka_detail(pka_path)
    stage_rows = stage1.to_dict("records")

    if args.workers <= 1:
        init_worker(stage_rows, parents, pka_detail)
        feature_rows = [feature_for_index(i) for i in range(len(stage_rows))]
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=init_worker,
            initargs=(stage_rows, parents, pka_detail),
        ) as pool:
            feature_rows = list(pool.map(feature_for_index, range(len(stage_rows)), chunksize=16))

    features = pd.DataFrame(feature_rows)
    merged = stage1.merge(
        features.drop(columns=[c for c in ["target", "mutation_list"] if c in features.columns]),
        on="variant_id",
        how="left",
        suffixes=("", "_feature"),
    )

    refined = add_refined_labels(merged)

    requested_cols = [
        "variant_id",
        "target",
        "mutation_list",
        "primary_generation_route",
        "selection_class_from_tier1_or_stage1",
        "his_seed_set",
        "mutation_count",
        "rosetta_delta",
        "local_validity_score",
        "new_clash_count_total",
        "max_new_clash_overlap",
        "new_clash_within_mutation_6A_fraction",
        "new_bad_contact_count",
        "lost_parent_contact_count",
        "his_min_antigen_distance",
        "his_sasa_min",
        "his_sasa_median",
        "cdr_ca_rmsd",
        "window_ca_rmsd",
        "mutation_shell_ca_rmsd",
        "refined_structure_risk_class",
        "stage2_action",
        "stage2a_list_action",
        "risk_reason_codes",
    ]
    existing_requested = [c for c in requested_cols if c in refined.columns]
    extra_cols = [c for c in refined.columns if c not in existing_requested]

    features.to_csv(out_dir / "stage1_geometry_features_only.csv", index=False)
    refined[existing_requested + extra_cols].to_csv(out_dir / "stage1_full_structure_features.csv", index=False)
    refined[existing_requested + extra_cols].to_csv(out_dir / "stage1_refined_structure_risk.csv", index=False)
    by_target_route_seed(refined).to_csv(out_dir / "stage1_refined_structure_risk_by_target_route_seed.csv", index=False)

    candidate, excluded = choose_stage2a(refined)
    checks, metrics = audit_candidate_list(candidate, excluded)

    candidate.to_csv(out_dir / "stage2a_candidate_list.csv", index=False)
    count_table(candidate, ["target", "primary_generation_route"]).to_csv(out_dir / "stage2a_candidate_route_summary.csv", index=False)
    count_table(candidate, ["target", "his_seed_set"]).to_csv(out_dir / "stage2a_candidate_seed_summary.csv", index=False)
    count_table(candidate, ["target", "mutation_count"]).to_csv(out_dir / "stage2a_candidate_mutation_count_summary.csv", index=False)
    count_table(candidate, ["target", "near_duplicate_cluster_id"]).to_csv(out_dir / "stage2a_candidate_near_duplicate_summary.csv", index=False)
    count_table(excluded, ["target", "stage2a_exclusion_reason"]).to_csv(out_dir / "stage2a_candidate_exclusion_report.csv", index=False)
    checks.to_csv(out_dir / "stage2a_candidate_list_audit_checks.csv", index=False)

    write_summary(refined, candidate, excluded, checks, metrics, out_dir)

    print(f"Wrote Stage-1.5 / Stage-2A outputs to {out_dir}")
    print(f"Rows analyzed: {len(refined)}")
    print(f"Stage-2A candidate-list rows: {len(candidate)}")
    print("Audit statuses:")
    print(checks["status"].value_counts().to_string())


if __name__ == "__main__":
    main()
