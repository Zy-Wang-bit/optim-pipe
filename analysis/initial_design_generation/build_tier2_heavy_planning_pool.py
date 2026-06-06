#!/usr/bin/env python3
"""Build capped targeted Tier2-heavy planning pools from Stage-2A results."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(".")
INPUT_CSV = ROOT / "results/initial_design_generation/stage2a_compute/stage2a_compute_results.csv"
OUT_DIR = ROOT / "results/initial_design_generation/tier2_heavy_planning"
CURRENT_STAGE_REPORT = ROOT / ".tasks/active/initial-design-generation/current_stage_report.md"
PROGRESS = ROOT / ".tasks/active/initial-design-generation/progress.md"


MAIN_CLASSES = {
    "T2A_high_confidence": 0,
    "T2A_good": 1,
    "T2A_supported_boundary_confirmed": 2,
}


@dataclass(frozen=True)
class TargetPolicy:
    target: str
    output_name: str
    default_size: int
    clean_min: int
    clean_max: int
    lower_clean_min: int
    lower_clean_max: int
    cluster_default_cap: int
    cluster_high_conf_exception: int
    top_cluster_absolute_cap: int
    top5_cluster_fraction_cap: float
    boundary_fraction_cap: float
    control_target: int


POLICIES = {
    "Ab_1E62": TargetPolicy(
        target="Ab_1E62",
        output_name="1E62",
        default_size=100,
        clean_min=80,
        clean_max=100,
        lower_clean_min=60,
        lower_clean_max=80,
        cluster_default_cap=4,
        cluster_high_conf_exception=5,
        top_cluster_absolute_cap=12,
        top5_cluster_fraction_cap=0.60,
        boundary_fraction_cap=0.45,
        control_target=20,
    ),
    "Ab_sdAb": TargetPolicy(
        target="Ab_sdAb",
        output_name="sdAb",
        default_size=120,
        clean_min=100,
        clean_max=120,
        lower_clean_min=80,
        lower_clean_max=100,
        cluster_default_cap=4,
        cluster_high_conf_exception=5,
        top_cluster_absolute_cap=15,
        top5_cluster_fraction_cap=0.60,
        boundary_fraction_cap=0.50,
        control_target=20,
    ),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def contains_seed(seed: object, token: str) -> bool:
    return token in str(seed or "")


def contains_mutation(mutation_list: object, token: str) -> bool:
    return token in str(mutation_list or "")


def class_rank(value: object) -> int:
    return MAIN_CLASSES.get(str(value), 99)


def source_rank(value: object) -> int:
    return 0 if str(value) == "promote_to_heavy_candidate_pool" else 1


def prepare_rank(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_class_rank"] = out["stage2a_final_class"].map(class_rank)
    out["_source_rank"] = out["stage2a_action"].map(source_rank)
    for col in [
        "stage2a_structure_score",
        "stage2a_pH_mechanism_support",
        "stage2a_local_validity",
        "stage2a_new_clash_count",
        "rosetta_delta_score",
    ]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.sort_values(
        [
            "_class_rank",
            "_source_rank",
            "stage2a_structure_score",
            "stage2a_pH_mechanism_support",
            "stage2a_local_validity",
            "stage2a_new_clash_count",
            "rosetta_delta_score",
            "variant_id",
        ],
        ascending=[True, True, False, False, False, True, True, True],
        kind="mergesort",
    )
    out["_would_have_ranked_position"] = range(1, len(out) + 1)
    return out


def cluster_cap_for(row: pd.Series, policy: TargetPolicy) -> int:
    if str(row["stage2a_final_class"]) in {"T2A_high_confidence", "T2A_good"}:
        return policy.cluster_high_conf_exception
    return policy.cluster_default_cap


def seed_key(row: pd.Series) -> str:
    return str(row.get("his_seed_set") or "")


def is_sdabhard_excluded(row: pd.Series) -> bool:
    mutation_list = row.get("mutation_list")
    return contains_mutation(mutation_list, "AG102H") or contains_mutation(mutation_list, "AV105H")


def can_add(
    row: pd.Series,
    selected: list[pd.Series],
    policy: TargetPolicy,
    target_size: int,
    cluster_counts: dict[str, int],
    seed_counts: dict[str, int],
    four_mut_count: int,
) -> tuple[bool, str]:
    cluster = str(row.get("near_duplicate_cluster_id") or "")
    current_cluster_count = cluster_counts.get(cluster, 0)
    if current_cluster_count >= cluster_cap_for(row, policy):
        return False, "cluster_cap"

    seed = seed_key(row)
    if policy.target == "Ab_sdAb":
        if is_sdabhard_excluded(row):
            return False, "sdab_AG102H_or_AV105H_excluded"
        if seed and seed_counts.get(seed, 0) >= 21:
            return False, "sdab_exact_seed_absolute_cap"
        if contains_seed(seed, "AY111H") and seed_counts.get("__contains_AY111H__", 0) >= 25:
            return False, "sdab_AY111H_selection_cap"
        if seed == "AD110H;AY111H" and seed_counts.get(seed, 0) >= 15:
            return False, "sdab_AD110H_AY111H_cap"
        if seed == "AQ100H;AD110H" and seed_counts.get(seed, 0) >= 10:
            return False, "sdab_AQ100H_AD110H_cap"
        if seed == "AQ100H;AY111H" and seed_counts.get(seed, 0) >= 3:
            return False, "sdab_AQ100H_AY111H_cap"
        mut_count = pd.to_numeric(row.get("mutation_count"), errors="coerce")
        if pd.notna(mut_count) and mut_count >= 4:
            max_four_mut = max(1, int(target_size * 0.10))
            if four_mut_count >= max_four_mut:
                return False, "sdab_4mut_cap"

    if str(row.get("stage2a_final_class")) == "T2A_supported_boundary_confirmed":
        boundary_count = sum(r.get("stage2a_final_class") == "T2A_supported_boundary_confirmed" for r in selected)
        if boundary_count >= int(target_size * policy.boundary_fraction_cap):
            return False, "confirmed_supported_boundary_cap"

    if policy.target == "Ab_1E62":
        if contains_seed(seed, "LK24H") and seed_counts.get("__contains_LK24H__", 0) >= 70:
            return False, "1E62_LK24H_containing_cap"

    return True, ""


def add_row(
    row: pd.Series,
    selected: list[pd.Series],
    cluster_counts: dict[str, int],
    seed_counts: dict[str, int],
) -> None:
    selected.append(row)
    cluster = str(row.get("near_duplicate_cluster_id") or "")
    seed = seed_key(row)
    cluster_counts[cluster] = cluster_counts.get(cluster, 0) + 1
    seed_counts[seed] = seed_counts.get(seed, 0) + 1
    if contains_seed(seed, "AY111H"):
        seed_counts["__contains_AY111H__"] = seed_counts.get("__contains_AY111H__", 0) + 1
    if contains_seed(seed, "LK24H"):
        seed_counts["__contains_LK24H__"] = seed_counts.get("__contains_LK24H__", 0) + 1


def greedy_add(
    rows: Iterable[pd.Series],
    selected: list[pd.Series],
    selected_ids: set[str],
    policy: TargetPolicy,
    target_size: int,
    cluster_counts: dict[str, int],
    seed_counts: dict[str, int],
) -> None:
    for _, row in rows:
        if len(selected) >= target_size:
            return
        variant_id = str(row["variant_id"])
        if variant_id in selected_ids:
            continue
        four_mut_count = sum(pd.to_numeric(r.get("mutation_count"), errors="coerce") >= 4 for r in selected)
        ok, _ = can_add(row, selected, policy, target_size, cluster_counts, seed_counts, int(four_mut_count))
        if ok:
            add_row(row, selected, cluster_counts, seed_counts)
            selected_ids.add(variant_id)


def select_main_pool(candidates: pd.DataFrame, policy: TargetPolicy) -> pd.DataFrame:
    selected: list[pd.Series] = []
    selected_ids: set[str] = set()
    cluster_counts: dict[str, int] = {}
    seed_counts: dict[str, int] = {}
    target_size = policy.default_size

    if policy.target == "Ab_1E62":
        high_good = candidates[candidates["stage2a_final_class"].isin(["T2A_high_confidence", "T2A_good"])]
        boundary = candidates[candidates["stage2a_final_class"].eq("T2A_supported_boundary_confirmed")]
        non_lk_high_good = high_good[~high_good["his_seed_set"].fillna("").str.contains("LK24H")]
        non_lk_boundary = boundary[~boundary["his_seed_set"].fillna("").str.contains("LK24H")]
        greedy_add(non_lk_high_good.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(high_good.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(non_lk_boundary.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(boundary.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(candidates.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)

    elif policy.target == "Ab_sdAb":
        high_good = candidates[candidates["stage2a_final_class"].isin(["T2A_high_confidence", "T2A_good"])]
        boundary = candidates[candidates["stage2a_final_class"].eq("T2A_supported_boundary_confirmed")]
        ad110_high_good = high_good[high_good["his_seed_set"].fillna("").eq("AD110H")]
        aq100_high_good = high_good[high_good["his_seed_set"].fillna("").eq("AQ100H")]
        non_ay_high_good = high_good[~high_good["his_seed_set"].fillna("").str.contains("AY111H")]
        ad110_boundary = boundary[boundary["his_seed_set"].fillna("").eq("AD110H")]
        aq100_boundary = boundary[boundary["his_seed_set"].fillna("").eq("AQ100H")]
        non_ay_boundary = boundary[~boundary["his_seed_set"].fillna("").str.contains("AY111H")]
        greedy_add(ad110_high_good.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(aq100_high_good.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(non_ay_high_good.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(high_good.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(ad110_boundary.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(aq100_boundary.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(non_ay_boundary.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(boundary.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(candidates.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)

    else:
        non_lk = candidates[~candidates["his_seed_set"].fillna("").str.contains("LK24H")]
        greedy_add(non_lk.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)
        greedy_add(candidates.iterrows(), selected, selected_ids, policy, target_size, cluster_counts, seed_counts)

    if not selected:
        return candidates.head(0).copy()

    out = pd.DataFrame([r.to_dict() for r in selected])
    out["planning_pool_role"] = "main_heavy_pool"
    out["planning_selected_rank"] = range(1, len(out) + 1)
    out["planning_pool_source"] = out["stage2a_action"].map(
        lambda x: "promoted" if x == "promote_to_heavy_candidate_pool" else "diverse_supported_boundary_reserve"
    )
    return out


def assign_routes(df: pd.DataFrame) -> pd.Series:
    routes = []
    cluster_seen: set[str] = set()
    seed_seen: set[str] = set()
    for _, row in df.iterrows():
        row_routes: list[str] = []
        seed = seed_key(row)
        cluster = str(row.get("near_duplicate_cluster_id") or "")
        pka_support = pd.to_numeric(row.get("stage2a_pH_mechanism_support"), errors="coerce")
        mutation_count = pd.to_numeric(row.get("mutation_count"), errors="coerce")

        if seed or (pd.notna(pka_support) and pka_support >= 0.55):
            row_routes.append("electrostatics_pKa_review")
        if cluster not in cluster_seen or seed not in seed_seen or row.get("stage2a_final_class") == "T2A_high_confidence":
            row_routes.append("AF3_complex_check")
        if row.get("target") == "Ab_sdAb" or row.get("stage2a_final_class") == "T2A_supported_boundary_confirmed" or (
            pd.notna(mutation_count) and mutation_count >= 3
        ):
            row_routes.append("SimpleFold_antibody_sampling")

        # No explicit glycan-risk field is present in the Stage-2A result table. Keep
        # glycan_check absent unless a future upstream feature marks it.
        if not row_routes:
            row_routes.append("AF3_complex_check")
        routes.append(";".join(dict.fromkeys(row_routes)))
        cluster_seen.add(cluster)
        seed_seen.add(seed)
    return pd.Series(routes, index=df.index)


def choose_controls(df: pd.DataFrame, policy: TargetPolicy) -> pd.DataFrame:
    controls = df[
        (df["target"].eq(policy.target))
        & (
            df["stage2a_final_class"].eq("T2A_control_or_audit")
            | df["control_anchor"].astype(str).eq("yes")
            | df["audit_only"].astype(str).eq("yes")
        )
    ].copy()
    if controls.empty:
        return controls
    controls["_seed"] = controls["his_seed_set"].fillna("WT")
    controls["_has_risk_control"] = controls["mutation_list"].fillna("").str.contains("AG102H|AV105H|V105H|D110H")
    controls["_mutation_count_num"] = pd.to_numeric(controls["mutation_count"], errors="coerce").fillna(0)
    controls = controls.sort_values(
        ["_has_risk_control", "_seed", "_mutation_count_num", "stage2a_structure_score", "variant_id"],
        ascending=[False, True, True, False, True],
        kind="mergesort",
    )
    picked = []
    seed_counts: dict[str, int] = {}
    # First retain all distinct controls/risk anchors at low per-seed cap.
    for _, row in controls.iterrows():
        if len(picked) >= policy.control_target:
            break
        seed = str(row["_seed"])
        cap = 3 if seed not in {"AD110H", "LD34H;LQ35H"} else 6
        if seed_counts.get(seed, 0) >= cap:
            continue
        picked.append(row)
        seed_counts[seed] = seed_counts.get(seed, 0) + 1
    # Fill any remaining slots with the best available controls, still avoiding extreme dominance.
    for _, row in controls.iterrows():
        if len(picked) >= policy.control_target:
            break
        if str(row["variant_id"]) in {str(r["variant_id"]) for r in picked}:
            continue
        seed = str(row["_seed"])
        cap = 10 if policy.target == "Ab_1E62" else 14
        if seed_counts.get(seed, 0) >= cap:
            continue
        picked.append(row)
        seed_counts[seed] = seed_counts.get(seed, 0) + 1
    controls = pd.DataFrame([r.to_dict() for r in picked]).head(policy.control_target)
    controls["planning_pool_role"] = "heavy_control_anchor_panel"
    controls["planning_selected_rank"] = range(1, len(controls) + 1)
    controls["planning_pool_source"] = "control_anchor"
    controls["heavy_route"] = "control_or_anchor_review"
    return controls


def candidate_universe(df: pd.DataFrame, policy: TargetPolicy) -> pd.DataFrame:
    base = df[
        (df["target"].eq(policy.target))
        & (df["stage2a_final_class"].isin(MAIN_CLASSES))
        & (~df["stage2a_final_class"].eq("T2A_control_or_audit"))
        & (~df["control_anchor"].astype(str).eq("yes"))
        & (~df["audit_only"].astype(str).eq("yes"))
        & (df["stage2a_action"].isin(["promote_to_heavy_candidate_pool", "retain_for_stage2b"]))
    ].copy()
    if policy.target == "Ab_sdAb":
        base = base[~base["mutation_list"].fillna("").str.contains("AG102H|AV105H")].copy()
    return prepare_rank(base)


def cap_audit(selected: pd.DataFrame, policy: TargetPolicy) -> dict[str, object]:
    total = len(selected)
    if total == 0:
        return {
            "target": policy.target,
            "main_count": 0,
            "verdict": "FAIL",
            "reasons": "empty_pool",
        }

    cluster_counts = selected["near_duplicate_cluster_id"].fillna("").value_counts()
    seed_counts = selected["his_seed_set"].fillna("").value_counts()
    class_counts = selected["stage2a_final_class"].value_counts()
    top_seed = seed_counts.index[0] if len(seed_counts) else ""
    top_cluster = cluster_counts.index[0] if len(cluster_counts) else ""
    boundary_count = int(class_counts.get("T2A_supported_boundary_confirmed", 0))
    boundary_fraction = boundary_count / total
    top_cluster_count = int(cluster_counts.iloc[0]) if len(cluster_counts) else 0
    top5_fraction = float(cluster_counts.head(5).sum() / total) if total else 0.0

    reasons: list[str] = []
    verdict = "PASS"
    if total < policy.clean_min:
        verdict = "PATCH"
        reasons.append("main_count_below_preferred_clean_range")
    if total < policy.lower_clean_min:
        verdict = "PATCH"
        reasons.append("main_count_below_lower_clean_range_due_diversity_caps")
    if top_cluster_count > policy.top_cluster_absolute_cap:
        verdict = "FAIL"
        reasons.append("top_cluster_absolute_cap_exceeded")
    if top5_fraction > policy.top5_cluster_fraction_cap:
        verdict = "FAIL"
        reasons.append("top5_cluster_fraction_cap_exceeded")
    if boundary_fraction > policy.boundary_fraction_cap:
        verdict = "PATCH" if verdict != "FAIL" else verdict
        reasons.append("confirmed_supported_boundary_fraction_above_cap")
    if selected["stage2a_final_class"].isin(["T2A_reject", "T2A_structure_risk", "T2A_mechanism_weak", "T2A_boundary_unresolved"]).any():
        verdict = "FAIL"
        reasons.append("ineligible_final_class_in_main_pool")
    if selected["heavy_route"].isna().any() or selected["heavy_route"].astype(str).eq("").any():
        verdict = "PATCH" if verdict != "FAIL" else verdict
        reasons.append("missing_heavy_route")

    if policy.target == "Ab_1E62":
        non_lk = int((~selected["his_seed_set"].fillna("").str.contains("LK24H")).sum())
        if non_lk < min(20, total):
            verdict = "PATCH" if verdict != "FAIL" else verdict
            reasons.append("non_LK24H_representatives_below_minimum")
        seed_specific = {
            "non_LK24H_count": non_lk,
            "LK24H_containing_count": int(selected["his_seed_set"].fillna("").str.contains("LK24H").sum()),
        }
    else:
        non_ay = int((~selected["his_seed_set"].fillna("").str.contains("AY111H")).sum())
        ay = int(selected["his_seed_set"].fillna("").str.contains("AY111H").sum())
        ag102 = int(selected["mutation_list"].fillna("").str.contains("AG102H").sum())
        av105 = int(selected["mutation_list"].fillna("").str.contains("AV105H").sum())
        if ay > 25:
            verdict = "FAIL"
            reasons.append("AY111H_selection_cap_exceeded")
        if seed_counts.iloc[0] > 21:
            verdict = "FAIL"
            reasons.append("top_seed_absolute_cap_exceeded")
        if seed_counts.iloc[0] / total > 0.28:
            verdict = "PATCH" if verdict != "FAIL" else verdict
            reasons.append("top_seed_fraction_above_28pct_manual_waiver_line")
        if ay > 22:
            verdict = "PATCH" if verdict != "FAIL" else verdict
            reasons.append("AY111H_containing_above_preferred_22")
        if non_ay / total < 0.70:
            verdict = "PATCH" if verdict != "FAIL" else verdict
            reasons.append("non_AY111H_fraction_below_minimum")
        if ag102 or av105:
            verdict = "FAIL"
            reasons.append("AG102H_or_AV105H_in_main_pool")
        seed_specific = {
            "non_AY111H_count": non_ay,
            "AY111H_containing_count": ay,
            "AG102H_count": ag102,
            "AV105H_count": av105,
            "AD110H_count": int(selected["his_seed_set"].fillna("").eq("AD110H").sum()),
            "AQ100H_count": int(selected["his_seed_set"].fillna("").eq("AQ100H").sum()),
            "AD110H_AY111H_count": int(selected["his_seed_set"].fillna("").eq("AD110H;AY111H").sum()),
            "AQ100H_AD110H_count": int(selected["his_seed_set"].fillna("").eq("AQ100H;AD110H").sum()),
        }

    return {
        "target": policy.target,
        "main_count": total,
        "preferred_min": policy.clean_min,
        "preferred_max": policy.clean_max,
        "lower_clean_min": policy.lower_clean_min,
        "boundary_count": boundary_count,
        "boundary_fraction": round(boundary_fraction, 4),
        "top_seed": top_seed,
        "top_seed_count": int(seed_counts.iloc[0]) if len(seed_counts) else 0,
        "top_seed_fraction": round(float(seed_counts.iloc[0] / total), 4) if len(seed_counts) else 0.0,
        "top_cluster": top_cluster,
        "top_cluster_count": top_cluster_count,
        "top_cluster_fraction": round(float(top_cluster_count / total), 4),
        "top5_cluster_fraction": round(top5_fraction, 4),
        "four_mut_count": int(pd.to_numeric(selected["mutation_count"], errors="coerce").ge(4).sum()),
        "four_mut_fraction": round(float(pd.to_numeric(selected["mutation_count"], errors="coerce").ge(4).sum() / total), 4),
        "verdict": verdict,
        "reasons": ";".join(reasons) if reasons else "none",
        **seed_specific,
    }


def exclusion_reason(row: pd.Series, selected_ids: set[str], policy: TargetPolicy) -> tuple[str, str]:
    if str(row["variant_id"]) in selected_ids:
        return "", ""
    if policy.target == "Ab_sdAb" and is_sdabhard_excluded(row):
        return "sdab_AG102H_or_AV105H_excluded", "seed_risk"
    if row["stage2a_final_class"] not in MAIN_CLASSES:
        return "ineligible_final_class", "quality"
    if row["stage2a_action"] != "promote_to_heavy_candidate_pool":
        return "not_promoted_or_reserve_not_selected", "rank_or_cap"
    return "diversity_cap_or_rank_limit", "seed_or_cluster_cap"


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    lines = [
        "| " + " | ".join(map(str, df.columns)) + " |",
        "| " + " | ".join("---" for _ in df.columns) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in df.columns) + " |")
    return "\n".join(lines)


def write_report(
    audits: pd.DataFrame,
    class_summary: pd.DataFrame,
    seed_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    route_summary: pd.DataFrame,
    controls_summary: pd.DataFrame,
) -> str:
    overall = "PASS"
    if audits["verdict"].eq("FAIL").any():
        overall = "FAIL"
    elif audits["verdict"].eq("PATCH").any():
        overall = "PATCH"

    report = [
        "# Targeted Tier2-heavy Planning Pool Audit",
        "",
        "## Executive Summary",
        "",
        f"Overall verdict: `{overall}`.",
        "",
        "This stage constructed capped targeted Tier2-heavy planning pools from Stage-2A results. It did not run Tier2-heavy compute, AF3, SimpleFold, MD, or final 10K selection.",
        "",
        "The main decision point is diversity capping. A smaller clean pool is preferred over refilling from dominant seed or near-duplicate clusters.",
        "",
        "## Input",
        "",
        f"Stage-2A compute results: `{INPUT_CSV.as_posix()}`",
        f"Input SHA256: `{sha256_file(INPUT_CSV)}`",
        "",
        "Stage-2A promoted candidates before capping:",
        "",
        "```text",
        "1E62: 183",
        "sdAb: 214",
        "```",
        "",
        "## Audit Verdict By Target",
        "",
        md_table(audits),
        "",
        "## Main Pool Final Class Summary",
        "",
        md_table(class_summary),
        "",
        "## Main Pool Seed Summary",
        "",
        md_table(seed_summary),
        "",
        "## Main Pool Cluster Summary",
        "",
        md_table(cluster_summary),
        "",
        "## Heavy Route Summary",
        "",
        md_table(route_summary),
        "",
        "## Control / Anchor Panel Summary",
        "",
        md_table(controls_summary),
        "",
        "Control / anchor rationale:",
        "",
        "`tier2_heavy_control_anchor_panel_rationale.md` explains why skewed control seeds are retained and how risk controls are labeled.",
        "",
        "## Interpretation",
        "",
        "A `PASS` target can proceed to manual review for Tier2-heavy-lite unlock. A `PATCH` target should not run heavy compute yet; it needs either acceptance of a smaller clean pool or targeted additions that do not violate seed / cluster caps.",
        "",
        "Tier2-heavy execution remains locked until the user explicitly unlocks it after reviewing this audit.",
        "",
    ]
    text = "\n".join(report)
    (OUT_DIR / "tier2_heavy_planning_pool_audit.md").write_text(text, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(text, encoding="utf-8")
    return overall


def write_control_rationale(controls_all: pd.DataFrame) -> None:
    lines = [
        "# Tier2-heavy Control / Anchor Panel Rationale",
        "",
        "This file explains the separate control / anchor panel used with the targeted Tier2-heavy planning pool. Controls are not counted toward the main heavy-pool target size.",
        "",
        "## Overall Role",
        "",
        "Controls / anchors are included for calibration, risk sanity checks, and mechanism interpretation. They are not main candidates and should not be promoted by main-pool ranking rules.",
        "",
    ]
    if controls_all.empty:
        lines.append("No controls were selected.")
    for target, group in controls_all.groupby("target"):
        seed_counts = group["his_seed_set"].fillna("WT").value_counts()
        lines.extend(
            [
                f"## {target}",
                "",
                f"Control count: {len(group)}",
                "",
                "Seed composition:",
                "",
                "```text",
                seed_counts.to_string(),
                "```",
                "",
            ]
        )
        if target == "Ab_1E62":
            lines.extend(
                [
                    "Rationale:",
                    "",
                    "- Parent / WT and single-His anchors are retained where available.",
                    "- `LD34H;LQ35H` remains enriched because the available 1E62 control source is strongly enriched for this calibration seed.",
                    "- The selector caps this enrichment where possible; remaining enrichment should be interpreted as route/control calibration, not main-pool diversity.",
                    "",
                ]
            )
        elif target == "Ab_sdAb":
            lines.extend(
                [
                    "Rationale:",
                    "",
                    "- `AG102H` and `AV105H` are retained only as explicit risk/audit controls.",
                    "- Parent / WT, `AY111H`, `AQ100H`, and `AE108H` anchors are retained where available.",
                    "- `AD110H` remains enriched because the available sdAb control candidates are strongly AD110H-heavy; this is calibration context, not a main-pool recommendation.",
                    "",
                ]
            )
    (OUT_DIR / "tier2_heavy_control_anchor_panel_rationale.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_CSV)
    main_pools = []
    control_panels = []
    excluded_tables = []
    audit_rows = []

    for target, policy in POLICIES.items():
        universe = candidate_universe(df, policy)
        selected = select_main_pool(universe, policy)
        if not selected.empty:
            selected["heavy_route"] = assign_routes(selected)
        selected_ids = set(selected["variant_id"].astype(str))

        controls = choose_controls(df, policy)
        main_pools.append(selected)
        control_panels.append(controls)
        audit_rows.append(cap_audit(selected, policy))

        promoted = prepare_rank(
            df[(df["target"].eq(target)) & (df["stage2a_action"].eq("promote_to_heavy_candidate_pool"))].copy()
        )
        excluded = promoted[~promoted["variant_id"].astype(str).isin(selected_ids)].copy()
        if not excluded.empty:
            reasons = excluded.apply(lambda r: exclusion_reason(r, selected_ids, policy), axis=1, result_type="expand")
            excluded["excluded_reason"] = reasons[0]
            excluded["excluded_by_cap_type"] = reasons[1]
            excluded["would_have_ranked_position"] = excluded["_would_have_ranked_position"]
        excluded_tables.append(excluded)

    main = pd.concat(main_pools, ignore_index=True) if main_pools else pd.DataFrame()
    controls_all = pd.concat(control_panels, ignore_index=True) if control_panels else pd.DataFrame()
    excluded_all = pd.concat(excluded_tables, ignore_index=True) if excluded_tables else pd.DataFrame()
    audit = pd.DataFrame(audit_rows)

    output_cols = [
        "variant_id",
        "target",
        "canonical_sequence_hash",
        "mutation_list",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "mutation_count",
        "stage2a_final_class",
        "stage2a_action",
        "input_class",
        "boundary_support_level",
        "stage2a_structure_score",
        "stage2a_local_validity",
        "stage2a_new_clash_count",
        "stage2a_pH_mechanism_support",
        "rosetta_delta_score",
        "planning_pool_role",
        "planning_pool_source",
        "planning_selected_rank",
        "heavy_route",
        "pdb_path",
        "source_bank",
        "source_round",
    ]
    for target, policy in POLICIES.items():
        target_main = main[main["target"].eq(target)].copy()
        target_main[output_cols].to_csv(OUT_DIR / f"tier2_heavy_planning_pool_{policy.output_name}.csv", index=False)

    controls_all[output_cols].to_csv(OUT_DIR / "tier2_heavy_control_anchor_panel.csv", index=False)

    excluded_cols = [
        "variant_id",
        "target",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "stage2a_final_class",
        "stage2a_action",
        "mutation_list",
        "excluded_reason",
        "excluded_by_cap_type",
        "would_have_ranked_position",
    ]
    excluded_all[excluded_cols].to_csv(OUT_DIR / "tier2_heavy_planning_excluded_promoted_candidates.csv", index=False)
    audit.to_csv(OUT_DIR / "tier2_heavy_planning_audit_metrics.csv", index=False)

    class_summary = (
        main.groupby(["target", "stage2a_final_class"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "stage2a_final_class"])
    )
    seed_summary = (
        main.groupby(["target", "his_seed_set"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )
    cluster_summary = (
        main.groupby(["target", "near_duplicate_cluster_id"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
        .groupby("target")
        .head(12)
        .reset_index(drop=True)
    )
    route_summary = (
        main.assign(heavy_route=main["heavy_route"].str.split(";"))
        .explode("heavy_route")
        .groupby(["target", "heavy_route"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )
    controls_summary = (
        controls_all.groupby(["target", "his_seed_set"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "count"], ascending=[True, False])
    )

    class_summary.to_csv(OUT_DIR / "tier2_heavy_planning_class_summary.csv", index=False)
    seed_summary.to_csv(OUT_DIR / "tier2_heavy_planning_seed_summary.csv", index=False)
    cluster_summary.to_csv(OUT_DIR / "tier2_heavy_planning_top_cluster_summary.csv", index=False)
    route_summary.to_csv(OUT_DIR / "tier2_heavy_planning_route_summary.csv", index=False)
    controls_summary.to_csv(OUT_DIR / "tier2_heavy_control_anchor_summary.csv", index=False)

    overall = write_report(audit, class_summary, seed_summary.groupby("target").head(12).reset_index(drop=True), cluster_summary, route_summary, controls_summary)
    write_control_rationale(controls_all)

    with PROGRESS.open("a", encoding="utf-8") as handle:
        handle.write(
            "\n## 2026-05-29 - Targeted Tier2-heavy planning pool construction\n\n"
            f"- Built capped planning pools in `{OUT_DIR.as_posix()}`.\n"
            f"- Overall audit verdict: `{overall}`.\n"
            f"- 1E62 main pool count: {len(main[main['target'].eq('Ab_1E62')])}.\n"
            f"- sdAb main pool count: {len(main[main['target'].eq('Ab_sdAb')])}.\n"
            f"- Control / anchor panel count: {len(controls_all)}.\n"
            "- Tier2-heavy compute remains locked pending manual review.\n"
        )

    print(f"overall={overall}")
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
