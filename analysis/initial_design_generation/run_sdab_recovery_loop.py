#!/usr/bin/env python3
"""Prepare and run sdAb recovery-loop round inputs.

This script implements the compute-light portion of
`.tasks/active/initial-design-generation/sdab_recovery_loop_execution_plan.md`.
It builds a recovery-round failure map, negative mask, positive seed bank, recovery
config, targeted generation pool, Tier1-style cheap gate, 1K Stage-1 mini-batch
list, and preflight report. It deliberately does not launch Stage-2A compute,
Tier2-heavy, AF3, SimpleFold, MD, or final 10K selection.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import re
import socket
import subprocess
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT
STAGE15 = ROOT / "results/initial_design_generation/stage1_5_stage2a/stage1_refined_structure_risk.csv"
STAGE2A = ROOT / "results/initial_design_generation/stage1_5_stage2a/stage2a_candidate_list.csv"
OUT_ROOT = ROOT / "results/initial_design_generation/sdab_recovery_loop"
ROUND_DIR = OUT_ROOT / "round_01"
ONEE62_DIR = OUT_ROOT / "1e62_sync_waiting_bank"

TARGET = "Ab_sdAb"
DRY_TARGET = "sdAb"
CHAIN = "A"
WINDOW_ID = "sdAb_VHH_072_111"

WETLAB_SUPPORTED_SEEDS = {
    "AQ100H",
    "AG102H",
    "AV105H",
    "AE108H",
    "AD110H",
    "AY111H",
    "AQ100H;AD110H",
    "AQ100H;AY111H",
    "AD110H;AY111H",
}
PRIORITY_SEED_CAPS = {
    "AD110H;AY111H": 0.20,
    "AQ100H;AD110H": 0.15,
    "AQ100H;AY111H": 0.15,
}
DEFAULT_SEED_CAP = 0.25
SEED_HARD_LINE_FIRST_ROUND = 0.35
CLUSTER_HARD_LINE_FIRST_ROUND = 0.20
FOUR_MUT_HARD_LINE_FIRST_ROUND = 0.20

SUPPORTED_ROUTES = {
    "His_rule",
    "His_plus_rescue",
    "ProteinMPNN_seeded_rescue",
    "wetlab_informed_expansion",
    "structure_or_interface_guided",
    "repaired_sparse_constrained_MPNN",
}
RECOVERY_ROUTES = [
    "His_plus_rescue",
    "His_rule",
    "ProteinMPNN_seeded_rescue",
    "wetlab_informed_expansion",
    "structure_or_interface_guided",
]
ROUTE_WEIGHTS = {
    "His_plus_rescue": 0.36,
    "His_rule": 0.24,
    "ProteinMPNN_seeded_rescue": 0.20,
    "wetlab_informed_expansion": 0.15,
    "structure_or_interface_guided": 0.05,
}
CONSERVATIVE_RESCUE_AA = ["A", "S", "T", "N", "Q", "Y", "D", "E", "K", "R", "G", "P"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage15", default=str(STAGE15))
    parser.add_argument("--stage2a", default=str(STAGE2A))
    parser.add_argument("--out-root", default=str(OUT_ROOT))
    parser.add_argument("--round", default="round_01")
    parser.add_argument("--raw-target", type=int, default=40_000)
    parser.add_argument("--post-hard-unique-target", type=int, default=18_000)
    parser.add_argument("--tier1-proposal-target", type=int, default=3_500)
    parser.add_argument("--stage1-mini-batch-target", type=int, default=1_000)
    parser.add_argument("--seed", type=int, default=20260528)
    parser.add_argument("--mode", choices=["prep", "generate", "all"], default="all")
    parser.add_argument("--current-bank", default="")
    return parser.parse_args()


def ensure_dirs(out_root: Path, round_name: str) -> tuple[Path, Path]:
    round_dir = out_root / round_name
    onee62_dir = out_root / "1e62_sync_waiting_bank"
    round_dir.mkdir(parents=True, exist_ok=True)
    onee62_dir.mkdir(parents=True, exist_ok=True)
    return round_dir, onee62_dir


def csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy() if max_rows else df.copy()
    lines = [
        "| " + " | ".join(str(c) for c in view.columns) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        values = []
        for col in view.columns:
            value = row[col]
            if isinstance(value, float):
                values.append(f"{value:.4f}")
            else:
                values.append(str(value).replace("\n", " "))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def norm_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    value = str(value).strip()
    return "" if value.lower() == "nan" else value


def numeric(series: pd.Series | Any, default: float = math.nan) -> pd.Series:
    if isinstance(series, pd.Series):
        return pd.to_numeric(series, errors="coerce")
    return pd.to_numeric(pd.Series(series), errors="coerce").fillna(default)


def mutation_labels(value: Any) -> list[str]:
    text = norm_text(value)
    if not text:
        return []
    return [x.strip() for x in text.replace(",", ";").split(";") if x.strip() and x.strip().lower() != "nan"]


def mutation_positions(value: Any) -> list[int]:
    positions: list[int] = []
    for label in mutation_labels(value):
        m = re.match(r"^[A-Za-z][A-Za-z](\d+)[A-Za-z]$", label)
        if m:
            positions.append(int(m.group(1)))
    return sorted(set(positions))


def mutation_position_pairs(value: Any) -> list[str]:
    positions = mutation_positions(value)
    return [f"{a};{b}" for a, b in itertools.combinations(positions, 2)]


def class_counts(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        total = len(sub)
        severe = int(sub["refined_structure_risk_class"].eq("T2_severe_structure_risk").sum())
        boundary = int(sub["refined_structure_risk_class"].eq("T2_boundary_structure_risk").sum())
        pass_like = int(sub["refined_structure_risk_class"].eq("T2_pass_like_structure").sum())
        manual = int(sub["refined_structure_risk_class"].eq("T2_manual_review").sum())
        row = {col: key for col, key in zip(group_cols, keys)}
        row.update(
            {
                "total_count": total,
                "severe_count": severe,
                "boundary_count": boundary,
                "pass_like_count": pass_like,
                "manual_count": manual,
                "severe_fraction": severe / max(1, total),
                "boundary_fraction": boundary / max(1, total),
                "pass_like_fraction": pass_like / max(1, total),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["severe_count", "total_count"], ascending=False)


def seed_support_type(seed: str, sub: pd.DataFrame) -> str:
    support: list[str] = []
    pass_like = int(sub["refined_structure_risk_class"].eq("T2_pass_like_structure").sum())
    boundary = sub[sub["refined_structure_risk_class"].eq("T2_boundary_structure_risk")].copy()
    if pass_like:
        support.append("pass_like_seed")
    if not boundary.empty:
        low_clash = (
            pd.to_numeric(boundary.get("new_clash_count_total", 0), errors="coerce").fillna(0).le(2)
            & pd.to_numeric(boundary.get("max_new_clash_overlap", 0), errors="coerce").fillna(0).lt(1.1)
        )
        if int(low_clash.sum()) > 0:
            support.append("low_clash_boundary")
    if seed in WETLAB_SUPPORTED_SEEDS:
        support.append("wetlab_supported_seed")
    route_text = ";".join(sub.get("primary_generation_route", pd.Series(dtype=str)).fillna("").astype(str).unique())
    route_text += ";" + ";".join(sub.get("all_source_routes", pd.Series(dtype=str)).fillna("").astype(str).unique())
    if "MPNN" in route_text:
        support.append("mpnn_rescue_provenance")
    neutral = pd.to_numeric(sub.get("neutral_retention_t2_score", pd.Series(dtype=float)), errors="coerce")
    if neutral.notna().any() and float(neutral.quantile(0.75)) >= 0.55:
        support.append("neutral_retention_support")
    return ";".join(sorted(set(support))) if support else "insufficient_support"


def build_failure_maps(stage15: pd.DataFrame, round_dir: Path) -> dict[str, pd.DataFrame]:
    sdab = stage15[stage15["target"].eq(TARGET)].copy()
    maps: dict[str, pd.DataFrame] = {}
    maps["seed"] = class_counts(sdab.assign(his_seed_set=sdab["his_seed_set"].fillna("")), ["his_seed_set"])
    maps["route"] = class_counts(sdab, ["primary_generation_route"])
    maps["mutation_count"] = class_counts(sdab, ["mutation_count"])
    if "rescue_objective_list" in sdab.columns:
        rescue_obj = sdab.assign(rescue_objective=sdab["rescue_objective_list"].fillna(""))
    else:
        rescue_obj = sdab.assign(rescue_objective=sdab.get("rescue_signature", "").fillna(""))
    maps["rescue_objective"] = class_counts(rescue_obj, ["rescue_objective"])

    pair_rows = []
    for _, row in sdab.iterrows():
        pairs = mutation_position_pairs(row.get("mutation_list", ""))
        if not pairs:
            pairs = ["single_or_no_pair"]
        for pair in pairs:
            pair_rows.append(
                {
                    "position_pair": pair,
                    "refined_structure_risk_class": row["refined_structure_risk_class"],
                    "variant_id": row["variant_id"],
                }
            )
    pair_df = pd.DataFrame(pair_rows)
    maps["position_pair"] = class_counts(pair_df, ["position_pair"])

    cluster_col = "near_duplicate_cluster_id"
    if cluster_col not in sdab.columns:
        cluster_col = "tier1_near_duplicate_cluster_id"
    maps["near_duplicate_cluster"] = class_counts(sdab.assign(cluster=sdab.get(cluster_col, "")), ["cluster"])

    severe_signature = (
        sdab[sdab["refined_structure_risk_class"].eq("T2_severe_structure_risk")]
        .groupby(["risk_reason_codes", "his_seed_set", "primary_generation_route"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )
    maps["severe_signature"] = severe_signature

    csv(maps["seed"], round_dir / "sdab_failure_by_seed.csv")
    csv(maps["route"], round_dir / "sdab_failure_by_route.csv")
    csv(maps["mutation_count"], round_dir / "sdab_failure_by_mutation_count.csv")
    csv(maps["rescue_objective"], round_dir / "sdab_failure_by_rescue_objective.csv")
    csv(maps["position_pair"], round_dir / "sdab_failure_by_position_pair.csv")
    csv(maps["near_duplicate_cluster"], round_dir / "sdab_failure_by_near_duplicate_cluster.csv")
    csv(maps["severe_signature"], round_dir / "sdab_severe_signature_table.csv")
    return maps


def build_negative_mask(maps: dict[str, pd.DataFrame], round_dir: Path, round_name: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def add(mask_type: str, seed_set: str = "", position_set: str = "", mutation_pattern: str = "", reason: str = "", severity: str = "hard") -> None:
        rows.append(
            {
                "target": TARGET,
                "mask_type": mask_type,
                "seed_set": seed_set,
                "position_set": position_set,
                "mutation_pattern": mutation_pattern,
                "reason": reason,
                "severity": severity,
                "source_iteration": f"{round_name}_prep",
            }
        )

    add("forbid_route_pattern", mutation_pattern="broad_AE108H_centered_expansion", reason="external_policy_hold", severity="hard")
    add("forbid_route_pattern", mutation_pattern="weak_filler", reason="external_policy_hold", severity="hard")
    add("forbid_route_pattern", mutation_pattern="broad_F_boundary_expansion", reason="external_policy_hold", severity="hard")
    add("forbid_pair", seed_set="AV105H;AD110H", position_set="105;110", reason="known_negative_pair", severity="hard")

    seed_map = maps["seed"].copy()
    seed_map["his_seed_set"] = seed_map["his_seed_set"].fillna("").astype(str)
    for _, row in seed_map.iterrows():
        seed = str(row["his_seed_set"])
        if not seed:
            continue
        severe_fraction = float(row["severe_fraction"])
        pass_like_fraction = float(row["pass_like_fraction"])
        total = int(row["total_count"])
        if "AE108H" in seed:
            add("seed_set_hold", seed_set=seed, reason="AE108H_broad_expansion_hold", severity="hard")
        elif total >= 5 and severe_fraction >= 0.75 and pass_like_fraction < 0.05:
            add(
                "seed_set_cap_or_hold",
                seed_set=seed,
                reason=f"high_severe_seed severe_fraction={severe_fraction:.3f} pass_like_fraction={pass_like_fraction:.3f}",
                severity="cap",
            )

    pair_map = maps["position_pair"].copy()
    for _, row in pair_map.iterrows():
        pair = str(row["position_pair"])
        total = int(row["total_count"])
        severe_fraction = float(row["severe_fraction"])
        pass_like_fraction = float(row["pass_like_fraction"])
        if pair == "single_or_no_pair":
            continue
        if total >= 5 and severe_fraction >= 0.65 and pass_like_fraction < 0.05:
            add(
                "position_pair_cap_or_hold",
                position_set=pair,
                reason=f"high_clash_or_severe_position_pair severe_fraction={severe_fraction:.3f}",
                severity="cap",
            )

    sig = maps["severe_signature"].head(25).copy()
    for _, row in sig.iterrows():
        add(
            "severe_signature_repeat",
            seed_set=norm_text(row.get("his_seed_set", "")),
            mutation_pattern=norm_text(row.get("risk_reason_codes", "")),
            reason=f"top_severe_signature_count={int(row.get('count', 0))}",
            severity="cap",
        )

    out = pd.DataFrame(rows)
    csv(out, round_dir / "sdab_negative_mask.csv")
    return out


def build_positive_seed_bank(stage15: pd.DataFrame, round_dir: Path) -> pd.DataFrame:
    sdab = stage15[stage15["target"].eq(TARGET)].copy()
    rows: list[dict[str, Any]] = []
    for seed, sub in sdab.groupby(sdab["his_seed_set"].fillna("").astype(str), dropna=False):
        if not seed:
            continue
        support = seed_support_type(seed, sub)
        severe = int(sub["refined_structure_risk_class"].eq("T2_severe_structure_risk").sum())
        boundary = int(sub["refined_structure_risk_class"].eq("T2_boundary_structure_risk").sum())
        pass_like = int(sub["refined_structure_risk_class"].eq("T2_pass_like_structure").sum())
        total = len(sub)
        seed_positions = [int(x) for x in re.findall(r"\d+", str(seed))]
        if support == "insufficient_support" and pass_like == 0:
            continue
        if pass_like == 0 and boundary == 0 and len(seed_positions) > 1:
            continue
        cap_fraction = PRIORITY_SEED_CAPS.get(seed, DEFAULT_SEED_CAP)
        if "AE108H" in seed and pass_like == 0:
            cap_fraction = 0.0
        rows.append(
            {
                "seed_set": seed,
                "support_type": support,
                "total_count": total,
                "pass_like_count": pass_like,
                "boundary_count": boundary,
                "severe_count": severe,
                "pass_like_fraction": pass_like / max(1, total),
                "severe_fraction": severe / max(1, total),
                "neutral_retention_support": "neutral_retention_support" in support,
                "near_duplicate_risk": "needs_cluster_cap",
                "allowed_max_fraction": cap_fraction,
                "allowed_max_count": int(math.floor(cap_fraction * 1000)),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(
            ["pass_like_count", "boundary_count", "severe_fraction"],
            ascending=[False, False, True],
        )
    csv(out, round_dir / "sdab_positive_seed_bank.csv")
    return out


def supported_boundary_status(row: pd.Series) -> str:
    cls = norm_text(row.get("refined_structure_risk_class", ""))
    if cls == "T2_pass_like_structure":
        return "pass_like"
    if cls != "T2_boundary_structure_risk":
        return "not_boundary"
    supports: list[str] = []
    new_clash = pd.to_numeric(row.get("new_clash_count_total"), errors="coerce")
    max_overlap = pd.to_numeric(row.get("max_new_clash_overlap"), errors="coerce")
    if pd.notna(new_clash) and pd.notna(max_overlap) and float(new_clash) <= 2 and float(max_overlap) < 1.1:
        supports.append("low_clash_boundary")
    seed = norm_text(row.get("his_seed_set", ""))
    if seed in WETLAB_SUPPORTED_SEEDS:
        supports.append("wetlab_seed")
    routes = norm_text(row.get("primary_generation_route", "")) + ";" + norm_text(row.get("all_source_routes", ""))
    if "MPNN" in routes:
        supports.append("mpnn_provenance")
    neutral = pd.to_numeric(row.get("neutral_retention_t2_score"), errors="coerce")
    if pd.notna(neutral) and float(neutral) >= 0.55:
        supports.append("neutral_retention")
    if norm_text(row.get("tier1_review_class", "")).startswith(("A_", "B_", "C_")):
        supports.append("tier1_support")
    return "supported_boundary:" + ";".join(sorted(set(supports))) if supports else "unsupported_boundary"


def build_current_banks(stage15: pd.DataFrame, stage2a: pd.DataFrame, round_dir: Path, onee62_dir: Path) -> None:
    candidate = stage2a.copy()
    if "supported_boundary_status" not in candidate.columns:
        candidate["supported_boundary_status"] = candidate.apply(supported_boundary_status, axis=1)

    sdab_bank = candidate[
        candidate["target"].eq(TARGET)
        & candidate["refined_structure_risk_class"].isin(["T2_pass_like_structure", "T2_boundary_structure_risk"])
        & ~candidate["supported_boundary_status"].eq("unsupported_boundary")
    ].copy()
    csv(sdab_bank, round_dir / "sdab_candidate_bank.csv")

    onee62 = candidate[candidate["target"].eq("Ab_1E62")].copy()
    onee62["supported_boundary_status"] = onee62.apply(supported_boundary_status, axis=1)
    csv(
        onee62[
            [
                "variant_id",
                "target",
                "mutation_list",
                "his_seed_set",
                "primary_generation_route",
                "refined_structure_risk_class",
                "supported_boundary_status",
                "stage2_action",
                "stage2a_list_action",
                "stage2a_rank_score",
            ]
            if "stage2a_rank_score" in onee62.columns
            else [
                "variant_id",
                "target",
                "mutation_list",
                "his_seed_set",
                "primary_generation_route",
                "refined_structure_risk_class",
                "supported_boundary_status",
                "stage2_action",
                "stage2a_list_action",
            ]
        ],
        onee62_dir / "1e62_supported_boundary_strata.csv",
    )
    seed_whitelist = class_counts(stage15[stage15["target"].eq("Ab_1E62")], ["his_seed_set"])
    seed_whitelist = seed_whitelist[
        (seed_whitelist["pass_like_count"].gt(0))
        | ((seed_whitelist["boundary_count"].gt(0)) & (seed_whitelist["severe_fraction"].lt(0.75)))
    ].copy()
    seed_whitelist["status"] = "sync_waiting_seed_whitelist"
    csv(seed_whitelist, onee62_dir / "1e62_seed_whitelist.csv")

    controls = stage15[(stage15["target"].eq("Ab_1E62")) & stage15.get("tier2_class", "").eq("T2_control_anchor")].copy()
    if controls.empty:
        controls = stage15[stage15["target"].eq("Ab_1E62")].head(0).copy()
    csv(controls, onee62_dir / "1e62_controls_anchors_draft.csv")

    report = [
        "# 1E62 Sync Waiting Bank Audit",
        "",
        "Status: `bank freeze / wait-for-sync`.",
        "",
        "No Tier2-heavy, final 10K, or broad new Tier2 compute is unlocked.",
        "",
        "## Current Stage-2A Candidate Rows",
        "",
        markdown_table(onee62.groupby(["refined_structure_risk_class", "supported_boundary_status"], dropna=False).size().reset_index(name="count")),
        "",
        "## Seed Whitelist Summary",
        "",
        markdown_table(seed_whitelist.head(30)),
    ]
    (onee62_dir / "1e62_sync_waiting_bank_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def recovery_config(seed_bank: pd.DataFrame, negative_mask: pd.DataFrame, round_dir: Path, args: argparse.Namespace) -> dict[str, Any]:
    seed_caps = {
        str(row["seed_set"]): float(row["allowed_max_fraction"])
        for _, row in seed_bank.iterrows()
        if float(row.get("allowed_max_fraction", 0)) > 0
    }
    cfg = {
        "run_id": f"sdab_recovery_{args.round}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "target": TARGET,
        "window": "VHH:72-111",
        "stage2a_compute_unlocked": False,
        "tier2_heavy_unlocked": False,
        "final_10k_selection_unlocked": False,
        args.round: {
            "raw_generated_target": args.raw_target,
            "post_hard_exact_unique_target": args.post_hard_unique_target,
            "tier1_proposal_target": args.tier1_proposal_target,
            "stage1_mini_batch_target": args.stage1_mini_batch_target,
        },
        "generation": {
            "max_his_count": 2,
            "default_max_mutations": 3,
            "four_mut_policy": {
                "allowed": True,
                "max_fraction": 0.05,
                "require": ["positive_seed_bank", "no_new_clash_signature", "local_validity_proxy_pass"],
            },
            "seed_caps": seed_caps,
            "top_seed_max_fraction": 0.25,
            "first_round_top_seed_hard_line": SEED_HARD_LINE_FIRST_ROUND,
            "forbidden": [
                "broad_AE108H_centered_expansion",
                "weak_filler",
                "broad_F_boundary_expansion",
                "severe_signature_repeat",
                "sdAb_V105H_D110H",
            ],
        },
        "bank_gate": {
            "count_scope": "cumulative_canonical_unique_candidate_bank",
            "minimal_sync_gate": {
                "unique_candidate_count_min": 600,
                "pass_like_fraction_min": 0.15,
                "severe_risk_count": 0,
                "unsupported_boundary_count": 0,
                "top_seed_fraction_max": 0.30,
                "top_near_duplicate_cluster_fraction_max": 0.15,
                "four_mut_fraction_max": 0.15,
            },
            "full_stage2a_gate": {
                "unique_candidate_count_min": 1000,
                "pass_like_count_min": 200,
                "pass_like_fraction_min": 0.20,
                "severe_risk_count": 0,
                "unsupported_boundary_count": 0,
                "boundary_fraction_max": 0.70,
                "top_seed_fraction_max": 0.25,
                "top_near_duplicate_cluster_fraction_max": 0.12,
                "four_mut_fraction_max": 0.10,
            },
        },
        "negative_mask_count": int(len(negative_mask)),
        "positive_seed_bank_count": int(len(seed_bank)),
    }
    with (round_dir / "sdab_recovery_config.yaml").open("w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)
    return cfg


def write_failure_summary(stage15: pd.DataFrame, maps: dict[str, pd.DataFrame], seed_bank: pd.DataFrame, negative_mask: pd.DataFrame, round_dir: Path) -> None:
    sdab = stage15[stage15["target"].eq(TARGET)].copy()
    class_summary = sdab.groupby("refined_structure_risk_class", dropna=False).size().reset_index(name="count")
    class_summary["fraction"] = class_summary["count"] / max(1, len(sdab))
    lines = [
        "# sdAb Failure Map Summary",
        "",
        "## Baseline",
        "",
        markdown_table(class_summary),
        "",
        "## Top Failure By Seed",
        "",
        markdown_table(maps["seed"].head(20)),
        "",
        "## Top Failure By Route",
        "",
        markdown_table(maps["route"].head(20)),
        "",
        "## Top Failure By Position Pair",
        "",
        markdown_table(maps["position_pair"].head(20)),
        "",
        "## Positive Seed Bank",
        "",
        markdown_table(seed_bank.head(30)),
        "",
        "## Negative Mask",
        "",
        markdown_table(negative_mask.head(40)),
    ]
    (round_dir / "sdab_failure_map_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_generation_context(config: dict) -> dict[str, Any]:
    inputs = dry.load_inputs()
    target_cfg = config["targets"][DRY_TARGET]
    reference = inputs["reference_sequence_map"]
    mask = inputs["window_mutation_mask"]
    parent_seq = dry.chain_sequence(reference, DRY_TARGET, CHAIN)
    parent_motifs = dry.nxs_motifs(parent_seq)
    window = mask[mask.window_id == WINDOW_ID].copy().sort_values("pos")
    positions = {int(r.pos): str(r.aa) for _, r in window.iterrows()}
    regions = {int(r.pos): str(r.region) for _, r in window.iterrows()}
    return {
        "inputs": inputs,
        "target_cfg": target_cfg,
        "parent_seq": parent_seq,
        "parent_motifs": parent_motifs,
        "positions": positions,
        "regions": regions,
    }


def choose_seed_set(seed_bank: pd.DataFrame, rng: random.Random) -> tuple[int, ...]:
    allowed = seed_bank[seed_bank["allowed_max_fraction"].astype(float).gt(0)].copy()
    allowed = allowed[~allowed["seed_set"].astype(str).str.contains("AE108H", regex=False, na=False)].copy()
    if allowed.empty:
        return (rng.choice([100, 102, 110, 111]),)
    seed_rows = []
    for _, row in allowed.iterrows():
        seed = str(row["seed_set"])
        positions = tuple(int(x) for x in re.findall(r"\d+", seed))
        if not positions:
            continue
        if set(positions) == {105, 110}:
            continue
        weight = 1.0 + float(row.get("pass_like_count", 0)) * 1.5 + float(row.get("boundary_count", 0)) * 0.15
        if len(positions) == 1:
            weight *= 1.6
        if len(positions) > 2:
            continue
        seed_rows.append((positions, weight))
    if not seed_rows:
        return (rng.choice([100, 102, 110, 111]),)
    total = sum(w for _, w in seed_rows)
    x = rng.random() * total
    acc = 0.0
    for positions, weight in seed_rows:
        acc += weight
        if x <= acc:
            return positions
    return seed_rows[-1][0]


def risky_position_pairs(negative_mask: pd.DataFrame) -> set[tuple[int, int]]:
    pairs = set()
    for _, row in negative_mask.iterrows():
        if str(row.get("mask_type", "")) != "position_pair_cap_or_hold":
            continue
        nums = [int(x) for x in re.findall(r"\d+", str(row.get("position_set", "")))]
        if len(nums) == 2:
            pairs.add(tuple(sorted(nums)))
    return pairs


def mutation_has_risky_pair(muts: list[dry.Mutation], pairs: set[tuple[int, int]]) -> bool:
    positions = sorted(m.pos for m in muts)
    return any(tuple(sorted(pair)) in pairs for pair in itertools.combinations(positions, 2))


def hard_filter_recovery(target_cfg: dict, muts: list[dry.Mutation], risky_pairs: set[tuple[int, int]]) -> str | None:
    if any(m.pos in set(target_cfg["hard_protect_positions"]) for m in muts):
        return "hard_protect_mutation"
    if dry.forbidden_pair_fail(DRY_TARGET, muts):
        return "forbidden_pair"
    if sum(m.is_his() for m in muts) > 2:
        return "his_count"
    if len(muts) > 4:
        return "mutation_count"
    if any(m.mutant == "H" and m.pos == 108 for m in muts):
        return "broad_AE108H_hold"
    if len(muts) == 1 and not muts[0].is_his():
        return "single_non_his_rescue_only"
    return None


def choose_recovery_mutations(
    route: str,
    seed_bank: pd.DataFrame,
    negative_mask: pd.DataFrame,
    context: dict[str, Any],
    rng: random.Random,
    four_mut_current: int,
    raw_generated: int,
) -> list[dry.Mutation]:
    target_cfg = context["target_cfg"]
    positions = context["positions"]
    regions = context["regions"]
    hard = set(target_cfg["hard_protect_positions"])
    seed_positions = list(choose_seed_set(seed_bank, rng))
    if set(seed_positions) == {105, 110}:
        seed_positions = [105, rng.choice([100, 102, 111])]
    seed_positions = [p for p in seed_positions if p not in hard and positions.get(p) != "H" and p != 108]
    if not seed_positions:
        seed_positions = [rng.choice([100, 102, 110, 111])]

    muts = [dry.Mutation(CHAIN, p, positions[p], "H") for p in sorted(seed_positions)]

    rescue_positions = [
        p
        for p in sorted(positions)
        if p not in hard and p not in seed_positions and positions[p] not in {"C", "H"}
    ]
    fr_positions = [p for p in rescue_positions if not str(regions.get(p, "")).startswith("CDR")]
    cdr_positions = [p for p in rescue_positions if str(regions.get(p, "")).startswith("CDR")]
    local_positions = [
        p
        for p in rescue_positions
        if any(2 <= abs(p - seed) <= 9 for seed in seed_positions)
    ]

    if route == "ProteinMPNN_seeded_rescue":
        pool = sorted(dict.fromkeys(local_positions + fr_positions + rescue_positions))
    elif route == "wetlab_informed_expansion":
        pool = sorted(dict.fromkeys(local_positions + rescue_positions))
    elif route == "structure_or_interface_guided":
        pool = sorted(dict.fromkeys(fr_positions + local_positions + cdr_positions + rescue_positions))
    else:
        pool = sorted(dict.fromkeys(local_positions + fr_positions + rescue_positions))
    pool = [p for p in pool if p not in seed_positions]

    desired_total = 2 if len(seed_positions) == 1 else 3
    if route == "His_plus_rescue":
        desired_total = 3
    if route == "ProteinMPNN_seeded_rescue" and rng.random() < 0.35:
        desired_total = 3
    if raw_generated > 0 and four_mut_current / max(1, raw_generated) < 0.05 and rng.random() < 0.05:
        desired_total = 4
    desired_total = max(len(muts), min(4, desired_total))
    rescue_count = max(0, desired_total - len(muts))
    rescue_count = min(rescue_count, len(pool))
    if rescue_count <= 0:
        return sorted(muts, key=lambda m: m.pos)

    rescue_sites = sorted(rng.sample(pool, rescue_count))
    for pos in rescue_sites:
        alphabet = [aa for aa in CONSERVATIVE_RESCUE_AA if aa != positions[pos] and aa not in {"C", "H"}]
        if not alphabet:
            continue
        aa = rng.choice(alphabet)
        muts.append(dry.Mutation(CHAIN, pos, positions[pos], aa))
    return sorted(muts, key=lambda m: m.pos)


def weighted_route(rng: random.Random) -> str:
    x = rng.random()
    acc = 0.0
    for route in RECOVERY_ROUTES:
        acc += ROUTE_WEIGHTS[route]
        if x <= acc:
            return route
    return RECOVERY_ROUTES[-1]


def select_without_cluster_caps(df: pd.DataFrame, target_n: int) -> pd.DataFrame:
    selected: list[pd.Series] = []
    seed_counts: Counter[str] = Counter()
    four_count = 0
    seed_default_cap = max(1, int(math.floor(DEFAULT_SEED_CAP * target_n)))
    four_cap = max(1, int(math.floor(0.05 * target_n)))
    for _, row in df.iterrows():
        if len(selected) >= target_n:
            break
        seed = norm_text(row.get("his_seed_set", "none")) or "none"
        is_four = int(row.get("mutation_count", 0)) >= 4
        seed_cap = int(math.floor(PRIORITY_SEED_CAPS.get(seed, DEFAULT_SEED_CAP) * target_n))
        seed_cap = max(1, min(seed_cap, seed_default_cap))
        if seed_counts[seed] >= seed_cap:
            continue
        if is_four and four_count >= four_cap:
            continue
        selected.append(row)
        seed_counts[seed] += 1
        four_count += int(is_four)
    out = pd.DataFrame(selected) if selected else df.head(0).copy()
    if len(out) < target_n:
        remainder = df[~df["variant_id"].isin(set(out["variant_id"]))].copy() if not out.empty else df.copy()
        filler = remainder.head(target_n - len(out)).copy()
        out = pd.concat([out, filler], ignore_index=True)
    return out.head(target_n).copy()


def seed_options_for_generation(seed_bank: pd.DataFrame) -> list[tuple[int, ...]]:
    options: list[tuple[int, ...]] = []
    for _, row in seed_bank.iterrows():
        seed = str(row.get("seed_set", ""))
        if not seed or "AE108H" in seed:
            continue
        if float(row.get("allowed_max_fraction", 0)) <= 0:
            continue
        positions = tuple(sorted(int(x) for x in re.findall(r"\d+", seed)))
        if not positions or len(positions) > 2:
            continue
        if set(positions) == {105, 110}:
            continue
        if len(positions) > 1 and int(row.get("pass_like_count", 0)) == 0 and int(row.get("boundary_count", 0)) == 0:
            continue
        if positions not in options:
            options.append(positions)
    fallback = [(110, 111), (110,), (111,), (100, 111), (100, 110), (100,), (102,), (105,)]
    for positions in fallback:
        if positions not in options and 108 not in positions and set(positions) != {105, 110}:
            options.append(positions)
    return options


def rescue_pool_for_route(route: str, seed_positions: tuple[int, ...], context: dict[str, Any]) -> list[int]:
    target_cfg = context["target_cfg"]
    positions = context["positions"]
    regions = context["regions"]
    hard = set(target_cfg["hard_protect_positions"])
    rescue_positions = [
        p
        for p in sorted(positions)
        if p not in hard and p not in seed_positions and positions[p] not in {"C", "H"}
    ]
    fr_positions = [p for p in rescue_positions if not str(regions.get(p, "")).startswith("CDR")]
    cdr_positions = [p for p in rescue_positions if str(regions.get(p, "")).startswith("CDR")]
    local_positions = [
        p
        for p in rescue_positions
        if any(2 <= abs(p - seed) <= 9 for seed in seed_positions)
    ]
    if route == "ProteinMPNN_seeded_rescue":
        return sorted(dict.fromkeys(local_positions + fr_positions + rescue_positions))
    if route == "wetlab_informed_expansion":
        return sorted(dict.fromkeys(local_positions + rescue_positions))
    if route == "structure_or_interface_guided":
        return sorted(dict.fromkeys(fr_positions + local_positions + cdr_positions + rescue_positions))
    return sorted(dict.fromkeys(local_positions + fr_positions + rescue_positions))


def deterministic_candidate_stream(
    route: str,
    seed_positions: tuple[int, ...],
    context: dict[str, Any],
    route_index: int,
):
    positions = context["positions"]
    his_muts = [dry.Mutation(CHAIN, p, positions[p], "H") for p in seed_positions]
    pool = rescue_pool_for_route(route, seed_positions, context)
    total_orders = [2, 3, 4] if len(seed_positions) == 1 else [3, 4]
    for total_order in total_orders:
        rescue_count = total_order - len(seed_positions)
        if rescue_count < 0 or rescue_count > len(pool):
            continue
        for combo in itertools.combinations(pool, rescue_count):
            aa_offsets = range(3 if rescue_count <= 2 else 2)
            for aa_offset in aa_offsets:
                muts = list(his_muts)
                for j, pos in enumerate(combo):
                    alphabet = [aa for aa in CONSERVATIVE_RESCUE_AA if aa != positions[pos] and aa not in {"C", "H"}]
                    if not alphabet:
                        continue
                    aa = alphabet[(pos + aa_offset * 5 + j * 7 + route_index * 3) % len(alphabet)]
                    muts.append(dry.Mutation(CHAIN, pos, positions[pos], aa))
                yield sorted(muts, key=lambda m: m.pos)


def generate_recovery_pool(
    seed_bank: pd.DataFrame,
    negative_mask: pd.DataFrame,
    round_dir: Path,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base_config = dry.read_config()
    cfg_text = yaml.safe_dump(base_config, sort_keys=True, allow_unicode=True) + f"_sdab_recovery_{args.round}"
    cfg_hash = dry.sha(cfg_text, 12)
    context = load_generation_context(base_config)
    target_cfg = context["target_cfg"]
    risky_pairs_set = risky_position_pairs(negative_mask)
    by_sequence: dict[str, dict[str, Any]] = {}
    rescue_by_variant: dict[str, list[dict[str, Any]]] = defaultdict(list)
    audit = Counter()
    audit_by_reason = Counter()
    audit_by_route = Counter()
    proposed_mutation_keys: set[str] = set()
    i = 0
    seed_options = seed_options_for_generation(seed_bank)
    route_targets = {
        route: int(round(args.raw_target * ROUTE_WEIGHTS[route]))
        for route in RECOVERY_ROUTES
    }
    route_targets[RECOVERY_ROUTES[0]] += args.raw_target - sum(route_targets.values())
    raw_4mut_cap = int(math.floor(0.05 * args.raw_target))
    seed_raw_caps = {
        seed: max(1, int(math.floor(1.50 * PRIORITY_SEED_CAPS.get(";".join(f"A{context['positions'][p]}{p}H" for p in seed), DEFAULT_SEED_CAP) * args.raw_target)))
        for seed in seed_options
    }
    global_seed_counts: Counter[tuple[int, ...]] = Counter()

    for route_index, route in enumerate(RECOVERY_ROUTES):
        route_target = route_targets[route]
        route_generated = 0
        while route_generated < route_target:
            made_progress = False
            for seed_positions in seed_options:
                if route_generated >= route_target:
                    break
                if global_seed_counts[seed_positions] >= seed_raw_caps.get(seed_positions, int(DEFAULT_SEED_CAP * args.raw_target)):
                    continue
                stream = deterministic_candidate_stream(route, seed_positions, context, route_index)
                for muts in stream:
                    if route_generated >= route_target:
                        break
                    if global_seed_counts[seed_positions] >= seed_raw_caps.get(seed_positions, int(DEFAULT_SEED_CAP * args.raw_target)):
                        break
                    if len(muts) >= 4 and audit["raw_4mut_count"] >= raw_4mut_cap:
                        continue
                    mut_key = dry.mutation_key(muts)
                    if mut_key in proposed_mutation_keys:
                        audit_by_reason["duplicate_mutation_key_pre_generation"] += 1
                        continue
                    proposed_mutation_keys.add(mut_key)
                    i += 1
                    audit["raw_generated_count"] += 1
                    route_generated += 1
                    made_progress = True
                    global_seed_counts[seed_positions] += 1
                    audit_by_route[(route, "raw")] += 1
                    if len(muts) >= 4:
                        audit["raw_4mut_count"] += 1
                    reason = hard_filter_recovery(target_cfg, muts, risky_pairs_set)
                    if reason:
                        audit_by_reason[reason] += 1
                        audit_by_route[(route, f"hard_fail:{reason}")] += 1
                        continue
                    row, rescue_rows, failure = dry.make_candidate(
                        DRY_TARGET,
                        CHAIN,
                        WINDOW_ID,
                        route,
                        context["parent_seq"],
                        muts,
                        cfg_hash,
                        i,
                        context["parent_motifs"],
                    )
                    if failure:
                        audit_by_reason[failure] += 1
                        audit_by_route[(route, f"make_candidate_fail:{failure}")] += 1
                        continue
                    audit["post_hard_count"] += 1
                    seq_hash = row["sequence_hash"]
                    if seq_hash in by_sequence:
                        existing = by_sequence[seq_hash]
                        routes = set(str(existing["all_source_routes"]).split(";")) | {route}
                        records = set(str(existing["all_generation_records"]).split(";")) | {row["generation_record_id"]}
                        existing["all_source_routes"] = ";".join(sorted(routes))
                        existing["all_generation_records"] = ";".join(sorted(records))
                        existing["duplicate_source_count"] = int(existing["duplicate_source_count"]) + 1
                        existing["exact_duplicate_status"] = "merged_duplicate"
                        audit_by_reason["exact_duplicate_merged"] += 1
                        continue
                    by_sequence[seq_hash] = row
                    rescue_by_variant[row["variant_id"]].extend(rescue_rows)
            if not made_progress:
                break
        audit_by_route[(route, "route_target")] = route_target
        audit_by_route[(route, "route_generated")] = route_generated

    proposals_all = pd.DataFrame(by_sequence.values())
    if proposals_all.empty:
        raise RuntimeError("No sdAb recovery proposals passed generation hard filters")
    if len(proposals_all) < args.post_hard_unique_target:
        raise RuntimeError(
            f"Only {len(proposals_all)} exact-unique legal recovery proposals; "
            f"expected at least {args.post_hard_unique_target}"
        )
    proposals_all = proposals_all.sort_values(
        ["hit_likelihood_score_v0", "neutral_retention_score", "acidic_release_support_score", "sequence_hash"],
        ascending=[False, False, False, True],
    )
    proposals = select_without_cluster_caps(proposals_all, args.post_hard_unique_target)
    proposals = dry.assign_near_duplicate_clusters(proposals)
    proposals["source_universe"] = f"sdAb_recovery_{args.round}_targeted_generation"
    proposals["recovery_round"] = args.round
    proposals["recovery_policy_status"] = "pre_stage1_candidate"
    proposals["generation_policy_tags"] = "seed_capped;geometry_aware;no_broad_AE108H;no_weak_filler"
    rescue_rows_all = []
    for variant_id in proposals["variant_id"]:
        rescue_rows_all.extend(rescue_by_variant.get(variant_id, []))
    rescue = pd.DataFrame(rescue_rows_all)

    route_rows = []
    for (route, event), count in sorted(audit_by_route.items()):
        route_rows.append({"primary_generation_route": route, "event": event, "count": count})
    audit["post_exact_unique_before_selection_count"] = len(proposals_all)
    audit["post_exact_unique_selected_count"] = len(proposals)
    audit_rows = [{"metric": k, "value": v} for k, v in audit.items()]
    audit_rows.extend({"metric": f"hard_or_generation_filter:{k}", "value": v} for k, v in audit_by_reason.items())
    audit_df = pd.DataFrame(audit_rows)
    route_audit = pd.DataFrame(route_rows)

    dry.write_csv(proposals, round_dir / "sdab_targeted_generation_pool.csv")
    dry.write_csv(rescue, round_dir / "sdab_targeted_generation_rescue_mutations.csv")
    csv(audit_df, round_dir / "sdab_targeted_generation_audit.csv")
    csv(route_audit, round_dir / "sdab_targeted_generation_route_audit.csv")
    return proposals, rescue, audit_df


def cheap_gate_status(row: pd.Series, negative_mask: pd.DataFrame) -> tuple[str, str, float]:
    reasons: list[str] = []
    score = float(row.get("hit_likelihood_score_v0", 0.0))
    seed = norm_text(row.get("his_seed_set", ""))
    muts = mutation_labels(row.get("mutation_list", ""))
    positions = mutation_positions(row.get("mutation_list", ""))
    if "AE108H" in seed:
        reasons.append("AE108H_hold")
        score -= 1.0
    if "AV105H" in seed and "AD110H" in seed:
        reasons.append("forbidden_pair")
        score -= 5.0
    if int(row.get("mutation_count", 0)) >= 4:
        reasons.append("4mut_requires_low_quota")
        score -= 0.2
    if int(row.get("His_count", 0)) == 1:
        reasons.append("one_His_rescue_preferred")
        score += 0.4
    if int(row.get("His_count", 0)) == 2:
        reasons.append("two_His_capped")
        score += 0.1
    if int(row.get("rescue_count", 0)) > 0:
        reasons.append("has_rescue")
        score += 0.25
    if "ProteinMPNN" in norm_text(row.get("all_source_routes", "")):
        reasons.append("mpnn_rescue_provenance")
        score += 0.15
    if norm_text(row.get("primary_generation_route", "")) in {"His_plus_rescue", "ProteinMPNN_seeded_rescue"}:
        reasons.append("recovery_priority_route")
        score += 0.2
    risky_pairs = risky_position_pairs(negative_mask)
    for pair in itertools.combinations(positions, 2):
        if tuple(sorted(pair)) in risky_pairs:
            reasons.append("negative_position_pair_signature")
            score -= 0.65
            break
    if "forbidden_pair" in reasons or "AE108H_hold" in reasons:
        return "reject_by_cheap_gate", ";".join(sorted(set(reasons))), score
    if "one_His_rescue_preferred" in reasons or "mpnn_rescue_provenance" in reasons or "recovery_priority_route" in reasons:
        return "stage1_minibatch_eligible", ";".join(sorted(set(reasons))), score
    return "tier1_recovery_candidate", ";".join(sorted(set(reasons))) if reasons else "basic_recovery_candidate", score


def select_with_caps(
    df: pd.DataFrame,
    target_n: int,
    seed_cap_fraction: float,
    cluster_cap_fraction: float,
    four_mut_fraction: float,
    absolute_seed_caps: dict[str, int] | None = None,
) -> pd.DataFrame:
    selected: list[pd.Series] = []
    seed_counts: Counter[str] = Counter()
    cluster_counts: Counter[str] = Counter()
    four_mut_count = 0
    seed_cap = max(1, int(math.floor(seed_cap_fraction * target_n)))
    cluster_cap = max(1, int(math.floor(cluster_cap_fraction * target_n)))
    four_cap = max(1, int(math.floor(four_mut_fraction * target_n)))
    for _, row in df.iterrows():
        if len(selected) >= target_n:
            break
        seed = norm_text(row.get("his_seed_set", "none")) or "none"
        cluster = norm_text(row.get("near_duplicate_cluster_id", "none")) or "none"
        is_four = int(row.get("mutation_count", 0)) >= 4
        cap = int(math.floor(PRIORITY_SEED_CAPS.get(seed, DEFAULT_SEED_CAP) * target_n))
        cap = max(1, min(seed_cap, cap))
        if absolute_seed_caps is not None and seed in absolute_seed_caps:
            cap = max(0, min(cap, int(absolute_seed_caps[seed])))
        if seed_counts[seed] >= cap:
            continue
        if cluster_counts[cluster] >= cluster_cap:
            continue
        if is_four and four_mut_count >= four_cap:
            continue
        selected.append(row)
        seed_counts[seed] += 1
        cluster_counts[cluster] += 1
        four_mut_count += int(is_four)
    return pd.DataFrame(selected) if selected else df.head(0).copy()


def current_bank_seed_pressure(current_bank_path: str) -> dict[str, int]:
    if not current_bank_path:
        return {}
    path = Path(current_bank_path)
    if not path.exists():
        return {}
    bank = pd.read_csv(path, low_memory=False)
    if bank.empty or "his_seed_set" not in bank.columns:
        return {}
    total = len(bank)
    counts = bank.groupby("his_seed_set", dropna=False).size().to_dict()
    caps: dict[str, int] = {}
    for seed, count in counts.items():
        seed = norm_text(seed) or "none"
        fraction = int(count) / max(1, total)
        if fraction > 0.24:
            caps[seed] = 50
        elif fraction > 0.22:
            caps[seed] = 100
    return caps


def build_tier1_and_minibatch(proposals: pd.DataFrame, negative_mask: pd.DataFrame, round_dir: Path, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = proposals.copy()
    statuses = df.apply(lambda row: cheap_gate_status(row, negative_mask), axis=1)
    df["cheap_gate_status"] = [x[0] for x in statuses]
    df["cheap_gate_reason"] = [x[1] for x in statuses]
    df["recovery_rank_score"] = [x[2] for x in statuses]
    df["canonical_recovery_sequence_hash"] = df["sequence"].astype(str).map(lambda x: dry.sha(x))
    df["canonical_recovery_mutation_key"] = df["mutation_list"].fillna("").astype(str)
    eligible = df[~df["cheap_gate_status"].eq("reject_by_cheap_gate")].copy()
    eligible = eligible.sort_values(
        ["cheap_gate_status", "recovery_rank_score", "neutral_retention_score", "acidic_release_support_score"],
        ascending=[True, False, False, False],
    )
    tier1 = select_with_caps(
        eligible,
        args.tier1_proposal_target,
        seed_cap_fraction=0.30,
        cluster_cap_fraction=0.08,
        four_mut_fraction=0.15,
        absolute_seed_caps=current_bank_seed_pressure(args.current_bank),
    )
    if len(tier1) < args.tier1_proposal_target:
        remainder = eligible[~eligible["variant_id"].isin(set(tier1["variant_id"]))].copy()
        filler = select_with_caps(
            remainder,
            args.tier1_proposal_target - len(tier1),
            seed_cap_fraction=0.35,
            cluster_cap_fraction=0.10,
            four_mut_fraction=0.20,
            absolute_seed_caps=current_bank_seed_pressure(args.current_bank),
        )
        tier1 = pd.concat([tier1, filler], ignore_index=True)
    tier1["tier1_recovery_role"] = f"{args.round}_tier1_proposal"

    minibatch_parts: list[pd.DataFrame] = []
    remaining = tier1.copy()

    def take(label: str, predicate: pd.Series, n: int) -> None:
        nonlocal remaining
        pool = remaining[predicate.reindex(remaining.index).fillna(False)].copy()
        pool = pool.sort_values(["recovery_rank_score", "neutral_retention_score"], ascending=[False, False])
        picked = select_with_caps(pool, n, 0.30, 0.05, 0.15, absolute_seed_caps=current_bank_seed_pressure(args.current_bank))
        if not picked.empty:
            picked = picked.copy()
            picked["stage1_minibatch_role"] = label
            minibatch_parts.append(picked)
            remaining = remaining[~remaining["variant_id"].isin(set(picked["variant_id"]))].copy()

    take("pass_like_prioritized_proxy", tier1["cheap_gate_status"].eq("stage1_minibatch_eligible"), 500)
    take("supported_boundary_proxy", tier1["cheap_gate_reason"].str.contains("has_rescue|recovery_priority_route|mpnn", regex=True, na=False), 200)
    take("one_His_rescue", tier1["His_count"].astype(int).eq(1) & tier1["rescue_count"].astype(int).gt(0), 150)
    take("seed_diversity_representatives", pd.Series(True, index=tier1.index), 100)
    take("controls_or_prior_positive_proxy", tier1["his_seed_set"].isin(WETLAB_SUPPORTED_SEEDS), 50)
    mini = pd.concat(minibatch_parts, ignore_index=True) if minibatch_parts else tier1.head(0).copy()
    if len(mini) < args.stage1_mini_batch_target:
        filler = remaining.sort_values(["recovery_rank_score"], ascending=False).head(args.stage1_mini_batch_target - len(mini)).copy()
        if not filler.empty:
            filler["stage1_minibatch_role"] = "quota_shortfall_best_remaining"
            mini = pd.concat([mini, filler], ignore_index=True)
    mini = mini.head(args.stage1_mini_batch_target).copy()
    mini["normalized_mutation_list"] = mini["mutation_list"]
    tier1["normalized_mutation_list"] = tier1["mutation_list"]
    dry.write_csv(tier1, round_dir / "sdab_tier1_recovery_results.csv")
    dry.write_csv(mini, round_dir / "sdab_stage1_minibatch_list.csv")
    return tier1, mini


def summarize_selection(df: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    total = len(df)
    if total == 0:
        return pd.DataFrame([{"scope": label, "metric": "total", "value": 0}])
    rows.extend(
        [
            {"scope": label, "metric": "total", "value": total},
            {"scope": label, "metric": "top_seed_fraction", "value": df.groupby("his_seed_set", dropna=False).size().max() / total},
            {"scope": label, "metric": "top_near_duplicate_cluster_fraction", "value": df.groupby("near_duplicate_cluster_id", dropna=False).size().max() / total},
            {"scope": label, "metric": "four_mut_fraction", "value": pd.to_numeric(df["mutation_count"], errors="coerce").ge(4).sum() / total},
            {"scope": label, "metric": "one_his_fraction", "value": pd.to_numeric(df["His_count"], errors="coerce").eq(1).sum() / total},
            {"scope": label, "metric": "two_his_fraction", "value": pd.to_numeric(df["His_count"], errors="coerce").eq(2).sum() / total},
        ]
    )
    return pd.DataFrame(rows)


def write_preflight(
    stage15: pd.DataFrame,
    seed_bank: pd.DataFrame,
    negative_mask: pd.DataFrame,
    proposals: pd.DataFrame | None,
    tier1: pd.DataFrame | None,
    mini: pd.DataFrame | None,
    round_dir: Path,
    args: argparse.Namespace,
) -> None:
    sdab = stage15[stage15["target"].eq(TARGET)].copy()
    baseline = sdab.groupby("refined_structure_risk_class", dropna=False).size().reset_index(name="count")
    baseline["fraction"] = baseline["count"] / max(1, len(sdab))
    checks = [
        {"check": "stage2a_compute_locked", "status": "PASS", "details": "no Stage-2A compute launched by this script"},
        {"check": "tier2_heavy_locked", "status": "PASS", "details": "no Tier2-heavy launched"},
        {"check": "final_10k_locked", "status": "PASS", "details": "no final selection launched"},
        {"check": "positive_seed_bank_nonempty", "status": "PASS" if len(seed_bank) > 0 else "FAIL", "details": len(seed_bank)},
        {"check": "negative_mask_nonempty", "status": "PASS" if len(negative_mask) > 0 else "FAIL", "details": len(negative_mask)},
    ]
    if proposals is not None:
        audit_path = round_dir / "sdab_targeted_generation_audit.csv"
        raw_count = len(proposals)
        unique_before_selection = len(proposals)
        if audit_path.exists():
            audit_table = pd.read_csv(audit_path)
            audit_map = dict(zip(audit_table["metric"].astype(str), audit_table["value"]))
            raw_count = int(float(audit_map.get("raw_generated_count", raw_count)))
            unique_before_selection = int(float(audit_map.get("post_exact_unique_before_selection_count", unique_before_selection)))
        checks.extend(
            [
                {"check": "raw_generated_30k_to_50k", "status": "PASS" if 30_000 <= raw_count <= 50_000 else "WARN", "details": raw_count},
                {
                    "check": "post_hard_unique_before_selection_ge_15k",
                    "status": "PASS" if unique_before_selection >= 15_000 else "WARN",
                    "details": unique_before_selection,
                },
                {
                    "check": "post_hard_unique_selected_15k_to_20k",
                    "status": "PASS" if 15_000 <= len(proposals) <= 20_000 else "WARN",
                    "details": len(proposals),
                },
                {
                    "check": "no_AE108H_generation",
                    "status": "PASS" if not proposals["his_seed_set"].astype(str).str.contains("AE108H", regex=False, na=False).any() else "FAIL",
                    "details": int(proposals["his_seed_set"].astype(str).str.contains("AE108H", regex=False, na=False).sum()),
                },
            ]
        )
    if tier1 is not None:
        checks.append(
            {
                "check": "tier1_proposal_3k_to_4k",
                "status": "PASS" if 3_000 <= len(tier1) <= 4_000 else "WARN",
                "details": len(tier1),
            }
        )
    if mini is not None:
        seed_frac = mini.groupby("his_seed_set", dropna=False).size().max() / max(1, len(mini))
        cluster_frac = mini.groupby("near_duplicate_cluster_id", dropna=False).size().max() / max(1, len(mini))
        four_frac = pd.to_numeric(mini["mutation_count"], errors="coerce").ge(4).sum() / max(1, len(mini))
        checks.extend(
            [
                {"check": "stage1_minibatch_eq_1000", "status": "PASS" if len(mini) == args.stage1_mini_batch_target else "WARN", "details": len(mini)},
                {"check": "stage1_minibatch_top_seed_le_35pct", "status": "PASS" if seed_frac <= 0.35 else "WARN", "details": round(seed_frac, 4)},
                {"check": "stage1_minibatch_top_cluster_le_20pct", "status": "PASS" if cluster_frac <= 0.20 else "WARN", "details": round(cluster_frac, 4)},
                {"check": "stage1_minibatch_4mut_le_20pct", "status": "PASS" if four_frac <= 0.20 else "WARN", "details": round(four_frac, 4)},
            ]
        )
    check_df = pd.DataFrame(checks)
    csv(check_df, round_dir / "sdab_recovery_preflight_checks.csv")

    metric_parts = []
    if proposals is not None:
        metric_parts.append(summarize_selection(proposals, "targeted_generation_pool"))
    if tier1 is not None:
        metric_parts.append(summarize_selection(tier1, "tier1_recovery_results"))
    if mini is not None:
        metric_parts.append(summarize_selection(mini, "stage1_minibatch_list"))
    metrics = pd.concat(metric_parts, ignore_index=True) if metric_parts else pd.DataFrame()
    if not metrics.empty:
        csv(metrics, round_dir / "sdab_recovery_preflight_metrics.csv")

    status_counts = check_df["status"].value_counts().to_dict()
    gate = "PASS" if status_counts.get("FAIL", 0) == 0 else "FAIL"
    lines = [
        f"# sdAb Recovery {args.round} Preflight Report",
        "",
        f"Generated at: `{datetime.now().isoformat(timespec='seconds')}`",
        f"Host: `{socket.gethostname()}`",
        "",
        "## Gate",
        "",
        f"Preflight status: `{gate}`",
        "",
        "This preflight does not unlock Stage-2A compute, Tier2-heavy, or final 10K selection.",
        "",
        "## Baseline sdAb Stage-1.5 Status",
        "",
        markdown_table(baseline),
        "",
        "## Checks",
        "",
        markdown_table(check_df),
    ]
    if not metrics.empty:
        lines.extend(["", "## Metrics", "", markdown_table(metrics)])
    if mini is not None:
        lines.extend(
            [
                "",
                "## Mini-batch Role Distribution",
                "",
                markdown_table(mini.groupby("stage1_minibatch_role", dropna=False).size().reset_index(name="count")),
                "",
                "## Mini-batch Seed Distribution Top 20",
                "",
                markdown_table(mini.groupby("his_seed_set", dropna=False).size().reset_index(name="count").sort_values("count", ascending=False).head(20)),
            ]
        )
    (round_dir / "sdab_recovery_preflight_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(round_dir: Path, args: argparse.Namespace) -> None:
    manifest = {
        "run_id": f"sdab_recovery_{args.round}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "mode": args.mode,
        "stage2a_compute_started": False,
        "tier2_heavy_started": False,
        "final_10k_selection_started": False,
        "git_commit": git_commit(),
        "host": socket.gethostname(),
        "operator": subprocess.getoutput("whoami").strip(),
        "raw_target": args.raw_target,
        "post_hard_unique_target": args.post_hard_unique_target,
        "tier1_proposal_target": args.tier1_proposal_target,
        "stage1_mini_batch_target": args.stage1_mini_batch_target,
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    with (round_dir / "sdab_recovery_round_manifest.yaml").open("w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "unavailable"


def main() -> None:
    args = parse_args()
    out_root = Path(args.out_root)
    round_dir, onee62_dir = ensure_dirs(out_root, args.round)
    stage15 = pd.read_csv(args.stage15, low_memory=False)
    stage2a = pd.read_csv(args.stage2a, low_memory=False)
    maps = build_failure_maps(stage15, round_dir)
    negative_mask = build_negative_mask(maps, round_dir, args.round)
    seed_bank = build_positive_seed_bank(stage15, round_dir)
    recovery_config(seed_bank, negative_mask, round_dir, args)
    build_current_banks(stage15, stage2a, round_dir, onee62_dir)
    write_failure_summary(stage15, maps, seed_bank, negative_mask, round_dir)
    proposals = tier1 = mini = None
    if args.mode in {"generate", "all"}:
        proposals, _, _ = generate_recovery_pool(seed_bank, negative_mask, round_dir, args)
        tier1, mini = build_tier1_and_minibatch(proposals, negative_mask, round_dir, args)
    write_preflight(stage15, seed_bank, negative_mask, proposals, tier1, mini, round_dir, args)
    write_manifest(round_dir, args)
    print(
        json.dumps(
            {
                "status": "complete",
                "round_dir": str(round_dir),
                "positive_seed_bank": len(seed_bank),
                "negative_mask": len(negative_mask),
                "generated_pool": 0 if proposals is None else len(proposals),
                "tier1_recovery_results": 0 if tier1 is None else len(tier1),
                "stage1_minibatch": 0 if mini is None else len(mini),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
