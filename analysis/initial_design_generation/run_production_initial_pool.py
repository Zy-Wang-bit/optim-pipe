#!/usr/bin/env python3
"""Generate the production initial pool without unlocking Tier 1 or Tier 2.

This script reuses the dry-run rule engine for legal proposal generation, but
adds the production-only layer required by the stress-run review:

oversample legal proposals -> hard filter -> exact dedup -> Hamming/Jaccard
clustering -> route/seed/rescue/cluster capped diversity selection -> audit.

It does not run ESM, FoldX actual, PyRosetta, AF3, SimpleFold, explicit glycan
modeling, Rosetta dddG_elec, or MD.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import socket
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT
OUT = dry.OUT
PRODUCTION = OUT / "production_initial_pool"
REPORTS = OUT / "reports"
REPORT_NAME = "production_pool_validation_report.md"


def production_config(config: dict) -> dict:
    cfg = yaml.safe_load(yaml.safe_dump(config, sort_keys=False, allow_unicode=True))
    cfg["run_mode"] = "production_initial_pool"
    cfg["dry_run_gate_thresholds"] = cfg.get("production_pool_guardrails", cfg["dry_run_gate_thresholds"])
    status = cfg.setdefault("status", {})
    status["approved_for_production_initial_pool_generation"] = "conditional_go"
    status["approved_for_tier1"] = False
    status["approved_for_tier2"] = False
    for target_cfg in cfg["targets"].values():
        target_cfg["dry_run_effective_design_pool"] = int(target_cfg["production_effective_design_pool"])
    return cfg


def route_targets_for(target_cfg: dict, config: dict, target_size: int) -> dict[str, int]:
    quotas = target_cfg.get("route_soft_quota", config["routes"]["route_soft_quota"])
    route_targets = {
        route: int(round(target_size * float(quotas[route])))
        for route in dry.PRIMARY_ROUTES
    }
    route_targets[dry.PRIMARY_ROUTES[0]] += target_size - sum(route_targets.values())
    return route_targets


def hard_filter_reason(target: str, target_cfg: dict, muts: list[dry.Mutation]) -> str | None:
    if not muts:
        return "empty_mutation"
    if any(m.pos in set(target_cfg["hard_protect_positions"]) for m in muts):
        return "hard_protect_mutation"
    if dry.forbidden_pair_fail(target, muts):
        return "forbidden_pair"
    if sum(m.is_his() for m in muts) > target_cfg["main_max_his"]:
        return "his_count"
    max_mutations = int(target_cfg.get("production_max_mutations", target_cfg["main_max_mutations"]))
    if len(muts) > max_mutations:
        return "mutation_count"
    if len(muts) == 1 and not muts[0].is_his():
        return "single_non_his_rescue_only"
    return None


def choose_production_mutations(
    target: str,
    chain: str,
    route: str,
    target_cfg: dict,
    positions: dict[int, str],
    regions: dict[int, str],
    rng: random.Random,
) -> list[dry.Mutation]:
    hard = set(target_cfg["hard_protect_positions"])
    his_seeds = [p for p in target_cfg["his_seed_positions"] if p not in hard and positions[p] != "H"]
    rescue_positions = [
        p for p in sorted(positions)
        if p not in hard and p not in his_seeds and positions[p] != "C"
    ]
    cdr_positions = [p for p in rescue_positions if str(regions[p]).startswith("CDR")]
    fr_positions = [p for p in rescue_positions if not str(regions[p]).startswith("CDR")]

    if target == "sdAb":
        single_probability = {
            "His_rule": 0.72,
            "His_plus_rescue": 0.62,
            "ProteinMPNN_seeded_rescue": 0.62,
            "wetlab_informed_expansion": 0.68,
            "structure_or_interface_guided": 0.68,
        }.get(route, 0.66)
    else:
        single_probability = {
            "His_rule": 0.42,
            "His_plus_rescue": 0.36,
            "ProteinMPNN_seeded_rescue": 0.34,
            "wetlab_informed_expansion": 0.40,
            "structure_or_interface_guided": 0.38,
        }.get(route, 0.38)

    use_single = rng.random() < single_probability or len(his_seeds) < 2
    if use_single:
        seed_positions = [rng.choice(his_seeds)]
    else:
        seed_positions = sorted(rng.sample(his_seeds, 2))
        if target == "sdAb" and set(seed_positions) == {105, 110}:
            alternatives = [p for p in his_seeds if p not in {105, 110}]
            seed_positions = sorted([105, rng.choice(alternatives)])

    muts = [dry.Mutation(chain, p, positions[p], "H") for p in seed_positions]

    if route == "ProteinMPNN_seeded_rescue":
        local = [
            p
            for p in rescue_positions
            if any(3 <= abs(p - seed) <= 8 for seed in seed_positions)
        ]
        pool = sorted(dict.fromkeys(local + fr_positions + rescue_positions))
    elif route == "structure_or_interface_guided":
        if target == "1E62":
            local = [
                p
                for p in rescue_positions
                if any(abs(p - seed) <= 10 for seed in seed_positions)
            ]
            pool = sorted(dict.fromkeys(cdr_positions + local + rescue_positions))
        else:
            pool = sorted(dict.fromkeys(fr_positions + cdr_positions + rescue_positions))
    elif route == "wetlab_informed_expansion":
        local = [
            p
            for p in rescue_positions
            if any(abs(p - seed) <= 12 for seed in seed_positions)
        ]
        pool = sorted(dict.fromkeys(local + rescue_positions))
    else:
        pool = rescue_positions

    if target == "sdAb":
        if len(seed_positions) == 1:
            rescue_count = 1 if rng.random() < 0.20 else 2
        else:
            rescue_count = 1 if rng.random() < 0.15 else 2
    else:
        if len(seed_positions) == 1:
            rescue_count = rng.choice([1, 2, 3])
        else:
            rescue_count = rng.choice([1, 2])

    max_mutations = int(target_cfg.get("production_max_mutations", target_cfg["main_max_mutations"]))
    rescue_count = min(rescue_count, max_mutations - len(muts), len(pool))
    if rescue_count <= 0:
        return sorted(muts, key=lambda m: m.pos)

    rescue_sites = sorted(rng.sample(pool, rescue_count))
    for pos in rescue_sites:
        residues = dry.ordered_rescue_residues(positions[pos])
        if not residues:
            continue
        aa = rng.choice(residues)
        muts.append(dry.Mutation(chain, pos, positions[pos], aa))
    return sorted(muts, key=lambda m: m.pos)


def generate_proposals(
    target: str,
    config: dict,
    inputs: dict[str, pd.DataFrame],
    cfg_hash: str,
) -> tuple[pd.DataFrame, dict[str, list[dict]], dict]:
    target_cfg = config["targets"][target]
    chain = target_cfg["chain_id"]
    window_id = target_cfg["window_id"]
    target_size = int(target_cfg["production_effective_design_pool"])
    multiplier = float(
        config.get("production_initial_pool", {})
        .get("oversample_post_hard_multiplier", {})
        .get(target, 1.5)
    )
    route_targets = route_targets_for(target_cfg, config, target_size)
    route_proposal_targets = {
        route: max(route_targets[route], int(math.ceil(route_targets[route] * multiplier)))
        for route in dry.PRIMARY_ROUTES
    }
    reference = inputs["reference_sequence_map"]
    mask = inputs["window_mutation_mask"]
    parent_seq = dry.chain_sequence(reference, target, chain)
    parent_motifs = dry.nxs_motifs(parent_seq)
    window = mask[mask.window_id == window_id].copy().sort_values("pos")
    positions = {int(r.pos): str(r.aa) for _, r in window.iterrows()}
    regions = {int(r.pos): str(r.region) for _, r in window.iterrows()}

    by_sequence: dict[str, dict] = {}
    rescue_by_variant: dict[str, list[dict]] = defaultdict(list)
    audit = {
        "raw_generated_count": 0,
        "post_generation_hard_filter_count": 0,
        "post_exact_dedup_count": 0,
        "effective_design_pool_count": 0,
        "control_anchor_pool_count": 0,
        "combined_pool_count": 0,
        "hard_filter_failures": Counter(),
        "raw_by_route": Counter(),
        "hard_fail_by_route_reason": Counter(),
        "accepted_by_route": Counter(),
        "dedup_by_route": Counter(),
        "diversity_cap_by_route": Counter(),
        "diversity_cap_failures": Counter(),
        "proposal_by_route": Counter(),
        "proposal_shortfall_by_route": Counter(),
        "route_targets": route_targets,
        "route_proposal_targets": route_proposal_targets,
    }

    i = 0
    rng = random.Random(int(config["random_seed"]) + (101 if target == "1E62" else 202))
    target_proposal_goal = int(math.ceil(target_size * multiplier))
    proposed_mutation_keys: dict[str, set[str]] = {route: set() for route in dry.PRIMARY_ROUTES}
    route_attempts: Counter[str] = Counter()
    stall_attempts: Counter[str] = Counter()
    active_routes = set(dry.PRIMARY_ROUTES)
    max_stall_attempts = {
        route: max(50000, route_proposal_targets[route] * 10)
        for route in dry.PRIMARY_ROUTES
    }
    max_total_attempts = max(100000, target_proposal_goal * 250)
    total_attempts = 0
    route_cursor = 0

    while len(by_sequence) < target_proposal_goal and active_routes and total_attempts < max_total_attempts:
        route = dry.PRIMARY_ROUTES[route_cursor % len(dry.PRIMARY_ROUTES)]
        route_cursor += 1
        if route not in active_routes:
            continue
        total_attempts += 1
        route_attempts[route] += 1
        muts = choose_production_mutations(target, chain, route, target_cfg, positions, regions, rng)
        i += 1
        proposed_key = dry.mutation_key(muts)
        if proposed_key in proposed_mutation_keys[route]:
            stall_attempts[route] += 1
        else:
            proposed_mutation_keys[route].add(proposed_key)
            audit["raw_generated_count"] += 1
            audit["raw_by_route"][route] += 1
            reason = hard_filter_reason(target, target_cfg, muts)
            if reason:
                stall_attempts[route] += 1
                audit["hard_filter_failures"][reason] += 1
                audit["hard_fail_by_route_reason"][(route, reason)] += 1
            else:
                row, rescue_rows, failure = dry.make_candidate(
                    target, chain, window_id, route, parent_seq, muts, cfg_hash, i, parent_motifs
                )
                if failure:
                    stall_attempts[route] += 1
                    audit["hard_filter_failures"][failure] += 1
                    audit["hard_fail_by_route_reason"][(route, failure)] += 1
                else:
                    audit["post_generation_hard_filter_count"] += 1
                    seq_hash = row["sequence_hash"]
                    if seq_hash in by_sequence:
                        stall_attempts[route] += 1
                        existing = by_sequence[seq_hash]
                        routes = set(existing["all_source_routes"].split(";")) | {route}
                        records = set(existing["all_generation_records"].split(";")) | {row["generation_record_id"]}
                        existing["all_source_routes"] = ";".join(sorted(routes))
                        existing["all_generation_records"] = ";".join(sorted(records))
                        existing["duplicate_source_count"] = int(existing["duplicate_source_count"]) + 1
                        existing["exact_duplicate_status"] = "merged_duplicate"
                        audit["dedup_by_route"][route] += 1
                    else:
                        by_sequence[seq_hash] = row
                        rescue_by_variant[row["variant_id"]].extend(rescue_rows)
                        audit["proposal_by_route"][route] += 1
                        stall_attempts[route] = 0

        if stall_attempts[route] >= max_stall_attempts[route]:
            active_routes.remove(route)
            print(
                json.dumps(
                    {
                        "stage": "proposal_route_deactivated",
                        "target": target,
                        "route": route,
                        "unique_legal_proposals": int(audit["proposal_by_route"].get(route, 0)),
                        "requested_route_proposals": route_proposal_targets[route],
                        "attempts": int(route_attempts[route]),
                        "stall_attempts": int(stall_attempts[route]),
                    },
                    ensure_ascii=False,
                ),
                flush=True,
            )

    for route in dry.PRIMARY_ROUTES:
        unique_count = int(audit["proposal_by_route"].get(route, 0))
        if unique_count < route_proposal_targets[route]:
            audit["proposal_shortfall_by_route"][route] = route_proposal_targets[route] - unique_count
        print(
            json.dumps(
                {
                    "stage": "proposal_route_complete",
                    "target": target,
                    "route": route,
                    "unique_legal_proposals": unique_count,
                    "requested_route_proposals": route_proposal_targets[route],
                    "attempts": int(route_attempts[route]),
                    "shortfall": int(audit["proposal_shortfall_by_route"].get(route, 0)),
                },
                ensure_ascii=False,
            ),
            flush=True,
        )

    if len(by_sequence) < target_size:
        raise RuntimeError(
            f"{target} produced {len(by_sequence)} exact-unique legal proposals; "
            f"expected at least {target_size}"
        )

    proposals = pd.DataFrame(by_sequence.values())
    print(
        json.dumps(
            {
                "stage": "assign_near_duplicate_clusters_start",
                "target": target,
                "proposal_count": len(proposals),
            },
            ensure_ascii=False,
        ),
        flush=True,
    )
    proposals = dry.assign_near_duplicate_clusters(proposals)
    audit["post_exact_dedup_count"] = len(proposals)
    return proposals, rescue_by_variant, audit


def adjust_route_targets_for_capacity(
    desired: dict[str, int],
    available: dict[str, int],
    target_size: int,
) -> dict[str, int]:
    adjusted = {
        route: min(int(desired.get(route, 0)), int(available.get(route, 0)))
        for route in dry.PRIMARY_ROUTES
    }
    deficit = target_size - sum(adjusted.values())
    while deficit > 0:
        candidates = [
            route
            for route in dry.PRIMARY_ROUTES
            if adjusted[route] < int(available.get(route, 0))
        ]
        if not candidates:
            break
        total_room = sum(int(available.get(route, 0)) - adjusted[route] for route in candidates)
        if total_room <= 0:
            break
        for route in sorted(candidates, key=lambda r: int(available.get(r, 0)) - adjusted[r], reverse=True):
            room = int(available.get(route, 0)) - adjusted[route]
            share = max(1, int(math.ceil(deficit * room / total_room)))
            add = min(room, share, deficit)
            adjusted[route] += add
            deficit -= add
            if deficit <= 0:
                break
    if sum(adjusted.values()) != target_size:
        raise RuntimeError(
            f"Only {sum(adjusted.values())} proposals can be allocated across routes; expected {target_size}"
        )
    return adjusted


def rescue_labels(row: pd.Series) -> list[str]:
    return [x for x in str(row.get("rescue_mutation_list", "")).split(";") if x and x != "nan"]


def select_pool(target: str, proposals: pd.DataFrame, config: dict, audit: dict) -> pd.DataFrame:
    target_size = int(config["targets"][target]["production_effective_design_pool"])
    available_by_route = proposals.groupby("primary_generation_route").size().to_dict()
    route_targets = adjust_route_targets_for_capacity(
        audit["route_targets"], available_by_route, target_size
    )
    audit["route_targets_soft"] = audit["route_targets"]
    audit["route_targets"] = route_targets
    caps = config.get("production_initial_pool", {}).get("diversity_selection_caps", {})
    cluster_cap = int(
        math.floor(
            target_size
            * float(caps.get("max_single_near_duplicate_cluster_fraction", {}).get(target, 0.03))
        )
    )
    seed_cap = int(
        math.floor(
            target_size
            * float(caps.get("max_single_his_seed_set_fraction", {}).get(target, 0.15))
        )
    )
    rescue_signature_cap = int(
        math.floor(target_size * float(caps.get("max_single_rescue_signature_fraction", 0.05)))
    )
    rescue_mutation_cap = int(
        math.floor(target_size * float(caps.get("max_single_rescue_mutation_fraction", 0.02)))
    )

    selected: list[int] = []
    selected_set: set[int] = set()
    route_counts: Counter[str] = Counter()
    cluster_counts: Counter[str] = Counter()
    seed_counts: Counter[str] = Counter()
    rescue_signature_counts: Counter[str] = Counter()
    rescue_mutation_counts: Counter[str] = Counter()
    his_count_counts: Counter[int] = Counter()
    mutation_count_counts: Counter[int] = Counter()

    sort_cols = ["hit_likelihood_score_v0", "duplicate_source_count", "sequence_hash"]
    sorted_proposals = proposals.sort_values(sort_cols, ascending=[False, False, True])
    records_by_route: dict[str, list[dict]] = {route: [] for route in dry.PRIMARY_ROUTES}
    for row in sorted_proposals.itertuples(index=True):
        rescue_list = [
            x
            for x in str(getattr(row, "rescue_mutation_list", "")).split(";")
            if x and x != "nan"
        ]
        record = {
            "index": row.Index,
            "route": str(row.primary_generation_route),
            "cluster": str(row.near_duplicate_cluster_id),
            "seed": str(row.his_seed_set),
            "rescue_signature": str(row.rescue_signature),
            "rescue_labels": rescue_list,
            "his_count": int(row.His_count),
            "mutation_count": int(row.mutation_count),
        }
        records_by_route.setdefault(record["route"], []).append(record)

    def can_select(record: dict, enforce_route_targets: bool = True) -> bool:
        route = record["route"]
        if enforce_route_targets and route_counts[route] >= route_targets[route]:
            return False
        if cluster_counts[record["cluster"]] >= max(1, cluster_cap):
            return False
        if seed_counts[record["seed"]] >= max(1, seed_cap):
            return False
        if rescue_signature_counts[record["rescue_signature"]] >= max(1, rescue_signature_cap):
            return False
        if any(rescue_mutation_counts[label] >= max(1, rescue_mutation_cap) for label in record["rescue_labels"]):
            return False
        if target == "sdAb":
            if record["his_count"] == 2 and his_count_counts[2] >= int(0.40 * target_size):
                return False
            if record["mutation_count"] == 3 and mutation_count_counts[3] >= int(0.85 * target_size):
                return False
        return True

    def take(record: dict) -> None:
        idx = record["index"]
        selected.append(idx)
        selected_set.add(idx)
        route = record["route"]
        route_counts[route] += 1
        cluster_counts[record["cluster"]] += 1
        seed_counts[record["seed"]] += 1
        rescue_signature_counts[record["rescue_signature"]] += 1
        for label in record["rescue_labels"]:
            rescue_mutation_counts[label] += 1
        his_count_counts[record["his_count"]] += 1
        mutation_count_counts[record["mutation_count"]] += 1

    def fill(
        predicate,
        desired_count: int | None = None,
        current_count=None,
        enforce_route_targets: bool = True,
    ) -> None:
        pointers = {route: 0 for route in dry.PRIMARY_ROUTES}
        while len(selected) < target_size:
            if desired_count is not None and current_count is not None and current_count() >= desired_count:
                return
            made_progress = False
            routes = sorted(
                dry.PRIMARY_ROUTES,
                key=lambda route: route_targets[route] - route_counts[route],
                reverse=True,
            )
            for route in routes:
                if enforce_route_targets and route_counts[route] >= route_targets[route]:
                    continue
                route_pool = records_by_route.get(route, [])
                while pointers[route] < len(route_pool):
                    record = route_pool[pointers[route]]
                    pointers[route] += 1
                    if record["index"] in selected_set:
                        continue
                    if not predicate(record):
                        continue
                    if not can_select(record, enforce_route_targets=enforce_route_targets):
                        continue
                    take(record)
                    made_progress = True
                    if len(selected) >= target_size:
                        return
                    if desired_count is not None and current_count is not None and current_count() >= desired_count:
                        return
                    break
            if not made_progress:
                return

    if target == "sdAb":
        min_2mut = int(math.ceil(0.15 * target_size))
        min_2his = int(math.ceil(0.25 * target_size))
        fill(lambda r: r["mutation_count"] == 2, min_2mut, lambda: mutation_count_counts[2])
        fill(lambda r: r["his_count"] == 2, min_2his, lambda: his_count_counts[2])

    fill(lambda r: True, target_size)
    if len(selected) < target_size:
        fill(lambda r: True, target_size, enforce_route_targets=False)

    if len(selected) != target_size:
        raise RuntimeError(
            f"{target} selected {len(selected)} candidates, expected {target_size}; "
            f"route_counts={dict(route_counts)}; route_targets={route_targets}; "
            f"his_count_counts={dict(his_count_counts)}; mutation_count_counts={dict(mutation_count_counts)}; "
            f"max_cluster_count={max(cluster_counts.values()) if cluster_counts else 0}; "
            f"max_seed_count={max(seed_counts.values()) if seed_counts else 0}; "
            f"max_rescue_signature_count={max(rescue_signature_counts.values()) if rescue_signature_counts else 0}; "
            f"max_rescue_mutation_count={max(rescue_mutation_counts.values()) if rescue_mutation_counts else 0}"
        )

    selected_df = proposals.loc[selected].copy().reset_index(drop=True)
    print(
        json.dumps(
            {
                "stage": "selection_complete",
                "target": target,
                "selected_count": len(selected_df),
                "route_targets": route_targets,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )
    selected_df = dry.assign_near_duplicate_clusters(selected_df)

    selected_by_route = selected_df.groupby("primary_generation_route").size().to_dict()
    proposal_by_route = proposals.groupby("primary_generation_route").size().to_dict()
    for route in dry.PRIMARY_ROUTES:
        selected_count = int(selected_by_route.get(route, 0))
        audit["accepted_by_route"][route] = selected_count
        audit["diversity_cap_by_route"][route] = max(0, int(proposal_by_route.get(route, 0)) - selected_count)
    audit["effective_design_pool_count"] = len(selected_df)
    audit["combined_pool_count"] = len(selected_df)
    return selected_df


def write_route_largest_cluster_summary(candidates: pd.DataFrame) -> None:
    rows = []
    for (target, route), sub in candidates.groupby(["target", "primary_generation_route"], sort=False):
        route_total = len(sub)
        cluster_counts = sub.groupby("near_duplicate_cluster_id").size().sort_values(ascending=False)
        largest = int(cluster_counts.iloc[0]) if not cluster_counts.empty else 0
        rows.append(
            {
                "target": target,
                "primary_generation_route": route,
                "route_total": route_total,
                "largest_near_duplicate_cluster_size": largest,
                "largest_near_duplicate_cluster_fraction": largest / max(1, route_total),
                "largest_near_duplicate_cluster_id": cluster_counts.index[0] if not cluster_counts.empty else "",
            }
        )
    pd.DataFrame(rows).pipe(dry.write_csv, PRODUCTION / "route_largest_near_duplicate_cluster_summary.csv")


def write_manifest(config: dict, cfg_hash: str, input_paths: list[Path], candidates: pd.DataFrame) -> None:
    manifest = {
        "run_id": f"production_initial_pool_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "run_mode": "production_initial_pool_generation_only",
        "tier1_unlocked": False,
        "tier2_unlocked": False,
        "final_10k_selection_unlocked": False,
        "git_commit": dry.git_commit(),
        "config_hash": cfg_hash,
        "input_file_hashes": {str(p.relative_to(ROOT)): dry.file_sha(p) for p in input_paths},
        "target_counts": candidates.groupby("target").size().to_dict(),
        "random_seed": config["random_seed"],
        "python_env": os.environ.get("CONDA_DEFAULT_ENV", "optim-pipe"),
        "tool_versions": {"production_initial_pool_generator": "rule_based_oversample_diversity_select_v1"},
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "completed_at": datetime.now().isoformat(timespec="seconds"),
        "operator": os.environ.get("USER", "unknown"),
        "hostname": socket.gethostname(),
    }
    with (PRODUCTION / "production_pool_manifest.yaml").open("w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)


def markdown_table(frame: pd.DataFrame) -> str:
    cols = list(frame.columns)
    out = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in frame.iterrows():
        vals = []
        for col in cols:
            value = row[col]
            vals.append(f"{value:.4f}" if isinstance(value, float) else str(value))
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)


def write_report(statuses: dict[str, str], audits: dict[str, dict], candidates: pd.DataFrame) -> None:
    fail = [key for key, value in statuses.items() if str(value).startswith("FAIL")]
    warn = [key for key, value in statuses.items() if str(value).startswith("WARN")]
    verdict = "FAIL" if fail else "WARN" if warn else "PASS"
    if fail:
        recommendation = "Do not start Tier 1. Fix production-pool failure modes and regenerate the initial pool."
    elif warn:
        recommendation = "Tier 1 remains locked pending manual review of bounded WARN diversity guardrails."
    else:
        recommendation = "Production initial pool passed validation; Tier 1 can be considered only after explicit manual approval."

    lines = [
        "# Production Initial Pool Validation Report",
        "",
        f"Overall verdict: **{verdict}**",
        "",
        f"Recommendation: {recommendation}",
        "",
        "Tier 1 / Tier 2 / final 10K selection were not run by this command.",
        "",
        "## Counts",
        "",
        "| target | selected effective candidates | raw generated | post hard filter | post exact dedup proposals | controls |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for target, audit in audits.items():
        lines.append(
            f"| {target} | {audit['effective_design_pool_count']} | {audit['raw_generated_count']} | "
            f"{audit['post_generation_hard_filter_count']} | {audit['post_exact_dedup_count']} | "
            f"{audit['control_anchor_pool_count']} |"
        )

    lines += ["", "## Validation Checks", "", "| check | status |", "|---|---|"]
    for key, value in statuses.items():
        lines.append(f"| {key} | {value} |")

    for title, filename in [
        ("Route Distribution", "route_summary.csv"),
        ("His Count Distribution", "his_distribution.csv"),
        ("Mutation Order Distribution", "mutation_order_distribution.csv"),
        ("Top Near-Duplicate Clusters", "top_near_duplicate_clusters.csv"),
    ]:
        path = PRODUCTION / filename
        if path.exists():
            lines += ["", f"## {title}", ""]
            frame = pd.read_csv(path)
            lines.append(markdown_table(frame.head(30)))

    lines += [
        "",
        "## Notes",
        "",
        "- This is production initial-pool generation only.",
        "- It used oversample + hard filter + exact dedup + Hamming/Jaccard diversity selection.",
        "- It did not run ESM, FoldX actual, PyRosetta, AF3, explicit glycan modeling, SimpleFold, dddG_elec, or MD.",
        "- `final_selection_bucket` remains unset.",
    ]
    text = "\n".join(lines) + "\n"
    (PRODUCTION / REPORT_NAME).write_text(text)
    (REPORTS / REPORT_NAME).write_text(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate production initial pool only.")
    parser.add_argument(
        "--target",
        choices=["1E62", "sdAb"],
        action="append",
        help="Optional target subset. Defaults to both targets.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = production_config(dry.read_config())
    selected_targets = args.target or list(config["targets"])
    for target in list(config["targets"]):
        if target not in selected_targets:
            config["targets"].pop(target)

    PRODUCTION.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    dry.TABLES.mkdir(parents=True, exist_ok=True)
    dry.DRY = PRODUCTION
    dry.REPORT_NAME = REPORT_NAME

    cfg_hash = dry.config_hash(config)
    config_text = yaml.safe_dump(config, sort_keys=False, allow_unicode=True)
    (PRODUCTION / "design_config.yaml").write_text(config_text)
    inputs = dry.load_inputs()
    input_paths = [dry.UPSTREAM / f"{name}.csv" for name in inputs]
    dry.build_evidence_ledger(config, inputs)

    all_candidates = []
    all_rescue = []
    all_controls = []
    audits: dict[str, dict] = {}
    for target in config["targets"]:
        proposals, rescue_by_variant, audit = generate_proposals(target, config, inputs, cfg_hash)
        proposals.pipe(dry.write_csv, PRODUCTION / f"production_initial_pool_proposals_{target}.csv")
        selected = select_pool(target, proposals, config, audit)
        selected.pipe(dry.write_csv, PRODUCTION / f"production_initial_pool_candidates_{target}.csv")
        rescue_rows = []
        for variant_id in selected["variant_id"]:
            rescue_rows.extend(rescue_by_variant.get(variant_id, []))
        if rescue_rows:
            all_rescue.append(pd.DataFrame(rescue_rows))

        target_cfg = config["targets"][target]
        chain = target_cfg["chain_id"]
        parent_seq = dry.chain_sequence(inputs["reference_sequence_map"], target, chain)
        parent_motifs = dry.nxs_motifs(parent_seq)
        window = inputs["window_mutation_mask"]
        window = window[window.window_id == target_cfg["window_id"]].copy().sort_values("pos")
        positions = {int(r.pos): str(r.aa) for _, r in window.iterrows()}
        controls = dry.make_control_panel(
            target,
            target_cfg,
            positions,
            chain,
            parent_seq,
            parent_motifs,
            cfg_hash,
            set(selected["sequence_hash"]),
        )
        audit["control_anchor_pool_count"] = len(controls)
        all_controls.extend(controls)
        all_candidates.append(selected)
        audits[target] = audit

    combined = pd.concat(all_candidates, ignore_index=True)
    combined.pipe(dry.write_csv, PRODUCTION / "production_initial_pool_candidates_all.csv")
    combined.pipe(dry.write_csv, dry.TABLES / "production_initial_pool_candidates_all.csv")
    if all_rescue:
        rescue_df = pd.concat(all_rescue, ignore_index=True)
    else:
        rescue_df = pd.DataFrame()
    rescue_df.pipe(dry.write_csv, PRODUCTION / "candidate_rescue_mutations.csv")
    rescue_df.pipe(dry.write_csv, dry.TABLES / "production_candidate_rescue_mutations.csv")
    pd.DataFrame(all_controls).pipe(dry.write_csv, PRODUCTION / "control_anchor_panel.csv")

    statuses = dry.write_audits(combined, audits, all_controls, config)
    write_route_largest_cluster_summary(combined)
    write_manifest(config, cfg_hash, input_paths, combined)
    write_report(statuses, audits, combined)
    print(json.dumps({"statuses": statuses, "counts": combined.groupby("target").size().to_dict()}, indent=2))


if __name__ == "__main__":
    main()
