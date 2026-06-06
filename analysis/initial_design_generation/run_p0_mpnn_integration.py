#!/usr/bin/env python3
"""Prepare P0 ProteinMPNN integration artifacts.

This command does not run ProteinMPNN. It materializes the P0 workflow changes
requested in the ProteinMPNN integration recommendation:

- status-coded scoring table for the current production initial pool;
- constrained-MPNN rescue generation plan;
- relaxed-MPNN / MPNN-only counterfactual plan;
- provenance fields that distinguish rule-generated, MPNN-generated, and
  MPNN-scored-only candidates;
- P0 comparison report scaffold.

Heavy ProteinMPNN jobs can be launched later from these plans once the target
backbone inputs are selected and the proteinmpnn environment is active.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import os
import socket
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT
OUT = dry.OUT
CONFIG = dry.CONFIG
PROD = OUT / "production_initial_pool"
P0 = OUT / "p0_mpnn"
TABLES = OUT / "tables"

TARGET_REVERSE = {"Ab_1E62": "1E62", "Ab_sdAb": "sdAb"}

MPNN_SCORE_COLUMNS = [
    "target",
    "variant_id",
    "sequence_hash",
    "mutation_list",
    "primary_generation_route",
    "his_seed_set",
    "rescue_signature",
    "mpnn_generation_status",
    "mpnn_score_status",
    "mpnn_failure_reason",
    "mpnn_model_version",
    "mpnn_input_backbone_id",
    "mpnn_fixed_positions",
    "mpnn_designed_positions",
    "mpnn_temperature",
    "mpnn_sample_id",
    "mpnn_random_seed",
    "mpnn_total_score",
    "mpnn_total_score_per_residue",
    "mpnn_window_score",
    "mpnn_mutated_position_score",
    "mpnn_his_seed_score",
    "mpnn_rescue_score",
    "mpnn_parent_delta",
    "mpnn_route_percentile",
    "mpnn_his_seed_percentile",
    "mpnn_rescue_type_percentile",
    "mpnn_score_confidence",
]

POOL_COLUMNS = [
    "target",
    "variant_id",
    "generation_record_id",
    "sequence",
    "sequence_hash",
    "mutation_list",
    "mutation_count",
    "His_count",
    "his_seed_set",
    "primary_generation_route",
    "rescue_mutation_list",
    "rescue_signature",
    "rescue_count",
    "mpnn_generation_status",
    "mpnn_model_version",
    "mpnn_input_backbone_id",
    "mpnn_fixed_positions",
    "mpnn_designed_positions",
    "mpnn_temperature",
    "mpnn_sample_id",
    "mpnn_random_seed",
    "hard_filter_status",
    "liability_flags",
    "selection_reason",
]


def read_config() -> dict:
    with CONFIG.open() as fh:
        return yaml.safe_load(fh)


def sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_target_labels(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "target" in out:
        out["target"] = out["target"].replace(TARGET_REVERSE)
    return out


def label_for_pos(reference: pd.DataFrame, target: str, chain: str, pos: int, mutant: str = "H") -> str:
    hit = reference[
        (reference["target"].replace(TARGET_REVERSE) == target)
        & (reference["chain"] == chain)
        & (reference["local_pos"].astype(int) == int(pos))
    ]
    if hit.empty:
        parent = "X"
    else:
        parent = str(hit.iloc[0]["aa"])
    return f"{chain}{parent}{int(pos)}{mutant}"


def parse_range_midpoint(value: int | str) -> int:
    if isinstance(value, int):
        return value
    text = str(value)
    if "-" in text:
        left, right = text.split("-", 1)
        return int((int(left) + int(right)) / 2)
    return int(text)


def as_semicolon(values: Iterable[int | str]) -> str:
    return ";".join(str(v) for v in values)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    # Reuse the project writer so CSV target labels are exported as Ab_1E62 /
    # Ab_sdAb. A bare "1E62" is otherwise parsed as scientific notation by
    # default CSV readers.
    dry.write_csv(df, path)


def current_pool_path() -> Path:
    return PROD / "production_initial_pool_candidates_all.csv"


def load_current_pool() -> pd.DataFrame:
    path = current_pool_path()
    if not path.exists():
        raise FileNotFoundError(f"Missing production pool: {path}")
    return normalize_target_labels(pd.read_csv(path))


def build_mpnn_score_overlay(pool: pd.DataFrame) -> pd.DataFrame:
    score = pool[
        [
            "target",
            "variant_id",
            "sequence_hash",
            "mutation_list",
            "primary_generation_route",
            "his_seed_set",
            "rescue_signature",
        ]
    ].copy()
    score["mpnn_generation_status"] = "not_generated_by_mpnn"
    score["mpnn_score_status"] = "pending_mpnn_run"
    score["mpnn_failure_reason"] = ""
    score["mpnn_model_version"] = ""
    score["mpnn_input_backbone_id"] = ""
    score["mpnn_fixed_positions"] = ""
    score["mpnn_designed_positions"] = ""
    score["mpnn_temperature"] = ""
    score["mpnn_sample_id"] = ""
    score["mpnn_random_seed"] = ""
    for col in [
        "mpnn_total_score",
        "mpnn_total_score_per_residue",
        "mpnn_window_score",
        "mpnn_mutated_position_score",
        "mpnn_his_seed_score",
        "mpnn_rescue_score",
        "mpnn_parent_delta",
        "mpnn_route_percentile",
        "mpnn_his_seed_percentile",
        "mpnn_rescue_type_percentile",
    ]:
        score[col] = pd.NA
    score["mpnn_score_confidence"] = "not_scored"
    return score[MPNN_SCORE_COLUMNS]


def summarize_pending_scores(score: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    summary = score.groupby(group_cols, dropna=False).size().reset_index(name="candidate_count")
    summary["mpnn_score_status"] = "pending_mpnn_run"
    summary["scored_count"] = 0
    summary["pending_count"] = summary["candidate_count"]
    summary["mpnn_score_mean"] = pd.NA
    summary["mpnn_parent_delta_mean"] = pd.NA
    return summary


def build_constrained_plan(config: dict, reference: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    mpnn_cfg = config["proteinmpnn_integration"]["constrained_rescue_generation"]
    rows: list[dict] = []
    for target, target_cfg in config["targets"].items():
        chain = target_cfg["chain_id"]
        his_positions = [int(p) for p in target_cfg["his_seed_positions"]]
        hard_positions = [int(p) for p in target_cfg["hard_protect_positions"]]
        rescue_positions = (
            evidence[
                (evidence["target"] == target)
                & (evidence["rescue_eligible"].astype(str).str.lower().isin(["true", "1"]))
            ]["position"]
            .astype(int)
            .sort_values()
            .tolist()
        )
        seed_sets: list[tuple[int, ...]] = [(p,) for p in his_positions]
        for pair in itertools.combinations(his_positions, 2):
            if target == "sdAb" and set(pair) == {105, 110}:
                continue
            seed_sets.append(pair)
        raw_samples = parse_range_midpoint(mpnn_cfg["raw_samples_per_target_initial"][target])
        retained = str(mpnn_cfg["retained_unique_legal_target"][target])
        per_seed_temp = max(1, raw_samples // max(1, len(seed_sets) * len(mpnn_cfg["temperatures"])))
        for seed_index, seed in enumerate(seed_sets, start=1):
            seed_labels = [
                label_for_pos(reference, target, chain, pos, "H")
                for pos in seed
            ]
            for temp in mpnn_cfg["temperatures"]:
                rows.append(
                    {
                        "target": target,
                        "branch": "constrained_mpnn_rescue",
                        "mpnn_generation_status": "planned_constrained_mpnn",
                        "his_seed_set": ";".join(seed_labels),
                        "fixed_his_positions": as_semicolon(seed),
                        "rescue_design_positions": as_semicolon(rescue_positions),
                        "hard_protect_positions": as_semicolon(hard_positions),
                        "forbidden_pair_policy": "sdAb_V105H_D110H" if target == "sdAb" else "none",
                        "disallowed_global_residues": as_semicolon(mpnn_cfg["disallowed_global"]),
                        "disallow_new_canonical_nxs_t": bool(mpnn_cfg["disallow_new_canonical_nxs_t"]),
                        "mpnn_fixed_positions": as_semicolon(sorted(set(hard_positions) | set(seed))),
                        "mpnn_designed_positions": as_semicolon(rescue_positions),
                        "mpnn_fixed_positions_semantics": "seed_and_hard_protect_only_not_runner_fixed_jsonl",
                        "full_fixed_position_complement_required": True,
                        "allowed_alphabet_status": "pending_position_specific_allowed_alphabet",
                        "backbone_index_map_status": "pending_backbone_and_index_mapping",
                        "mpnn_temperature": temp,
                        "planned_raw_samples": per_seed_temp,
                        "expected_retained_unique_legal": retained,
                        "mpnn_model_version": "",
                        "mpnn_input_backbone_id": "pending_backbone_selection",
                        "mpnn_random_seed": int(config["random_seed"]) + seed_index,
                        "runner_ready": False,
                        "plan_status": "scaffold_ready_not_direct_runner_input",
                    }
                )
    return pd.DataFrame(rows)


def build_relaxed_plan(config: dict, evidence: pd.DataFrame) -> pd.DataFrame:
    relax_cfg = config["proteinmpnn_integration"]["relaxed_counterfactual_generation"]
    rows: list[dict] = []
    for target, target_cfg in config["targets"].items():
        hard_positions = [int(p) for p in target_cfg["hard_protect_positions"]]
        design_positions = (
            evidence[
                (evidence["target"] == target)
                & (~evidence["protected_from_mutation"].astype(str).str.lower().isin(["true", "1"]))
            ]["position"]
            .astype(int)
            .sort_values()
            .tolist()
        )
        for mode in relax_cfg["modes"]:
            rows.append(
                {
                    "target": target,
                    "branch": mode,
                    "mpnn_generation_status": "planned_relaxed_mpnn",
                    "his_policy": "His optional" if mode.endswith("His_optional") else "at least one configured His seed",
                    "hard_protect_positions": as_semicolon(hard_positions),
                    "mpnn_fixed_positions": as_semicolon(hard_positions),
                    "mpnn_designed_positions": as_semicolon(design_positions),
                    "mpnn_fixed_positions_semantics": "hard_protect_only_not_runner_fixed_jsonl",
                    "full_fixed_position_complement_required": True,
                    "allowed_alphabet_status": "pending_position_specific_allowed_alphabet",
                    "backbone_index_map_status": "pending_backbone_and_index_mapping",
                    "seed_enumeration_status": "not_applicable" if mode.endswith("His_optional") else "pending_his_seed_set_enumeration",
                    "effective_candidate_target": f"per_target_total_{relax_cfg['effective_candidates_per_target'][target]}_not_per_branch_quota",
                    "planned_raw_samples_per_target_for_branch": relax_cfg.get("raw_samples_per_target_by_mode", {}).get(mode, ""),
                    "final_10k_allocation_policy": yaml.safe_dump(
                        relax_cfg["final_10k_allocation_policy"],
                        sort_keys=False,
                    ).strip().replace("\n", " | "),
                    "cannot_enter_final_without": as_semicolon(relax_cfg["cannot_enter_final_without"]),
                    "runner_ready": False,
                    "plan_status": "audit_only_until_p0_or_tier2_support",
                }
            )
    return pd.DataFrame(rows)


def empty_pool(target: str, status: str) -> pd.DataFrame:
    df = pd.DataFrame(columns=POOL_COLUMNS)
    df.loc[0, "target"] = target
    df.loc[0, "mpnn_generation_status"] = status
    df.loc[0, "hard_filter_status"] = "not_run"
    df.loc[0, "selection_reason"] = "placeholder_schema_row_not_a_candidate"
    return df


def build_constrained_audit(plan: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target, sub in plan.groupby("target", sort=False):
        rows.append(
            {
                "target": target,
                "branch": "constrained_mpnn_rescue",
                "planned_seed_temperature_jobs": len(sub),
                "unique_his_seed_sets": sub["his_seed_set"].nunique(),
                "planned_raw_samples": int(sub["planned_raw_samples"].sum()),
                "mpnn_generation_status": "pending_mpnn_run",
                "hard_filter_pass_rate": pd.NA,
                "exact_duplicate_rate": pd.NA,
                "near_duplicate_cluster_status": "pending",
                "overlap_with_existing_pool_status": "pending",
            }
        )
    return pd.DataFrame(rows)


def empty_rescue_mutations() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "variant_id",
            "generation_record_id",
            "target",
            "rescue_index",
            "rescue_position_scheme",
            "rescue_position",
            "rescue_window_local_position",
            "parent_aa",
            "mutant_aa",
            "rescue_type",
            "rescue_objective",
            "linked_his_seed",
            "linked_pH_sensitive_module",
            "source_evidence",
            "source_route",
            "risk_flags",
            "allowed_by_evidence_ledger",
            "mpnn_generation_status",
            "mpnn_model_version",
            "mpnn_sample_id",
            "mpnn_random_seed",
        ]
    )


def write_counterfactual_audit(plan: pd.DataFrame) -> None:
    lines = [
        "# Relaxed-MPNN / MPNN-only Counterfactual Audit",
        "",
        "Status: `pending_mpnn_run`",
        "",
        "This is an audit branch only before P0 evidence. It does not reserve a fixed final-library allocation.",
        "",
        "| target | mode | effective target | final allocation policy |",
        "|---|---|---:|---|",
    ]
    for _, row in plan.iterrows():
        policy = str(row["final_10k_allocation_policy"]).replace(" | ", "<br>")
        lines.append(
            f"| {row['target']} | {row['branch']} | {row['effective_candidate_target']} | "
            f"{policy} |"
        )
    lines.append("")
    lines.append("Candidates from this branch cannot enter final 10K without pH mechanism, neutral-retention, liability, and diversity support.")
    (P0 / "relaxed_mpnn_counterfactual_audit.md").write_text("\n".join(lines) + "\n")


def write_runner_preflight(constrained: pd.DataFrame, relaxed: pd.DataFrame) -> None:
    requirements = {
        "status": "blocked_pending_backbone_index_map",
        "runner_ready": False,
        "safe_to_submit_to_proteinmpnn": False,
        "reason": (
            "P0 plan rows are scaffold rows. They are not direct ProteinMPNN "
            "fixed_positions_jsonl or bias/omit inputs."
        ),
        "required_before_true_mpnn_run": {
            "backbone_inputs": [
                "target",
                "mpnn_input_backbone_id",
                "pdb_or_cif_path",
                "design_chain_id",
                "antigen_chain_ids",
                "chain_length_map",
                "local_position_to_mpnn_index_map",
            ],
            "runner_jsonl_inputs": [
                "chain_id_jsonl",
                "full_fixed_positions_jsonl",
                "per_position_allowed_or_omit_AA_jsonl",
                "seed_specific_bias_or_fixed_His_constraints",
            ],
            "constrained_generation_rules": [
                "fix all antigen positions",
                "fix all non-target antibody chains",
                "fix all target-chain positions outside selected 40 aa window",
                "fix His seed positions to H",
                "open only evidence-ledger rescue positions",
                "apply per-position allowed amino acid rules",
            ],
            "relaxed_counterfactual_rules": [
                "His_optional mode must still keep hard-protect and liability constraints",
                "His_seeded mode must enumerate or otherwise constrain His seed sets",
                "relaxed candidates remain audit-only until manual unlock",
            ],
            "post_filters": [
                "no hard-protect mutation",
                "no sdAb V105H-D110H",
                "no new Cys",
                "no new canonical N-X-S/T",
                "max His count",
                "max non-His rescue count",
                "exact duplicate and near-duplicate audit",
                "overlap with existing production pool",
            ],
        },
        "current_plan_summary": {
            "constrained_plan_rows": int(len(constrained)),
            "relaxed_plan_rows": int(len(relaxed)),
            "constrained_runner_ready_values": sorted(
                str(x) for x in constrained["runner_ready"].dropna().unique()
            ),
            "relaxed_runner_ready_values": sorted(
                str(x) for x in relaxed["runner_ready"].dropna().unique()
            ),
        },
    }
    with (P0 / "mpnn_runner_input_requirements.yaml").open("w") as fh:
        yaml.safe_dump(requirements, fh, sort_keys=False, allow_unicode=True)

    lines = [
        "# ProteinMPNN Runner Input Preflight",
        "",
        "Status: `blocked_pending_backbone_index_map`",
        "",
        "The P0 plan tables are scaffold/audit artifacts. They are not direct ProteinMPNN runner inputs.",
        "",
        "## Why Runner Is Not Ready",
        "",
        "- Full chain-aware fixed-position complements have not been generated.",
        "- Per-position allowed/omitted amino acid rules have not been generated.",
        "- Backbone IDs, chain IDs, and local-position to ProteinMPNN-index maps have not been selected.",
        "- Relaxed His-seeded mode still needs explicit seed-set enumeration or an equivalent constraint.",
        "",
        "## Current Plan Rows",
        "",
        f"- Constrained plan rows: {len(constrained)}",
        f"- Relaxed plan rows: {len(relaxed)}",
        "- `runner_ready`: false for all current plan rows.",
        "",
        "## Required Before True MPNN Run",
        "",
        "1. Select backbone inputs and chain mappings for each target.",
        "2. Generate full `fixed_positions_jsonl`, not just seed/hard-protect summaries.",
        "3. Generate per-position allowed/omit/bias JSONL files.",
        "4. Encode fixed His seed constraints for constrained generation.",
        "5. Add post-filters for hard-protect, forbidden pair, Cys, N-X-S/T, His count, rescue count, duplicates, and overlap.",
    ]
    (P0 / "mpnn_runner_preflight_report.md").write_text("\n".join(lines) + "\n")


def write_comparison_report(pool: pd.DataFrame, score: pd.DataFrame, constrained: pd.DataFrame, relaxed: pd.DataFrame) -> None:
    counts = pool.groupby("target").size().to_dict()
    route = pool.groupby(["target", "primary_generation_route"]).size().reset_index(name="count")
    route["fraction"] = route.groupby("target")["count"].transform(lambda s: s / s.sum())
    lines = [
        "# P0 ProteinMPNN Comparison Report",
        "",
        "Status: `scaffold_ready_pending_mpnn_run`",
        "",
        "This report verifies the workflow change and required tables. It does not claim ProteinMPNN scores or MPNN-generated candidates have already been produced.",
        "",
        "## Current Production Pool",
        "",
        "| target | candidates |",
        "|---|---:|",
    ]
    for target, count in counts.items():
        lines.append(f"| {target} | {count} |")
    lines += [
        "",
        "## MPNN Scoring Status",
        "",
        f"- Score rows prepared: {len(score)}",
        "- `mpnn_score_status`: `pending_mpnn_run`",
        "- Current production pool provenance: `not_generated_by_mpnn`",
        "",
        "## Constrained-MPNN Rescue Plan",
        "",
        f"- Planned seed/temperature jobs: {len(constrained)}",
        "- This branch is supplemental and may compete in His-rescue buckets after P0 evidence.",
        "- Scaffold plan rows are not direct ProteinMPNN runner inputs; runner-specific JSONL must be generated separately with reviewed backbones before compute.",
        "",
        "## Relaxed-MPNN / MPNN-only Counterfactual Plan",
        "",
        f"- Planned mode rows: {len(relaxed)}",
        "- Pre-Tier1 status: audit only.",
        "- No fixed 3-5% final allocation is reserved.",
        "- His-seeded relaxed mode requires explicit seed-set enumeration in runner inputs before execution.",
        "",
        "## Route Snapshot",
        "",
        "| target | route | count | fraction |",
        "|---|---|---:|---:|",
    ]
    for _, row in route.iterrows():
        lines.append(
            f"| {row['target']} | {row['primary_generation_route']} | {int(row['count'])} | {row['fraction']:.4f} |"
        )
    lines += [
        "",
        "## Next Decision",
        "",
        "Run ProteinMPNN scoring and constrained/relaxed generation, then replace pending status-coded tables with computed outputs before approving Tier 2-core.",
    ]
    (P0 / "p0_mpnn_comparison_report.md").write_text("\n".join(lines) + "\n")


def write_manifest(config: dict, output_files: list[Path]) -> None:
    manifest = {
        "run_id": f"p0_mpnn_scaffold_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "run_mode": "p0_mpnn_integration_scaffold",
        "proteinmpnn_ran": False,
        "tier2_unlocked": False,
        "final_10k_selection_unlocked": False,
        "config_hash": dry.config_hash(config),
        "input_files": {
            str(current_pool_path().relative_to(ROOT)): sha_file(current_pool_path()),
            str(CONFIG.relative_to(ROOT)): sha_file(CONFIG),
        },
        "output_files": {str(p.relative_to(ROOT)): sha_file(p) for p in output_files if p.exists()},
        "operator": os.environ.get("USER", "unknown"),
        "hostname": socket.gethostname(),
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    with (P0 / "p0_mpnn_manifest.yaml").open("w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, allow_unicode=True)


def validate_outputs() -> list[str]:
    required = [
        "mpnn_scores_current_pool.csv",
        "mpnn_score_distribution_by_route.csv",
        "mpnn_score_distribution_by_his_seed.csv",
        "mpnn_score_distribution_by_rescue_type.csv",
        "mpnn_parent_delta_by_route.csv",
        "mpnn_low_score_high_pH_support_candidates.csv",
        "mpnn_high_score_no_pH_mechanism_candidates.csv",
        "constrained_mpnn_rescue_plan.csv",
        "constrained_mpnn_rescue_pool_1E62.csv",
        "constrained_mpnn_rescue_pool_sdAb.csv",
        "constrained_mpnn_rescue_audit.csv",
        "candidate_rescue_mutations_mpnn.csv",
        "relaxed_mpnn_counterfactual_plan.csv",
        "relaxed_mpnn_counterfactual_pool_1E62.csv",
        "relaxed_mpnn_counterfactual_pool_sdAb.csv",
        "relaxed_mpnn_counterfactual_audit.md",
        "mpnn_runner_input_requirements.yaml",
        "mpnn_runner_preflight_report.md",
        "p0_mpnn_comparison_report.md",
        "p0_mpnn_manifest.yaml",
    ]
    missing = [name for name in required if not (P0 / name).exists()]
    if missing:
        return [f"missing output: {name}" for name in missing]
    score_cols = set(pd.read_csv(P0 / "mpnn_scores_current_pool.csv", nrows=1).columns)
    missing_cols = [col for col in MPNN_SCORE_COLUMNS if col not in score_cols]
    if missing_cols:
        return [f"mpnn_scores_current_pool.csv missing columns: {missing_cols}"]
    return []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare P0 ProteinMPNN integration scaffold.")
    parser.add_argument("--output-dir", default=None, help="Override P0 output directory.")
    return parser.parse_args()


def main() -> None:
    global P0
    args = parse_args()
    config = read_config()
    if args.output_dir:
        P0 = ROOT / args.output_dir
    else:
        P0 = ROOT / config["proteinmpnn_integration"]["output_dir"]
    P0.mkdir(parents=True, exist_ok=True)

    pool = load_current_pool()
    score = build_mpnn_score_overlay(pool)
    write_csv(score, P0 / "mpnn_scores_current_pool.csv")

    write_csv(summarize_pending_scores(score, ["target", "primary_generation_route"]), P0 / "mpnn_score_distribution_by_route.csv")
    write_csv(summarize_pending_scores(score, ["target", "his_seed_set"]), P0 / "mpnn_score_distribution_by_his_seed.csv")
    write_csv(summarize_pending_scores(score, ["target", "primary_generation_route", "rescue_signature"]), P0 / "mpnn_score_distribution_by_rescue_type.csv")
    write_csv(summarize_pending_scores(score, ["target", "primary_generation_route"]), P0 / "mpnn_parent_delta_by_route.csv")

    empty_flag_cols = list(score.columns) + ["reason"]
    write_csv(pd.DataFrame(columns=empty_flag_cols), P0 / "mpnn_low_score_high_pH_support_candidates.csv")
    write_csv(pd.DataFrame(columns=empty_flag_cols), P0 / "mpnn_high_score_no_pH_mechanism_candidates.csv")

    evidence = normalize_target_labels(pd.read_csv(TABLES / "evidence_ledger.csv"))
    reference = normalize_target_labels(pd.read_csv(dry.UPSTREAM / "reference_sequence_map.csv"))
    constrained = build_constrained_plan(config, reference, evidence)
    relaxed = build_relaxed_plan(config, evidence)
    write_csv(constrained, P0 / "constrained_mpnn_rescue_plan.csv")
    write_csv(relaxed, P0 / "relaxed_mpnn_counterfactual_plan.csv")
    for target in config["targets"]:
        write_csv(empty_pool(target, "pending_constrained_mpnn_run"), P0 / f"constrained_mpnn_rescue_pool_{target}.csv")
        write_csv(empty_pool(target, "pending_relaxed_mpnn_run"), P0 / f"relaxed_mpnn_counterfactual_pool_{target}.csv")
    write_csv(build_constrained_audit(constrained), P0 / "constrained_mpnn_rescue_audit.csv")
    write_csv(empty_rescue_mutations(), P0 / "candidate_rescue_mutations_mpnn.csv")
    write_counterfactual_audit(relaxed)
    write_runner_preflight(constrained, relaxed)
    write_comparison_report(pool, score, constrained, relaxed)

    output_files = [
        p
        for p in P0.iterdir()
        if p.is_file() and p.name != "p0_mpnn_manifest.yaml"
    ]
    write_manifest(config, output_files)
    errors = validate_outputs()
    if errors:
        raise RuntimeError("; ".join(errors))
    print(
        {
            "status": "p0_mpnn_scaffold_ready",
            "output_dir": str(P0.relative_to(ROOT)),
            "score_rows": len(score),
            "constrained_plan_rows": len(constrained),
            "relaxed_plan_rows": len(relaxed),
        }
    )


if __name__ == "__main__":
    main()
