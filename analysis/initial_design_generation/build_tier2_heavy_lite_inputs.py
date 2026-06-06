#!/usr/bin/env python3
"""Prepare manual-waived Tier2-heavy-lite pilot input package."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(".")
PLANNING_DIR = ROOT / "results/initial_design_generation/tier2_heavy_planning"
OUT_DIR = ROOT / "results/initial_design_generation/tier2_heavy_lite"
CURRENT_STAGE_REPORT = ROOT / ".tasks/active/initial-design-generation/current_stage_report.md"


MAIN_1E62 = PLANNING_DIR / "tier2_heavy_planning_pool_1E62.csv"
MAIN_SDAB = PLANNING_DIR / "tier2_heavy_planning_pool_sdAb.csv"
CONTROL = PLANNING_DIR / "tier2_heavy_control_anchor_panel.csv"
AUDIT_REPORT = PLANNING_DIR / "tier2_heavy_planning_pool_audit.md"
WAIVER_NOTE = ROOT / ".tasks/active/initial-design-generation/tier2_heavy_lite_manual_waiver_note.md"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


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


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    main_1e62 = pd.read_csv(MAIN_1E62)
    main_sdab = pd.read_csv(MAIN_SDAB)
    control = pd.read_csv(CONTROL)

    main = pd.concat([main_1e62, main_sdab], ignore_index=True)
    control = control.copy()
    main["tier2_heavy_lite_role"] = "main_candidate"
    control["tier2_heavy_lite_role"] = "control_anchor"
    combined = pd.concat([main, control], ignore_index=True)
    combined["tier2_heavy_lite_input_rank"] = range(1, len(combined) + 1)

    input_path = OUT_DIR / "tier2_heavy_lite_input_list.csv"
    combined.to_csv(input_path, index=False)

    by_role = combined.groupby(["target", "tier2_heavy_lite_role"]).size().reset_index(name="count")
    route = (
        combined.assign(heavy_route=combined["heavy_route"].fillna("").str.split(";"))
        .explode("heavy_route")
        .groupby(["target", "tier2_heavy_lite_role", "heavy_route"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "tier2_heavy_lite_role", "count"], ascending=[True, True, False])
    )
    seed = (
        combined.groupby(["target", "tier2_heavy_lite_role", "his_seed_set"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "tier2_heavy_lite_role", "count"], ascending=[True, True, False])
    )
    by_role.to_csv(OUT_DIR / "tier2_heavy_lite_input_by_role.csv", index=False)
    route.to_csv(OUT_DIR / "tier2_heavy_lite_input_by_route.csv", index=False)
    seed.to_csv(OUT_DIR / "tier2_heavy_lite_input_by_seed.csv", index=False)

    manifest = [
        "stage: Tier2-heavy-lite",
        "status: input_prepared_manual_waiver",
        "manual_waiver:",
        f"  path: {WAIVER_NOTE.as_posix()}",
        f"  sha256: {sha256_file(WAIVER_NOTE)}",
        "planning_audit:",
        f"  path: {AUDIT_REPORT.as_posix()}",
        f"  sha256: {sha256_file(AUDIT_REPORT)}",
        "inputs:",
        f"  main_1E62: {MAIN_1E62.as_posix()}",
        f"  main_1E62_sha256: {sha256_file(MAIN_1E62)}",
        f"  main_sdAb: {MAIN_SDAB.as_posix()}",
        f"  main_sdAb_sha256: {sha256_file(MAIN_SDAB)}",
        f"  controls: {CONTROL.as_posix()}",
        f"  controls_sha256: {sha256_file(CONTROL)}",
        f"  combined_input: {input_path.as_posix()}",
        f"  combined_input_sha256: {sha256_file(input_path)}",
        "counts:",
        f"  main_1E62: {len(main_1e62)}",
        f"  main_sdAb: {len(main_sdab)}",
        f"  controls: {len(control)}",
        f"  total: {len(combined)}",
        "locked:",
        "  broad_tier2_heavy: true",
        "  broad_af3_simplefold_review: true",
        "  md: true",
        "  final_10k: true",
        "",
    ]
    (OUT_DIR / "tier2_heavy_lite_input_manifest.yaml").write_text("\n".join(manifest), encoding="utf-8")

    config = [
        "stage: Tier2-heavy-lite",
        "status: draft_until_explicit_run_command",
        f"input_list: {input_path.as_posix()}",
        f"input_list_sha256: {sha256_file(input_path)}",
        "manual_waiver_unlock: true",
        "scope:",
        "  main_candidates: 127",
        "  controls_anchors: 40",
        "  total_rows: 167",
        "heavy_routes:",
        "  electrostatics_pKa_review:",
        "    candidates: all main candidates",
        "  SimpleFold_antibody_sampling:",
        "    candidates: all sdAb main candidates; all 1E62 main candidates if resources allow",
        "  AF3_complex_check:",
        "    candidates: assigned representatives only",
        "locked:",
        "  broad_tier2_heavy: true",
        "  md: true",
        "  final_10k: true",
        "",
    ]
    (OUT_DIR / "tier2_heavy_lite_compute_config_draft.yaml").write_text("\n".join(config), encoding="utf-8")

    command_plan = [
        "# Tier2-heavy-lite Command Plan",
        "",
        "This is a command plan only. No heavy-lite compute has been started by this preparation step.",
        "",
        "## Input",
        "",
        f"`{input_path.as_posix()}`",
        "",
        "## Scope",
        "",
        "```text",
        "1E62 main = 52",
        "sdAb main = 75",
        "controls / anchors = 40",
        "total = 167",
        "```",
        "",
        "## Route execution intent",
        "",
        "```text",
        "electrostatics_pKa_review: all main candidates",
        "SimpleFold_antibody_sampling: all sdAb main candidates; all 1E62 main candidates if resources allow",
        "AF3_complex_check: assigned representatives only",
        "```",
        "",
        "## Required reports after compute",
        "",
        "```text",
        "tier2_heavy_lite_results.csv",
        "tier2_heavy_lite_by_target.md",
        "tier2_heavy_lite_by_seed.md",
        "tier2_heavy_lite_by_cluster.md",
        "tier2_heavy_lite_by_mutation_count.md",
        "tier2_heavy_lite_control_anchor_report.md",
        "tier2_heavy_lite_next_gate.md",
        "```",
        "",
        "Broad Tier2-heavy, MD, and final 10K remain locked.",
        "",
    ]
    (OUT_DIR / "tier2_heavy_lite_command_plan.md").write_text("\n".join(command_plan), encoding="utf-8")

    report = [
        "# Manual-waived Tier2-heavy-lite Pilot Preparation",
        "",
        "## Executive Summary",
        "",
        "The secondary PATCH review accepts the current small clean planning pool by manual waiver for a targeted Tier2-heavy-lite pilot. This is not an automatic PASS and does not unlock broad Tier2-heavy, broad AF3 / SimpleFold review, MD, or final 10K.",
        "",
        "## Pilot Input",
        "",
        md_table(by_role),
        "",
        "## Heavy Route Assignment",
        "",
        md_table(route),
        "",
        "## Seed Summary",
        "",
        md_table(seed.groupby(["target", "tier2_heavy_lite_role"]).head(12).reset_index(drop=True)),
        "",
        "## Prepared Files",
        "",
        f"- Input list: `{input_path.as_posix()}`",
        f"- Manifest: `{(OUT_DIR / 'tier2_heavy_lite_input_manifest.yaml').as_posix()}`",
        f"- Draft config: `{(OUT_DIR / 'tier2_heavy_lite_compute_config_draft.yaml').as_posix()}`",
        f"- Command plan: `{(OUT_DIR / 'tier2_heavy_lite_command_plan.md').as_posix()}`",
        "",
        "## Boundary",
        "",
        "No compute has been started by this preparation step. Broad Tier2-heavy, broad AF3 / SimpleFold review, MD, and final 10K remain locked.",
        "",
    ]
    current_text = "\n".join(report)
    (OUT_DIR / "tier2_heavy_lite_preparation_report.md").write_text(current_text, encoding="utf-8")
    CURRENT_STAGE_REPORT.write_text(current_text, encoding="utf-8")

    print(f"input_rows={len(combined)}")
    print(by_role.to_string(index=False))


if __name__ == "__main__":
    main()
