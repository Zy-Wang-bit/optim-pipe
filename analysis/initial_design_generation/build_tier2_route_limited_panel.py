#!/usr/bin/env python3
"""Build Tier2-heavy-lite Stage B route-limited review panel.

This step only assigns next-review routes from the local heavy-lite results.
It does not run AF3, SimpleFold, glycan review, MD, or final selection.
"""

from __future__ import annotations

import hashlib
import math
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "results/initial_design_generation/tier2_heavy_lite"
OUT_DIR = ROOT / "results/initial_design_generation/tier2_route_limited_stage_b"
TASK_DIR = ROOT / ".tasks/active/initial-design-generation"

SOURCE_RESULTS = SOURCE_DIR / "tier2_heavy_lite_results.csv"
SOURCE_GATE = SOURCE_DIR / "tier2_heavy_lite_next_gate.md"
CURRENT_STAGE_REPORT = TASK_DIR / "current_stage_report.md"

SUPPORTED_CLASSES = {"HL_promotable", "HL_boundary_supported"}
TARGET_ROUTE_COUNTS = {
    "Ab_1E62": {
        "af3": 20,
        "simplefold": 20,
        "af3_range": "15-25",
        "simplefold_range": "15-25",
        "af3_seed_cap": 5,
        "af3_cluster_cap": 2,
        "simplefold_seed_cap": 5,
        "simplefold_cluster_cap": 2,
    },
    "Ab_sdAb": {
        "af3": 25,
        "simplefold": 40,
        "af3_range": "20-30",
        "simplefold_range": "30-45",
        "af3_seed_cap": 7,
        "af3_cluster_cap": 2,
        "simplefold_seed_cap": 10,
        "simplefold_cluster_cap": 3,
    },
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    view = df.copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
    lines = [
        "| " + " | ".join(map(str, view.columns)) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in view.columns) + " |")
    return "\n".join(lines)


def numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def route_rank(df: pd.DataFrame, route: str) -> pd.DataFrame:
    """Return a ranked view for a specific route assignment."""
    out = df.copy()
    class_priority = {"HL_promotable": 0, "HL_boundary_supported": 1}
    stage_priority = {
        "T2A_high_confidence": 0,
        "T2A_good": 1,
        "T2A_supported_boundary_confirmed": 2,
    }
    out["_class_priority"] = out["heavy_lite_final_class"].map(class_priority).fillna(9)
    out["_stage_priority"] = out["stage2a_final_class"].map(stage_priority).fillna(9)
    out["_validity_sort"] = -out["pyrosetta_local_structure_validity_t2_score"].fillna(-1)
    out["_delta_sort"] = out["pyrosetta_rosetta_delta_score"].fillna(math.inf)
    out["_pka_sort"] = -out["his_pka_support_t2_score"].fillna(-1)
    out["_clash_sort"] = out["pyrosetta_local_clash_count"].fillna(math.inf)

    # SimpleFold is the route where loop/window sensitivity matters most, so
    # let high mutation-order and boundary-confirmed representatives surface.
    if route == "simplefold":
        out["_mutation_sort"] = -out["mutation_count"].fillna(0)
        out["_boundary_sort"] = out["heavy_lite_final_class"].ne("HL_boundary_supported").astype(int)
        keys = [
            "_class_priority",
            "_boundary_sort",
            "_mutation_sort",
            "_stage_priority",
            "_validity_sort",
            "_delta_sort",
            "_pka_sort",
            "_clash_sort",
            "variant_id",
        ]
    else:
        keys = [
            "_class_priority",
            "_stage_priority",
            "_validity_sort",
            "_delta_sort",
            "_pka_sort",
            "_clash_sort",
            "variant_id",
        ]
    return out.sort_values(keys, ascending=True)


def diverse_select(
    sub: pd.DataFrame,
    target_count: int,
    route: str,
    seed_cap: int,
    cluster_cap: int,
) -> set[str]:
    """Select a capped representative subset with seed and cluster coverage."""
    ranked = route_rank(sub, route)
    selected: list[str] = []
    selected_set: set[str] = set()
    seed_counts: Counter[str] = Counter()
    cluster_counts: Counter[str] = Counter()

    def can_take(row: pd.Series, seed_limit: int, cluster_limit: int) -> bool:
        seed = str(row.get("his_seed_set", ""))
        cluster = str(row.get("near_duplicate_cluster_id", ""))
        return seed_counts[seed] < seed_limit and cluster_counts[cluster] < cluster_limit

    def take(row: pd.Series) -> None:
        variant_id = str(row["variant_id"])
        if variant_id in selected_set:
            return
        selected.append(variant_id)
        selected_set.add(variant_id)
        seed_counts[str(row.get("his_seed_set", ""))] += 1
        cluster_counts[str(row.get("near_duplicate_cluster_id", ""))] += 1

    # First pass: guarantee at least one representative per seed where possible.
    for _, seed_rows in ranked.groupby("his_seed_set", dropna=False, sort=False):
        for _, row in seed_rows.iterrows():
            if can_take(row, seed_cap, cluster_cap):
                take(row)
                break
        if len(selected) >= target_count:
            return selected_set

    # Main capped fill.
    for _, row in ranked.iterrows():
        if len(selected) >= target_count:
            break
        if can_take(row, seed_cap, cluster_cap):
            take(row)

    # Conservative relaxation only if target count cannot be reached. This
    # keeps route panels in the requested size range without reintroducing a
    # dominant seed or cluster.
    for relaxed_cluster_cap in (cluster_cap + 1, cluster_cap + 2):
        if len(selected) >= target_count:
            break
        for _, row in ranked.iterrows():
            if len(selected) >= target_count:
                break
            if can_take(row, seed_cap, relaxed_cluster_cap):
                take(row)

    return selected_set


def split_routes(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["af3_selected"] = False
    out["simplefold_selected"] = False
    out["electrostatics_selected"] = False
    out["glycan_selected"] = False
    out["control_or_anchor"] = out["tier2_heavy_lite_role"].eq("control_anchor")
    out["selection_reason"] = ""

    supported_mask = (
        out["tier2_heavy_lite_role"].eq("main_candidate")
        & out["heavy_lite_final_class"].isin(SUPPORTED_CLASSES)
    )
    out.loc[supported_mask, "electrostatics_selected"] = True
    out.loc[supported_mask, "selection_reason"] = "local_heavy_lite_supported;electrostatics_all_supported"

    for target, cfg in TARGET_ROUTE_COUNTS.items():
        sub = out[supported_mask & out["target"].eq(target)].copy()
        af3 = diverse_select(
            sub,
            target_count=cfg["af3"],
            route="af3",
            seed_cap=cfg["af3_seed_cap"],
            cluster_cap=cfg["af3_cluster_cap"],
        )
        simplefold = diverse_select(
            sub,
            target_count=cfg["simplefold"],
            route="simplefold",
            seed_cap=cfg["simplefold_seed_cap"],
            cluster_cap=cfg["simplefold_cluster_cap"],
        )
        out.loc[out["variant_id"].isin(af3), "af3_selected"] = True
        out.loc[out["variant_id"].isin(simplefold), "simplefold_selected"] = True

    out.loc[out["af3_selected"], "selection_reason"] = (
        out.loc[out["af3_selected"], "selection_reason"] + ";AF3_seed_cluster_representative"
    ).str.strip(";")
    out.loc[out["simplefold_selected"], "selection_reason"] = (
        out.loc[out["simplefold_selected"], "selection_reason"] + ";SimpleFold_loop_or_boundary_representative"
    ).str.strip(";")
    out.loc[out["control_or_anchor"], "selection_reason"] = "control_or_anchor_retained_for_calibration"

    # No explicit glycan-risk flag is present in the heavy-lite results table.
    # Keep this route empty until a glycan-near flag is available.
    out["glycan_selection_note"] = "not_selected_no_glycan_risk_flag_available"
    out.loc[out["glycan_selected"], "glycan_selection_note"] = "selected_by_glycan_risk_flag"

    def join_routes(row: pd.Series) -> str:
        routes: list[str] = []
        if row["af3_selected"]:
            routes.append("AF3_complex_check")
        if row["simplefold_selected"]:
            routes.append("SimpleFold_antibody_sampling")
        if row["electrostatics_selected"]:
            routes.append("electrostatics_pKa_review")
        if row["glycan_selected"]:
            routes.append("glycan_risk_review")
        if row["control_or_anchor"]:
            routes.append("control_or_anchor_review")
        return ";".join(routes) if routes else "not_selected"

    out["route_assignment"] = out.apply(join_routes, axis=1)
    return out


def summarize_route_by_target(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target, sub in panel.groupby("target"):
        main_supported = sub[
            sub["tier2_heavy_lite_role"].eq("main_candidate")
            & sub["heavy_lite_final_class"].isin(SUPPORTED_CLASSES)
        ]
        rows.append(
            {
                "target": target,
                "supported_main_count": len(main_supported),
                "af3_selected": int(main_supported["af3_selected"].sum()),
                "simplefold_selected": int(main_supported["simplefold_selected"].sum()),
                "electrostatics_selected": int(main_supported["electrostatics_selected"].sum()),
                "glycan_selected": int(main_supported["glycan_selected"].sum()),
                "control_anchor_count": int(sub["control_or_anchor"].sum()),
                "top_af3_seed_fraction": top_fraction(main_supported[main_supported["af3_selected"]], "his_seed_set"),
                "top_simplefold_seed_fraction": top_fraction(main_supported[main_supported["simplefold_selected"]], "his_seed_set"),
                "top_af3_cluster_fraction": top_fraction(main_supported[main_supported["af3_selected"]], "near_duplicate_cluster_id"),
                "top_simplefold_cluster_fraction": top_fraction(main_supported[main_supported["simplefold_selected"]], "near_duplicate_cluster_id"),
            }
        )
    return pd.DataFrame(rows)


def top_fraction(df: pd.DataFrame, col: str) -> float:
    if df.empty:
        return 0.0
    return float(df[col].value_counts(normalize=True, dropna=False).max())


def route_explode(panel: pd.DataFrame) -> pd.DataFrame:
    exploded = panel.copy()
    exploded["route"] = exploded["route_assignment"].str.split(";")
    exploded = exploded.explode("route")
    return exploded[exploded["route"].ne("not_selected")]


def write_outputs(panel: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    supported_or_control = (
        panel["control_or_anchor"]
        | (
            panel["tier2_heavy_lite_role"].eq("main_candidate")
            & panel["heavy_lite_final_class"].isin(SUPPORTED_CLASSES)
        )
    )
    panel_out = panel[supported_or_control].copy()

    panel_cols = [
        "variant_id",
        "target",
        "mutation_list",
        "heavy_lite_final_class",
        "heavy_lite_action",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "mutation_count",
        "route_assignment",
        "af3_selected",
        "simplefold_selected",
        "electrostatics_selected",
        "glycan_selected",
        "control_or_anchor",
        "selection_reason",
        "glycan_selection_note",
        "stage2a_final_class",
        "boundary_support_level",
        "pyrosetta_rosetta_delta_score",
        "pyrosetta_local_structure_validity_t2_score",
        "pyrosetta_local_clash_count",
        "pka_pka_combined_mean",
        "his_pka_support_t2_score",
        "pyrosetta_pdb_path",
    ]
    panel_path = OUT_DIR / "tier2_heavy_route_limited_panel.csv"
    panel_out[panel_cols].to_csv(panel_path, index=False)
    shutil.copyfile(panel_path, OUT_DIR / "tier2_route_limited_panel.csv")

    supported = panel_out[
        panel_out["tier2_heavy_lite_role"].eq("main_candidate")
        & panel_out["heavy_lite_final_class"].isin(SUPPORTED_CLASSES)
    ].copy()
    supported_evidence = supported.rename(
        columns={
            "pyrosetta_local_structure_validity_t2_score": "local_validity",
            "pyrosetta_rosetta_delta_score": "rosetta_delta",
            "his_pka_support_t2_score": "his_pka_support",
        }
    )
    supported_evidence["selection_status"] = "local_heavy_lite_supported"
    evidence_cols = [
        "variant_id",
        "target",
        "mutation_list",
        "his_seed_set",
        "near_duplicate_cluster_id",
        "mutation_count",
        "heavy_lite_final_class",
        "local_validity",
        "rosetta_delta",
        "his_pka_support",
        "control_or_anchor",
        "selection_status",
    ]
    supported_evidence[evidence_cols].to_csv(
        OUT_DIR / "tier2_heavy_lite_supported_evidence_table.csv",
        index=False,
    )

    controls = panel_out[panel_out["control_or_anchor"]].copy()
    controls[panel_cols].to_csv(OUT_DIR / "tier2_heavy_route_controls.csv", index=False)

    panel_out[panel_out["af3_selected"]][panel_cols].to_csv(
        OUT_DIR / "af3_complex_candidate_list.csv",
        index=False,
    )
    panel_out[panel_out["simplefold_selected"]][panel_cols].to_csv(
        OUT_DIR / "simplefold_candidate_list.csv",
        index=False,
    )
    panel_out[panel_out["electrostatics_selected"]][panel_cols].to_csv(
        OUT_DIR / "electrostatics_pka_review_list.csv",
        index=False,
    )
    glycan_rows = panel_out[panel_out["glycan_selected"]][panel_cols].copy()
    glycan_rows.to_csv(OUT_DIR / "glycan_review_candidate_list.csv", index=False)

    by_target = summarize_route_by_target(panel_out)
    by_target.to_csv(OUT_DIR / "tier2_heavy_route_by_target.csv", index=False)

    route_rows = route_explode(panel_out)
    by_seed = (
        route_rows.groupby(["target", "route", "his_seed_set", "heavy_lite_final_class"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "route", "count"], ascending=[True, True, False])
    )
    by_cluster = (
        route_rows.groupby(["target", "route", "near_duplicate_cluster_id"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "route", "count"], ascending=[True, True, False])
    )
    by_seed.to_csv(OUT_DIR / "tier2_heavy_route_by_seed.csv", index=False)
    by_cluster.to_csv(OUT_DIR / "tier2_heavy_route_by_cluster.csv", index=False)

    write_manifest(panel_path, by_target, panel_out)
    write_audit(panel_path, panel_out, by_target, by_seed, by_cluster)


def write_manifest(panel_path: Path, by_target: pd.DataFrame, panel: pd.DataFrame) -> None:
    lines = [
        "stage: Tier2-heavy-lite Stage B",
        "status: route_limited_panel_ready_manual_unlock_required",
        "source:",
        f"  local_heavy_lite_results: {SOURCE_RESULTS.as_posix()}",
        f"  local_heavy_lite_results_sha256: {sha256_file(SOURCE_RESULTS)}",
        f"  local_gate_report: {SOURCE_GATE.as_posix()}",
        f"  local_gate_report_sha256: {sha256_file(SOURCE_GATE)}",
        "outputs:",
        f"  route_limited_panel: {panel_path.as_posix()}",
        f"  route_limited_panel_sha256: {sha256_file(panel_path)}",
        f"  supported_evidence_table: {(OUT_DIR / 'tier2_heavy_lite_supported_evidence_table.csv').as_posix()}",
        f"  supported_evidence_table_sha256: {sha256_file(OUT_DIR / 'tier2_heavy_lite_supported_evidence_table.csv')}",
        f"  af3_candidate_list: {(OUT_DIR / 'af3_complex_candidate_list.csv').as_posix()}",
        f"  simplefold_candidate_list: {(OUT_DIR / 'simplefold_candidate_list.csv').as_posix()}",
        f"  electrostatics_pka_review_list: {(OUT_DIR / 'electrostatics_pka_review_list.csv').as_posix()}",
        f"  glycan_review_candidate_list: {(OUT_DIR / 'glycan_review_candidate_list.csv').as_posix()}",
        "counts:",
    ]
    for _, row in by_target.iterrows():
        target = row["target"]
        lines.extend(
            [
                f"  {target}:",
                f"    supported_main_count: {int(row['supported_main_count'])}",
                f"    af3_selected: {int(row['af3_selected'])}",
                f"    simplefold_selected: {int(row['simplefold_selected'])}",
                f"    electrostatics_selected: {int(row['electrostatics_selected'])}",
                f"    glycan_selected: {int(row['glycan_selected'])}",
                f"    controls: {int(row['control_anchor_count'])}",
            ]
        )
    lines.extend(
        [
            f"  total_panel_rows_including_controls: {len(panel)}",
            "manual_unlock:",
            "  required: true",
            "  unlocked_by: null",
            "  unlocked_at: null",
            "locked:",
            "  broad_tier2_heavy: true",
            "  md: true",
            "  final_10k: true",
            "",
        ]
    )
    (OUT_DIR / "tier2_heavy_route_input_manifest.yaml").write_text("\n".join(lines), encoding="utf-8")

    cfg = [
        "stage: Tier2-heavy-lite Stage B",
        "status: draft_only_until_manual_unlock",
        "purpose: route_limited_global_ensemble_review",
        "input_list:",
        f"  path: {panel_path.as_posix()}",
        f"  sha256: {sha256_file(panel_path)}",
        "manual_unlock:",
        "  required: true",
        "  unlocked_by: null",
        "  unlocked_at: null",
        "routes:",
        "  AF3_complex_check:",
        "    scope: selected representatives only",
        "  SimpleFold_antibody_sampling:",
        "    scope: selected representatives only",
        "  electrostatics_pKa_review:",
        "    scope: all locally supported main candidates",
        "  glycan_risk_review:",
        "    scope: glycan-risk flagged representatives only",
        "    current_selection: none_no_glycan_flag_available",
        "locked:",
        "  broad_tier2_heavy: true",
        "  md: true",
        "  final_10k: true",
        "",
    ]
    (OUT_DIR / "tier2_heavy_route_compute_config.yaml").write_text("\n".join(cfg), encoding="utf-8")
    shutil.copyfile(
        OUT_DIR / "tier2_heavy_route_compute_config.yaml",
        OUT_DIR / "route_limited_compute_config.yaml",
    )


def write_audit(
    panel_path: Path,
    panel: pd.DataFrame,
    by_target: pd.DataFrame,
    by_seed: pd.DataFrame,
    by_cluster: pd.DataFrame,
) -> None:
    class_summary = (
        panel[panel["tier2_heavy_lite_role"].eq("main_candidate")]
        .groupby(["target", "heavy_lite_final_class"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    route_count = (
        route_explode(panel)
        .groupby(["target", "route"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "route"])
    )
    seed_view = (
        by_seed.groupby(["target", "route"], group_keys=False)
        .head(12)
        .reset_index(drop=True)
    )
    cluster_view = (
        by_cluster.groupby(["target", "route"], group_keys=False)
        .head(12)
        .reset_index(drop=True)
    )

    report = [
        "# Tier2-heavy-lite Stage B Route-limited Panel Audit",
        "",
        "## Executive Summary",
        "",
        "Verdict: `ROUTE_LIMITED_PANEL_READY_FOR_MANUAL_UNLOCK`.",
        "",
        "This step converts the local Tier2-heavy-lite PASS into a route-limited Stage B review panel. It does not run AF3, SimpleFold, MD, glycan ensemble review, or final 10K selection.",
        "",
        "The panel keeps all locally supported main candidates as the evidence table and assigns only representative subsets to AF3 and SimpleFold. Electrostatics / pKa review remains assigned to every locally supported main candidate.",
        "",
        "## Source",
        "",
        f"- Local heavy-lite results: `{SOURCE_RESULTS.as_posix()}`",
        f"- Source SHA256: `{sha256_file(SOURCE_RESULTS)}`",
        f"- Local next-gate report: `{SOURCE_GATE.as_posix()}`",
        "",
        "## Route Counts By Target",
        "",
        md_table(by_target),
        "",
        "## Main Candidate Class Summary",
        "",
        md_table(class_summary),
        "",
        "## Route Assignment Counts",
        "",
        md_table(route_count),
        "",
        "## Seed Coverage Snapshot",
        "",
        md_table(seed_view),
        "",
        "## Cluster Coverage Snapshot",
        "",
        md_table(cluster_view),
        "",
        "## Interpretation",
        "",
        "- 1E62 keeps 46 locally supported main candidates in the evidence table; 20 are selected for AF3 complex check and 20 for SimpleFold antibody-only sampling.",
        "- sdAb keeps 66 locally supported main candidates in the evidence table; 25 are selected for AF3 complex check and 40 for SimpleFold antibody-only sampling.",
        "- All 112 supported main candidates are assigned to electrostatics / pKa review.",
        "- No glycan route is selected because the local heavy-lite table does not contain a glycan-risk flag. This should remain skipped unless a glycan-near flag is added.",
        "- Controls / anchors are retained in a separate route label for calibration and risk sanity checks; they are not promoted main candidates.",
        "",
        "## Locked Items",
        "",
        "- Broad Tier2-heavy execution remains locked.",
        "- AF3 and SimpleFold compute require explicit manual unlock from this route-limited panel.",
        "- MD remains locked.",
        "- Final 10K remains locked.",
        "",
        "## Output Files",
        "",
        f"- `{(OUT_DIR / 'tier2_heavy_lite_supported_evidence_table.csv').as_posix()}`",
        f"- `{panel_path.as_posix()}`",
        f"- `{(OUT_DIR / 'tier2_route_limited_panel.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'af3_complex_candidate_list.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'simplefold_candidate_list.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'electrostatics_pka_review_list.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'glycan_review_candidate_list.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'tier2_heavy_route_by_target.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'tier2_heavy_route_by_seed.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'tier2_heavy_route_by_cluster.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'tier2_heavy_route_controls.csv').as_posix()}`",
        f"- `{(OUT_DIR / 'tier2_heavy_route_compute_config.yaml').as_posix()}`",
        f"- `{(OUT_DIR / 'route_limited_compute_config.yaml').as_posix()}`",
        "",
    ]
    audit_text = "\n".join(report)
    (OUT_DIR / "tier2_heavy_route_assignment_audit.md").write_text(audit_text, encoding="utf-8")
    shutil.copyfile(
        OUT_DIR / "tier2_heavy_route_assignment_audit.md",
        OUT_DIR / "tier2_route_assignment_audit.md",
    )
    CURRENT_STAGE_REPORT.write_text(audit_text, encoding="utf-8")


def main() -> None:
    df = pd.read_csv(SOURCE_RESULTS)
    df = numeric(
        df,
        [
            "mutation_count",
            "pyrosetta_rosetta_delta_score",
            "pyrosetta_local_structure_validity_t2_score",
            "pyrosetta_local_clash_count",
            "pka_pka_combined_mean",
            "his_pka_support_t2_score",
        ],
    )
    panel = split_routes(df)
    write_outputs(panel)

    by_target = summarize_route_by_target(
        panel[
            panel["control_or_anchor"]
            | (
                panel["tier2_heavy_lite_role"].eq("main_candidate")
                & panel["heavy_lite_final_class"].isin(SUPPORTED_CLASSES)
            )
        ]
    )
    print("route_limited_stage_b_ready")
    print(by_target.to_string(index=False))


if __name__ == "__main__":
    main()
