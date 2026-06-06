#!/usr/bin/env python3
"""Review Tier2 Stage-1.5 calibration labels and freeze the next decision.

This is table-only. It does not build a Tier2-C candidate list and does not
launch structural compute.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


DEFAULT_DIAG = Path("results/initial_design_generation/tier2_stage1_diagnostics")
DEFAULT_OUT = Path("results/initial_design_generation/tier2_stage1_calibration_review")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def markdown_table(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy()
    lines = [
        "| " + " | ".join(view.columns.astype(str)) + " |",
        "| " + " | ".join("---" for _ in view.columns) + " |",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]).replace("\n", " ") for c in view.columns) + " |")
    return "\n".join(lines)


def rate_summary(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for keys, sub in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row["count"] = len(sub)
        row["strong_good_count"] = int(sub["tier2_class"].isin({"T2_strong_candidate", "T2_good_candidate"}).sum())
        row["strong_good_rate"] = round(row["strong_good_count"] / len(sub), 6) if len(sub) else 0.0
        row["structure_risk_count"] = int(sub["tier2_class"].eq("T2_structure_risk").sum())
        row["structure_risk_rate"] = round(row["structure_risk_count"] / len(sub), 6) if len(sub) else 0.0
        for col in [
            "neutral_retention_t2_score",
            "acidic_release_mechanism_t2_score",
            "local_structure_validity_t2_score",
            "rosetta_delta_score",
            "stage1_5_review_support_score",
        ]:
            if col in sub.columns:
                vals = pd.to_numeric(sub[col], errors="coerce")
                row[f"{col}_median"] = round(float(vals.median()), 6) if vals.notna().any() else pd.NA
        rows.append(row)
    return pd.DataFrame(rows).sort_values(group_cols + ["count"], ascending=[True] * len(group_cols) + [False])


def consistency_checks(cal: pd.DataFrame, boundary: pd.DataFrame, targeted: pd.DataFrame) -> pd.DataFrame:
    checks: list[dict[str, object]] = []

    def add(name: str, status: str, details: object = "") -> None:
        checks.append({"check": name, "status": status, "details": details})

    add("calibrated_row_count_7000", "PASS" if len(cal) == 7000 else "FAIL", len(cal))
    add("variant_id_unique", "PASS" if cal["variant_id"].is_unique else "FAIL", int(cal["variant_id"].duplicated().sum()))

    required = [
        "stage1_5_structure_risk_band",
        "stage1_5_recommended_action",
        "stage1_5_seed_policy",
        "stage1_5_review_support_score",
    ]
    missing = [c for c in required if c not in cal.columns]
    add("required_stage1_5_columns", "PASS" if not missing else "FAIL", ";".join(missing))

    severe = cal["stage1_5_structure_risk_band"].eq("severe_structure_risk")
    bad_severe = int((severe & ~cal["stage1_5_recommended_action"].eq("hold_severe_structure_risk")).sum())
    add("severe_structure_risk_held", "PASS" if bad_severe == 0 else "FAIL", bad_severe)

    controls = cal["tier2_class"].eq("T2_control_anchor")
    bad_controls = int((controls & ~cal["stage1_5_recommended_action"].eq("keep_control_anchor")).sum())
    add("control_anchors_preserved", "PASS" if bad_controls == 0 else "FAIL", bad_controls)

    boundary_severe = 0
    if not boundary.empty and "stage1_5_structure_risk_band" in boundary.columns:
        boundary_severe = int(boundary["stage1_5_structure_risk_band"].eq("severe_structure_risk").sum())
    add("boundary_queue_excludes_severe", "PASS" if boundary_severe == 0 else "FAIL", boundary_severe)

    targeted_risk = 0
    if not targeted.empty and "tier2_class" in targeted.columns:
        targeted_risk = int(targeted["tier2_class"].eq("T2_structure_risk").sum())
    add("targeted_queue_excludes_structure_risk", "PASS" if targeted_risk == 0 else "FAIL", targeted_risk)

    ae108h_targeted = 0
    if not targeted.empty and "stage1_5_seed_policy" in targeted.columns:
        ae108h_targeted = int(targeted["stage1_5_seed_policy"].eq("deprioritize_ae108h_seed").sum())
    add("targeted_queue_excludes_ae108h_hold", "PASS" if ae108h_targeted == 0 else "FAIL", ae108h_targeted)
    return pd.DataFrame(checks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics-dir", default=str(DEFAULT_DIAG))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    diag = Path(args.diagnostics_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cal = pd.read_csv(diag / "tier2_stage1_calibrated_results.csv", low_memory=False)
    boundary = pd.read_csv(diag / "tier2_boundary_review_candidates.csv", low_memory=False)
    targeted = pd.read_csv(diag / "tier2_targeted_pool_review_candidates.csv", low_memory=False)

    checks = consistency_checks(cal, boundary, targeted)
    action_summary = cal.groupby(["target", "stage1_5_recommended_action"], dropna=False).size().reset_index(name="count")
    action_summary["fraction_of_target"] = (
        action_summary["count"] / action_summary.groupby("target")["count"].transform("sum")
    ).round(6)
    boundary_summary = rate_summary(boundary, ["target", "stage1_5_structure_risk_band", "stage1_5_seed_policy"])
    targeted_route_summary = rate_summary(targeted, ["target", "stage1_5_seed_policy", "primary_generation_route"])
    targeted_seed_summary = rate_summary(targeted, ["target", "his_seed_set"])

    write_csv(checks, out_dir / "calibration_review_consistency_checks.csv")
    write_csv(action_summary, out_dir / "calibration_review_action_summary.csv")
    write_csv(boundary_summary, out_dir / "boundary_review_summary.csv")
    write_csv(targeted_route_summary, out_dir / "targeted_pool_review_route_summary.csv")
    write_csv(targeted_seed_summary, out_dir / "targeted_pool_review_seed_summary.csv")

    pass_all = checks["status"].eq("PASS").all()
    unlock_status = {
        "tier2c_candidate_list_preparation": "conditional_go" if pass_all else "hold",
        "tier2c_structural_compute": "hold",
        "tier2_heavy": "hold",
        "final_10k_selection": "no_go",
        "conditions": [
            "Use Stage-1.5 labels as candidate-list policy only; do not treat them as final biological evidence.",
            "Expand only priority or secondary seed/route regions supported by Stage-1.",
            "Keep severe_structure_risk held.",
            "Keep sdAb AE108H broad expansion held.",
            "Run a Stage-2 candidate-list audit before any additional PyRosetta/pKa/FoldX compute.",
        ],
        "target_policy": {
            "Ab_1E62": {
                "status": "candidate_list_preparation_allowed",
                "priority_seeds": ["LQ35H;LY38H", "LK24H;LY38H", "LK24H;LQ35H", "LY31H;LQ35H"],
                "mutation_order": "<=3 preferred; 4-mut only for priority seed/route support",
            },
            "Ab_sdAb": {
                "status": "narrow_candidate_list_preparation_allowed",
                "priority_seeds": ["AD110H;AY111H"],
                "secondary_seeds": ["AQ100H;AD110H", "AQ100H;AY111H"],
                "holds": ["broad AE108H-centered expansion", "severe_structure_risk"],
                "mutation_order": "<=3 preferred; 4-mut cautious only with priority seed support",
            },
        },
    }
    (out_dir / "tier2c_unlock_recommendation.yaml").write_text(yaml.safe_dump(unlock_status, sort_keys=False), encoding="utf-8")

    report = [
        "# Tier2 Stage-1.5 Calibration Review",
        "",
        "Status: `calibration_review_complete`",
        "",
        "This review uses Stage-1.5 diagnostic tables only. It does not run Tier2-C, Tier2-heavy, AF3, SimpleFold, MD, or final selection.",
        "",
        "## Verdict",
        "",
        "```text",
        f"Tier2-C candidate-list preparation: {unlock_status['tier2c_candidate_list_preparation']}",
        "Tier2-C structural compute: HOLD",
        "Tier2-heavy: HOLD",
        "Final 10K selection: NO-GO",
        "```",
        "",
        "Interpretation: the calibration labels are internally consistent enough to prepare a targeted Stage-2 candidate list, but not enough to launch additional structure compute without a new candidate-list audit.",
        "",
        "## Consistency Checks",
        markdown_table(checks, 20),
        "",
        "## Action Summary",
        markdown_table(action_summary, 40),
        "",
        "## Boundary Review Queue",
        "",
        f"Boundary queue rows: {len(boundary)}. These rows are reviewable exceptions, not automatic promotions.",
        "",
        markdown_table(boundary_summary, 30),
        "",
        "## Targeted Pool Review Queue",
        "",
        f"Targeted review rows: {len(targeted)}. These rows define the initial policy for Stage-2 candidate-list preparation.",
        "",
        "### By Route",
        markdown_table(targeted_route_summary, 40),
        "",
        "### By Seed",
        markdown_table(targeted_seed_summary, 40),
        "",
        "## Recommended Stage-2 Policy",
        "",
        "1E62:",
        "",
        "- Prepare a targeted Tier2-C candidate list around `LQ35H;LY38H`, `LK24H;LY38H`, `LK24H;LQ35H`, and `LY31H;LQ35H`.",
        "- Prefer variants with no more than 3 mutations; allow 4-mut only when route/seed evidence is strong.",
        "- Do not broadly rescue all F-class rows.",
        "",
        "sdAb:",
        "",
        "- Prepare a narrow Tier2-C candidate list around `AD110H;AY111H` first.",
        "- Keep `AQ100H;AD110H` and `AQ100H;AY111H` as secondary seeds.",
        "- Keep broad `AE108H`-centered expansion held.",
        "- Keep severe structure-risk rows held.",
        "",
        "Global:",
        "",
        "- Do not globally relax Rosetta/local thresholds.",
        "- Use `manual_boundary_review` only as a limited exception source.",
        "- Build and audit the Stage-2 list before running additional PyRosetta/pKa/FoldX.",
    ]
    (out_dir / "tier2_stage1_calibration_review_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
