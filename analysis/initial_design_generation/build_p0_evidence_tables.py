#!/usr/bin/env python3
"""Build P0 evidence tables after ProteinMPNN compute has completed."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT
P0 = ROOT / "results/initial_design_generation/p0_mpnn"
PRODUCTION = ROOT / "results/initial_design_generation/production_initial_pool/production_initial_pool_candidates_all.csv"


def normalize_flags(value: object) -> list[str]:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return ["none"]
    return [item for item in str(value).split(";") if item]


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy() if max_rows else df.copy()
    cols = [str(c) for c in view.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in view.iterrows():
        values = [str(row[col]).replace("\n", " ") for col in view.columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_production_features(p0_dir: Path) -> pd.DataFrame:
    pool = pd.read_csv(PRODUCTION)
    score = pd.read_csv(p0_dir / "mpnn_scores_current_pool.csv")
    score_cols = [
        "variant_id",
        "mpnn_score_status",
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
    merged = pool.merge(score[score_cols], on="variant_id", how="left", validate="one_to_one")
    dry.write_csv(merged, p0_dir / "production_pool_features.csv")
    return merged


def write_score_audits(features: pd.DataFrame, p0_dir: Path) -> None:
    target_summary = (
        features.groupby("target", dropna=False)
        .agg(
            candidate_count=("variant_id", "count"),
            scored_count=("mpnn_score_status", lambda s: int((s == "scored_by_mpnn").sum())),
            mpnn_score_mean=("mpnn_total_score_per_residue", "mean"),
            mpnn_score_median=("mpnn_total_score_per_residue", "median"),
            mpnn_score_p05=("mpnn_total_score_per_residue", lambda s: s.quantile(0.05)),
            mpnn_score_p95=("mpnn_total_score_per_residue", lambda s: s.quantile(0.95)),
            mpnn_parent_delta_mean=("mpnn_parent_delta", "mean"),
            mpnn_parent_delta_median=("mpnn_parent_delta", "median"),
        )
        .reset_index()
    )
    dry.write_csv(target_summary, p0_dir / "mpnn_score_distribution_by_target.csv")

    high_ph = (
        (features["acidic_release_support_score"].fillna(0) >= 0.65)
        & (features["neutral_retention_score"].fillna(0) >= 0.70)
        & (features["His_count"].fillna(0) >= 1)
    )
    favorable_mpnn = features["mpnn_route_percentile"].fillna(1) <= 0.10
    weak_ph = (
        (features["acidic_release_support_score"].fillna(0) < 0.50)
        | (features["His_count"].fillna(0) < 1)
    )
    unfavorable_mpnn = features["mpnn_route_percentile"].fillna(0) >= 0.90

    low_score_high_support = features[favorable_mpnn & high_ph].copy()
    low_score_high_support["reason"] = "favorable_mpnn_percentile_and_high_pH_support"
    low_score_high_support = low_score_high_support.sort_values(
        ["target", "mpnn_route_percentile", "hit_likelihood_score_v0"],
        ascending=[True, True, False],
    )
    dry.write_csv(low_score_high_support.head(1000), p0_dir / "mpnn_low_score_high_pH_support_candidates.csv")

    high_score_no_mechanism = features[unfavorable_mpnn & weak_ph].copy()
    high_score_no_mechanism["reason"] = "unfavorable_mpnn_percentile_and_weak_pH_mechanism"
    high_score_no_mechanism = high_score_no_mechanism.sort_values(
        ["target", "mpnn_route_percentile", "hit_likelihood_score_v0"],
        ascending=[True, False, True],
    )
    dry.write_csv(high_score_no_mechanism.head(1000), p0_dir / "mpnn_high_score_no_pH_mechanism_candidates.csv")


def write_relaxed_failure_summary(p0_dir: Path) -> pd.DataFrame:
    pools = []
    for target in ["1E62", "sdAb"]:
        path = p0_dir / f"relaxed_mpnn_counterfactual_pool_{target}.csv"
        if path.exists():
            pools.append(pd.read_csv(path))
    if not pools:
        summary = pd.DataFrame()
        dry.write_csv(summary, p0_dir / "relaxed_counterfactual_failure_summary.csv")
        return summary
    relaxed = pd.concat(pools, ignore_index=True)
    rows: list[dict] = []
    for row in relaxed.itertuples(index=False):
        flags = normalize_flags(getattr(row, "liability_flags", ""))
        for flag in flags:
            rows.append(
                {
                    "target": row.target,
                    "failure_reason": flag,
                    "route_or_mode": row.primary_generation_route,
                    "his_seed_set": getattr(row, "his_seed_set", ""),
                    "His_count": int(getattr(row, "His_count", 0)),
                    "mutation_count": int(getattr(row, "mutation_count", 0)),
                    "new_canonical_nxs_t_flag": flag == "new_canonical_nxs_t" or "new_canonical_nxs_t" in flags,
                    "non_his_rescue_count": int(getattr(row, "rescue_count", 0)),
                    "forbidden_pair_flag": "sdAb_V105H_D110H" in flags,
                    "hard_protect_violation": any(str(x).startswith("hard_protect_") for x in flags),
                }
            )
    exploded = pd.DataFrame(rows)
    summary = (
        exploded.groupby(
            [
                "target",
                "route_or_mode",
                "failure_reason",
                "His_count",
                "mutation_count",
                "new_canonical_nxs_t_flag",
                "non_his_rescue_count",
                "forbidden_pair_flag",
                "hard_protect_violation",
                "his_seed_set",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="count")
    )
    totals = relaxed.groupby(["target", "primary_generation_route"], dropna=False).size().rename("total").reset_index()
    summary = summary.merge(
        totals,
        left_on=["target", "route_or_mode"],
        right_on=["target", "primary_generation_route"],
        how="left",
    ).drop(columns=["primary_generation_route"])
    summary["fraction"] = summary["count"] / summary["total"]
    summary = summary.sort_values(["target", "route_or_mode", "count"], ascending=[True, True, False])
    dry.write_csv(summary, p0_dir / "relaxed_counterfactual_failure_summary.csv")
    return summary


def write_reports(features: pd.DataFrame, failure_summary: pd.DataFrame, p0_dir: Path) -> None:
    score = pd.read_csv(p0_dir / "mpnn_scores_current_pool.csv")
    con_audit = pd.read_csv(p0_dir / "constrained_mpnn_generation_collection_audit.csv")
    rel_audit = pd.read_csv(p0_dir / "relaxed_mpnn_generation_collection_audit.csv")
    route = pd.read_csv(p0_dir / "mpnn_score_distribution_by_route.csv")
    constrained_generated = int(con_audit["generated_count"].sum())
    constrained_pass = int(con_audit["passed_hard_filter_count"].sum())
    relaxed_generated = int(rel_audit["generated_count"].sum())
    relaxed_pass = int(rel_audit["passed_hard_filter_count"].sum())

    lines = [
        "# P0 ProteinMPNN Evidence Report",
        "",
        "Status: `p0_compute_complete_evidence_partial`",
        "",
        "## Compute Completeness",
        "",
        f"- Score-only rows: {int((score['mpnn_score_status'] == 'scored_by_mpnn').sum())} / {len(score)}",
        f"- Constrained generated: {constrained_generated}; hard-filter pass: {constrained_pass}",
        f"- Relaxed generated: {relaxed_generated}; hard-filter pass: {relaxed_pass}",
        "",
        "## Score-Only Decision",
        "",
        "Score-only ProteinMPNN features are complete and may be used for P0 evidence review. Use target/route/seed/rescue-specific ranks, not a global cutoff.",
        "",
        "## Route Score Summary",
        "",
        markdown_table(route),
        "",
        "## Generation Branch Decisions",
        "",
        "- Current constrained-MPNN generated candidates are rejected for downstream use because zero candidates passed hard filters.",
        "- Relaxed-MPNN / MPNN-only candidates remain audit-only and have no preallocated final-library quota.",
        "",
        "## Required Next Step",
        "",
        "Repair constrained-MPNN input generation to enforce sparse rescue subsets, run a constrained smoke test, then decide whether to rerun the full constrained branch or exclude it from this round.",
    ]
    (p0_dir / "p0_mpnn_evidence_report.md").write_text("\n".join(lines) + "\n")

    comparison = [
        "# P0 ProteinMPNN Comparison Report",
        "",
        "Status: `p0_compute_complete_evidence_partial`",
        "",
        "This report reflects true ProteinMPNN compute results. It does not unlock Tier 2 or final 10K selection.",
        "",
        "## Candidate Universes",
        "",
        "| universe | status | count | downstream decision |",
        "|---|---|---:|---|",
        f"| current production pool + score-only features | scored | {len(features)} | use for P0 evidence review |",
        f"| constrained-MPNN generated pool | generated, hard-filter failed | {constrained_generated} | exclude until sparse rescue repair |",
        f"| relaxed counterfactual pool | generated, audit-only | {relaxed_generated} | audit only |",
        "",
        "## Key Metrics",
        "",
        f"- Score-only scored rows: {int((score['mpnn_score_status'] == 'scored_by_mpnn').sum())} / {len(score)}",
        f"- Constrained hard-filter pass: {constrained_pass} / {constrained_generated}",
        f"- Relaxed hard-filter pass: {relaxed_pass} / {relaxed_generated}",
        "",
        "## Tier Unlock",
        "",
        "Tier 2-core remains locked until P0 evidence review is complete and constrained-MPNN is either repaired or explicitly excluded from the current round.",
    ]
    (p0_dir / "p0_mpnn_comparison_report.md").write_text("\n".join(comparison) + "\n")

    relaxed_lines = [
        "# Relaxed-MPNN / MPNN-only Counterfactual Audit",
        "",
        "Status: `compute_complete_audit_only`",
        "",
        f"Generated candidates: {relaxed_generated}",
        f"Hard-filter pass: {relaxed_pass}",
        "",
        "This branch remains audit-only. It does not reserve or unlock any final 10K allocation.",
        "",
        "## Failure Summary",
        "",
        markdown_table(failure_summary, max_rows=40) if not failure_summary.empty else "No relaxed candidates were collected.",
    ]
    (p0_dir / "relaxed_mpnn_counterfactual_audit.md").write_text("\n".join(relaxed_lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build P0 evidence tables after ProteinMPNN compute.")
    parser.add_argument("--p0-dir", default=str(P0.relative_to(ROOT)))
    args = parser.parse_args()
    p0_dir = ROOT / args.p0_dir
    features = build_production_features(p0_dir)
    write_score_audits(features, p0_dir)
    failure_summary = write_relaxed_failure_summary(p0_dir)
    write_reports(features, failure_summary, p0_dir)
    print(
        {
            "status": "p0_evidence_tables_built",
            "production_feature_rows": len(features),
            "relaxed_failure_summary_rows": len(failure_summary),
        }
    )


if __name__ == "__main__":
    main()
