#!/usr/bin/env python3
"""Collect constrained-MPNN smoke results without requiring score-only outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from analysis.initial_design_generation import collect_p0_mpnn_results as collect
from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT


def markdown_table(df):
    if df.empty:
        return "_No rows._"
    cols = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in df.iterrows():
        values = [str(row[col]).replace("\n", " ") for col in df.columns]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Collect sparse constrained-MPNN smoke outputs.")
    parser.add_argument("--runner-dir", required=True)
    parser.add_argument("--p0-dir", required=True)
    args = parser.parse_args()
    runner_dir = ROOT / args.runner_dir if not Path(args.runner_dir).is_absolute() else Path(args.runner_dir)
    p0_dir = ROOT / args.p0_dir if not Path(args.p0_dir).is_absolute() else Path(args.p0_dir)
    p0_dir.mkdir(parents=True, exist_ok=True)
    pool, audit = collect.collect_generation(runner_dir, p0_dir)
    generated = int(audit["generated_count"].sum()) if not audit.empty else 0
    passed = int(audit["passed_hard_filter_count"].sum()) if not audit.empty else 0
    pass_rate = passed / generated if generated else 0.0
    if pool.empty:
        target_summary = audit.groupby("target", dropna=False).agg(
            generated_count=("generated_count", "sum"),
            passed_hard_filter_count=("passed_hard_filter_count", "sum"),
        ).reset_index()
        target_summary["pass_rate"] = target_summary["passed_hard_filter_count"] / target_summary["generated_count"]
    else:
        pool_flags = pool["liability_flags"].fillna("").astype(str)
        pool = pool.assign(
            new_canonical_nxs_t_failure=pool_flags.str.contains("new_canonical_nxs_t", regex=False),
            non_his_rescue_count_failure=pool_flags.str.contains(
                "non_his_rescue_count_exceeds_max", regex=False
            ),
        )
        target_summary = pool.groupby("target", dropna=False).agg(
            generated_count=("variant_id", "count"),
            passed_hard_filter_count=("hard_filter_status", lambda s: int((s == "pass").sum())),
            unique_sequences=("sequence_hash", "nunique"),
            unique_mutation_lists=("mutation_list", "nunique"),
            max_rescue_count=("rescue_count", "max"),
            new_canonical_nxs_t_failure_count=("new_canonical_nxs_t_failure", "sum"),
            non_his_rescue_count_failure_count=("non_his_rescue_count_failure", "sum"),
        ).reset_index()
        target_summary["pass_rate"] = (
            target_summary["passed_hard_filter_count"] / target_summary["generated_count"]
        )
        target_summary["new_canonical_nxs_t_failure_fraction"] = (
            target_summary["new_canonical_nxs_t_failure_count"] / target_summary["generated_count"]
        )
        target_summary["non_his_rescue_count_failure_fraction"] = (
            target_summary["non_his_rescue_count_failure_count"] / target_summary["generated_count"]
        )
    lines = [
        "# Sparse Constrained-MPNN Smoke Report",
        "",
        f"Generated candidates: {generated}",
        f"Hard-filter pass: {passed}",
        f"Hard-filter pass rate: {pass_rate:.4f}",
        "",
        "Acceptance thresholds from review:",
        "",
        "- 1E62 pass rate: preferred >= 30%, minimum >= 15%.",
        "- sdAb pass rate: preferred >= 20%, minimum >= 10%.",
        "- New canonical N-X-S/T failure: < 5%.",
        "- Non-His rescue-count failure: < 10%.",
        "",
        "## Target Summary",
        "",
        markdown_table(target_summary),
    ]
    (p0_dir / "sparse_constrained_mpnn_smoke_report.md").write_text("\n".join(lines) + "\n")
    print(
        {
            "status": "sparse_constrained_smoke_collected",
            "generated_candidates": generated,
            "passed_hard_filter": passed,
            "pass_rate": pass_rate,
        }
    )


if __name__ == "__main__":
    main()
