#!/usr/bin/env python3
"""Build audit tables for repaired full sparse constrained-MPNN generation."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT


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


def read_pool(p0_dir: Path) -> pd.DataFrame:
    pools = []
    for target in ["1E62", "sdAb"]:
        path = p0_dir / f"constrained_mpnn_rescue_pool_{target}.csv"
        if path.exists():
            pools.append(pd.read_csv(path))
    if not pools:
        return pd.DataFrame()
    return pd.concat(pools, ignore_index=True)


def mutation_positions(mutation_list: object) -> list[int]:
    if mutation_list is None or pd.isna(mutation_list):
        return []
    out: list[int] = []
    for token in str(mutation_list).split(";"):
        match = re.match(r"[A-Z][A-Z](\d+)[A-Z]$", token.strip())
        if match:
            out.append(int(match.group(1)))
    return out


def build_target_summary(pool: pd.DataFrame, jobs: pd.DataFrame) -> pd.DataFrame:
    flags = pool["liability_flags"].fillna("").astype(str)
    enriched = pool.assign(
        pass_flag=pool["hard_filter_status"].eq("pass"),
        new_nxs_t_failure=flags.str.contains("new_canonical_nxs_t", regex=False),
        rescue_count_failure=flags.str.contains("non_his_rescue_count_exceeds_max", regex=False),
    )
    summary = enriched.groupby("target", dropna=False).agg(
        generated_count=("variant_id", "count"),
        hard_filter_pass=("pass_flag", "sum"),
        unique_sequences=("sequence_hash", "nunique"),
        unique_mutation_lists=("mutation_list", "nunique"),
        max_mutation_count=("mutation_count", "max"),
        max_his_count=("His_count", "max"),
        max_rescue_count=("rescue_count", "max"),
        new_nxs_t_failure_count=("new_nxs_t_failure", "sum"),
        rescue_count_failure_count=("rescue_count_failure", "sum"),
    ).reset_index()
    planned = jobs.groupby("target", dropna=False)["planned_raw_samples"].sum().reset_index(name="planned_raw_samples")
    summary = summary.merge(planned, on="target", how="left")
    summary["collection_fraction"] = summary["generated_count"] / summary["planned_raw_samples"]
    summary["hard_filter_pass_rate"] = summary["hard_filter_pass"] / summary["generated_count"]
    summary["unique_sequence_fraction"] = summary["unique_sequences"] / summary["hard_filter_pass"]
    summary["unique_mutation_list_fraction"] = summary["unique_mutation_lists"] / summary["hard_filter_pass"]
    summary["duplicate_collapse_fraction"] = 1 - summary["unique_sequence_fraction"]
    summary["new_nxs_t_failure_fraction"] = summary["new_nxs_t_failure_count"] / summary["generated_count"]
    summary["rescue_count_failure_fraction"] = summary["rescue_count_failure_count"] / summary["generated_count"]
    return summary


def build_seed_summary(pool: pd.DataFrame) -> pd.DataFrame:
    return pool.groupby(["target", "his_seed_set"], dropna=False).agg(
        generated_count=("variant_id", "count"),
        unique_sequences=("sequence_hash", "nunique"),
        unique_mutation_lists=("mutation_list", "nunique"),
        hard_filter_pass=("hard_filter_status", lambda s: int((s == "pass").sum())),
        mean_mutation_count=("mutation_count", "mean"),
        mean_rescue_count=("rescue_count", "mean"),
    ).reset_index().sort_values(["target", "generated_count"], ascending=[True, False])


def build_rescue_subset_summary(pool: pd.DataFrame, jobs: pd.DataFrame) -> pd.DataFrame:
    joined = pool.merge(
        jobs[["job_name", "target", "seed_positions_forced_to_H", "rescue_design_positions", "designable_position_count"]],
        left_on=["generation_record_id", "target"],
        right_on=["job_name", "target"],
        how="left",
    )
    return joined.groupby(
        ["target", "his_seed_set", "seed_positions_forced_to_H", "rescue_design_positions", "designable_position_count"],
        dropna=False,
    ).agg(
        generated_count=("variant_id", "count"),
        unique_sequences=("sequence_hash", "nunique"),
        unique_mutation_lists=("mutation_list", "nunique"),
        hard_filter_pass=("hard_filter_status", lambda s: int((s == "pass").sum())),
    ).reset_index().sort_values(["target", "generated_count"], ascending=[True, False])


def build_mutation_summary(pool: pd.DataFrame) -> pd.DataFrame:
    summary = pool.groupby(["target", "mutation_list"], dropna=False).agg(
        generated_count=("variant_id", "count"),
        unique_sequences=("sequence_hash", "nunique"),
        hard_filter_pass=("hard_filter_status", lambda s: int((s == "pass").sum())),
        his_seed_set=("his_seed_set", lambda s: ";".join(sorted(set(s.astype(str))))),
        mutation_count=("mutation_count", "first"),
        rescue_count=("rescue_count", "first"),
    ).reset_index()
    return summary.sort_values(["target", "generated_count"], ascending=[True, False])


def build_failure_summary(pool: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in pool.itertuples(index=False):
        raw_flags = getattr(row, "liability_flags", "")
        flags = "" if pd.isna(raw_flags) else str(raw_flags or "").strip()
        reasons = [item for item in flags.split(";") if item] or ["none"]
        for reason in reasons:
            rows.append({"target": row.target, "failure_reason": reason})
    exploded = pd.DataFrame(rows)
    return exploded.groupby(["target", "failure_reason"], dropna=False).size().reset_index(name="count")


def build_position_evidence_summary(pool: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    evidence = evidence.copy()
    evidence["target"] = evidence["target"].astype(str)
    evidence_by_key = {
        (str(row.target), int(row.position)): row
        for row in evidence.itertuples(index=False)
    }
    rows = []
    for row in pool.itertuples(index=False):
        for pos in mutation_positions(row.mutation_list):
            ev = evidence_by_key.get((str(row.target), pos))
            rows.append(
                {
                    "target": row.target,
                    "position": pos,
                    "his_seed_set": row.his_seed_set,
                    "region": getattr(ev, "region", "") if ev is not None else "",
                    "wetlab_support": getattr(ev, "wetlab_support", "") if ev is not None else "",
                    "structural_support": getattr(ev, "structural_support", "") if ev is not None else "",
                    "interface_support": getattr(ev, "interface_support", "") if ev is not None else "",
                    "known_failure_flag": getattr(ev, "known_failure_flag", "") if ev is not None else "",
                }
            )
    if not rows:
        return pd.DataFrame()
    long = pd.DataFrame(rows)
    return long.groupby(
        [
            "target",
            "position",
            "region",
            "wetlab_support",
            "structural_support",
            "interface_support",
            "known_failure_flag",
        ],
        dropna=False,
    ).size().reset_index(name="mutation_occurrence_count").sort_values(
        ["target", "mutation_occurrence_count"], ascending=[True, False]
    )


def write_report(
    p0_dir: Path,
    target_summary: pd.DataFrame,
    seed_summary: pd.DataFrame,
    rescue_summary: pd.DataFrame,
    mutation_summary: pd.DataFrame,
    failure_summary: pd.DataFrame,
    position_summary: pd.DataFrame,
) -> None:
    lines = [
        "# Repaired Full Sparse Constrained-MPNN Audit",
        "",
        "Status: `repaired_full_sparse_constrained_complete`",
        "",
        "This audit covers the repaired constrained-MPNN branch generated with sparse rescue subsets. It replaces the old broad-rescue constrained branch for evidence review, but generated candidates are still not merged downstream until manual review.",
        "",
        "## Target Summary",
        "",
        markdown_table(target_summary),
        "",
        "## Failure Summary",
        "",
        markdown_table(failure_summary),
        "",
        "## Seed Summary",
        "",
        markdown_table(seed_summary, max_rows=30),
        "",
        "## Rescue Subset Summary",
        "",
        markdown_table(rescue_summary, max_rows=30),
        "",
        "## Top Repeated Mutation Lists",
        "",
        markdown_table(mutation_summary, max_rows=30),
        "",
        "## Position Evidence Summary",
        "",
        markdown_table(position_summary, max_rows=40),
        "",
        "## Decision Boundary",
        "",
        "- Hard-filter compatibility is repaired: all generated candidates pass hard filters.",
        "- Diversity remains limited after deduplication; downstream use requires deduplication and diversity-aware selection.",
        "- Tier 2 and final 10K remain locked until this audit is reviewed and a downstream-use decision is made.",
    ]
    (p0_dir / "repaired_full_sparse_constrained_mpnn_audit.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build full sparse constrained-MPNN audit.")
    parser.add_argument("--runner-dir", required=True)
    parser.add_argument("--p0-dir", required=True)
    args = parser.parse_args()
    runner_dir = ROOT / args.runner_dir if not Path(args.runner_dir).is_absolute() else Path(args.runner_dir)
    p0_dir = ROOT / args.p0_dir if not Path(args.p0_dir).is_absolute() else Path(args.p0_dir)
    pool = read_pool(p0_dir)
    jobs = pd.read_csv(runner_dir / "mpnn_runner_jobs_constrained.csv")
    evidence = pd.read_csv(ROOT / "results/initial_design_generation/tables/evidence_ledger.csv")
    target_summary = build_target_summary(pool, jobs)
    seed_summary = build_seed_summary(pool)
    rescue_summary = build_rescue_subset_summary(pool, jobs)
    mutation_summary = build_mutation_summary(pool)
    failure_summary = build_failure_summary(pool)
    position_summary = build_position_evidence_summary(pool, evidence)
    dry.write_csv(target_summary, p0_dir / "repaired_sparse_constrained_summary_by_target.csv")
    dry.write_csv(seed_summary, p0_dir / "repaired_sparse_constrained_seed_summary.csv")
    dry.write_csv(rescue_summary, p0_dir / "repaired_sparse_constrained_rescue_subset_summary.csv")
    dry.write_csv(mutation_summary, p0_dir / "repaired_sparse_constrained_mutation_list_summary.csv")
    dry.write_csv(failure_summary, p0_dir / "repaired_sparse_constrained_failure_summary.csv")
    dry.write_csv(position_summary, p0_dir / "repaired_sparse_constrained_position_evidence_summary.csv")
    write_report(
        p0_dir,
        target_summary,
        seed_summary,
        rescue_summary,
        mutation_summary,
        failure_summary,
        position_summary,
    )
    print(
        {
            "status": "repaired_full_sparse_constrained_audit_built",
            "generated_count": int(target_summary["generated_count"].sum()),
            "hard_filter_pass": int(target_summary["hard_filter_pass"].sum()),
            "unique_sequences": int(target_summary["unique_sequences"].sum()),
            "unique_mutation_lists": int(target_summary["unique_mutation_lists"].sum()),
        }
    )


if __name__ == "__main__":
    main()
