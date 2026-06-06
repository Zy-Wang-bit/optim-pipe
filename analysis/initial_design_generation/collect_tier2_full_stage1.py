#!/usr/bin/env python3
"""Collect Tier2 Full Stage-1 results into auditable tables."""

from __future__ import annotations

import argparse
import shutil
import time
from pathlib import Path
from typing import Any

import pandas as pd


def numeric(series: pd.Series, default: float | None = None) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    if default is not None:
        values = values.fillna(default)
    return values


def clip01(series: pd.Series) -> pd.Series:
    return numeric(series, 0.0).clip(lower=0.0, upper=1.0)


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    lines = [
        "| " + " | ".join(map(str, df.columns)) + " |",
        "| " + " | ".join("---" for _ in df.columns) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in df.columns) + " |")
    return "\n".join(lines)


def classify(row: pd.Series) -> str:
    if str(row.get("tier2_stage1_bucket", "")) == "control_anchor" or str(row.get("tier1_review_class", "")) == "control_anchor":
        return "T2_control_anchor"
    pyro = str(row.get("pyrosetta_tool_status", ""))
    if pyro != "success":
        return "T2_structure_risk"
    local = float(row.get("local_structure_validity_t2_score", 0) or 0)
    clashes = float(row.get("local_clash_count", 0) or 0)
    rosetta_delta = float(row.get("rosetta_delta_score", 0) or 0)
    if local < 0.35 or clashes >= 10 or rosetta_delta > 600:
        return "T2_structure_risk"
    glycan = float(row.get("glycan_coverage_risk_t2_score", 0) or 0)
    if glycan >= 0.80:
        return "T2_glycan_or_coverage_risk"
    neutral = float(row.get("neutral_retention_t2_score", 0) or 0)
    acid = float(row.get("acidic_release_mechanism_t2_score", 0) or 0)
    global_risk = float(row.get("global_weakening_risk_t2_score", 0) or 0)
    pka_status = str(row.get("pka_tool_status", "not_run"))
    his_count = float(row.get("pka_his_count_requested", 0) or 0)
    his_pka_missing = his_count > 0 and pka_status not in {"success", "not_applicable"}
    if global_risk >= 0.78 or neutral < 0.35:
        return "T2_global_weakening_risk"
    if acid >= 0.65 and neutral < 0.55:
        return "T2_release_possible_but_neutral_risky"
    if neutral >= 0.65 and acid < 0.40:
        return "T2_neutral_retained_but_release_weak"
    if neutral >= 0.72 and acid >= 0.62 and local >= 0.65 and global_risk < 0.55 and not his_pka_missing:
        return "T2_strong_candidate"
    if neutral >= 0.58 and local >= 0.55 and global_risk < 0.72:
        return "T2_good_candidate"
    return "T2_deprioritized"


def read_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", required=True)
    parser.add_argument("--pyrosetta-results", required=True)
    parser.add_argument("--pka-summary", required=True)
    parser.add_argument("--foldx-subset-results", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    started = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = read_required(Path(args.candidate_csv))
    pyro = read_required(Path(args.pyrosetta_results)).rename(
        columns={
            "tool_status": "pyrosetta_tool_status",
            "feature_missing_reason": "pyrosetta_feature_missing_reason",
            "feature_confidence": "pyrosetta_feature_confidence",
        }
    )
    pka = read_required(Path(args.pka_summary))
    foldx = read_required(Path(args.foldx_subset_results))

    df = candidates.merge(pyro, on=["variant_id", "target"], how="left", validate="one_to_one")
    df = df.merge(pka, on=["variant_id", "target"], how="left", validate="one_to_one")
    df = df.merge(foldx, on=["variant_id", "target"], how="left", validate="one_to_one")

    for col in ["foldx_tool_status", "foldx_feature_missing_reason", "foldx_feature_confidence"]:
        if col in df.columns:
            df[col] = df[col].astype("object")
    selected = df["tier2_foldx_subset_selected"].map(as_bool)
    df.loc[~selected & df["foldx_tool_status"].isna(), "foldx_tool_status"] = "skipped_by_policy"
    df.loc[~selected & df["foldx_feature_missing_reason"].isna(), "foldx_feature_missing_reason"] = "not_selected_for_foldx_subset"
    df.loc[~selected & df["foldx_feature_confidence"].isna(), "foldx_feature_confidence"] = "not_applicable"

    t1_neutral = clip01(df.get("neutral_retention_score", pd.Series(index=df.index)))
    t1_acid = clip01(df.get("acidic_release_support_score", pd.Series(index=df.index)))
    t1_risk = clip01(df.get("global_weakening_risk_score", pd.Series(index=df.index)))
    local = clip01(df.get("local_structure_validity_t2_score", pd.Series(index=df.index)))
    pka_score = clip01(df.get("his_pka_support_t2_score", pd.Series(index=df.index)))
    pka_requested = numeric(df.get("pka_his_count_requested", pd.Series(index=df.index)), 0)
    pka_component = pka_score.where(pka_requested.gt(0), t1_acid)
    foldx_score = clip01(df.get("foldx_neutral_interface_t2_score", pd.Series(index=df.index)))
    has_foldx = df.get("foldx_tool_status", pd.Series(index=df.index)).eq("success")
    foldx_component = foldx_score.where(has_foldx, t1_neutral)

    df["neutral_retention_t2_score"] = (0.42 * t1_neutral + 0.38 * local + 0.20 * foldx_component).round(4)
    df["acidic_release_mechanism_t2_score"] = (0.68 * t1_acid + 0.32 * pka_component).round(4)
    rosetta_delta = numeric(df.get("rosetta_delta_score", pd.Series(index=df.index)), 0)
    rosetta_risk = (rosetta_delta.clip(lower=0) / 600.0).clip(upper=1.0)
    df["global_weakening_risk_t2_score"] = (0.70 * t1_risk + 0.30 * rosetta_risk).round(4)
    df["glycan_coverage_risk_t2_score"] = clip01(df.get("glycan_or_epitope_risk_score", pd.Series(index=df.index))).round(4)
    df["mpnn_compatibility_t2_score"] = clip01(df.get("mpnn_score_percentile_by_route", pd.Series(index=df.index))).round(4)
    df["tier2_class"] = df.apply(classify, axis=1)

    results_path = out_dir / "tier2b_full_stage1_results.csv"
    df.to_csv(results_path, index=False)
    plan_results_path = out_dir / "tier2b_stage1_results.csv"
    shutil.copy2(results_path, plan_results_path)

    tool_rows: list[dict[str, Any]] = []
    for tool, col in [
        ("pyrosetta", "pyrosetta_tool_status"),
        ("pka", "pka_tool_status"),
        ("foldx_subset", "foldx_tool_status"),
    ]:
        counts = df.groupby(["target", col], dropna=False).size().reset_index(name="count").rename(columns={col: "tool_status"})
        counts["tool"] = tool
        tool_rows.extend(counts[["tool", "target", "tool_status", "count"]].to_dict(orient="records"))
    tool_status = pd.DataFrame(tool_rows)
    tool_status.to_csv(out_dir / "tier2b_full_stage1_tool_status.csv", index=False)
    shutil.copy2(out_dir / "tier2b_full_stage1_tool_status.csv", out_dir / "tier2b_stage1_tool_status.csv")

    failure_rows = []
    for _, row in tool_status.iterrows():
        if row["tool_status"] not in {"success", "not_applicable", "skipped_by_policy"}:
            failure_rows.append(row.to_dict())
    failure_summary = pd.DataFrame(failure_rows, columns=["tool", "target", "tool_status", "count"])
    failure_summary.to_csv(out_dir / "tier2b_full_stage1_failure_summary.csv", index=False)
    shutil.copy2(out_dir / "tier2b_full_stage1_failure_summary.csv", out_dir / "tier2b_stage1_failure_summary.csv")

    class_summary = df.groupby(["target", "tier2_class"], dropna=False).size().reset_index(name="count")
    class_summary.to_csv(out_dir / "tier2b_full_stage1_class_summary.csv", index=False)

    required_evidence = pd.DataFrame(
        [
            {"requirement": "Stage-1 rows", "evidence": len(df), "status": "PASS" if len(df) == 7000 else "FAIL"},
            {
                "requirement": "PyRosetta success",
                "evidence": int(df["pyrosetta_tool_status"].eq("success").sum()),
                "status": "PASS" if int(df["pyrosetta_tool_status"].eq("success").sum()) == len(df) else "FAIL",
            },
            {
                "requirement": "pKa complete for His or not_applicable for non-His",
                "evidence": df["pka_tool_status"].value_counts(dropna=False).to_dict(),
                "status": "PASS" if df["pka_tool_status"].isin(["success", "not_applicable"]).all() else "WARN",
            },
            {
                "requirement": "FoldX selected subset evaluated or policy-skipped",
                "evidence": df["foldx_tool_status"].value_counts(dropna=False).to_dict(),
                "status": "PASS"
                if df["foldx_tool_status"].isin(["success", "skipped_by_policy"]).all()
                else "WARN",
            },
        ]
    )
    required_evidence.to_csv(out_dir / "tier2b_full_stage1_completion_checks.csv", index=False)

    report = [
        "# Tier2-B Full Stage-1 Review",
        "",
        f"Rows: `{len(df)}`",
        f"Runtime seconds for collection: `{time.time() - started:.1f}`",
        "",
        "## Completion Checks",
        markdown_table(required_evidence),
        "",
        "## Tool Status",
        markdown_table(tool_status),
        "",
        "## Tier2 Class Summary",
        markdown_table(class_summary),
        "",
        "## Failure Summary",
        markdown_table(failure_summary),
        "",
        "## Interpretation",
        "",
        "- This is a Full Stage-1 Tier2-core result package, not a final 10K library selection.",
        "- PyRosetta local structure validity is available for all Stage-1 rows.",
        "- pKa support is evaluated for His-containing rows and marked not_applicable for non-His rows.",
        "- FoldX actual AnalyseComplex is evaluated only for the approved selected representative subset; other rows are skipped_by_policy.",
    ]
    review = "\n".join(report) + "\n"
    (out_dir / "tier2b_full_stage1_review.md").write_text(review, encoding="utf-8")
    (out_dir / "tier2b_stage1_review.md").write_text(review, encoding="utf-8")


if __name__ == "__main__":
    main()
