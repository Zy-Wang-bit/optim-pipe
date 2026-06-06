#!/usr/bin/env python3
"""Validate P0 ProteinMPNN integration artifacts.

This validator checks the engineering boundary before true ProteinMPNN compute:
scaffold outputs, AF3-derived runner inputs, score-only shards, and result
collection tables. By default pending ProteinMPNN results are allowed. Use
``--require-mpnn-results`` when the real score/generation jobs should already
have completed.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "analysis/initial_design_generation/schemas/audit_schema.yaml"
P0 = ROOT / "results/initial_design_generation/p0_mpnn"
RUNNER = ROOT / "results/initial_design_generation/p0_mpnn_runner_inputs"
BACKBONES = ROOT / "results/initial_design_generation/p0_mpnn_backbones"
PRODUCTION = ROOT / "results/initial_design_generation/production_initial_pool"


@dataclass
class Check:
    status: str
    name: str
    detail: str


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def count_jsonl(path: Path) -> int:
    count = 0
    with path.open() as fh:
        for line in fh:
            if line.strip():
                json.loads(line)
                count += 1
    return count


def add_file_checks(checks: list[Check], base: Path, names: Iterable[str], group: str) -> None:
    missing = []
    for name in names:
        path = base / name
        if not path.exists():
            missing.append(name)
    if missing:
        checks.append(Check("FAIL", f"{group} required files", f"missing: {', '.join(missing)}"))
    else:
        checks.append(Check("PASS", f"{group} required files", f"{len(list(names))} files present under {rel(base)}"))


def validate(require_mpnn_results: bool) -> list[Check]:
    checks: list[Check] = []
    schema = yaml.safe_load(SCHEMA.read_text())

    add_file_checks(checks, P0, schema["p0_mpnn_required_outputs"], "P0")
    add_file_checks(checks, RUNNER, schema["p0_mpnn_runner_required_outputs"], "runner input")
    add_file_checks(checks, BACKBONES, schema["p0_mpnn_backbone_required_outputs"], "backbone")

    production_pool = pd.read_csv(PRODUCTION / "production_initial_pool_candidates_all.csv")
    score = pd.read_csv(P0 / "mpnn_scores_current_pool.csv")
    expected_rows = len(production_pool)
    if len(score) == expected_rows:
        checks.append(Check("PASS", "score table row count", f"{len(score)} rows match production pool"))
    else:
        checks.append(Check("FAIL", "score table row count", f"{len(score)} rows, expected {expected_rows}"))

    score_status_counts = score["mpnn_score_status"].value_counts(dropna=False).to_dict()
    scored_rows = int((score["mpnn_score_status"] == "scored_by_mpnn").sum())
    pending_rows = int((score["mpnn_score_status"] == "pending_mpnn_run").sum())
    if require_mpnn_results and scored_rows != len(score):
        checks.append(Check("FAIL", "ProteinMPNN score results", f"required complete scoring but only {scored_rows}/{len(score)} rows are scored"))
    elif scored_rows == 0 and pending_rows == len(score):
        checks.append(Check("PASS", "ProteinMPNN score results", "all rows explicitly pending; no fabricated scores"))
    elif scored_rows > 0:
        checks.append(Check("PASS", "ProteinMPNN score results", f"{scored_rows} scored rows; statuses={score_status_counts}"))
    else:
        checks.append(Check("WARN", "ProteinMPNN score results", f"mixed non-scored statuses={score_status_counts}"))

    score_audit = pd.read_csv(P0 / "mpnn_score_only_collection_audit.csv")
    manifest = pd.read_csv(RUNNER / "mpnn_score_only_input_manifest.csv")
    manifest_candidates = int(manifest["candidate_count"].sum())
    audit_candidates = int(score_audit["candidate_count"].sum())
    if manifest_candidates == expected_rows == audit_candidates:
        checks.append(Check("PASS", "score-only manifest coverage", f"{manifest_candidates} candidates across {len(manifest)} shards"))
    else:
        checks.append(
            Check(
                "FAIL",
                "score-only manifest coverage",
                f"manifest={manifest_candidates}, audit={audit_candidates}, expected={expected_rows}",
            )
        )

    fasta_records = sum(count_fasta_records(Path(path)) for path in manifest["fasta_path"])
    if fasta_records == expected_rows:
        checks.append(Check("PASS", "score-only FASTA records", f"{fasta_records} FASTA records"))
    else:
        checks.append(Check("FAIL", "score-only FASTA records", f"{fasta_records} records, expected {expected_rows}"))

    jobs = pd.read_csv(RUNNER / "mpnn_runner_jobs_constrained.csv")
    if bool(jobs["runner_ready"].all()) and len(jobs) > 0:
        checks.append(Check("PASS", "constrained runner jobs", f"{len(jobs)} runner-ready jobs"))
    else:
        checks.append(Check("FAIL", "constrained runner jobs", "runner_ready is not true for all jobs"))

    parsed_count = count_jsonl(RUNNER / "jsonl/parsed_pdbs.jsonl")
    if parsed_count == len(jobs):
        checks.append(Check("PASS", "parsed PDB JSONL", f"{parsed_count} entries match constrained jobs"))
    else:
        checks.append(Check("FAIL", "parsed PDB JSONL", f"{parsed_count} entries, expected {len(jobs)}"))

    command_plan = pd.read_csv(RUNNER / "mpnn_constrained_command_plan.csv")
    planned_samples = int(jobs["planned_raw_samples"].sum())
    command_planned_samples = int(command_plan["planned_raw_samples"].sum())
    command_job_count = int(command_plan["job_count"].sum())
    if command_planned_samples == planned_samples and command_job_count == len(jobs):
        checks.append(
            Check(
                "PASS",
                "constrained command scale",
                f"{len(command_plan)} commands encode {planned_samples} planned raw samples",
            )
        )
    else:
        checks.append(
            Check(
                "FAIL",
                "constrained command scale",
                f"command samples={command_planned_samples}, planned={planned_samples}; command jobs={command_job_count}, jobs={len(jobs)}",
            )
        )

    gen_audit = pd.read_csv(P0 / "constrained_mpnn_generation_collection_audit.csv")
    if len(gen_audit) == len(jobs):
        checks.append(Check("PASS", "generation audit coverage", f"{len(gen_audit)} jobs audited"))
    else:
        checks.append(Check("FAIL", "generation audit coverage", f"{len(gen_audit)} audited jobs, expected {len(jobs)}"))

    generated = int(gen_audit["generated_count"].sum()) if "generated_count" in gen_audit else 0
    gen_with_plan = gen_audit.merge(jobs[["job_name", "planned_raw_samples"]], on="job_name", how="left")
    incomplete_jobs = int((gen_with_plan["generated_count"] < gen_with_plan["planned_raw_samples"]).sum()) if not gen_with_plan.empty else len(jobs)
    if require_mpnn_results and (generated < planned_samples or incomplete_jobs):
        checks.append(
            Check(
                "FAIL",
                "constrained MPNN generation results",
                f"required complete generation but generated {generated}/{planned_samples}; incomplete jobs={incomplete_jobs}",
            )
        )
    elif generated == 0:
        checks.append(Check("PASS", "constrained MPNN generation results", "all jobs explicitly pending; no fabricated candidates"))
    else:
        passed = int(gen_audit.get("passed_hard_filter_count", pd.Series(dtype=int)).sum())
        checks.append(Check("PASS", "constrained MPNN generation results", f"{generated} generated, {passed} passed hard filter"))

    relaxed_jobs = pd.read_csv(RUNNER / "mpnn_runner_jobs_relaxed.csv")
    if bool(relaxed_jobs["runner_ready"].all()) and len(relaxed_jobs) > 0:
        checks.append(Check("PASS", "relaxed runner jobs", f"{len(relaxed_jobs)} runner-ready jobs"))
    else:
        checks.append(Check("FAIL", "relaxed runner jobs", "runner_ready is not true for all relaxed jobs"))

    relaxed_parsed_count = count_jsonl(RUNNER / "jsonl/relaxed/parsed_pdbs.jsonl")
    if relaxed_parsed_count == len(relaxed_jobs):
        checks.append(Check("PASS", "relaxed parsed PDB JSONL", f"{relaxed_parsed_count} entries match relaxed jobs"))
    else:
        checks.append(Check("FAIL", "relaxed parsed PDB JSONL", f"{relaxed_parsed_count} entries, expected {len(relaxed_jobs)}"))

    relaxed_command_plan = pd.read_csv(RUNNER / "mpnn_relaxed_command_plan.csv")
    planned_relaxed_samples = int(relaxed_jobs["planned_raw_samples"].sum())
    relaxed_command_planned_samples = int(relaxed_command_plan["planned_raw_samples"].sum())
    relaxed_command_job_count = int(relaxed_command_plan["job_count"].sum())
    if relaxed_command_planned_samples == planned_relaxed_samples and relaxed_command_job_count == len(relaxed_jobs):
        checks.append(
            Check(
                "PASS",
                "relaxed command scale",
                f"{len(relaxed_command_plan)} commands encode {planned_relaxed_samples} audit-only raw samples",
            )
        )
    else:
        checks.append(
            Check(
                "FAIL",
                "relaxed command scale",
                f"command samples={relaxed_command_planned_samples}, planned={planned_relaxed_samples}; command jobs={relaxed_command_job_count}, jobs={len(relaxed_jobs)}",
            )
        )

    relaxed_audit = pd.read_csv(P0 / "relaxed_mpnn_generation_collection_audit.csv")
    if len(relaxed_audit) == len(relaxed_jobs):
        checks.append(Check("PASS", "relaxed generation audit coverage", f"{len(relaxed_audit)} jobs audited"))
    else:
        checks.append(Check("FAIL", "relaxed generation audit coverage", f"{len(relaxed_audit)} audited jobs, expected {len(relaxed_jobs)}"))

    relaxed_generated = int(relaxed_audit["generated_count"].sum()) if "generated_count" in relaxed_audit else 0
    relaxed_with_plan = relaxed_audit.merge(relaxed_jobs[["job_name", "planned_raw_samples"]], on="job_name", how="left")
    relaxed_incomplete_jobs = int((relaxed_with_plan["generated_count"] < relaxed_with_plan["planned_raw_samples"]).sum()) if not relaxed_with_plan.empty else len(relaxed_jobs)
    if require_mpnn_results and (relaxed_generated < planned_relaxed_samples or relaxed_incomplete_jobs):
        checks.append(
            Check(
                "FAIL",
                "relaxed MPNN generation results",
                f"required complete relaxed generation but generated {relaxed_generated}/{planned_relaxed_samples}; incomplete jobs={relaxed_incomplete_jobs}",
            )
        )
    elif relaxed_generated == 0:
        checks.append(Check("PASS", "relaxed MPNN generation results", "all jobs explicitly pending; no fabricated audit candidates"))
    else:
        relaxed_passed = int(relaxed_audit.get("passed_hard_filter_count", pd.Series(dtype=int)).sum())
        checks.append(Check("PASS", "relaxed MPNN generation results", f"{relaxed_generated} generated, {relaxed_passed} passed hard filter"))

    backbone_audit = pd.read_csv(P0 / "mpnn_backbone_manifest_audit.csv")
    sequence_match_col = "sequence_match_reference" if "sequence_match_reference" in backbone_audit else "sequence_match"
    if len(backbone_audit) == 2 and bool(backbone_audit[sequence_match_col].astype(bool).all()):
        checks.append(Check("PASS", "AF3 backbone audit", "2 selected backbones with matching design-chain sequences"))
    else:
        checks.append(Check("FAIL", "AF3 backbone audit", "expected 2 sequence-matched selected backbones"))

    return checks


def write_report(checks: list[Check], output: Path) -> None:
    counts = pd.Series([c.status for c in checks]).value_counts().to_dict()
    lines = [
        "# P0 ProteinMPNN Integration Validation Report",
        "",
        f"PASS: {counts.get('PASS', 0)}",
        f"WARN: {counts.get('WARN', 0)}",
        f"FAIL: {counts.get('FAIL', 0)}",
        "",
        "| status | check | detail |",
        "|---|---|---|",
    ]
    for check in checks:
        detail = check.detail.replace("|", "\\|")
        lines.append(f"| {check.status} | {check.name} | {detail} |")
    output.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate P0 ProteinMPNN integration artifacts.")
    parser.add_argument("--require-mpnn-results", action="store_true")
    default_report = "results/initial_design_generation/p0_mpnn/p0_mpnn_validation_report.md"
    parser.add_argument(
        "--report",
        default=default_report,
    )
    args = parser.parse_args()

    checks = validate(args.require_mpnn_results)
    report_name = args.report
    if args.require_mpnn_results and args.report == default_report:
        report_name = "results/initial_design_generation/p0_mpnn/p0_mpnn_require_results_validation_report.md"
    report = ROOT / report_name
    write_report(checks, report)
    counts = pd.Series([c.status for c in checks]).value_counts().to_dict()
    print(
        {
            "status": "p0_mpnn_validation_complete",
            "report": rel(report),
            "pass": int(counts.get("PASS", 0)),
            "warn": int(counts.get("WARN", 0)),
            "fail": int(counts.get("FAIL", 0)),
        }
    )
    if counts.get("FAIL", 0):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
