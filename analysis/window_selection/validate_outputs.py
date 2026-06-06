#!/usr/bin/env python
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.window_selection.common import (
    load_capacity_config,
    load_decision_rules,
    load_inputs,
    load_source_policy,
    load_yaml,
    read_fasta,
    repo_path,
)


@dataclass
class Check:
    category: str
    name: str
    status: str
    detail: str


def _tables_dir(output_root: str | Path) -> Path:
    return repo_path(output_root) / "tables"


def _truthy(value: Any) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _numeric(value: Any, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _read_csv(path: Path) -> tuple[pd.DataFrame, str | None]:
    if not path.exists():
        return pd.DataFrame(), "missing"
    try:
        return pd.read_csv(path, dtype=str, keep_default_na=False), None
    except Exception as exc:  # noqa: BLE001 - reported verbatim.
        return pd.DataFrame(), f"unparseable:{exc}"


def _status(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def schema_checks(output_root: str | Path) -> list[Check]:
    schema = load_yaml("analysis/window_selection/schemas/output_schemas.yaml")
    required = schema.get("required_columns", {})
    td = _tables_dir(output_root)
    checks: list[Check] = []
    for filename, columns in required.items():
        df, error = _read_csv(td / filename)
        if error:
            checks.append(Check("schema", filename, "FAIL", error))
            continue
        missing = [column for column in columns if column not in df.columns]
        status = _status(not missing)
        detail = "all required columns present" if not missing else f"missing columns: {', '.join(missing)}"
        checks.append(Check("schema", filename, status, detail))
    return checks


def fasta_length_checks() -> list[Check]:
    cfg = load_inputs()
    expected = {
        "1E62:heavy": ("experiments/1E62/data/heavy.fasta", None, 115),
        "1E62:light": ("experiments/1E62/data/light.fasta", None, 113),
        "1E62:AeS": ("experiments/1E62/data/antigen_genotypes.fasta", "AeS", 226),
        "sdAb:sdab": ("experiments/sdab/R2/data/sdab.fasta", None, 122),
        "sdAb:AeS": ("experiments/sdab/R2/data/antigen.fasta", None, 226),
    }
    checks: list[Check] = []
    configured_files = str(cfg)
    for name, (path, record_id, length) in expected.items():
        try:
            _, _, seq = read_fasta(path, record_id)
            detail = f"length={len(seq)}, expected={length}"
            checks.append(Check("fasta_length", name, _status(len(seq) == length), detail))
        except Exception as exc:  # noqa: BLE001
            checks.append(Check("fasta_length", name, "FAIL", str(exc)))
        if path not in configured_files:
            checks.append(Check("fasta_config", name, "WARN", f"{path} not referenced in inputs.yaml text"))
    return checks


def wet_source_count_checks() -> list[Check]:
    expected = {
        "1E62 R1 KD": ("experiments/1E62/R1/1e62_R1_kd.csv", 52),
        "1E62 R2 ELISA": ("experiments/1E62/R2/wet_lab/elisa_summary.csv", 20),
        "1E62 R2 expression": ("experiments/1E62/R2/wet_lab/expression.csv", 20),
        "sdAb R2 single": ("experiments/sdab/R2/data/sdab_variants.csv", 31),
        "sdAb R2 raw ELISA": ("experiments/sdab/R2/data/sdab_elisa_raw.csv", 576),
        "sdAb R2 combo": ("experiments/sdab/R2/data/hs32-92_8ng_ml_Elisa_results.csv", 61),
        "sdAb R4 training": ("experiments/sdab/R4/data/training_data.csv", 93),
    }
    checks: list[Check] = []
    for name, (path_text, expected_rows) in expected.items():
        path = repo_path(path_text)
        df, error = _read_csv(path)
        if error:
            checks.append(Check("wet_source_rows", name, "FAIL", f"{error}:{path_text}"))
            continue
        detail = f"rows={len(df)}, expected={expected_rows}"
        checks.append(Check("wet_source_rows", name, _status(len(df) == expected_rows), detail))
    return checks


def source_policy_checks(output_root: str | Path) -> list[Check]:
    policy = load_source_policy()
    td = _tables_dir(output_root)
    checks: list[Check] = []
    excluded = [str(pattern).lower() for pattern in policy.get("excluded_patterns", [])]

    for target, files in policy.get("scoring_allowed_files", {}).items():
        for file_name in files:
            lower = str(file_name).lower()
            if any(pattern in lower for pattern in excluded):
                checks.append(Check("source_policy", f"{target}:{file_name}", "FAIL", "scoring file matches excluded pattern"))
            elif not repo_path(file_name).exists():
                checks.append(Check("source_policy", f"{target}:{file_name}", "FAIL", "scoring file missing"))
            else:
                checks.append(Check("source_policy", f"{target}:{file_name}", "PASS", "allowed scoring source exists"))

    for filename in ["wet_observation_table_1e62.csv", "wet_observation_table_sdab.csv", "prior_constraints_table.csv"]:
        df, error = _read_csv(td / filename)
        if error:
            checks.append(Check("source_policy", filename, "WARN", f"cannot inspect scoring sources: {error}"))
            continue
        source_cols = [col for col in ["source_file", "evidence_source"] if col in df.columns]
        violations: list[str] = []
        for col in source_cols:
            for value in df[col].dropna().astype(str).unique():
                lower = value.lower()
                if any(pattern in lower for pattern in excluded):
                    violations.append(f"{col}={value}")
        checks.append(
            Check(
                "source_policy",
                filename,
                _status(not violations),
                "no excluded sources in scoring inputs" if not violations else "; ".join(violations[:10]),
            )
        )

    constraints, error = _read_csv(td / "prior_constraints_table.csv")
    if error:
        checks.append(Check("source_policy", "prior_constraints_non_AeS_hard_filters", "WARN", error))
    else:
        violations = []
        if not constraints.empty and {"allowed_usage", "source_antigen", "applies_to_AeS_window_selection"}.issubset(constraints.columns):
            hard = constraints[constraints["allowed_usage"].astype(str).eq("hard_filter")]
            for _, row in hard.iterrows():
                antigen = str(row.get("source_antigen", "")).strip()
                if antigen and antigen != "AeS" and not _truthy(row.get("applies_to_AeS_window_selection")):
                    violations.append(str(row.get("constraint_id", "")))
        checks.append(
            Check(
                "source_policy",
                "prior_constraints_non_AeS_hard_filters",
                _status(not violations),
                "all hard filters are AeS-applicable" if not violations else ", ".join(violations),
            )
        )
    return checks


def context_only_checks(output_root: str | Path) -> list[Check]:
    td = _tables_dir(output_root)
    checks: list[Check] = []
    constraints, constraints_error = _read_csv(td / "prior_constraints_table.csv")
    scores, scores_error = _read_csv(td / "window_scores.csv")

    context_ids: list[str] = []
    if constraints_error:
        checks.append(Check("context_only", "prior_constraints_table", "WARN", constraints_error))
    elif "allowed_usage" in constraints.columns:
        context = constraints[constraints["allowed_usage"].astype(str).eq("context_only")]
        context_ids = context.get("constraint_id", pd.Series(dtype=str)).dropna().astype(str).tolist()
        illegal = context[
            context.get("constraint_type", pd.Series(dtype=str)).astype(str).eq("hard_protect")
            | context.get("applies_to_AeS_window_selection", pd.Series(dtype=bool)).map(_truthy)
        ]
        checks.append(
            Check(
                "context_only",
                "context_constraints_are_non_scoring",
                _status(illegal.empty),
                "context_only constraints do not request scoring use" if illegal.empty else f"illegal rows={len(illegal)}",
            )
        )

    if scores_error:
        checks.append(Check("context_only", "window_scores", "WARN", f"cannot inspect score influence: {scores_error}"))
    else:
        scoring_cols = [
            col
            for col in [
                "grade",
                "recommendation_status",
                "decision_path",
                "triggered_rules",
                "failed_rules",
                "hard_filter_reasons",
                "soft_downgrade_rules",
            ]
            if col in scores.columns
        ]
        combined = scores[scoring_cols].fillna("").astype(str).agg("|".join, axis=1) if scoring_cols else pd.Series(dtype=str)
        leaked_ids = [cid for cid in context_ids if cid and combined.str.contains(re.escape(cid), regex=True).any()]
        legacy_terms = combined.str.lower().str.contains("context_only|legacy_md|md_metrics|mechanism_context_only").any()
        ok = not leaked_ids and not legacy_terms
        detail = "context_only/legacy terms absent from scoring fields"
        if not ok:
            detail = f"context ids leaked={leaked_ids}; legacy_terms={legacy_terms}"
        checks.append(Check("context_only", "score_fields_ignore_context_only", _status(ok), detail))

    for filename in ["legacy_md_inventory_1e62.csv", "legacy_md_inventory_sdab.csv"]:
        df, error = _read_csv(td / filename)
        if error:
            checks.append(Check("context_only", filename, "WARN", error))
            continue
        allowed = set(df.get("allowed_usage", pd.Series(dtype=str)).dropna().astype(str))
        ok = allowed.issubset({"context_only", "mechanism_context_only"})
        checks.append(Check("context_only", filename, _status(ok), f"allowed_usage_values={sorted(allowed)}"))
    return checks


def blocker_dry_run_checks(output_root: str | Path) -> list[Check]:
    cfg = load_inputs()
    rules = load_decision_rules()
    capacity_cfg = load_capacity_config()
    td = _tables_dir(output_root)
    targets = [target["name"] for target in cfg.get("targets", [])]
    threshold = int(rules.get("thresholds", {}).get("min_qualified_structures_for_main_scoring", 15))
    min_capacity = int(capacity_cfg.get("min_usable_noncontrol_variant_count", 10000))
    checks: list[Check] = []

    af3, af3_error = _read_csv(td / "af3_model_manifest.csv")
    scores, scores_error = _read_csv(td / "window_scores.csv")
    candidates, candidates_error = _read_csv(td / "candidate_windows.csv")
    capacity, capacity_error = _read_csv(td / "window_design_capacity_table.csv")

    for target in targets:
        if af3_error:
            qualified = 0
            reason = f"af3_model_manifest {af3_error}"
        else:
            target_af3 = af3[af3.get("target", pd.Series(dtype=str)).astype(str).eq(target)]
            quality = target_af3.get("quality_status", pd.Series(dtype=str)).astype(str).str.lower()
            blocker = target_af3.get("blocker_flag", pd.Series([False] * len(target_af3))).map(_truthy)
            qualified = int((quality.eq("pass") & ~blocker).sum())
            reason = f"qualified_AF3_models={qualified}, threshold={threshold}"
        expected_af3_blocker = qualified < threshold
        if scores_error:
            checks.append(Check("blocker_dry_run", f"{target}:insufficient_AF3_models", "WARN", f"{reason}; window_scores {scores_error}"))
        else:
            target_scores = scores[scores.get("target", pd.Series(dtype=str)).astype(str).eq(target)]
            has_blocker = (
                not target_scores.empty
                and target_scores.get("recommendation_status", pd.Series(dtype=str)).astype(str).eq("blocked").any()
                and target_scores.get("blocker_type", pd.Series(dtype=str)).astype(str).eq("insufficient_AF3_models").any()
            )
            ok = has_blocker if expected_af3_blocker else not has_blocker
            checks.append(Check("blocker_dry_run", f"{target}:insufficient_AF3_models", _status(ok), reason))

        if candidates_error or capacity_error:
            checks.append(Check("blocker_dry_run", f"{target}:capacity_failure", "WARN", f"candidate/capacity unavailable: {candidates_error or capacity_error}"))
            continue
        merged = candidates.merge(capacity, on="window_id", how="left", suffixes=("", "_capacity"))
        group = merged[merged.get("target", pd.Series(dtype=str)).astype(str).eq(target)]
        legal = group[
            group.get("length", pd.Series(dtype=float)).map(_numeric).eq(40)
            & ~group.get("hard_excluded", pd.Series([False] * len(group))).map(_truthy)
            & ~group.get("crosses_linker", pd.Series([False] * len(group))).map(_truthy)
        ]
        expected_capacity_blocker = not legal.empty and not legal.get("usable_noncontrol_variant_count", pd.Series(dtype=float)).map(_numeric).ge(min_capacity).any()
        if scores_error:
            checks.append(Check("blocker_dry_run", f"{target}:capacity_failure", "WARN", f"expected={expected_capacity_blocker}; window_scores {scores_error}"))
        else:
            target_scores = scores[scores.get("target", pd.Series(dtype=str)).astype(str).eq(target)]
            has_capacity_blocker = (
                not target_scores.empty
                and target_scores.get("recommendation_status", pd.Series(dtype=str)).astype(str).eq("blocked").any()
                and target_scores.get("blocker_type", pd.Series(dtype=str)).astype(str).eq("capacity_failure").any()
            )
            ok = has_capacity_blocker if expected_capacity_blocker else True
            detail = f"legal_windows={len(legal)}, min_capacity={min_capacity}, expected_target_blocker={expected_capacity_blocker}"
            checks.append(Check("blocker_dry_run", f"{target}:capacity_failure", _status(ok), detail))

    if not scores_error and "wet_prior_summary" in scores.columns and "blocker_type" in scores.columns:
        no_prior = scores["wet_prior_summary"].fillna("").astype(str).str.contains("unknown_or_neutral")
        bad = scores[no_prior & scores["blocker_type"].fillna("").astype(str).eq("missing_wet_data")]
        checks.append(
            Check(
                "blocker_dry_run",
                "window_no_prior_is_not_missing_wet_data",
                _status(bad.empty),
                "no window-level unknown prior triggers missing_wet_data" if bad.empty else f"bad_rows={len(bad)}",
            )
        )
    return checks


def render_report(checks: list[Check]) -> str:
    status_order = {"FAIL": 0, "WARN": 1, "PASS": 2}
    counts = {status: sum(1 for check in checks if check.status == status) for status in ["PASS", "WARN", "FAIL"]}
    lines = [
        "# Validation Report",
        "",
        f"- PASS: {counts['PASS']}",
        f"- WARN: {counts['WARN']}",
        f"- FAIL: {counts['FAIL']}",
        "",
        "| Category | Check | Status | Detail |",
        "|---|---|---:|---|",
    ]
    for check in sorted(checks, key=lambda item: (item.category, status_order.get(item.status, 9), item.name)):
        detail = str(check.detail).replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {check.category} | {check.name} | {check.status} | {detail} |")
    lines.append("")
    return "\n".join(lines)


def run_validation(output_root: str | Path) -> list[Check]:
    checks: list[Check] = []
    checks.extend(schema_checks(output_root))
    checks.extend(fasta_length_checks())
    checks.extend(wet_source_count_checks())
    checks.extend(source_policy_checks(output_root))
    checks.extend(context_only_checks(output_root))
    checks.extend(blocker_dry_run_checks(output_root))
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate pH-sensitive 40 aa window-selection outputs.")
    parser.add_argument("--output-root", default=load_inputs()["output_root"], help="Output root containing tables/.")
    parser.add_argument("--strict", action="store_true", help="Exit non-zero if any FAIL checks are present.")
    args = parser.parse_args()

    checks = run_validation(args.output_root)
    report = render_report(checks)
    out_path = repo_path(args.output_root) / "validation_report.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report)
    fail_count = sum(1 for check in checks if check.status == "FAIL")
    warn_count = sum(1 for check in checks if check.status == "WARN")
    print(f"Wrote {out_path} with {fail_count} FAIL and {warn_count} WARN checks.")
    if args.strict and fail_count:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
