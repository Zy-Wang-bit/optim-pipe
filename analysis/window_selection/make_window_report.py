#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.window_selection.common import load_inputs, repo_path


GRADE_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "not_scored": 4}


def _tables_dir(output_root: str | Path) -> Path:
    return repo_path(output_root) / "tables"


def _read_optional_table(output_root: str | Path, filename: str) -> pd.DataFrame:
    path = _tables_dir(output_root) / filename
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _value(row: pd.Series, key: str, default: str = "") -> str:
    value = row.get(key, default)
    if pd.isna(value):
        return default
    return str(value)


def _target_names() -> list[str]:
    return [target["name"] for target in load_inputs().get("targets", [])]


def _window_label(row: pd.Series) -> str:
    segment = _value(row, "segment")
    start = _value(row, "start")
    end = _value(row, "end")
    window_id = _value(row, "window_id")
    coords = f"{segment}:{start}-{end}" if start or end else segment
    return f"{window_id} ({coords})" if coords.strip(":") else window_id


def _sort_scores(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["_grade_order"] = out.get("grade", pd.Series(dtype=str)).map(lambda grade: GRADE_ORDER.get(str(grade), 9))
    out["_capacity"] = pd.to_numeric(out.get("usable_noncontrol_variant_count", pd.Series(dtype=float)), errors="coerce").fillna(0)
    out = out.sort_values(["_grade_order", "_capacity", "window_id"], ascending=[True, False, True])
    return out.drop(columns=["_grade_order", "_capacity"], errors="ignore")


def _render_window(row: pd.Series, prefix: str) -> list[str]:
    lines = [
        f"{prefix} `{_window_label(row)}`",
        f"  - grade/status: `{_value(row, 'grade')}` / `{_value(row, 'recommendation_status')}`",
        f"  - evidence route: `{_value(row, 'evidence_route')}`",
        f"  - pH mechanism hypothesis: `{_value(row, 'pH_mechanism_hypothesis_strength')}`",
        f"  - capacity: raw `{_value(row, 'raw_feasible_variant_count', '0')}`, usable non-control `{_value(row, 'usable_noncontrol_variant_count', '0')}`",
        f"  - wet prior: {_value(row, 'wet_prior_summary')}",
        f"  - structure: {_value(row, 'structure_summary')}",
        f"  - glycan: {_value(row, 'glycan_risk_summary')}",
        f"  - decision path: `{_value(row, 'decision_path')}`",
    ]
    soft = _value(row, "soft_downgrade_rules")
    hard = _value(row, "hard_filter_reasons")
    if soft:
        lines.append(f"  - soft downgrade rules: `{soft}`")
    if hard:
        lines.append(f"  - hard filter reasons: `{hard}`")
    return lines


def _blocked_section(target: str, scores: pd.DataFrame) -> list[str]:
    blocker_rows = scores[scores.get("recommendation_status", pd.Series(dtype=str)).astype(str).eq("blocked")]
    blocker_type = ""
    blocker_reason = ""
    if not blocker_rows.empty:
        blocker_type = _value(blocker_rows.iloc[0], "blocker_type")
        blocker_reason = _value(blocker_rows.iloc[0], "blocker_reason")
    remediation = {
        "insufficient_AF3_models": "run or parse enough qualified new AF3 parent-complex models before final ranking",
        "missing_wet_data": "restore or fix required target wet-lab input files",
        "schema_error": "generate the missing upstream tables and re-run validation",
        "capacity_failure": "relax library design constraints or choose a different legal window definition",
        "source_policy_violation": "remove excluded/context-only sources from scoring inputs",
        "numbering_failure": "fix antibody/reference numbering and regenerate candidate windows",
    }.get(blocker_type, "resolve the blocker and re-run scoring")
    return [
        f"## {target}",
        "",
        "`no_recommendation_due_to_blocker`",
        "",
        f"- blocker_type: `{blocker_type or 'unknown'}`",
        f"- blocker_reason: {blocker_reason or 'not recorded'}",
        f"- remediation: {remediation}",
        "",
    ]


def _target_section(target: str, scores: pd.DataFrame) -> list[str]:
    target_scores = scores[scores.get("target", pd.Series(dtype=str)).astype(str).eq(target)].copy()
    if target_scores.empty:
        return [
            f"## {target}",
            "",
            "`no_recommendation_due_to_blocker`",
            "",
            "- blocker_type: `schema_error`",
            "- blocker_reason: no rows found in window_scores.csv for this target",
            "- remediation: regenerate candidate windows and scoring table",
            "",
        ]

    if target_scores.get("recommendation_status", pd.Series(dtype=str)).astype(str).eq("blocked").any():
        return _blocked_section(target, target_scores)

    recommended = target_scores[target_scores.get("recommendation_status", pd.Series(dtype=str)).astype(str).eq("recommended")]
    recommended = _sort_scores(recommended)
    exploratory = _sort_scores(
        target_scores[target_scores.get("recommendation_status", pd.Series(dtype=str)).astype(str).ne("recommended")]
    )

    lines = [f"## {target}", ""]
    if recommended.empty:
        lines.extend(
            [
                "`no_primary_recommendation`",
                "",
                "- blocker_type: `insufficient_evidence`",
                "- blocker_reason: scoring completed but produced no A/B recommended windows",
                "- remediation: inspect exploratory windows, structure evidence, and capacity constraints",
                "",
            ]
        )
    else:
        primary = recommended.iloc[0]
        lines.extend(_render_window(primary, "Primary:"))
        lines.append("")
        backups = recommended.iloc[1:3]
        if backups.empty:
            lines.append("Backups: none scored as recommended.")
        else:
            lines.append("Backups:")
            for _, row in backups.iterrows():
                lines.extend(_render_window(row, "-"))
        lines.append("")

    if not exploratory.empty:
        lines.append("Not Recommended / Exploratory:")
        for _, row in exploratory.head(5).iterrows():
            lines.append(
                f"- `{_window_label(row)}`: grade `{_value(row, 'grade')}`, "
                f"decision `{_value(row, 'decision_path')}`, hard filters `{_value(row, 'hard_filter_reasons')}`"
            )
        lines.append("")
    return lines


def _glycan_mode_note(output_root: str | Path, scores: pd.DataFrame) -> str:
    guardrails = _read_optional_table(output_root, "glycan_guardrail_table.csv")
    modes: set[str] = set()
    if not guardrails.empty and "glycan_guardrail_mode" in guardrails.columns:
        modes.update(guardrails["glycan_guardrail_mode"].dropna().astype(str).unique())
    if "glycan_risk_summary" in scores.columns:
        for value in scores["glycan_risk_summary"].dropna().astype(str):
            if "mode=heuristic" in value:
                modes.add("heuristic")
            if "mode=explicit" in value:
                modes.add("explicit")
    if "heuristic" in modes:
        return (
            "Glycan note: at least one window uses `glycan_guardrail_mode=heuristic`; "
            "these risk labels come from N146 geometry/radius rules, not explicit glycan conformer modeling."
        )
    if "explicit" in modes:
        return "Glycan note: glycan guardrails use explicit N146 glycan conformer evidence where available."
    return "Glycan note: no glycan guardrail mode was available in the current scoring outputs."


def build_report(output_root: str | Path) -> str:
    scores = _read_optional_table(output_root, "window_scores.csv")
    lines: list[str] = [
        "# Final 40 aa Window Recommendation",
        "",
    ]
    if scores.empty:
        lines.extend(
            [
                "`no_recommendation_due_to_blocker`",
                "",
                "- blocker_type: `schema_error`",
                "- blocker_reason: `window_scores.csv` is missing or empty",
                "- remediation: run `analysis/window_selection/score_windows.py` after upstream tables are generated",
                "",
            ]
        )
        return "\n".join(lines)

    lines.append(_glycan_mode_note(output_root, scores))
    lines.append("")
    for target in _target_names():
        lines.extend(_target_section(target, scores))
    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate final 40 aa window recommendation report.")
    parser.add_argument("--output-root", default=load_inputs()["output_root"], help="Output root containing tables/.")
    args = parser.parse_args()

    report = build_report(args.output_root)
    out_path = repo_path(args.output_root) / "final_recommendation.md"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
