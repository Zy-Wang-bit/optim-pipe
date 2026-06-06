#!/usr/bin/env python
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd

if __package__ is None or __package__ == "":
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from analysis.window_selection.common import (
    load_source_policy,
    make_arg_parser,
    repo_path,
    write_csv,
)


INVENTORY_COLUMNS = [
    "target",
    "batch",
    "variant_id",
    "pH_pair",
    "protocol_note",
    "available_metrics",
    "allowed_usage",
    "source_file",
]


def rel_path(path: Path) -> str:
    try:
        return path.relative_to(repo_path(".")).as_posix()
    except ValueError:
        return path.as_posix()


def read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        data = json.load(handle)
    return data if isinstance(data, dict) else {"json_type": type(data).__name__}


def metric_family_summary(keys: list[str]) -> str:
    families: dict[str, int] = {}
    for key in keys:
        family = key.split("_", 1)[0]
        families[family] = families.get(family, 0) + 1
    return "|".join(f"{name}:{count}" for name, count in sorted(families.items()))


def json_metric_summary(data: dict[str, Any]) -> str:
    parts = []
    for section in ["base_summary", "target_summary", "deltas"]:
        value = data.get(section)
        if isinstance(value, dict):
            parts.append(f"{section}({metric_family_summary(list(value))})")
    if parts:
        return ";".join(parts)
    metric_keys = [key for key, value in data.items() if isinstance(value, int | float | str)]
    return metric_family_summary(metric_keys) or "json_metadata"


def infer_batch(target: str, path: Path) -> str:
    parts = path.parts
    if target == "1E62":
        if "R2" in parts and "md_com8_vs_com18_ph_verification" in parts:
            return "R2_md_com8_vs_com18_ph_verification"
        if "R3" in parts and "md" in parts:
            return "R3_md_ph_comparison"
    if target == "sdAb":
        if "R2" in parts and "md" in parts:
            return "R2_md_ph_comparison"
        if "R3" in parts and "md" in parts:
            return "R3_md_ph_comparison"
        if "R4" in parts:
            return "R4_md_metrics"
    return path.parent.name


def infer_variant_from_path(path: Path) -> str:
    if path.name == "summary.json" and len(path.parts) >= 4:
        return path.parts[-4]
    if path.name == "ph_comparison.json":
        return path.parent.name
    if path.suffix.lower() == ".md":
        return path.stem
    return path.parent.name


def infer_ph_from_path(path: Path) -> str:
    matches = re.findall(r"pH[_-]?(\d+(?:\.\d+)?)", path.as_posix(), flags=re.IGNORECASE)
    if matches:
        return "|".join(dict.fromkeys(matches))
    return ""


def ph_pair_from_json(data: dict[str, Any], path: Path) -> str:
    base = data.get("base_ph")
    target = data.get("target_ph")
    if base is not None and target is not None:
        return f"{base}_vs_{target}"
    ph_from_path = infer_ph_from_path(path)
    return f"single_pH_{ph_from_path}" if ph_from_path else ""


def allowed_usage_for(path: Path) -> str:
    if path.suffix.lower() == ".md":
        return "context_only"
    return "mechanism_context_only"


def row_from_json(target: str, path: Path) -> dict[str, str]:
    data = read_json(path)
    variant_id = str(data.get("variant") or infer_variant_from_path(path))
    return {
        "target": target,
        "batch": infer_batch(target, path),
        "variant_id": variant_id,
        "pH_pair": ph_pair_from_json(data, path),
        "protocol_note": "legacy_md_json_or_summary_not_used_for_scoring",
        "available_metrics": json_metric_summary(data),
        "allowed_usage": allowed_usage_for(path),
        "source_file": rel_path(path),
    }


def row_from_markdown(target: str, path: Path) -> dict[str, str]:
    text = path.read_text(errors="ignore")
    phs = "|".join(dict.fromkeys(re.findall(r"pH\s*[_-]?(\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)))
    return {
        "target": target,
        "batch": infer_batch(target, path),
        "variant_id": infer_variant_from_path(path),
        "pH_pair": phs,
        "protocol_note": "legacy_md_report_not_used_for_scoring",
        "available_metrics": "report_markdown",
        "allowed_usage": "context_only",
        "source_file": rel_path(path),
    }


def rows_from_metrics_csv(target: str, path: Path) -> list[dict[str, str]]:
    df = pd.read_csv(path)
    metric_cols = [col for col in df.columns if col not in {"name", "variant", "variant_id", "status"}]
    rows = []
    for index, rec in df.iterrows():
        variant_id = rec.get("variant_id", rec.get("variant", rec.get("name", f"row_{index + 1}")))
        status = rec.get("status", "")
        rows.append(
            {
                "target": target,
                "batch": infer_batch(target, path),
                "variant_id": str(variant_id),
                "pH_pair": "7.4_vs_6.0" if any(str(col).startswith("delta_") for col in metric_cols) else "",
                "protocol_note": f"legacy_md_metrics_csv_status={status};not_used_for_scoring",
                "available_metrics": "|".join(metric_cols),
                "allowed_usage": "mechanism_context_only",
                "source_file": rel_path(path),
            }
        )
    return rows


def is_inventory_source(path: Path) -> bool:
    name = path.name.lower()
    if name == "ph_comparison.json":
        return True
    if name == "summary.json" and path.parent.name == "analysis":
        return True
    if path.suffix.lower() == ".md" and "report" in name:
        return True
    if path.suffix.lower() == ".csv" and ("metrics" in name or "report" in name):
        return True
    return False


def scan_context_path(target: str, configured_path: str) -> list[dict[str, str]]:
    path = repo_path(configured_path)
    if not path.exists():
        return [
            {
                "target": target,
                "batch": Path(configured_path).name,
                "variant_id": "",
                "pH_pair": "",
                "protocol_note": "configured_legacy_md_path_missing",
                "available_metrics": "",
                "allowed_usage": "context_only",
                "source_file": configured_path,
            }
        ]

    files = [path] if path.is_file() else sorted(p for p in path.rglob("*") if p.is_file() and is_inventory_source(p))
    rows: list[dict[str, str]] = []
    for file_path in files:
        suffix = file_path.suffix.lower()
        if suffix == ".json":
            rows.append(row_from_json(target, file_path))
        elif suffix == ".md":
            rows.append(row_from_markdown(target, file_path))
        elif suffix == ".csv":
            rows.extend(rows_from_metrics_csv(target, file_path))
    return rows


def build_inventory() -> tuple[pd.DataFrame, pd.DataFrame]:
    policy = load_source_policy()
    rows_by_target: dict[str, list[dict[str, str]]] = {"1E62": [], "sdAb": []}
    for target in rows_by_target:
        for configured_path in policy.get("context_only_files", {}).get(target, []):
            rows_by_target[target].extend(scan_context_path(target, configured_path))

    frames = []
    for target in ["1E62", "sdAb"]:
        df = pd.DataFrame(rows_by_target[target], columns=INVENTORY_COLUMNS)
        if not df.empty:
            df = df.sort_values(["batch", "variant_id", "source_file"]).reset_index(drop=True)
        frames.append(df)
    return frames[0], frames[1]


def main() -> None:
    make_arg_parser("Build legacy MD inventory tables for context-only evidence.").parse_args()
    one, sdab = build_inventory()
    write_csv(one, "legacy_md_inventory_1e62.csv")
    write_csv(sdab, "legacy_md_inventory_sdab.csv")
    print(f"Wrote legacy MD inventory rows: 1E62={len(one)}, sdAb={len(sdab)}.")


if __name__ == "__main__":
    main()
