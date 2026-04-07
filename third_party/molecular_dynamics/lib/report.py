# third_party/molecular_dynamics/lib/report.py
"""汇总所有变体的 MD 分析结果为 md_report.csv。"""

import json
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def generate_report(output_dir: Path, ph_values: list[float] | None = None) -> Path:
    """遍历 output_dir 下所有变体目录，汇总分析结果。"""
    output_dir = Path(output_dir)
    rows = []

    for variant_dir in sorted(output_dir.iterdir()):
        if not variant_dir.is_dir():
            continue
        variant_id = variant_dir.name

        ph_dirs = sorted(variant_dir.glob("pH_*"))
        if not ph_dirs:
            summary_path = variant_dir / "analysis" / "summary.json"
            if summary_path.exists():
                row = _load_summary_row(summary_path, variant_id, ph=None)
                rows.append(row)
            continue

        for ph_dir in ph_dirs:
            ph_str = ph_dir.name.replace("pH_", "")
            try:
                ph = float(ph_str)
            except ValueError:
                continue
            if ph_values and ph not in ph_values:
                continue

            summary_path = ph_dir / "analysis" / "summary.json"
            if summary_path.exists():
                row = _load_summary_row(summary_path, variant_id, ph=ph)
                rows.append(row)

        comparison_path = variant_dir / "ph_comparison.json"
        if comparison_path.exists():
            with open(comparison_path) as f:
                comp = json.load(f)
            deltas = comp.get("deltas", {})
            for row in rows:
                if row["variant_id"] == variant_id and row.get("ph") == comp.get("target_ph"):
                    row.update(deltas)

    if not rows:
        logger.warning("No analysis results found in %s", output_dir)
        return output_dir / "md_report.csv"

    df = pd.DataFrame(rows)
    out = output_dir / "md_report.csv"
    df.to_csv(out, index=False, float_format="%.4f")
    logger.info("Report saved: %s (%d variants)", out, len(df))
    return out


def _load_summary_row(summary_path: Path, variant_id: str, ph: float | None) -> dict:
    with open(summary_path) as f:
        summary = json.load(f)
    row = {"variant_id": variant_id, "ph": ph}
    row.update(summary)
    return row
