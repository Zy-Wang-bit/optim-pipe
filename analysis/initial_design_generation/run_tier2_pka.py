#!/usr/bin/env python3
"""Run Tier2 pKa support on PyRosetta-generated structures."""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import multiprocessing as mp
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.pka.run_pka import predict_pka


MUT_RE = re.compile(r"^([A-Za-z])([A-Za-z])(\d+)([A-Za-z])$")


def parse_mutations(value: object) -> list[tuple[str, str, int, str]]:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return []
    out = []
    for token in str(value).replace(",", ";").split(";"):
        token = token.strip()
        if not token or token.lower() == "nan":
            continue
        m = MUT_RE.match(token)
        if not m:
            raise ValueError(f"cannot_parse_mutation:{token}")
        chain, old, pos, new = m.groups()
        out.append((chain, old, int(pos), new))
    return out


def chunks(items: list[dict[str, Any]], n: int) -> list[list[dict[str, Any]]]:
    if not items:
        return []
    size = math.ceil(len(items) / max(1, n))
    return [items[i : i + size] for i in range(0, len(items), size)]


def summarize_pka(df: pd.DataFrame) -> dict[str, Any]:
    if df.empty:
        return {}
    propka = pd.to_numeric(df.get("pKa_propka"), errors="coerce")
    pkai = pd.to_numeric(df.get("pKa_pkai"), errors="coerce")
    combined = pd.concat([propka, pkai], ignore_index=True).dropna()
    if combined.empty:
        score = pd.NA
    else:
        # Values near/above pH 6 are useful for His protonation hypotheses.
        score = max(0.0, min(1.0, (float(combined.mean()) - 5.3) / 1.7))
    return {
        "pka_his_count_observed": int(len(df)),
        "pka_propka_mean": float(propka.mean()) if propka.notna().any() else pd.NA,
        "pka_pkai_mean": float(pkai.mean()) if pkai.notna().any() else pd.NA,
        "pka_combined_mean": float(combined.mean()) if not combined.empty else pd.NA,
        "his_pka_support_t2_score": score,
    }


def worker(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summary_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    for row in rows:
        started = time.time()
        variant_id = str(row["variant_id"])
        target = str(row["target"])
        pdb_path = str(row.get("pdb_path", ""))
        base = {
            "variant_id": variant_id,
            "target": target,
            "pka_tool_status": "not_run",
            "pka_feature_missing_reason": "",
            "pka_feature_confidence": "not_run",
            "pka_runtime_seconds": pd.NA,
            "pka_his_count_requested": 0,
            "pka_his_count_observed": 0,
            "pka_propka_mean": pd.NA,
            "pka_pkai_mean": pd.NA,
            "pka_combined_mean": pd.NA,
            "his_pka_support_t2_score": pd.NA,
        }
        try:
            muts = parse_mutations(row.get("normalized_mutation_list", row.get("mutation_list", "")))
            his_filter = [(chain, pos) for chain, _, pos, new in muts if new == "H"]
            base["pka_his_count_requested"] = len(his_filter)
            if not his_filter:
                base.update(
                    {
                        "pka_tool_status": "not_applicable",
                        "pka_feature_missing_reason": "non_his_candidate",
                        "pka_feature_confidence": "not_applicable",
                        "pka_runtime_seconds": round(time.time() - started, 3),
                    }
                )
                summary_rows.append(base)
                continue
            if not pdb_path or not Path(pdb_path).exists():
                raise FileNotFoundError(f"pdb_missing:{pdb_path}")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                detail = predict_pka(pdb_path, his_filter=his_filter)
            if detail.empty:
                base.update(
                    {
                        "pka_tool_status": "failed_runtime",
                        "pka_feature_missing_reason": "pka_no_his_results",
                        "pka_feature_confidence": "low",
                        "pka_runtime_seconds": round(time.time() - started, 3),
                    }
                )
                summary_rows.append(base)
                continue
            detail["variant_id"] = variant_id
            detail["target"] = target
            detail_rows.extend(detail.to_dict(orient="records"))
            base.update(summarize_pka(detail))
            base.update(
                {
                    "pka_tool_status": "success",
                    "pka_feature_confidence": "propka_pkai" if detail["pKa_pkai"].notna().any() else "propka_only",
                    "pka_runtime_seconds": round(time.time() - started, 3),
                }
            )
            summary_rows.append(base)
        except Exception as exc:
            base.update(
                {
                    "pka_tool_status": "failed_runtime",
                    "pka_feature_missing_reason": str(exc),
                    "pka_feature_confidence": "low",
                    "pka_runtime_seconds": round(time.time() - started, 3),
                }
            )
            summary_rows.append(base)
    return summary_rows, detail_rows


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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-csv", required=True)
    parser.add_argument("--pyrosetta-results", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--stage", choices=["smoke", "calibration", "full"], required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(args.candidate_csv)
    pyro = pd.read_csv(args.pyrosetta_results)
    df = candidates.merge(
        pyro[["variant_id", "tool_status", "pdb_path"]],
        on="variant_id",
        how="left",
        validate="one_to_one",
        suffixes=("_candidate", ""),
    )
    if args.limit:
        df = df.head(args.limit).copy()
    rows = df.to_dict(orient="records")
    workers = max(1, min(args.workers, len(rows)))
    started = time.time()
    if workers == 1:
        nested = [worker(rows)]
    else:
        with mp.Pool(processes=workers) as pool:
            nested = pool.map(worker, chunks(rows, workers))
    summary_rows = [r for chunk, _ in nested for r in chunk]
    detail_rows = [r for _, chunk in nested for r in chunk]
    summary = pd.DataFrame(summary_rows)
    detail = pd.DataFrame(detail_rows)
    summary.to_csv(out_dir / f"tier2b_{args.stage}_pka_summary.csv", index=False)
    detail.to_csv(out_dir / f"tier2b_{args.stage}_pka_detail.csv", index=False)
    status = summary.groupby(["target", "pka_tool_status"], dropna=False).size().reset_index(name="count")
    status.to_csv(out_dir / f"tier2b_{args.stage}_pka_status_summary.csv", index=False)
    report = [
        f"# Tier2-B {args.stage.title()} pKa Report",
        "",
        f"Rows: `{len(rows)}`",
        f"Workers: `{workers}`",
        f"Runtime seconds: `{time.time() - started:.1f}`",
        "",
        "## Status Counts",
        markdown_table(status),
    ]
    (out_dir / f"tier2b_{args.stage}_pka_review.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    his_fail = summary[summary["pka_his_count_requested"].fillna(0).astype(int).gt(0) & ~summary["pka_tool_status"].eq("success")]
    if len(his_fail) == len(summary[summary["pka_his_count_requested"].fillna(0).astype(int).gt(0)]) and len(his_fail) > 0:
        raise SystemExit("All His-containing candidates failed pKa")


if __name__ == "__main__":
    main()
