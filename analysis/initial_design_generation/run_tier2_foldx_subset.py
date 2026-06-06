#!/usr/bin/env python3
"""Run Tier2 selected FoldX AnalyseComplex on PyRosetta-generated structures."""

from __future__ import annotations

import argparse
import math
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FOLDX_BIN = ROOT / "third_party/foldx/foldx"
PH_POINTS = (7.4, 6.0)
CHAIN_GROUPS = {
    "Ab_1E62": "H,L;C",
    "Ab_sdAb": "A;B",
}


def chunks(items: list[dict[str, Any]], n: int) -> list[list[dict[str, Any]]]:
    if not items:
        return []
    size = math.ceil(len(items) / max(1, n))
    return [items[i : i + size] for i in range(0, len(items), size)]


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)[:180]


def run(cmd: list[str], cwd: Path, env: dict[str, str], timeout: int) -> subprocess.CompletedProcess[str]:
    env2 = dict(os.environ)
    env2.update(env)
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env2,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout,
        check=True,
    )


def find_interaction_file(work_dir: Path, pdb_noext: str, ph: float) -> Path | None:
    labels = []
    s = str(ph)
    labels.extend([s, f"{ph:.1f}", s.replace(".", "_"), str(int(round(ph))), s.replace(".", "")])
    labels = list(dict.fromkeys(labels))
    patterns: list[str] = []
    for label in labels:
        patterns.extend(
            [
                f"Interaction_AC_{pdb_noext}_{label}_AC.fxout",
                f"Interaction_AC_{pdb_noext}_{label}.fxout",
                f"Interaction_*{pdb_noext}*_{label}*_AC.fxout",
                f"Interaction_*{pdb_noext}*_{label}*.fxout",
            ]
        )
    numeric_label = re.sub(r"\D", "", str(ph))
    patterns.append(f"Interaction_*{pdb_noext}*{numeric_label}*.fxout")
    for pattern in patterns:
        hits = sorted(work_dir.glob(pattern))
        if hits:
            return hits[0]
    return None


def parse_interaction_energy(path: Path, pdb_base: str, groups: str) -> float | None:
    left, right = groups.split(";")
    left_set = set(x for x in left.split(",") if x)
    right_set = set(x for x in right.split(",") if x)
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    header_index = None
    for i, line in enumerate(lines):
        if "\t" in line and "Pdb" in line and "Group1" in line and "Group2" in line and "Interaction Energy" in line:
            header_index = i
            break
    if header_index is None:
        return None
    header = lines[header_index].split("\t")
    try:
        c_pdb = header.index("Pdb")
        c_g1 = header.index("Group1")
        c_g2 = header.index("Group2")
        c_ie = header.index("Interaction Energy")
    except ValueError:
        return None
    total = 0.0
    found = False
    for line in lines[header_index + 1 :]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) <= max(c_pdb, c_g1, c_g2, c_ie):
            continue
        if Path(parts[c_pdb].strip()).name.lstrip("./") != pdb_base:
            continue
        g1 = parts[c_g1].strip()
        g2 = parts[c_g2].strip()
        if not ((g1 in left_set and g2 in right_set) or (g2 in left_set and g1 in right_set)):
            continue
        try:
            total += float(parts[c_ie])
            found = True
        except ValueError:
            continue
    return total if found else None


def foldx_score_from_dg(dg_pH7_4: object) -> object:
    value = pd.to_numeric(pd.Series([dg_pH7_4]), errors="coerce").iloc[0]
    if pd.isna(value):
        return pd.NA
    # More negative interface energy is better. Keep this deliberately coarse:
    # it is a subset-level risk feature, not a final ranking truth.
    return max(0.0, min(1.0, (-float(value) - 2.0) / 18.0))


def worker(rows: list[dict[str, Any]], foldx_bin: str, work_root: str, timeout: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    env = {"OMP_NUM_THREADS": "1"}
    for row in rows:
        started = time.time()
        variant_id = str(row["variant_id"])
        target = str(row["target"])
        work_dir = Path(work_root) / safe_name(variant_id)
        result: dict[str, Any] = {
            "variant_id": variant_id,
            "target": target,
            "foldx_tool_status": "not_run",
            "foldx_feature_missing_reason": "",
            "foldx_feature_confidence": "foldx_analysecomplex_subset",
            "foldx_groups": CHAIN_GROUPS.get(target, ""),
            "foldx_work_dir": str(work_dir),
            "foldx_dG_pH7_4": pd.NA,
            "foldx_dG_pH6_0": pd.NA,
            "foldx_delta_pH6_minus_pH7_4": pd.NA,
            "foldx_neutral_interface_t2_score": pd.NA,
            "foldx_runtime_seconds": pd.NA,
        }
        try:
            if target not in CHAIN_GROUPS:
                raise ValueError(f"unknown_target:{target}")
            pdb_path = Path(str(row.get("pdb_path", "")))
            if not pdb_path.exists():
                raise FileNotFoundError(f"pdb_missing:{pdb_path}")
            if str(row.get("tool_status", "")) != "success":
                raise ValueError(f"pyrosetta_not_success:{row.get('tool_status')}")
            work_dir.mkdir(parents=True, exist_ok=True)
            input_pdb = work_dir / "input.pdb"
            if not input_pdb.exists() or input_pdb.stat().st_size != pdb_path.stat().st_size:
                shutil.copy2(pdb_path, input_pdb)
            values: dict[float, float] = {}
            groups = CHAIN_GROUPS[target]
            for ph in PH_POINTS:
                tag = f"AC_input_{str(ph).replace('.', '_')}"
                cmd = [
                    foldx_bin,
                    "--command=AnalyseComplex",
                    "--pdb=input.pdb",
                    f"--analyseComplexChains={groups}",
                    f"--pH={ph}",
                    f"--output-file={tag}",
                ]
                run(cmd, cwd=work_dir, env=env, timeout=timeout)
                fxout = find_interaction_file(work_dir, "input", ph)
                if fxout is None:
                    raise FileNotFoundError(f"foldx_interaction_missing:pH={ph}")
                value = parse_interaction_energy(fxout, "input.pdb", groups)
                if value is None:
                    raise ValueError(f"foldx_parse_failed:{fxout.name}")
                values[ph] = value
            d74 = values[7.4]
            d60 = values[6.0]
            result.update(
                {
                    "foldx_tool_status": "success",
                    "foldx_dG_pH7_4": d74,
                    "foldx_dG_pH6_0": d60,
                    "foldx_delta_pH6_minus_pH7_4": d60 - d74,
                    "foldx_neutral_interface_t2_score": foldx_score_from_dg(d74),
                }
            )
        except subprocess.TimeoutExpired as exc:
            result.update(
                {
                    "foldx_tool_status": "timeout",
                    "foldx_feature_missing_reason": str(exc),
                    "foldx_feature_confidence": "low",
                }
            )
        except subprocess.CalledProcessError as exc:
            reason = (exc.stderr or exc.stdout or str(exc)).strip().replace("\n", " ")[:500]
            result.update(
                {
                    "foldx_tool_status": "failed_runtime",
                    "foldx_feature_missing_reason": reason,
                    "foldx_feature_confidence": "low",
                }
            )
        except Exception as exc:
            reason = str(exc)
            status = "failed_input" if "pdb_missing" in reason or "unknown_target" in reason or "pyrosetta_not_success" in reason else "failed_runtime"
            result.update(
                {
                    "foldx_tool_status": status,
                    "foldx_feature_missing_reason": reason[:500],
                    "foldx_feature_confidence": "low",
                }
            )
        result["foldx_runtime_seconds"] = round(time.time() - started, 3)
        out.append(result)
    return out


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
    parser.add_argument("--stage", default="full")
    parser.add_argument("--foldx-bin", default=str(DEFAULT_FOLDX_BIN))
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work_root = out_dir / "structures/foldx_subset_work"
    work_root.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidate_csv)
    if "tier2_foldx_subset_selected" not in candidates.columns:
        raise SystemExit("tier2_foldx_subset_selected column missing")
    selected = candidates[candidates["tier2_foldx_subset_selected"].map(as_bool)].copy()
    if args.limit:
        selected = selected.head(args.limit).copy()
    pyro = pd.read_csv(args.pyrosetta_results)
    selected = selected.merge(pyro[["variant_id", "tool_status", "pdb_path"]], on="variant_id", how="left", validate="one_to_one")
    rows = selected.to_dict(orient="records")
    workers = max(1, min(args.workers, len(rows)))
    started = time.time()
    if workers == 1:
        nested = [worker(rows, args.foldx_bin, str(work_root), args.timeout)]
    else:
        with mp.Pool(processes=workers) as pool:
            nested = pool.starmap(
                worker,
                [(chunk, args.foldx_bin, str(work_root), args.timeout) for chunk in chunks(rows, workers)],
            )
    results = pd.DataFrame([item for chunk in nested for item in chunk])
    if not results.empty:
        for target, group in results.groupby("target"):
            ok = group[group["foldx_tool_status"].eq("success")]
            zero_mut_ids = set(
                selected[
                    selected["target"].eq(target)
                    & pd.to_numeric(selected.get("mutation_count", pd.Series(index=selected.index)), errors="coerce").fillna(1).eq(0)
                ]["variant_id"].astype(str)
            )
            parent = ok[ok["variant_id"].astype(str).isin(zero_mut_ids)]
            if parent.empty:
                parent = ok.head(1)
            for ph_col, out_col in [
                ("foldx_dG_pH7_4", "foldx_ddG_vs_parent_pH7_4"),
                ("foldx_dG_pH6_0", "foldx_ddG_vs_parent_pH6_0"),
            ]:
                baseline = pd.to_numeric(parent[ph_col], errors="coerce").dropna()
                if baseline.empty:
                    results.loc[results["target"].eq(target), out_col] = pd.NA
                else:
                    base_value = float(baseline.median())
                    values = pd.to_numeric(results.loc[results["target"].eq(target), ph_col], errors="coerce")
                    results.loc[results["target"].eq(target), out_col] = values - base_value
    results.to_csv(out_dir / f"tier2b_{args.stage}_foldx_subset_results.csv", index=False)
    status = results.groupby(["target", "foldx_tool_status"], dropna=False).size().reset_index(name="count")
    status.to_csv(out_dir / f"tier2b_{args.stage}_foldx_subset_status_summary.csv", index=False)
    report = [
        f"# Tier2-B {args.stage.title()} FoldX Subset Report",
        "",
        f"Selected rows: `{len(rows)}`",
        f"Workers: `{workers}`",
        f"Runtime seconds: `{time.time() - started:.1f}`",
        "",
        "## Status Counts",
        markdown_table(status),
    ]
    if not results.empty:
        report.extend(
            [
                "",
                "## Target Counts",
                markdown_table(results.groupby("target").size().reset_index(name="count")),
            ]
        )
    (out_dir / f"tier2b_{args.stage}_foldx_subset_review.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    selected_count = len(rows)
    success_count = int(results["foldx_tool_status"].eq("success").sum()) if not results.empty else 0
    if selected_count > 0 and success_count == 0:
        raise SystemExit("All selected FoldX subset rows failed")


if __name__ == "__main__":
    main()
