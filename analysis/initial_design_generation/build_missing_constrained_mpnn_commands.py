#!/usr/bin/env python3
"""Build ProteinMPNN commands for missing constrained-generation jobs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.initial_design_generation import run_dry_run as dry


ROOT = dry.ROOT


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def generated_sample_count(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for header, _ in parse_fasta(path) if header.startswith("T="))


def load_jsonl_dict(path: Path) -> dict:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"Expected one JSON object in {path}, found {len(lines)}")
    return json.loads(lines[0])


def load_parsed_by_name(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        out[str(row["name"])] = row
    return out


def write_filtered_group(
    source_group_dir: Path,
    filtered_group_dir: Path,
    names: list[str],
) -> None:
    filtered_group_dir.mkdir(parents=True, exist_ok=True)
    for filename in ["chain_id.jsonl", "fixed_positions.jsonl", "omit_AA.jsonl"]:
        data = load_jsonl_dict(source_group_dir / filename)
        with (filtered_group_dir / filename).open("w") as fh:
            fh.write(json.dumps({name: data[name] for name in names}) + "\n")
    parsed = load_parsed_by_name(source_group_dir / "parsed_pdbs.jsonl")
    with (filtered_group_dir / "parsed_pdbs.jsonl").open("w") as fh:
        for name in names:
            fh.write(json.dumps(parsed[name]) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build commands for missing constrained MPNN jobs.")
    parser.add_argument("--runner-dir", required=True)
    parser.add_argument("--output-command-file", default="missing_constrained_generation_runner_commands.sh")
    args = parser.parse_args()

    runner_dir = ROOT / args.runner_dir if not Path(args.runner_dir).is_absolute() else Path(args.runner_dir)
    jobs = pd.read_csv(runner_dir / "mpnn_runner_jobs_constrained.csv")
    missing_rows: list[dict] = []
    for job in jobs.itertuples(index=False):
        out_folder = runner_dir / "mpnn_out"
        fasta_matches = list(out_folder.glob(f"**/seqs/{job.job_name}.fa"))
        observed = max((generated_sample_count(path) for path in fasta_matches), default=0)
        expected = int(job.planned_raw_samples)
        if observed < expected:
            missing_rows.append(
                {
                    "job_name": job.job_name,
                    "target": job.target,
                    "mpnn_temperature": float(job.mpnn_temperature),
                    "planned_raw_samples": expected,
                    "observed_generated_samples": observed,
                }
            )
    missing_columns = [
        "job_name",
        "target",
        "mpnn_temperature",
        "planned_raw_samples",
        "observed_generated_samples",
    ]
    missing = pd.DataFrame(missing_rows, columns=missing_columns)
    dry.write_csv(missing, runner_dir / "missing_constrained_mpnn_jobs.csv")
    command_rows: list[dict] = []
    commands = [
        "# ProteinMPNN commands for missing constrained-generation jobs only.",
        "# Generated from current FASTA output completeness.",
        "",
    ]
    if not missing.empty:
        plan = pd.read_csv(runner_dir / "mpnn_constrained_command_plan.csv")
        for command_index, ((temp, samples), sub) in enumerate(
            missing.groupby(["mpnn_temperature", "planned_raw_samples"], sort=True)
        ):
            temp_tag = str(temp).replace(".", "p")
            group_tag = f"T_{temp_tag}_N{int(samples)}"
            source_group_dir = runner_dir / "jsonl" / "by_temperature" / group_tag
            filtered_group_dir = runner_dir / "jsonl" / "missing_by_temperature" / group_tag
            names = sub["job_name"].astype(str).tolist()
            write_filtered_group(source_group_dir, filtered_group_dir, names)
            out_folder = plan[
                (plan["mpnn_temperature"].astype(float) == float(temp))
                & (plan["num_seq_per_target"].astype(int) == int(samples))
            ]["out_folder"].iloc[0]
            commands.append(
                " ".join(
                    [
                        "/data/ziyang/mamba/envs/proteinmpnn/bin/python",
                        "third_party/ProteinMPNN/protein_mpnn_run.py",
                        f"--jsonl_path {filtered_group_dir / 'parsed_pdbs.jsonl'}",
                        f"--chain_id_jsonl {filtered_group_dir / 'chain_id.jsonl'}",
                        f"--fixed_positions_jsonl {filtered_group_dir / 'fixed_positions.jsonl'}",
                        f"--omit_AA_jsonl {filtered_group_dir / 'omit_AA.jsonl'}",
                        f"--out_folder {out_folder}",
                        f"--num_seq_per_target {int(samples)}",
                        f"--sampling_temp {temp}",
                        "--seed 20260523",
                        "--batch_size 1",
                    ]
                )
            )
            command_rows.append(
                {
                    "command_index": command_index,
                    "mpnn_temperature": temp,
                    "num_seq_per_target": int(samples),
                    "missing_job_count": int(len(sub)),
                    "planned_missing_raw_samples": int(sub["planned_raw_samples"].sum()),
                    "observed_generated_samples": int(sub["observed_generated_samples"].sum()),
                    "parsed_pdbs_jsonl": str(filtered_group_dir / "parsed_pdbs.jsonl"),
                    "out_folder": out_folder,
                }
            )
    (runner_dir / args.output_command_file).write_text("\n".join(commands) + "\n")
    command_columns = [
        "command_index",
        "mpnn_temperature",
        "num_seq_per_target",
        "missing_job_count",
        "planned_missing_raw_samples",
        "observed_generated_samples",
        "parsed_pdbs_jsonl",
        "out_folder",
    ]
    dry.write_csv(pd.DataFrame(command_rows, columns=command_columns), runner_dir / "missing_constrained_command_plan.csv")
    print(
        {
            "status": "missing_commands_built",
            "runner_dir": str(runner_dir),
            "missing_job_count": int(len(missing)),
            "missing_command_count": int(len(command_rows)),
            "planned_missing_raw_samples": int(missing["planned_raw_samples"].sum()) if not missing.empty else 0,
            "observed_generated_samples": int(missing["observed_generated_samples"].sum()) if not missing.empty else 0,
            "command_file": str(runner_dir / args.output_command_file),
        }
    )


if __name__ == "__main__":
    main()
