#!/usr/bin/env python3
"""Run generated ProteinMPNN command files with bounded parallelism.

Default mode is dry-run. Pass ``--execute`` only on an appropriate compute node.
The script is intentionally generic so the same runner can execute score-only
shards or constrained-generation temperature groups.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


@dataclass
class CommandResult:
    index: int
    returncode: int
    stdout_path: Path
    stderr_path: Path


def load_commands(path: Path) -> list[str]:
    commands: list[str] = []
    with path.open() as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            commands.append(stripped)
    return commands


def largest_divisor_at_most(value: int, max_value: int) -> int:
    for candidate in range(min(value, max_value), 0, -1):
        if value % candidate == 0:
            return candidate
    return 1


def apply_command_overrides(
    command: str,
    python_executable: str | None = None,
    auto_batch_size_max: int | None = None,
) -> list[str]:
    command_args = shlex.split(command)
    if python_executable:
        command_args[0] = python_executable
    if auto_batch_size_max:
        if auto_batch_size_max < 1:
            raise SystemExit("--auto-batch-size-max must be >= 1")
        try:
            n_ix = command_args.index("--num_seq_per_target")
            batch_ix = command_args.index("--batch_size")
        except ValueError:
            return command_args
        num_seq = int(command_args[n_ix + 1])
        command_args[batch_ix + 1] = str(largest_divisor_at_most(num_seq, auto_batch_size_max))
    return command_args


def run_one(
    index: int,
    command: str,
    log_dir: Path,
    cuda_devices: list[str] | None = None,
    python_executable: str | None = None,
    auto_batch_size_max: int | None = None,
    device_locks: dict[str, threading.Lock] | None = None,
) -> CommandResult:
    stdout_path = log_dir / f"cmd_{index:04d}.stdout.log"
    stderr_path = log_dir / f"cmd_{index:04d}.stderr.log"
    env = os.environ.copy()
    device = None
    if cuda_devices:
        device = cuda_devices[index % len(cuda_devices)]
        env["CUDA_VISIBLE_DEVICES"] = device
    command_args = apply_command_overrides(command, python_executable, auto_batch_size_max)
    lock = device_locks.get(device) if device and device_locks else None
    if lock is None:
        with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
            proc = subprocess.run(command_args, cwd=ROOT, env=env, stdout=stdout, stderr=stderr)
    else:
        with lock, stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
            proc = subprocess.run(command_args, cwd=ROOT, env=env, stdout=stdout, stderr=stderr)
    return CommandResult(index=index, returncode=proc.returncode, stdout_path=stdout_path, stderr_path=stderr_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a generated ProteinMPNN command file.")
    parser.add_argument("command_file", help="Path to a generated command file.")
    parser.add_argument("--max-workers", type=int, default=1, help="Maximum commands to run concurrently.")
    parser.add_argument("--execute", action="store_true", help="Actually run commands. Default is dry-run.")
    parser.add_argument("--start-index", type=int, default=0, help="Start at this zero-based command index.")
    parser.add_argument("--limit", type=int, help="Run only the first N commands after parsing.")
    parser.add_argument("--indices", help="Comma-separated zero-based original command indices to run.")
    parser.add_argument(
        "--cuda-devices",
        help="Comma-separated physical GPU IDs to assign round-robin via CUDA_VISIBLE_DEVICES.",
    )
    parser.add_argument(
        "--python-executable",
        help="Override the Python executable at the start of each generated command.",
    )
    parser.add_argument(
        "--auto-batch-size-max",
        type=int,
        help="Replace --batch_size with the largest divisor of --num_seq_per_target up to this value.",
    )
    parser.add_argument(
        "--log-dir",
        default="results/initial_design_generation/p0_mpnn_runner_inputs/mpnn_command_logs",
        help="Directory for stdout/stderr logs when --execute is used.",
    )
    args = parser.parse_args()

    command_file = ROOT / args.command_file if not Path(args.command_file).is_absolute() else Path(args.command_file)
    commands = load_commands(command_file)
    if args.indices:
        selected_indices = [int(item.strip()) for item in args.indices.split(",") if item.strip()]
        commands = [commands[index] for index in selected_indices]
    if args.start_index < 0:
        raise SystemExit("--start-index must be >= 0")
    if not args.indices:
        commands = commands[args.start_index :]
    if args.limit is not None:
        if args.limit < 1:
            raise SystemExit("--limit must be >= 1")
        commands = commands[: args.limit]
    if not commands:
        raise SystemExit(f"No runnable commands found in {command_file}")
    if args.max_workers < 1:
        raise SystemExit("--max-workers must be >= 1")

    if not args.execute:
        print(
            {
                "status": "dry_run",
                "command_file": str(command_file),
                "command_count": len(commands),
                "start_index": args.start_index,
                "max_workers": args.max_workers,
                "cuda_devices": args.cuda_devices,
                "python_executable": args.python_executable,
                "auto_batch_size_max": args.auto_batch_size_max,
                "first_command": commands[0],
                "first_command_args_after_overrides": apply_command_overrides(
                    commands[0],
                    args.python_executable,
                    args.auto_batch_size_max,
                ),
            }
        )
        return

    log_dir = ROOT / args.log_dir if not Path(args.log_dir).is_absolute() else Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    cuda_devices = None
    if args.cuda_devices:
        cuda_devices = [device.strip() for device in args.cuda_devices.split(",") if device.strip()]
        if not cuda_devices:
            raise SystemExit("--cuda-devices did not contain any usable device IDs")
    device_locks = {device: threading.Lock() for device in cuda_devices or []}
    results: list[CommandResult] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(
                run_one,
                index,
                command,
                log_dir,
                cuda_devices,
                args.python_executable,
                args.auto_batch_size_max,
                device_locks,
            ): index
            for index, command in enumerate(commands)
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                {
                    "command_index": result.index,
                    "returncode": result.returncode,
                    "stdout": str(result.stdout_path),
                    "stderr": str(result.stderr_path),
                }
            )

    failed = [r for r in results if r.returncode != 0]
    print(
        {
            "status": "completed_with_failures" if failed else "completed",
            "command_file": str(command_file),
            "command_count": len(commands),
            "failed_count": len(failed),
            "log_dir": str(log_dir),
        }
    )
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
