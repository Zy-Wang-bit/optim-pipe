#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run analyze_trajectory.py on all 32 r4 trajectories in parallel,
then compare_ph.py on all 16 variants, then aggregate deltas into md_metrics.csv.

Runs on h100. Expects:
- /home/ziyang/sdab_r4_md/runs/r4_NN/pH_{7.4,6.0}/production.xtc + .tpr
- /home/ziyang/sdab_r4_md/molecular_dynamics/configs/md_config_sdab.yaml

CPU strategy:
- 32 analyze_trajectory.py jobs in parallel (1 per trajectory)
- Each with OMP_NUM_THREADS=4, MKL_NUM_THREADS=4 → 32 × 4 = 128 cores used (h100: 192 total)
- Afterwards, 16 compare_ph.py jobs in parallel (lightweight, JSON reads)
- Final aggregation: single process reads 16 ph_comparison.json → md_metrics.csv
"""

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

ROOT = Path("/home/ziyang/sdab_r4_md")
MD_MODULE = ROOT / "molecular_dynamics"
ANALYZE = MD_MODULE / "analyze_trajectory.py"
COMPARE = MD_MODULE / "compare_ph.py"
CONFIG = MD_MODULE / "configs/md_config_sdab_r4.yaml"
RUNS = ROOT / "runs"
LOG_DIR = ROOT / "logs"
METRICS_CSV = ROOT / "md_metrics.csv"

VARIANTS = [f"r4_{i:02d}" for i in range(1, 17)]
PHS = [7.4, 6.0]

OMP_THREADS_PER_JOB = 4   # 32 jobs × 4 = 128 cores (fits in 192)
MAX_PARALLEL_ANALYZE = 32
MAX_PARALLEL_COMPARE = 16


def run_subprocess(cmd: list[str], log_path: Path, env: dict | None = None) -> int:
    """Run a subprocess with logs, return rc."""
    with open(log_path, "w") as lf:
        lf.write(f"# cmd: {' '.join(cmd)}\n")
        lf.flush()
        env2 = dict(os.environ)
        if env:
            env2.update(env)
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env2)
    return p.returncode


def analyze_one(traj_dir: str, log_path: str) -> tuple[str, int, float]:
    start = time.time()
    cmd = [
        sys.executable, str(ANALYZE),
        "--traj", traj_dir,
        "--config", str(CONFIG),
    ]
    env = {
        "OMP_NUM_THREADS": str(OMP_THREADS_PER_JOB),
        "MKL_NUM_THREADS": str(OMP_THREADS_PER_JOB),
        "OPENBLAS_NUM_THREADS": str(OMP_THREADS_PER_JOB),
    }
    rc = run_subprocess(cmd, Path(log_path), env)
    return (traj_dir, rc, time.time() - start)


def compare_one(variant_dir: str, log_path: str) -> tuple[str, int, float]:
    start = time.time()
    cmd = [
        sys.executable, str(COMPARE),
        "--variant-dir", variant_dir,
        "--base-ph", "7.4",
        "--target-ph", "6.0",
    ]
    rc = run_subprocess(cmd, Path(log_path))
    return (variant_dir, rc, time.time() - start)


def build_traj_jobs() -> list[tuple[str, str]]:
    jobs = []
    for v in VARIANTS:
        for ph in PHS:
            traj_dir = RUNS / v / f"pH_{ph}"
            log = LOG_DIR / f"analyze_{v}_pH{ph}.log"
            jobs.append((str(traj_dir), str(log)))
    return jobs


def build_compare_jobs() -> list[tuple[str, str]]:
    return [
        (str(RUNS / v), str(LOG_DIR / f"compare_{v}.log"))
        for v in VARIANTS
    ]


def aggregate_deltas() -> pd.DataFrame:
    rows = []
    for v in VARIANTS:
        cmp_path = RUNS / v / "ph_comparison.json"
        if not cmp_path.exists():
            print(f"  [missing] {cmp_path}")
            rows.append({"name": v, "status": "missing"})
            continue
        with open(cmp_path) as f:
            data = json.load(f)
        row = {"name": v, "status": "ok"}
        row.update(data.get("deltas", {}))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[analyze] launching {len(VARIANTS)*len(PHS)} analyze_trajectory jobs (max {MAX_PARALLEL_ANALYZE} parallel)")
    t0 = time.time()
    traj_jobs = build_traj_jobs()

    with ProcessPoolExecutor(max_workers=MAX_PARALLEL_ANALYZE) as pool:
        futures = [pool.submit(analyze_one, td, lp) for td, lp in traj_jobs]
        n_ok, n_fail = 0, 0
        for fut in as_completed(futures):
            td, rc, dur = fut.result()
            tag = "OK" if rc == 0 else "FAIL"
            print(f"  [{tag}] {td} rc={rc} ({dur:.1f}s)")
            if rc == 0:
                n_ok += 1
            else:
                n_fail += 1
    print(f"[analyze] done in {(time.time()-t0)/60:.1f} min | ok={n_ok} fail={n_fail}")

    print(f"\n[compare] launching {len(VARIANTS)} compare_ph jobs")
    t1 = time.time()
    cmp_jobs = build_compare_jobs()
    with ProcessPoolExecutor(max_workers=MAX_PARALLEL_COMPARE) as pool:
        futures = [pool.submit(compare_one, vd, lp) for vd, lp in cmp_jobs]
        for fut in as_completed(futures):
            vd, rc, dur = fut.result()
            tag = "OK" if rc == 0 else "FAIL"
            print(f"  [{tag}] {vd} rc={rc} ({dur:.1f}s)")
    print(f"[compare] done in {(time.time()-t1)/60:.1f} min")

    print(f"\n[aggregate] reading 16 ph_comparison.json -> {METRICS_CSV}")
    df = aggregate_deltas()
    df.to_csv(METRICS_CSV, index=False)
    print(df.to_string(index=False))
    print(f"\n[done] total elapsed: {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
