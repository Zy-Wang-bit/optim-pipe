#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resume the 32 sdab R4 MD simulations from production.cpt after the original
run was killed by a 24h Python subprocess timeout in run_md.py.

Bypasses run_md.py entirely — spawns gmx mdrun directly with -cpi/-append,
so no Python-side timeout can cut it short again.

Distributes sims across GPUs given on the command line (default 2,3,4,5,6,7
because another user is occupying GPU 0/1 on h100 at the time of this incident).

Usage:
    python3 resume_md_r4.py                     # uses GPUs 2-7
    python3 resume_md_r4.py --gpus 2,3,4,5,6,7  # explicit
"""

import argparse
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/root/sdab_r4_md")
RUNS = ROOT / "runs"
LOG_DIR = ROOT / "logs"

VARIANTS = [f"r4_{i:02d}" for i in range(1, 17)]
PHS = [7.4, 6.0]

GMX_BIN = "gmx"
NTOMP = 5
NTMPI = 1


def build_tasks(gpus):
    tasks = []
    idx = 0
    for variant in VARIANTS:
        for ph in PHS:
            workdir = RUNS / variant / f"pH_{ph}"
            cpt = workdir / "production.cpt"
            tpr = workdir / "production.tpr"
            log = LOG_DIR / f"{variant}_pH{ph}.log"
            gpu = gpus[idx % len(gpus)]
            idx += 1
            tasks.append({
                "variant": variant, "ph": ph,
                "workdir": workdir, "tpr": tpr, "cpt": cpt,
                "log": log, "gpu": gpu,
            })
    return tasks


def launch(task):
    if not task["cpt"].exists():
        print(f"  [skip] {task['variant']} pH{task['ph']} — no cpt")
        return None
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(task["gpu"])
    cmd = [
        GMX_BIN, "mdrun",
        "-s", str(task["tpr"]),
        "-deffnm", "production",
        "-cpi", str(task["cpt"]),
        "-append",
        "-ntomp", str(NTOMP),
        "-ntmpi", str(NTMPI),
    ]
    log_fp = open(task["log"], "a")
    log_fp.write(f"\n# === RESUME {time.strftime('%F %T')} ===\n")
    log_fp.write(f"# cmd: {' '.join(shlex.quote(c) for c in cmd)}\n")
    log_fp.write(f"# CUDA_VISIBLE_DEVICES={task['gpu']}\n")
    log_fp.flush()
    proc = subprocess.Popen(
        cmd, env=env,
        cwd=task["workdir"],
        stdout=log_fp, stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    task["proc"] = proc
    task["log_fp"] = log_fp
    return proc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpus", default="2,3,4,5,6,7",
                    help="comma-separated GPU indices to use")
    ap.add_argument("--stagger", type=float, default=3.0,
                    help="seconds between spawns (default 3)")
    args = ap.parse_args()

    gpus = [int(x) for x in args.gpus.split(",") if x.strip()]
    print(f"[resume] GPUs: {gpus}")

    tasks = build_tasks(gpus)
    print(f"[resume] {len(tasks)} sims to resume, {len(gpus)} GPUs, "
          f"{len(tasks)/len(gpus):.1f} sims/GPU avg")

    started = time.time()
    for t in tasks:
        p = launch(t)
        if p is not None:
            print(f"  [spawn] pid={p.pid} gpu={t['gpu']} {t['variant']} pH{t['ph']}")
            time.sleep(args.stagger)

    status_file = ROOT / "resume_status.txt"
    with open(status_file, "w") as sf:
        sf.write(f"resumed_at\t{time.strftime('%F %T')}\n")
        sf.write(f"gpus\t{gpus}\n")
        for t in tasks:
            p = t.get("proc")
            sf.write(f"{t['variant']}\tpH_{t['ph']}\tgpu={t['gpu']}\t"
                    f"pid={p.pid if p else 'none'}\tlog={t['log']}\n")
    print(f"[resume] status -> {status_file}")

    print(f"[resume] all spawned in {time.time()-started:.0f}s; waiting ...")
    for t in tasks:
        if t.get("proc") is not None:
            rc = t["proc"].wait()
            t["log_fp"].close()
            print(f"  [done] {t['variant']} pH{t['ph']} rc={rc}")
    print(f"[resume] elapsed: {(time.time()-started)/3600:.2f} h")


if __name__ == "__main__":
    main()
