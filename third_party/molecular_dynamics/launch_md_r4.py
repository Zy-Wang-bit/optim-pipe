#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Launch 32 parallel GROMACS MD simulations on 8 H100 GPUs (4 per GPU).

Runs on the h100 node; assumes:
- /root/sdab_r4_md/molecular_dynamics/  (this module, rsynced from optim-pipe)
- /root/sdab_r4_md/inputs/r4_01.pdb..r4_16.pdb
- /root/sdab_r4_md/runs/     (output trajectories)
- /root/sdab_r4_md/logs/     (per-sim stdout/stderr)
"""

import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/root/sdab_r4_md")
MD_MODULE = ROOT / "molecular_dynamics"
RUN_MD = MD_MODULE / "run_md.py"
CONFIG = MD_MODULE / "configs/md_config_sdab.yaml"
INPUTS = ROOT / "inputs"
OUT_DIR = ROOT / "runs"
LOG_DIR = ROOT / "logs"

VARIANTS = [f"r4_{i:02d}" for i in range(1, 17)]
PHS = [7.4, 6.0]
N_GPUS = 8


def build_tasks():
    tasks = []
    for vi, variant in enumerate(VARIANTS):
        for pi, ph in enumerate(PHS):
            gpu = (vi * len(PHS) + pi) % N_GPUS
            pdb_path = INPUTS / f"{variant}.pdb"
            log_path = LOG_DIR / f"{variant}_pH{ph}.log"
            tasks.append({
                "variant": variant,
                "ph": ph,
                "gpu": gpu,
                "pdb": pdb_path,
                "log": log_path,
            })
    return tasks


def launch(task):
    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(task["gpu"])
    cmd = [
        sys.executable, str(RUN_MD),
        "--pdb", str(task["pdb"]),
        "--ph", str(task["ph"]),
        "--output-dir", str(OUT_DIR),
        "--config", str(CONFIG),
        "--protonation", "fixed",
    ]
    log_fp = open(task["log"], "w")
    log_fp.write(f"# cmd: {' '.join(shlex.quote(c) for c in cmd)}\n")
    log_fp.write(f"# CUDA_VISIBLE_DEVICES={task['gpu']}\n")
    log_fp.flush()
    proc = subprocess.Popen(
        cmd, env=env, stdout=log_fp, stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    task["proc"] = proc
    task["log_fp"] = log_fp
    return proc


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    tasks = build_tasks()
    print(f"[launch] {len(tasks)} tasks, {N_GPUS} GPUs, {len(tasks)//N_GPUS} sims/GPU")
    print(f"[launch] MD_MODULE={MD_MODULE} CONFIG={CONFIG}")

    started = time.time()
    for t in tasks:
        launch(t)
        print(f"  [spawn] pid={t['proc'].pid} gpu={t['gpu']} {t['variant']} pH{t['ph']}")
        time.sleep(2.0)  # stagger pdb2gmx file I/O

    status_file = ROOT / "launch_status.txt"
    with open(status_file, "w") as sf:
        sf.write(f"launched_at\t{time.strftime('%F %T')}\n")
        for t in tasks:
            sf.write(f"{t['variant']}\tpH_{t['ph']}\tgpu={t['gpu']}\tpid={t['proc'].pid}\tlog={t['log']}\n")
    print(f"[launch] status -> {status_file}")

    print(f"[launch] all 32 MD spawned; waiting for completion ...")
    results = []
    for t in tasks:
        rc = t["proc"].wait()
        results.append((t["variant"], t["ph"], rc))
        t["log_fp"].close()
        print(f"  [done] {t['variant']} pH{t['ph']} rc={rc}")

    elapsed = time.time() - started
    print(f"[launch] all done in {elapsed/3600:.2f} h")
    ok = sum(1 for _, _, rc in results if rc == 0)
    print(f"[launch] success: {ok}/{len(results)}")


if __name__ == "__main__":
    main()
