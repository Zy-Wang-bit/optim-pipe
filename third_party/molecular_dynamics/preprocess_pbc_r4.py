#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Preprocess all 32 sdab R4 production.xtc files with PBC unwrap.

Follows the repo-standard recipe from analysis/md/unwrap_trajectories.py:
MDAnalysis transformations.unwrap(protein) + center_in_box(protein). Writes
production_unwrap.xtc next to production.xtc; BaseAnalyzer auto-prefers it.
"""

import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/home/ziyang/sdab_r4_md")
RUNS = ROOT / "runs"

VARIANTS = [f"r4_{i:02d}" for i in range(1, 17)]
PHS = [7.4, 6.0]

MAX_WORKERS = 10  # I/O-bound; matches analysis/md/unwrap_trajectories.py


def unwrap_one(traj_dir: str) -> str:
    import MDAnalysis as mda
    from MDAnalysis.transformations import unwrap, center_in_box

    traj_dir = Path(traj_dir)
    tpr = traj_dir / "production.tpr"
    xtc = traj_dir / "production.xtc"
    out = traj_dir / "production_unwrap.xtc"

    if out.exists() and out.stat().st_size > 1_000_000:
        return f"SKIP {traj_dir.parent.name}/{traj_dir.name}"

    t0 = time.time()
    u = mda.Universe(str(tpr), str(xtc))
    protein = u.select_atoms("protein")
    transforms = [unwrap(protein), center_in_box(protein, center="mass")]
    u.trajectory.add_transformations(*transforms)

    all_atoms = u.select_atoms("all")
    with mda.Writer(str(out), all_atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(all_atoms)

    size_mb = out.stat().st_size // 1024 // 1024
    return f"OK   {traj_dir.parent.name}/{traj_dir.name} → {size_mb}MB ({time.time()-t0:.1f}s)"


def main():
    traj_dirs = []
    for v in VARIANTS:
        for ph in PHS:
            d = RUNS / v / f"pH_{ph}"
            if (d / "production.xtc").exists():
                traj_dirs.append(str(d))

    print(f"[pbc] Unwrapping {len(traj_dirs)} trajectories (workers={MAX_WORKERS})")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(unwrap_one, d): d for d in traj_dirs}
        for fut in as_completed(futures):
            try:
                print(fut.result())
            except Exception as e:
                print(f"FAIL {futures[fut]}: {e}")
    print(f"[pbc] total: {(time.time()-t0)/60:.1f} min")

    n_done = sum(1 for d in traj_dirs
                 if (Path(d) / "production_unwrap.xtc").exists())
    print(f"[pbc] production_unwrap.xtc present: {n_done}/{len(traj_dirs)}")


if __name__ == "__main__":
    main()
