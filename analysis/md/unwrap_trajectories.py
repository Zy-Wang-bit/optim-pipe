#!/usr/bin/env python3
"""预处理: 对所有轨迹做 PBC unwrap, 保存为 production_unwrap.xtc"""

import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import MDAnalysis as mda
from MDAnalysis.transformations import unwrap, center_in_box


def unwrap_one(traj_dir: str) -> str:
    """Unwrap 单条轨迹, 保存 production_unwrap.xtc"""
    traj_dir = Path(traj_dir)
    tpr = traj_dir / "production.tpr"
    xtc = traj_dir / "production.xtc"
    out = traj_dir / "production_unwrap.xtc"

    if out.exists() and out.stat().st_size > 1_000_000:
        return f"SKIP {traj_dir.parent.name}/{traj_dir.name}"

    u = mda.Universe(str(tpr), str(xtc))
    protein = u.select_atoms("protein")
    transforms = [unwrap(protein), center_in_box(protein, center="mass")]
    u.trajectory.add_transformations(*transforms)

    all_atoms = u.select_atoms("all")
    with mda.Writer(str(out), all_atoms.n_atoms) as w:
        for ts in u.trajectory:
            w.write(all_atoms)

    return f"OK   {traj_dir.parent.name}/{traj_dir.name} → {out.stat().st_size // 1024 // 1024}MB"


def main():
    repo = Path(__file__).resolve().parent.parent.parent
    md_dir = repo / "experiments/sdab_variants/md"

    traj_dirs = []
    for i in range(15):
        vid = f"sdab_v2_{i:04d}"
        for ph in ["7.4", "6.0"]:
            d = md_dir / vid / f"pH_{ph}"
            if (d / "production.xtc").exists():
                traj_dirs.append(str(d))

    print(f"Unwrapping {len(traj_dirs)} trajectories...")

    # 并行, 但限制 worker 数 (每个 worker 内存密集)
    max_workers = min(10, len(traj_dirs))
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(unwrap_one, d): d for d in traj_dirs}
        for fut in as_completed(futures):
            try:
                print(fut.result())
            except Exception as e:
                print(f"FAIL {futures[fut]}: {e}")

    n_done = sum(1 for d in traj_dirs if (Path(d) / "production_unwrap.xtc").exists())
    print(f"\nDone: {n_done}/{len(traj_dirs)}")


if __name__ == "__main__":
    main()
