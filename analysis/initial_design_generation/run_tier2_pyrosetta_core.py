#!/usr/bin/env python3
"""Run Tier2 PyRosetta local mutation/repack scoring.

Run this script with:
  /data/ziyang/mamba/envs/pyrosetta/bin/python -m analysis.initial_design_generation.run_tier2_pyrosetta_core ...
"""

from __future__ import annotations

import argparse
import math
import multiprocessing as mp
import os
import re
import time
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BACKBONES = {
    "Ab_1E62": {
        "pdb": ROOT / "results/initial_design_generation/p0_mpnn_backbones/af3_1E62_seed7_model5_rank1.pdb",
        "design_chain": "L",
        "binder_chains": ["H", "L"],
        "target_chains": ["C"],
        "protected": {("L", 23)},
    },
    "Ab_sdAb": {
        "pdb": ROOT / "results/initial_design_generation/p0_mpnn_backbones/af3_sdAb_seed16_model5_rank1.pdb",
        "design_chain": "A",
        "binder_chains": ["A"],
        "target_chains": ["B"],
        "protected": {("A", 96)},
    },
}
MUT_RE = re.compile(r"^([A-Za-z])([A-Za-z])(\d+)([A-Za-z])$")
AA_1TO3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def parse_mutations(value: object) -> list[tuple[str, str, int, str]]:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return []
    muts: list[tuple[str, str, int, str]] = []
    for token in str(value).replace(",", ";").split(";"):
        token = token.strip()
        if not token or token.lower() == "nan":
            continue
        m = MUT_RE.match(token)
        if not m:
            raise ValueError(f"cannot_parse_mutation:{token}")
        chain, old, pos, new = m.groups()
        muts.append((chain, old, int(pos), new))
    return muts


def chunks(items: list[dict[str, Any]], n: int) -> list[list[dict[str, Any]]]:
    if not items:
        return []
    n = max(1, n)
    size = math.ceil(len(items) / n)
    return [items[i : i + size] for i in range(0, len(items), size)]


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


def atom_is_h(atom_name: str) -> bool:
    return atom_name.strip().upper().startswith("H")


def local_clash_count(pose, residue_indices: set[int], cutoff: float = 2.1) -> int:
    indices = sorted(residue_indices)
    count = 0
    for i_pos, i in enumerate(indices):
        ri = pose.residue(i)
        chain_i = pose.pdb_info().chain(i)
        num_i = pose.pdb_info().number(i)
        for j in indices[i_pos + 1 :]:
            rj = pose.residue(j)
            chain_j = pose.pdb_info().chain(j)
            num_j = pose.pdb_info().number(j)
            if chain_i == chain_j and abs(num_i - num_j) <= 1:
                continue
            for ai in range(1, ri.natoms() + 1):
                if atom_is_h(ri.atom_name(ai)):
                    continue
                xyz_i = ri.xyz(ai)
                for aj in range(1, rj.natoms() + 1):
                    if atom_is_h(rj.atom_name(aj)):
                        continue
                    if xyz_i.distance(rj.xyz(aj)) < cutoff:
                        count += 1
                        if count > 999:
                            return count
    return count


def worker_process(args: tuple[list[dict[str, Any]], str, float, bool]) -> list[dict[str, Any]]:
    rows, out_pdb_dir, repack_shell, minimize = args
    import pyrosetta
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.protocols.minimization_packing import MinMover, PackRotamersMover
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init("-mute all -ignore_unrecognized_res true")
    try:
        scorefxn = pyrosetta.create_score_function("ref2015")
    except Exception:
        scorefxn = pyrosetta.get_score_function()

    base_poses: dict[str, Any] = {}
    base_scores: dict[str, float] = {}
    results: list[dict[str, Any]] = []
    out_dir = Path(out_pdb_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for row in rows:
        start = time.time()
        target = str(row["target"])
        variant_id = str(row["variant_id"])
        result = {
            "variant_id": variant_id,
            "target": target,
            "tool_status": "not_run",
            "feature_missing_reason": "",
            "feature_confidence": "pyrosetta_local_repack",
            "pdb_path": "",
            "n_mutations_applied": 0,
            "rosetta_total_score": pd.NA,
            "rosetta_parent_score": pd.NA,
            "rosetta_delta_score": pd.NA,
            "local_repack_residue_count": pd.NA,
            "local_clash_count": pd.NA,
            "local_structure_validity_t2_score": pd.NA,
            "runtime_seconds": pd.NA,
        }
        try:
            if target not in BACKBONES:
                raise ValueError(f"unknown_target:{target}")
            info = BACKBONES[target]
            if target not in base_poses:
                base_poses[target] = pyrosetta.pose_from_pdb(str(info["pdb"]))
                base_scores[target] = float(scorefxn(base_poses[target]))
            pose = base_poses[target].clone()
            pdb_info = pose.pdb_info()
            mutations = parse_mutations(row.get("normalized_mutation_list", row.get("mutation_list", "")))
            mutated_indices: list[int] = []
            for chain, old, pos, new in mutations:
                if (chain, pos) in info["protected"]:
                    raise ValueError(f"protected_residue_mutated:{chain}{pos}")
                ros_idx = pdb_info.pdb2pose(chain, pos)
                if ros_idx == 0:
                    raise ValueError(f"pdb_mapping_missing:{chain}{pos}")
                actual = pose.residue(ros_idx).name1()
                if actual != old:
                    raise ValueError(f"parent_aa_mismatch:{chain}{old}{pos}{new}:pdb={actual}")
                if new not in AA_1TO3:
                    raise ValueError(f"unsupported_mutant_aa:{new}")
                mutator = MutateResidue(ros_idx, AA_1TO3[new])
                mutator.set_preserve_atom_coords(False)
                mutator.apply(pose)
                mutated_indices.append(ros_idx)

            repack_set: set[int] = set(mutated_indices)
            for mut_idx in mutated_indices:
                mut_xyz = pose.residue(mut_idx).nbr_atom_xyz()
                for i in range(1, pose.total_residue() + 1):
                    if i in repack_set:
                        continue
                    if mut_xyz.distance(pose.residue(i).nbr_atom_xyz()) <= repack_shell:
                        repack_set.add(i)

            if repack_set:
                task = TaskFactory().create_packer_task(pose)
                task.restrict_to_repacking()
                for i in range(1, pose.total_residue() + 1):
                    task.nonconst_residue_task(i).prevent_repacking()
                for i in repack_set:
                    task.nonconst_residue_task(i).restrict_to_repacking()
                PackRotamersMover(scorefxn, task).apply(pose)
                if minimize:
                    mm = MoveMap()
                    mm.set_bb(False)
                    mm.set_chi(False)
                    for i in repack_set:
                        mm.set_chi(i, True)
                    minimizer = MinMover()
                    minimizer.movemap(mm)
                    minimizer.score_function(scorefxn)
                    minimizer.apply(pose)

            score = float(scorefxn(pose))
            parent_score = base_scores[target]
            delta = score - parent_score
            clashes = local_clash_count(pose, repack_set) if repack_set else 0
            validity = 1.0 - min(0.65, clashes / 30.0) - min(0.35, max(delta, 0.0) / 80.0)
            validity = max(0.0, min(1.0, validity))
            pdb_path = out_dir / f"{variant_id}.pdb"
            pose.dump_pdb(str(pdb_path))
            result.update(
                {
                    "tool_status": "success",
                    "pdb_path": str(pdb_path),
                    "n_mutations_applied": len(mutated_indices),
                    "rosetta_total_score": score,
                    "rosetta_parent_score": parent_score,
                    "rosetta_delta_score": delta,
                    "local_repack_residue_count": len(repack_set),
                    "local_clash_count": clashes,
                    "local_structure_validity_t2_score": validity,
                    "runtime_seconds": round(time.time() - start, 3),
                }
            )
        except Exception as exc:
            result.update(
                {
                    "tool_status": "failed_input" if "mismatch" in str(exc) or "mapping" in str(exc) or "protected" in str(exc) else "failed_runtime",
                    "feature_missing_reason": str(exc),
                    "runtime_seconds": round(time.time() - start, 3),
                }
            )
        results.append(result)
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--stage", choices=["smoke", "calibration", "full"], required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--repack-shell", type=float, default=8.0)
    parser.add_argument("--minimize", action="store_true")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input_csv)
    if args.limit:
        df = df.head(args.limit).copy()
    rows = df.to_dict(orient="records")
    worker_count = max(1, min(args.workers, len(rows)))
    task_chunks = chunks(rows, worker_count)
    pdb_dir = out_dir / "structures" / "pyrosetta"
    started = time.time()
    if worker_count == 1:
        nested = [worker_process((task_chunks[0], str(pdb_dir), args.repack_shell, args.minimize))]
    else:
        with mp.get_context("spawn").Pool(processes=worker_count) as pool:
            nested = pool.map(worker_process, [(c, str(pdb_dir), args.repack_shell, args.minimize) for c in task_chunks])
    results = [row for chunk in nested for row in chunk]
    result_df = pd.DataFrame(results)
    result_path = out_dir / f"tier2b_{args.stage}_pyrosetta_results.csv"
    result_df.to_csv(result_path, index=False)
    summary = (
        result_df.groupby(["target", "tool_status"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["target", "tool_status"])
    )
    summary_path = out_dir / f"tier2b_{args.stage}_pyrosetta_failure_summary.csv"
    summary.to_csv(summary_path, index=False)
    report = [
        f"# Tier2-B {args.stage.title()} PyRosetta Core Report",
        "",
        f"Input: `{args.input_csv}`",
        f"Rows: `{len(rows)}`",
        f"Workers: `{worker_count}`",
        f"Runtime seconds: `{time.time() - started:.1f}`",
        "",
        "## Status Counts",
        markdown_table(summary),
    ]
    (out_dir / f"tier2b_{args.stage}_pyrosetta_review.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    if result_df["tool_status"].eq("success").sum() == 0:
        raise SystemExit("No PyRosetta candidates succeeded")


if __name__ == "__main__":
    main()
