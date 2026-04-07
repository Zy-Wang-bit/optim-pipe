#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 8 (Phase C): 突变体结构生成

从候选列表 CSV 读取突变，使用 PyRosetta 或 SimpleFold 生成突变体 PDB。
支持 FoldX 格式 (FA102H) 和统一格式 (HF102H) 的突变标注。

输出到 {phase_c_dir}/structures/{rosetta,simplefold}/
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

# 复用 analysis/naming/convert.py 的突变解析
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.naming.convert import foldx_to_unified, parse_unified, hl_to_pdb_chain

# ── 氨基酸映射 ─────────────────────────────────────────────────────────────

AA_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}
AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


# ── 突变解析 ───────────────────────────────────────────────────────────────

def parse_foldx_mutations(mut_str):
    """
    解析 FoldX 格式突变串。

    'FA102H,SA40A' → [('A', 102, 'F', 'H'), ('A', 40, 'S', 'A')]
    格式: {OrigAA}{PDB_Chain}{Pos}{NewAA}
    """
    if not mut_str or (isinstance(mut_str, float) and pd.isna(mut_str)):
        return []
    mutations = []
    for token in str(mut_str).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        unified = foldx_to_unified(token)
        hl, orig, pos, new = parse_unified(unified)
        pdb_chain = hl_to_pdb_chain(hl)
        mutations.append((pdb_chain, pos, orig, new))
    return mutations


def parse_unified_mutations(mut_str, chain_map):
    """
    解析统一格式突变串。

    'HF102H;LS42H' → [('A', 102, 'F', 'H'), ('B', 42, 'S', 'H')]
    格式: {H|L}{OrigAA}{Pos}{NewAA}，chain_map 将 H/L 映射到 PDB 链号
    """
    if not mut_str or (isinstance(mut_str, float) and pd.isna(mut_str)):
        return []
    mutations = []
    for token in str(mut_str).replace(",", ";").split(";"):
        token = token.strip()
        if not token or token.startswith("none"):
            continue
        hl, orig, pos, new = parse_unified(token)
        pdb_chain = chain_map[hl]
        mutations.append((pdb_chain, pos, orig, new))
    return mutations


def load_candidates(cfg):
    """
    从 Phase C 输入 CSV 加载候选变体。

    Returns: [{"variant_id": str, "mutations": [(chain, pos, wt, mut), ...]}]
    """
    pc = cfg["phase_c"]
    csv_path = pc["input"]["csv"]
    mut_col = pc["input"]["mutations_column"]
    fmt = pc["input"]["mutations_format"]

    df = pd.read_csv(csv_path)
    if mut_col not in df.columns:
        raise ValueError(f"候选 CSV 缺少 '{mut_col}' 列: {csv_path}")

    # 统一格式需要 H/L → PDB chain 映射
    chain_map = None
    if fmt == "unified":
        chains = pc["chains"]
        chain_map = {"H": chains["heavy"], "L": chains["light"]}

    candidates = []
    for idx, row in df.iterrows():
        # variant_id: 优先用 mpdb stem，否则用行号
        if "mpdb" in row and pd.notna(row["mpdb"]):
            vid = Path(str(row["mpdb"])).stem
        elif "variant_name" in row and pd.notna(row["variant_name"]):
            vid = str(row["variant_name"])
        else:
            vid = f"var_{idx:06d}"

        if fmt == "foldx":
            muts = parse_foldx_mutations(row[mut_col])
        elif fmt == "unified":
            muts = parse_unified_mutations(row[mut_col], chain_map)
        else:
            raise ValueError(f"未知突变格式: {fmt}")

        candidates.append({"variant_id": vid, "mutations": muts})

    return candidates


# ── PyRosetta 方法 ───────────────────────────────────────────────────────────


def build_with_rosetta(candidates, template_pdb, out_dir, repack_shell, minimize, scorefxn):
    """用 PyRosetta 生成突变体 PDB。"""
    import pyrosetta
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
    from pyrosetta.rosetta.core.kinematics import MoveMap

    pyrosetta.init("-mute all")

    os.makedirs(out_dir, exist_ok=True)
    sfxn = pyrosetta.create_score_function(scorefxn)

    template_pose = pyrosetta.pose_from_pdb(template_pdb)
    pdb_info = template_pose.pdb_info()

    results = []
    for cand in candidates:
        vid = cand["variant_id"]
        mutations = cand["mutations"]
        print(f"[Rosetta] {vid}: {len(mutations)} mutations...", end=" ")

        pose = template_pose.clone()

        mutated_indices = []
        for pdb_chain, pdb_resid, wt_aa, mut_aa in mutations:
            ros_idx = pdb_info.pdb2pose(pdb_chain, pdb_resid)
            if ros_idx == 0:
                print(f"\n  警告: 位点 {pdb_chain}:{pdb_resid} 在模板中不存在，跳过")
                continue

            actual_aa = pose.residue(ros_idx).name1()
            if actual_aa != wt_aa:
                print(f"\n  警告: {pdb_chain}:{pdb_resid} 期望 {wt_aa} 实际 {actual_aa}")

            new_type = AA_1TO3[mut_aa]
            mutator = MutateResidue(ros_idx, new_type)
            mutator.set_preserve_atom_coords(False)
            mutator.apply(pose)
            mutated_indices.append(ros_idx)

        if not mutated_indices:
            print("无突变，保存模板副本")
            out_path = os.path.join(out_dir, f"{vid}.pdb")
            pose.dump_pdb(out_path)
            results.append({"variant_id": vid, "pdb_file": f"{vid}.pdb", "n_mutations": 0})
            continue

        # Repack: 突变位点 + 周围 shell 内残基
        task = TaskFactory().create_packer_task(pose)
        task.restrict_to_repacking()

        for i in range(1, pose.total_residue() + 1):
            task.nonconst_residue_task(i).prevent_repacking()

        repack_set = set(mutated_indices)
        for mut_idx in mutated_indices:
            mut_xyz = pose.residue(mut_idx).nbr_atom_xyz()
            for i in range(1, pose.total_residue() + 1):
                if i in repack_set:
                    continue
                nbr_xyz = pose.residue(i).nbr_atom_xyz()
                dist = mut_xyz.distance(nbr_xyz)
                if dist <= repack_shell:
                    repack_set.add(i)

        for i in repack_set:
            task.nonconst_residue_task(i).restrict_to_repacking()

        packer = PackRotamersMover(sfxn, task)
        packer.apply(pose)

        if minimize:
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            for i in repack_set:
                mm.set_chi(i, True)
            minimizer = MinMover()
            minimizer.movemap(mm)
            minimizer.score_function(sfxn)
            minimizer.apply(pose)

        out_path = os.path.join(out_dir, f"{vid}.pdb")
        pose.dump_pdb(out_path)
        print(f"OK → {vid}.pdb")

        # 将突变记为 FoldX 格式用于追溯
        mut_str = ",".join(f"{wt}{ch}{pos}{mt}" for ch, pos, wt, mt in mutations)
        results.append({
            "variant_id": vid,
            "pdb_file": f"{vid}.pdb",
            "n_mutations": len(mutations),
            "mutations": mut_str,
        })

    # 保存映射表
    map_path = os.path.join(os.path.dirname(out_dir), "mutant_map.csv")
    with open(map_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["variant_id", "pdb_file", "n_mutations", "mutations"])
        writer.writeheader()
        writer.writerows(results)

    print(f"\n[Rosetta] 完成: {len(results)} 个变体 → {out_dir}")
    print(f"[Rosetta] 映射表: {map_path}")


# ── SimpleFold 方法 ──────────────────────────────────────────────────────────


def build_with_simplefold(candidates, template_pdb, out_dir, sf_cfg):
    """用 SimpleFold 从序列预测结构。注意: 需要候选中包含序列信息。"""
    os.makedirs(out_dir, exist_ok=True)

    fasta_dir = os.path.join(out_dir, "fasta_input")
    os.makedirs(fasta_dir, exist_ok=True)

    for cand in candidates:
        vid = cand["variant_id"]
        seq = cand.get("sequence", "")
        if not seq:
            print(f"[SimpleFold] 跳过 {vid}: 无序列信息")
            continue
        fasta_path = os.path.join(fasta_dir, f"{vid}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{vid}\n{seq}\n")

    cmd = [
        "simplefold",
        "--simplefold_model", sf_cfg["model"],
        "--num_steps", str(sf_cfg["num_steps"]),
        "--tau", str(sf_cfg["tau"]),
        "--nsample_per_protein", str(sf_cfg["nsample_per_protein"]),
        "--plddt",
        "--fasta_path", fasta_dir,
        "--output_dir", out_dir,
        "--output_format", "pdb",
    ]

    print(f"[SimpleFold] 预测 {len(candidates)} 个变体...")
    print(f"[SimpleFold] 命令: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"[SimpleFold] 完成 → {out_dir}")
    except FileNotFoundError:
        print("[SimpleFold] 错误: simplefold 命令未找到，请确认 simplefold conda 环境已激活")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"[SimpleFold] 错误:\n{e.stderr}")
        sys.exit(1)


# ── 主函数 ───────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="Phase C Step 8: 突变体结构生成")
    parser.add_argument("config", nargs="?", default="configs/config.yaml",
                        help="配置文件路径 (默认: configs/config.yaml)")
    parser.add_argument("--method", choices=["rosetta", "simplefold"], default=None,
                        help="结构生成方法 (默认: 读取 config)")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    pc = cfg["phase_c"]
    method = args.method or pc["structure_generation"]["primary"]
    pc_dir = pc["paths"]["phase_c_dir"]
    template_pdb = pc["paths"]["template_pdb"]

    candidates = load_candidates(cfg)
    print(f"加载 {len(candidates)} 个候选变体\n")

    if not candidates:
        print("无候选变体，退出")
        return

    if method == "rosetta":
        ros_cfg = pc["structure_generation"]["rosetta"]
        out_dir = os.path.join(pc_dir, "structures", "rosetta")
        build_with_rosetta(
            candidates, template_pdb, out_dir,
            repack_shell=ros_cfg["repack_shell"],
            minimize=ros_cfg["minimize"],
            scorefxn=ros_cfg["scorefxn"],
        )
    elif method == "simplefold":
        sf_cfg = pc["structure_generation"]["simplefold"]
        out_dir = os.path.join(pc_dir, "structures", "simplefold")
        build_with_simplefold(candidates, template_pdb, out_dir, sf_cfg)


if __name__ == "__main__":
    main()
