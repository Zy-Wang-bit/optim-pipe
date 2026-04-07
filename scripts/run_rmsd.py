#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 11 (Phase C): CDR RMSD 结构稳定性评估

计算突变体与 WT 模板的全局 RMSD 和 CDR 区域 RMSD。
支持 PyRosetta 和 SimpleFold 两种来源的结构。
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio.PDB import PDBParser, Superimposer

# 复用 analysis/ 中的 Ca 提取函数（返回 (chain, resseq, icode) 三元组 key）
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.structure_compare.rmsd_ca_global import extract_ca_atoms


def _align_and_compute_rmsds(wt_ca, mut_ca, cdr_keys):
    """一次全局叠合，同时计算全局 RMSD 和各 CDR 子集 RMSD。"""
    common_keys = sorted(set(wt_ca.keys()) & set(mut_ca.keys()))
    if len(common_keys) < 3:
        result = {"global_rmsd": float("nan")}
        for cdr_name in cdr_keys:
            result[f"{cdr_name.lower()}_rmsd"] = float("nan")
        return result

    # 一次全局叠合
    sup = Superimposer()
    wt_atoms = [wt_ca[k] for k in common_keys]
    mut_atoms = [mut_ca[k] for k in common_keys]
    sup.set_atoms(wt_atoms, mut_atoms)

    result = {"global_rmsd": sup.rms}

    # 应用变换到所有突变体原子（一次性）
    sup.apply([mut_ca[k] for k in mut_ca])

    # 从已对齐坐标计算各 CDR 子集 RMSD
    common_set = set(common_keys)
    for cdr_name, keys in cdr_keys.items():
        subset_keys = [k for k in keys if k in common_set]
        if not subset_keys:
            result[f"{cdr_name.lower()}_rmsd"] = float("nan")
            continue
        diffs = np.array([wt_ca[k].coord - mut_ca[k].coord for k in subset_keys])
        result[f"{cdr_name.lower()}_rmsd"] = np.sqrt(np.mean(np.sum(diffs ** 2, axis=1)))

    return result


def get_cdr_keys(cdr_regions):
    """从 config 的 CDR 定义生成残基 key 列表（三元组格式含 icode）。"""
    cdr_keys = {}
    for cdr_name, spec in cdr_regions.items():
        chain = spec["chain"]
        keys = [(chain, r, " ") for r in range(spec["start"], spec["end"] + 1)]
        cdr_keys[cdr_name] = keys
    return cdr_keys


def evaluate_structure_set(wt_pdb, pdb_paths, cdr_regions, source_label):
    """评估一组结构的 RMSD。"""
    parser = PDBParser(QUIET=True)
    wt_struct = parser.get_structure("wt", wt_pdb)
    wt_ca = extract_ca_atoms(wt_struct)
    cdr_keys = get_cdr_keys(cdr_regions)

    results = []
    for pdb_path in pdb_paths:
        name = Path(pdb_path).stem
        print(f"  [{source_label}] {name}...", end=" ")

        try:
            mut_struct = parser.get_structure("mut", pdb_path)
            mut_ca = extract_ca_atoms(mut_struct)

            rmsds = _align_and_compute_rmsds(wt_ca, mut_ca, cdr_keys)
            row = {"variant_id": name, "source": source_label}
            row.update(rmsds)

            print(f"global={row['global_rmsd']:.3f}")
        except Exception as e:
            print(f"失败: {e}")
            row = {
                "variant_id": name,
                "source": source_label,
                "global_rmsd": None,
            }
            for cdr_name in cdr_keys:
                row[f"{cdr_name.lower()}_rmsd"] = None

        results.append(row)

    return results


def main():
    parser = argparse.ArgumentParser(description="Phase C Step 11: CDR RMSD 评估")
    parser.add_argument("config", nargs="?", default="configs/config.yaml",
                        help="配置文件路径 (默认: configs/config.yaml)")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    pc = cfg["phase_c"]
    wt_pdb = pc["paths"]["template_pdb"]
    pc_dir = pc["paths"]["phase_c_dir"]
    out_dir = os.path.join(pc_dir, "rmsd")
    os.makedirs(out_dir, exist_ok=True)

    cdr_regions = pc["cdr_regions"]

    print("== Phase C: CDR RMSD 评估 ==")
    print(f"WT 参考: {wt_pdb}")
    print(f"CDR 区域: {list(cdr_regions.keys())}\n")

    all_results = []

    # 评估 primary 方法的结构
    methods = [pc["structure_generation"]["primary"]]
    if pc["structure_generation"].get("auxiliary"):
        methods.append(pc["structure_generation"]["auxiliary"])

    for method in methods:
        struct_dir = os.path.join(pc_dir, "structures", method)
        pdbs = sorted(glob.glob(os.path.join(struct_dir, "*.pdb")))
        if pdbs:
            print(f"[{method}] {len(pdbs)} 个结构")
            all_results.extend(
                evaluate_structure_set(wt_pdb, pdbs, cdr_regions, method)
            )
        else:
            print(f"[{method}] 无结构文件")

    if not all_results:
        print("错误: 没有结构文件可评估")
        sys.exit(1)

    df = pd.DataFrame(all_results)
    out_path = os.path.join(out_dir, "rmsd_summary.csv")
    df.to_csv(out_path, index=False)

    print(f"\n== 完成 ==")
    print(f"结果: {out_path}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
