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


def extract_ca_atoms(structure):
    """提取所有链的 Ca 原子，返回 {(chain_id, resseq): atom}"""
    ca_dict = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CA" in residue:
                    key = (chain.id, residue.id[1])
                    ca_dict[key] = residue["CA"]
        break
    return ca_dict


def calc_rmsd_subset(wt_ca, mut_ca, residue_keys):
    """计算指定残基子集的 RMSD（先全局对齐再计算子集 RMSD）。"""
    common_keys = sorted(set(wt_ca.keys()) & set(mut_ca.keys()))
    if len(common_keys) < 3:
        return float("nan")

    sup = Superimposer()
    wt_atoms = [wt_ca[k] for k in common_keys]
    mut_atoms = [mut_ca[k] for k in common_keys]
    sup.set_atoms(wt_atoms, mut_atoms)

    sup.apply([mut_ca[k] for k in mut_ca])

    subset_keys = [k for k in residue_keys if k in wt_ca and k in mut_ca]
    if not subset_keys:
        return float("nan")

    diffs = []
    for k in subset_keys:
        diff = wt_ca[k].coord - mut_ca[k].coord
        diffs.append(np.sum(diff ** 2))

    return np.sqrt(np.mean(diffs))


def calc_global_rmsd(wt_ca, mut_ca):
    """计算全局 Ca RMSD。"""
    common_keys = sorted(set(wt_ca.keys()) & set(mut_ca.keys()))
    if len(common_keys) < 3:
        return float("nan")

    sup = Superimposer()
    sup.set_atoms(
        [wt_ca[k] for k in common_keys],
        [mut_ca[k] for k in common_keys],
    )
    return sup.rms


def get_cdr_keys(cdr_regions):
    """从 config 的 CDR 定义生成残基 key 列表。"""
    cdr_keys = {}
    for cdr_name, spec in cdr_regions.items():
        chain = spec["chain"]
        keys = [(chain, r) for r in range(spec["start"], spec["end"] + 1)]
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

            row = {
                "variant_id": name,
                "source": source_label,
                "global_rmsd": calc_global_rmsd(wt_ca, mut_ca),
            }

            for cdr_name, keys in cdr_keys.items():
                row[f"{cdr_name.lower()}_rmsd"] = calc_rmsd_subset(wt_ca, mut_ca, keys)

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
