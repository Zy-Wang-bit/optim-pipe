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


def evaluate_simplefold_3x(wt_pdb, struct_dir, cdr_regions, outlier_threshold=2.0):
    """评估 SimpleFold 3x 采样的 CDR RMSD。

    每个变体有多个样本，剔除 global RMSD > outlier_threshold 的异常样本，
    取各 CDR RMSD 的中位数。
    """
    import re
    from collections import defaultdict

    all_pdbs = sorted(glob.glob(os.path.join(struct_dir, "*.pdb")))
    if not all_pdbs:
        return []

    # 按变体 ID 分组
    variant_pdbs = defaultdict(list)
    for p in all_pdbs:
        stem = Path(p).stem
        # SimpleFold 输出格式: {variant_id}_sampled_{N} 或 {variant_id}_sample_{N} 或 {variant_id}_{N}
        match = re.match(r"(.+?)(?:_sampled|_sample)?_(\d+)$", stem)
        if match:
            vid = match.group(1)
            variant_pdbs[vid].append(p)
        else:
            variant_pdbs[stem].append(p)

    print(f"[RMSD 3x] {len(variant_pdbs)} 个变体，共 {len(all_pdbs)} 个 PDB")

    results = []
    for vid, pdbs in variant_pdbs.items():
        sample_results = []
        for pdb_path in pdbs:
            try:
                r = evaluate_structure_set(wt_pdb, [pdb_path], cdr_regions, "simplefold_3x")
                if r:
                    sample_results.append(r[0])
            except Exception as e:
                print(f"  {vid} 样本 {pdb_path} 失败: {e}")

        if not sample_results:
            continue

        # 剔除 global RMSD 异常值
        valid = [r for r in sample_results
                 if r.get("global_rmsd") is not None and r["global_rmsd"] <= outlier_threshold]
        if not valid:
            print(f"  {vid}: 所有 {len(sample_results)} 样本均为异常值 (global RMSD > {outlier_threshold}Å)")
            continue

        # 取中位数
        median_result = {"variant_id": vid, "source": "simplefold_3x",
                         "n_valid_samples": len(valid)}
        rmsd_keys = [k for k in valid[0] if k.endswith("_rmsd")]
        for k in rmsd_keys:
            vals = [r[k] for r in valid if r.get(k) is not None]
            median_result[k] = float(np.median(vals)) if vals else None

        results.append(median_result)

    return results


def main():
    parser = argparse.ArgumentParser(description="Step 11: CDR RMSD 评估")
    parser.add_argument("config", nargs="?", default="configs/config.yaml",
                        help="配置文件路径 (默认: configs/config.yaml)")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    pc = cfg.get("tier2", cfg.get("phase_c", {}))
    wt_pdb = pc["paths"].get("wt_ab_pdb", pc["paths"].get("template_pdb"))
    pc_dir = pc["paths"].get("tier2_dir", pc["paths"].get("phase_c_dir", "tier2"))
    out_dir = os.path.join(pc_dir, "rmsd")
    os.makedirs(out_dir, exist_ok=True)

    cdr_regions = pc["cdr_regions"]

    print("== CDR RMSD 评估 ==")
    print(f"WT 参考: {wt_pdb}")
    print(f"CDR 区域: {list(cdr_regions.keys())}\n")

    all_results = []

    # PyRosetta 线结构
    rosetta_dir = os.path.join(pc_dir, "structures", "rosetta")
    rosetta_pdbs = sorted(glob.glob(os.path.join(rosetta_dir, "*.pdb")))
    if rosetta_pdbs:
        print(f"[rosetta] {len(rosetta_pdbs)} 个结构")
        all_results.extend(
            evaluate_structure_set(wt_pdb, rosetta_pdbs, cdr_regions, "rosetta")
        )

    # SimpleFold 3x 线结构
    sf_dir = os.path.join(pc_dir, "structures", "simplefold")
    if os.path.isdir(sf_dir) and glob.glob(os.path.join(sf_dir, "*.pdb")):
        outlier_th = pc.get("simplefold", {}).get("outlier_global_rmsd", 2.0)
        sf3x_results = evaluate_simplefold_3x(wt_pdb, sf_dir, cdr_regions, outlier_th)
        all_results.extend(sf3x_results)
        print(f"[SimpleFold 3x] {len(sf3x_results)} 变体 (中位数)")

    # 旧模式兼容: structure_generation.primary/auxiliary
    sg = pc.get("structure_generation", {})
    if sg and not rosetta_pdbs and not os.path.isdir(sf_dir):
        methods = [sg.get("primary")]
        if sg.get("auxiliary"):
            methods.append(sg["auxiliary"])
        for method in methods:
            if not method:
                continue
            struct_dir = os.path.join(pc_dir, "structures", method)
            pdbs = sorted(glob.glob(os.path.join(struct_dir, "*.pdb")))
            if pdbs:
                print(f"[{method}] {len(pdbs)} 个结构")
                all_results.extend(
                    evaluate_structure_set(wt_pdb, pdbs, cdr_regions, method)
                )

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
