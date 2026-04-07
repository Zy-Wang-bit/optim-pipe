#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 10 (Phase C): Rosetta 评分 (pH-score + dddG_elec)

整合 analysis/rosetta/ 下的 pHScoreCalculator 和 dddGElecCalculator，
批量计算突变体的 pH-score 和 dddG_elec。
需在 pyrosetta conda 环境中运行。
"""

import argparse
import glob
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

# 确保可以 import analysis 包
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)


def run_ph_score(pdb_paths, binder_chains, target_chains, out_path):
    """批量计算 pH-score（抗体内部稳定性版本）。"""
    from analysis.rosetta.calc_ph_score_antibody import pHScoreCalculatorInternal as pHScoreCalculator

    results = []

    for i, pdb_path in enumerate(pdb_paths, 1):
        pdb_name = os.path.basename(pdb_path)
        variant_id = Path(pdb_path).stem
        print(f"  [pH-score {i}/{len(pdb_paths)}] {pdb_name}...", end=" ")
        try:
            calc = pHScoreCalculator(
                pdb_path=pdb_path,
                binder_chains=binder_chains,
                target_chains=target_chains,
            )
            result = calc.calculate_ph_score()
            row = {
                "variant_id": variant_id,
                "ph_score": result["ph_score"],
                "n_his_binder": result["n_histidines_binder"],
                "status": "success",
            }
            for term, count in result["terms"].items():
                row[term] = count
            print(f"score={result['ph_score']:.1f}")
        except Exception as e:
            print(f"失败: {e}")
            row = {
                "variant_id": variant_id,
                "ph_score": None,
                "n_his_binder": None,
                "status": f"error: {e}",
            }
        results.append(row)

    df = pd.DataFrame(results)
    df.to_csv(out_path, index=False)
    return df


def run_dddg_elec(pdb_paths, binder_chains, target_chains, out_path):
    """批量计算 dddG_elec。"""
    from analysis.rosetta.batch_calc_dddg_elec import dddGElecCalculator

    import pyrosetta
    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init("-mute all")

    results = []
    for i, pdb_path in enumerate(pdb_paths, 1):
        pdb_name = os.path.basename(pdb_path)
        variant_id = Path(pdb_path).stem
        print(f"  [dddG_elec {i}/{len(pdb_paths)}] {pdb_name}...", end=" ")
        try:
            calc = dddGElecCalculator(
                pdb_path=pdb_path,
                binder_chains=binder_chains,
                target_chains=target_chains,
            )
            result = calc.calculate_dddg_elec()
            row = {
                "variant_id": variant_id,
                "dddG_elec": result["dddG_elec"],
                "ddG_elec_pH7": result["ddG_elec_pH7"],
                "ddG_elec_pH5": result["ddG_elec_pH5"],
                "n_his_binder": result.get("n_histidines_binder"),
                "status": "success" if result["error"] is None else f"error: {result['error']}",
            }
            if result["dddG_elec"] is not None:
                print(f"dddG_elec={result['dddG_elec']:.3f}")
            else:
                print(f"失败: {result['error']}")
        except Exception as e:
            print(f"失败: {e}")
            row = {
                "variant_id": variant_id,
                "dddG_elec": None,
                "ddG_elec_pH7": None,
                "ddG_elec_pH5": None,
                "n_his_binder": None,
                "status": f"error: {e}",
            }
        results.append(row)

    df = pd.DataFrame(results)
    df.to_csv(out_path, index=False)
    return df


def main():
    parser = argparse.ArgumentParser(description="Phase C Step 10: Rosetta 评分")
    parser.add_argument("config", nargs="?", default="configs/config.yaml",
                        help="配置文件路径 (默认: configs/config.yaml)")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    pc = cfg.get("tier2", cfg.get("phase_c", {}))
    pc_dir = pc["paths"].get("tier2_dir", pc["paths"].get("phase_c_dir", "tier2"))

    if "rosetta" in pc:
        primary = "rosetta"
    elif "structure_generation" in pc:
        primary = pc["structure_generation"]["primary"]
    else:
        primary = "rosetta"
    struct_dir = os.path.join(pc_dir, "structures", primary)
    out_dir = os.path.join(pc_dir, "rosetta")
    os.makedirs(out_dir, exist_ok=True)

    binder_chains = pc["chains"]["binder"]
    target_chains = pc["chains"]["target"]

    pdb_paths = sorted(glob.glob(os.path.join(struct_dir, "*.pdb")))
    if not pdb_paths:
        print(f"错误: {struct_dir} 下没有 PDB 文件")
        sys.exit(1)

    print(f"== Phase C: Rosetta 评分 ==")
    print(f"变体: {len(pdb_paths)} 个")
    print(f"Binder chains: {binder_chains}")
    print(f"Target chains: {target_chains}")
    print()

    # pH-score
    print("[1/2] pH-score 计算...")
    ph_path = os.path.join(out_dir, "ph_scores.csv")
    run_ph_score(pdb_paths, binder_chains, target_chains, ph_path)
    print(f"  → {ph_path}\n")

    # dddG_elec
    print("[2/2] dddG_elec 计算...")
    elec_path = os.path.join(out_dir, "dddg_elec.csv")
    run_dddg_elec(pdb_paths, binder_chains, target_chains, elec_path)
    print(f"  → {elec_path}\n")

    print("== Rosetta 评分完成 ==")


if __name__ == "__main__":
    main()
