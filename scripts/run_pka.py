#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 9 (Phase C): pKa 批量预测 (PROPKA3 + pKAI+)

调用 analysis/pka/run_pka.py 的 batch_predict，
然后将 per-His 长表转换为 per-variant 宽表。
"""

import argparse
import glob
import os
import sys

import pandas as pd
import yaml

# 确保可以 import analysis 包
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def main():
    parser = argparse.ArgumentParser(description="Phase C Step 9: pKa 批量预测")
    parser.add_argument("config", nargs="?", default="configs/config.yaml",
                        help="配置文件路径 (默认: configs/config.yaml)")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    pc = cfg.get("tier2", cfg.get("phase_c", {}))
    pc_dir = pc["paths"].get("tier2_dir", pc["paths"].get("phase_c_dir", "tier2"))
    wt_pdb = pc["paths"]["wt_pdb"]
    chains = pc["chains"]

    from analysis.pka.run_pka import batch_predict

    # 结构目录: tier2 模式用 rosetta/; 旧模式从 structure_generation.primary 读取
    if "rosetta" in pc:
        primary = "rosetta"
    elif "structure_generation" in pc:
        primary = pc["structure_generation"]["primary"]
    else:
        primary = "rosetta"
    struct_dir = os.path.join(pc_dir, "structures", primary)

    out_dir = os.path.join(pc_dir, "pka")
    os.makedirs(out_dir, exist_ok=True)

    # His 位点过滤
    his_filter = [
        (hp["chain"], hp["resid"]) for hp in pc["pka"]["his_positions"]
    ]

    # 收集突变体 PDB
    pdb_paths = sorted(glob.glob(os.path.join(struct_dir, "*.pdb")))
    if not pdb_paths:
        print(f"错误: {struct_dir} 下没有 PDB 文件，请先运行 build_structures.py")
        sys.exit(1)

    print(f"== Phase C: pKa 预测 ==")
    print(f"WT: {wt_pdb}")
    print(f"突变体: {len(pdb_paths)} 个")
    print(f"His 过滤: {his_filter}")
    print()

    # 批量预测
    long_df = batch_predict(pdb_paths, wt_pdb=wt_pdb, his_filter=his_filter)

    if long_df.empty:
        print("未检测到目标 His 位点，请检查结构文件和 his_positions 配置")
        sys.exit(1)

    # 保存长表
    long_path = os.path.join(out_dir, "pka_detail.csv")
    long_df.to_csv(long_path, index=False)
    print(f"\n长表已保存: {long_path}")

    # 转宽表: 每个 variant 一行
    _pdb_to_hl = {chains["heavy"]: "H", chains["light"]: "L"}

    wide_rows = []
    for pdb_name, group in long_df.groupby("pdb_name"):
        variant_id = pdb_name.replace(".pdb", "")
        row = {"variant_id": variant_id}

        for _, r in group.iterrows():
            resid = int(r["resid"])
            hl = _pdb_to_hl.get(r.get("chain", ""), r.get("chain", ""))
            prefix = f"{hl}{resid}"
            row[f"{prefix}_propka"] = r.get("pKa_propka")
            row[f"{prefix}_pkai"] = r.get("pKa_pkai")
            if "shift_propka" in r:
                row[f"{prefix}_shift_propka"] = r["shift_propka"]
                row[f"{prefix}_shift_pkai"] = r["shift_pkai"]
                row[f"{prefix}_consensus"] = r.get("consensus", "")

        # 均值 shift
        shift_cols_propka = [c for c in row if c.endswith("_shift_propka")]
        shift_cols_pkai = [c for c in row if c.endswith("_shift_pkai")]
        vals_propka = [row[c] for c in shift_cols_propka if row[c] is not None and not pd.isna(row[c])]
        vals_pkai = [row[c] for c in shift_cols_pkai if row[c] is not None and not pd.isna(row[c])]
        row["avg_shift_propka"] = sum(vals_propka) / len(vals_propka) if vals_propka else None
        row["avg_shift_pkai"] = sum(vals_pkai) / len(vals_pkai) if vals_pkai else None

        # 整体共识
        consensus_cols = [c for c in row if c.endswith("_consensus")]
        consensus_vals = [row[c] for c in consensus_cols if row[c]]
        if all(v == "agree" for v in consensus_vals):
            row["overall_consensus"] = "agree"
        elif any(v == "disagree" for v in consensus_vals):
            row["overall_consensus"] = "disagree"
        else:
            row["overall_consensus"] = "mixed"

        wide_rows.append(row)

    wide_df = pd.DataFrame(wide_rows)
    wide_path = os.path.join(out_dir, "pka_summary.csv")
    wide_df.to_csv(wide_path, index=False)
    print(f"宽表已保存: {wide_path}")
    print(f"\n{wide_df.to_string(index=False)}")


if __name__ == "__main__":
    main()
