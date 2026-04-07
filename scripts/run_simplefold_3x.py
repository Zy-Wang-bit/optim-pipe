#!/usr/bin/env python3
"""
Step 10 (Tier 2): SimpleFold 3x 采样

从 tier1_candidates.csv 读取序列，用 SimpleFold 3B 每变体 3 次采样。
输出到 {tier2_dir}/structures/simplefold/

独立于 PyRosetta 线运行：
- 使用单抗体序列（不含抗原）
- 3x 采样用于后续 RMSD 中位数计算
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


def _extract_antibody_sequence(row, cfg):
    """从候选行提取单抗体序列。

    候选 CSV 的 sequence 列可能是多链格式 (A/B/C)，
    需要去掉抗原链只保留抗体链。
    """
    seq = row.get("sequence", "")
    if not seq or (isinstance(seq, float) and pd.isna(seq)):
        return None

    t2 = cfg.get("tier2", cfg.get("phase_c", {}))
    chains_cfg = t2["chains"]
    ab_chains = [chains_cfg["heavy"], chains_cfg["light"]]

    # 多链格式: "HEAVY_SEQ/LIGHT_SEQ/ANTIGEN_SEQ"
    parts = str(seq).split("/")
    if len(parts) >= 2:
        return "/".join(parts[:len(ab_chains)])
    return seq


def main():
    parser = argparse.ArgumentParser(description="Tier 2 Step 10: SimpleFold 3x 采样")
    parser.add_argument("config", nargs="?", default="configs/config.yaml")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    t2 = cfg.get("tier2", cfg.get("phase_c", {}))
    sf_cfg = t2["simplefold"]
    t2_dir = t2["paths"].get("tier2_dir", t2["paths"].get("phase_c_dir", "tier2"))
    input_csv = t2["input"]["csv"]

    out_dir = os.path.join(t2_dir, "structures", "simplefold")
    fasta_dir = os.path.join(out_dir, "fasta_input")
    os.makedirs(fasta_dir, exist_ok=True)

    df = pd.read_csv(input_csv)
    print(f"[SimpleFold 3x] 读取 {len(df)} 个候选")

    # 生成 FASTA 文件
    n_fasta = 0
    for idx, row in df.iterrows():
        if "mpdb" in row and pd.notna(row["mpdb"]):
            vid = Path(str(row["mpdb"])).stem
        else:
            vid = f"var_{idx:06d}"

        seq = _extract_antibody_sequence(row, cfg)
        if not seq:
            print(f"  跳过 {vid}: 无序列")
            continue

        fasta_path = os.path.join(fasta_dir, f"{vid}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{vid}\n{seq}\n")
        n_fasta += 1

    print(f"[SimpleFold 3x] 生成 {n_fasta} 个 FASTA 文件")

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

    print(f"[SimpleFold 3x] 命令: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"[SimpleFold 3x] 完成 → {out_dir}")
    except FileNotFoundError:
        print("[SimpleFold 3x] 错误: simplefold 命令未找到，请确认 simplefold conda 环境已激活")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"[SimpleFold 3x] 错误:\n{e.stderr}")
        sys.exit(1)


if __name__ == "__main__":
    main()
