#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
后筛选：基于 R2 ELISA SAR 规则过滤 1E62 候选变体。

输入:  final_candidates.csv (pipeline 输出)
输出:  sar_filtered_candidates.csv

规则来源: experiments/1E62/R2/wet_lab/analysis_report.md
"""
import argparse
import os
import re
import sys

import pandas as pd

# ── SAR 规则定义 ────────────────────────────────────────────────────────────

# 硬排除：L41(W) 突变 → 100% 失活
HARD_EXCLUDE_POSITIONS = {
    ("B", 41),  # L41(W) 绝对保留
}

# 高风险位点（软标记）
HIGH_RISK_POSITIONS = {
    ("B", 22),   # L22 高频失活
    ("B", 48),   # L48 高频失活
    ("B", 49),   # L49 高频失活
    ("B", 50),   # L50 高频失活
    ("A", 101),  # H101 可能影响结合活性
    ("A", 32),   # H32 可能影响结合活性
    ("A", 102),  # H102 可能影响结合活性
    ("A", 35),   # H35 可能影响结合活性
    ("A", 103),  # H103 可能影响结合活性
    ("B", 40),   # L40 可能影响结合活性
    ("B", 61),   # L61 可能影响结合活性
}

# R2 已验证的安全/有利位点
VALIDATED_SAFE = {
    ("B", 42, "H"),  # S42H - com8/9/18 共享
    ("B", 43, "H"),  # Q43H - com8/9/18 共享
    ("B", 44, "H"),  # Q44H - com8/9/18 共享
    ("B", 51, "Q"),  # K51Q - com8/9/18 共享
}


def parse_foldx_mutation(token):
    """解析 FoldX 格式突变: FA102H → (chain='A', pos=102, wt='F', mut='H')"""
    token = token.strip()
    if len(token) < 4:
        return None
    wt_aa = token[0]
    mut_aa = token[-1]
    middle = token[1:-1]
    # chain letter + position number
    chain = middle[0]
    pos = int(middle[1:])
    return (chain, pos, wt_aa, mut_aa)


def apply_sar_filter(df, mutations_col="mutations"):
    """应用 SAR 规则，返回带标记的 DataFrame。"""
    results = []

    for idx, row in df.iterrows():
        mut_str = str(row.get(mutations_col, ""))
        if not mut_str or mut_str == "nan":
            results.append({
                "sar_excluded": False,
                "sar_risk_flags": "",
                "sar_validated_overlap": 0,
            })
            continue

        tokens = mut_str.replace(";", ",").split(",")
        mutations = []
        for t in tokens:
            parsed = parse_foldx_mutation(t)
            if parsed:
                mutations.append(parsed)

        # 硬排除检查
        excluded = False
        for chain, pos, wt, mut in mutations:
            if (chain, pos) in HARD_EXCLUDE_POSITIONS:
                excluded = True
                break

        # 高风险标记
        risk_flags = []
        for chain, pos, wt, mut in mutations:
            if (chain, pos) in HIGH_RISK_POSITIONS:
                risk_flags.append(f"{chain}{pos}")

        # 已验证位点重叠度
        validated_count = sum(
            1 for chain, pos, wt, mut in mutations
            if (chain, pos, mut) in VALIDATED_SAFE
        )

        results.append({
            "sar_excluded": excluded,
            "sar_risk_flags": ";".join(risk_flags),
            "sar_validated_overlap": validated_count,
        })

    sar_df = pd.DataFrame(results, index=df.index)
    return pd.concat([df, sar_df], axis=1)


def main():
    parser = argparse.ArgumentParser(description="1E62 SAR 后筛选")
    parser.add_argument("input_csv", help="pipeline 输出的候选 CSV")
    parser.add_argument("-o", "--output", default=None,
                        help="输出路径 (默认: 同目录下 sar_filtered_*.csv)")
    parser.add_argument("--mutations-col", default="mutations",
                        help="突变列名 (默认: mutations)")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    print(f"[SAR] 输入: {args.input_csv} ({len(df)} 行)")

    df = apply_sar_filter(df, args.mutations_col)

    # 统计
    n_excluded = df["sar_excluded"].sum()
    n_risky = (df["sar_risk_flags"] != "").sum()
    n_validated = (df["sar_validated_overlap"] > 0).sum()
    print(f"[SAR] 硬排除: {n_excluded}, 高风险标记: {n_risky}, 含已验证位点: {n_validated}")

    # 过滤掉硬排除
    passed = df[~df["sar_excluded"]].copy()
    print(f"[SAR] 通过: {len(passed)} / {len(df)}")

    # 输出
    if args.output:
        out_path = args.output
    else:
        base = os.path.basename(args.input_csv)
        out_path = os.path.join(os.path.dirname(args.input_csv),
                                f"sar_filtered_{base}")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    passed.to_csv(out_path, index=False)
    print(f"[SAR] 输出: {out_path}")


if __name__ == "__main__":
    main()
