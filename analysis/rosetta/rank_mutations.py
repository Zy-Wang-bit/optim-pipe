#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pH敏感型抗体His突变体排序工具

用法:
    python rank_mutations.py <input_csv> <output_csv> [选项]

选项:
    --w_ddd     delta_delta_dG权重 (默认0.25)
    --w_elec    dddG_elec权重 (默认0.25)
    --w_ph      ph_score权重 (默认0.25)
    --w_ddg     ddG_pH7.4权重 (默认0.25)
    
    --d_ddd     delta_delta_dG方向: higher/lower (理论预期: higher)
    --d_elec    dddG_elec方向: higher/lower (理论预期: higher)
    --d_ph      ph_score方向: higher/lower (理论预期: higher)
    --d_ddg     ddG_pH7.4方向: higher/lower (理论预期: lower)
    
    --exp_col   实验结果列名 (默认exp_ratio或Final ratio)

理论预期方向说明:
    - delta_delta_dG: 越高越好 (正值=酸性下结合弱)
    - dddG_elec: 越高越好 (正值=酸性下静电不利)
    - ph_score: 越高越好 (高分=pH敏感潜力强)
    - ddG_pH7.4: 越低越好 (低值=中性pH下结合稳定)

示例:
    python rank_mutations.py data.csv ranked.csv --w_elec 0.5 --w_ddd 0.2 --w_ph 0.2 --w_ddg 0.1
"""

import sys
import argparse
import pandas as pd
import numpy as np
from scipy import stats

def normalize_indicator(series, higher_better=True):
    """归一化指标到[0,1]"""
    min_val, max_val = series.min(), series.max()
    if max_val == min_val:
        return pd.Series([0.5] * len(series), index=series.index)
    norm = (series - min_val) / (max_val - min_val)
    return norm if higher_better else (1 - norm)

def compute_composite_score(df, weights, directions):
    """计算复合评分"""
    indicators = ['delta_delta_dG', 'dddG_elec', 'ph_score', 'ddG_pH7.4']
    score = np.zeros(len(df))
    
    for ind in indicators:
        if ind not in df.columns:
            print(f"[警告] 列 {ind} 不存在，跳过")
            continue
        w = weights.get(ind, 0.25)
        higher_better = (directions.get(ind, 'higher') == 'higher')
        norm = normalize_indicator(df[ind], higher_better)
        score += w * norm
    
    return score

def calculate_hit_rate(df, top_n, exp_col, threshold_percentile=0.3):
    """计算Top-N命中率（命中=实际也在前X%）"""
    n_positive = int(len(df) * threshold_percentile)
    positive_set = set(df.nlargest(n_positive, exp_col).index)
    
    top_n_set = set(df.nlargest(top_n, 'composite_score').index)
    hits = len(positive_set & top_n_set)
    
    return hits, top_n, hits / top_n * 100

def main():
    parser = argparse.ArgumentParser(description='pH敏感型抗体His突变体排序工具')
    parser.add_argument('input_csv', help='输入CSV文件路径')
    parser.add_argument('output_csv', help='输出CSV文件路径')
    
    # 权重参数
    parser.add_argument('--w_ddd', type=float, default=0.25, help='delta_delta_dG权重')
    parser.add_argument('--w_elec', type=float, default=0.25, help='dddG_elec权重')
    parser.add_argument('--w_ph', type=float, default=0.25, help='ph_score权重')
    parser.add_argument('--w_ddg', type=float, default=0.25, help='ddG_pH7.4权重')
    
    # 方向参数 - 默认使用理论预期方向
    parser.add_argument('--d_ddd', default='higher', choices=['higher', 'lower'], 
                        help='delta_delta_dG方向 (理论预期: higher)')
    parser.add_argument('--d_elec', default='higher', choices=['higher', 'lower'],
                        help='dddG_elec方向 (理论预期: higher)')
    parser.add_argument('--d_ph', default='higher', choices=['higher', 'lower'],
                        help='ph_score方向 (理论预期: higher)')
    parser.add_argument('--d_ddg', default='lower', choices=['higher', 'lower'],
                        help='ddG_pH7.4方向 (理论预期: lower)')
    
    # 其他参数
    parser.add_argument('--exp_col', default=None, help='实验结果列名')
    
    args = parser.parse_args()
    
    # 读取数据
    print(f"[读取] {args.input_csv}")
    df = pd.read_csv(args.input_csv)
    print(f"[数据] {len(df)} 行, {len(df.columns)} 列")
    
    # 确定实验结果列
    exp_col = args.exp_col
    if exp_col is None:
        if 'exp_ratio' in df.columns:
            exp_col = 'exp_ratio'
        elif 'Final ratio' in df.columns:
            exp_col = 'Final ratio'
        else:
            print("[警告] 未找到实验结果列，无法计算命中率和相关性")
            exp_col = None
    
    # 构建权重和方向字典
    weights = {
        'delta_delta_dG': args.w_ddd,
        'dddG_elec': args.w_elec,
        'ph_score': args.w_ph,
        'ddG_pH7.4': args.w_ddg
    }
    
    directions = {
        'delta_delta_dG': args.d_ddd,
        'dddG_elec': args.d_elec,
        'ph_score': args.d_ph,
        'ddG_pH7.4': args.d_ddg
    }
    
    # 归一化权重
    total_weight = sum(weights.values())
    if abs(total_weight - 1.0) > 0.01:
        print(f"[归一化] 权重总和={total_weight:.3f}，自动归一化到1.0")
        weights = {k: v/total_weight for k, v in weights.items()}
    
    print("\n[配置]")
    print(f"{'指标':<20} {'权重':>10} {'方向':>10}")
    print("-" * 42)
    for ind in ['delta_delta_dG', 'dddG_elec', 'ph_score', 'ddG_pH7.4']:
        d_str = '越高越好 ↑' if directions[ind] == 'higher' else '越低越好 ↓'
        print(f"{ind:<20} {weights[ind]:>10.4f} {d_str:>10}")
    
    # 计算复合评分
    df['composite_score'] = compute_composite_score(df, weights, directions)
    
    # 排序
    df_sorted = df.sort_values('composite_score', ascending=False).reset_index(drop=True)
    df_sorted['pred_rank'] = range(1, len(df_sorted) + 1)
    
    # 输出
    df_sorted.to_csv(args.output_csv, index=False)
    print(f"\n[输出] {args.output_csv}")
    
    # 评估（如果有实验数据）
    if exp_col and exp_col in df.columns:
        print("\n" + "=" * 60)
        print("[评估结果]")
        print("=" * 60)
        
        # 相关性
        rho, p = stats.spearmanr(df_sorted['composite_score'], df_sorted[exp_col])
        print(f"\n复合评分与{exp_col}的Spearman相关性:")
        print(f"  ρ = {rho:+.4f} (p = {p:.4f})")
        
        if rho > 0:
            print("  → 正相关：评分越高，实验效果越好 ✓")
        else:
            print("  → 负相关：评分越高，实验效果越差 ⚠️")
        
        # 命中率
        print(f"\nTop-N 命中率 (命中=实际也在前30%):")
        print(f"{'Top N':>10} {'命中数':>10} {'命中率':>12} {'随机期望':>12}")
        print("-" * 46)
        
        expected = 30.0  # 随机期望30%
        for top_n in [5, 10, 15, 20]:
            if top_n > len(df_sorted):
                continue
            hits, total, rate = calculate_hit_rate(df_sorted, top_n, exp_col, 0.3)
            better = "✓" if rate > expected else "✗"
            print(f"{top_n:>10} {hits:>10} {rate:>11.1f}% {expected:>11.1f}% {better}")
        
        # 各单一指标的相关性
        print(f"\n各单一指标与{exp_col}的相关性:")
        print(f"{'指标':<20} {'ρ':>10} {'p值':>12} {'配置方向':>12} {'实际一致':>10}")
        print("-" * 66)
        
        for ind in ['delta_delta_dG', 'dddG_elec', 'ph_score', 'ddG_pH7.4']:
            if ind not in df_sorted.columns:
                continue
            r, p = stats.spearmanr(df_sorted[ind], df_sorted[exp_col])
            config_dir = directions[ind]
            # 实际方向：正相关=higher好，负相关=lower好
            actual_dir = 'higher' if r > 0 else 'lower'
            match = "✓" if config_dir == actual_dir else "✗"
            sig = "*" if p < 0.05 else ""
            print(f"{ind:<20} {r:>+10.3f}{sig} {p:>12.4f} {config_dir:>12} {match:>10}")
    
    print("\n[完成]")

if __name__ == '__main__':
    main()
