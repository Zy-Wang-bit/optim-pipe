#!/usr/bin/env python3
"""
Wet-Lab ELISA 数据 vs 计算预测的全面关联分析
===========================================
输入：
  - 3 份 ELISA CSV (基因 Ae, B, D1)  → 20 个 variant × pH 6.0/7.4 × 6 浓度 × 3 重复
  - top20_by_score.csv               → 计算 composite_score, dddG_elec, ddG_pH7.4
  - 设计序列映射 CSV                  → com# ↔ mut_id

输出 (到 analysis/ 目录)：
  1. elisa_summary.csv         — 每个 variant 的实验汇总指标
  2. correlation_report.md     — 计算 vs 实验的相关性报告
  3. elisa_heatmap.png         — 热图可视化
  4. correlation_scatter.png   — 散点图
"""

import csv
import os
import sys
import math
from collections import defaultdict
from pathlib import Path

# ─── 路径 ─────────────────────────────────────────────
BASE = Path("/public/home/ziyang/code/optim-pipe")
WETLAB = BASE / "experiments" / "1E62_R2" / "wet_lab"
EXP_DIR = BASE / "experiments" / "1E62_R2"
OUT_DIR = BASE / "analysis"
OUT_DIR.mkdir(exist_ok=True)

ELISA_FILES = {
    "Ae": WETLAB / "20260313-1E62-系列改造-基因Ae.csv",
    "B":  WETLAB / "20260313-1E62-系列改造-基因B.csv",
    "D1": WETLAB / "20260313-1E62-系列改造-基因D1.csv",
}

DESIGN_FILE = WETLAB / "20260105 AI设计1E62变体序列.csv"
SCORE_FILE  = EXP_DIR / "1E62_R2_top20_by_score.csv"

CONCENTRATIONS = [500, 100, 20, 4, 0.8, 0.16]  # ng/mL


# ─── 1. 解析设计映射 ──────────────────────────────────
def parse_design_mapping():
    """返回 {com#: mut_id} 映射"""
    mapping = {}
    with open(DESIGN_FILE, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 7:
                mut_id = row[0].strip()
                com_label = row[6].strip()  # e.g. "com1"
                mapping[com_label] = mut_id
    return mapping


# ─── 2. 解析 ELISA CSV ───────────────────────────────
def parse_elisa_csv(filepath):
    """
    返回: {com#: {"pH6": [[rep1_vals], [rep2_vals], [rep3_vals]],
                  "pH7": [[rep1_vals], [rep2_vals], [rep3_vals]]}}
    每个 rep_vals 是长度 6 的 float list (对应 6 个浓度)
    """
    data = defaultdict(lambda: {"pH6": [], "pH7": []})
    
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        rows = list(reader)
    
    # 前 3 行是 header: [antigen], [pH], [concentrations]
    # 数据从第 4 行开始，每 3 行是一个 com 的 3 个重复
    for i in range(3, len(rows)):
        row = rows[i]
        if not row or not row[0].strip():
            continue
        
        com_label = row[0].strip()
        try:
            ph6_vals = [float(x) for x in row[1:7]]
            ph7_vals = [float(x) for x in row[7:13]]
        except (ValueError, IndexError):
            continue
        
        data[com_label]["pH6"].append(ph6_vals)
        data[com_label]["pH7"].append(ph7_vals)
    
    return dict(data)


# ─── 3. 计算实验指标 ──────────────────────────────────
def compute_auc(conc, od_vals):
    """梯形法计算 OD-浓度曲线下面积 (log-scale concentration)"""
    if len(conc) != len(od_vals):
        return 0.0
    # 使用 log10(conc) 作为 x 轴
    log_conc = [math.log10(c) for c in conc]
    auc = 0.0
    for i in range(len(log_conc) - 1):
        dx = log_conc[i] - log_conc[i+1]  # log_conc 是递减的
        avg_y = (od_vals[i] + od_vals[i+1]) / 2
        auc += dx * avg_y
    return auc


def compute_variant_metrics(elisa_data_per_antigen):
    """
    给定一个 antigen 的完整 ELISA data，计算每个 com 的关键指标：
    - auc_ph6: pH 6.0 下的 AUC (平均 3 重复)
    - auc_ph7: pH 7.4 下的 AUC (平均 3 重复)
    - ph_ratio: auc_ph6 / auc_ph7 (越小 = pH 切换越好)
    - od_max_ph7: pH 7.4, 500 ng/mL 下最大 OD (平均值)
    - binding_class: "strong" / "weak" / "none" 基于最大 OD
    """
    results = {}
    
    for com, ph_data in elisa_data_per_antigen.items():
        ph6_reps = ph_data["pH6"]
        ph7_reps = ph_data["pH7"]
        
        if not ph6_reps or not ph7_reps:
            continue
        
        # AUC for each replicate
        auc6_list = [compute_auc(CONCENTRATIONS, rep) for rep in ph6_reps]
        auc7_list = [compute_auc(CONCENTRATIONS, rep) for rep in ph7_reps]
        
        auc6_mean = sum(auc6_list) / len(auc6_list)
        auc7_mean = sum(auc7_list) / len(auc7_list)
        auc6_std = (sum((x - auc6_mean)**2 for x in auc6_list) / max(len(auc6_list)-1, 1))**0.5
        auc7_std = (sum((x - auc7_mean)**2 for x in auc7_list) / max(len(auc7_list)-1, 1))**0.5
        
        # pH ratio
        ph_ratio = auc6_mean / auc7_mean if auc7_mean > 0.05 else float('nan')
        
        # Max OD at pH 7.4, highest conc
        od_max_ph7_vals = [rep[0] for rep in ph7_reps]  # idx 0 = 500 ng/mL
        od_max_ph7 = sum(od_max_ph7_vals) / len(od_max_ph7_vals)
        
        # Max OD at pH 6.0, highest conc
        od_max_ph6_vals = [rep[0] for rep in ph6_reps]
        od_max_ph6 = sum(od_max_ph6_vals) / len(od_max_ph6_vals)
        
        # Binding classification
        if od_max_ph7 > 1.0:
            binding_class = "strong"
        elif od_max_ph7 > 0.1:
            binding_class = "weak"
        else:
            binding_class = "none"
        
        # pH-switch quality: ideal = strong at pH7.4, weak at pH6.0
        if binding_class == "strong" and not math.isnan(ph_ratio) and ph_ratio < 0.7:
            ph_switch = "GOOD"
        elif binding_class == "strong" and not math.isnan(ph_ratio) and ph_ratio < 0.9:
            ph_switch = "MODERATE"
        elif binding_class == "strong":
            ph_switch = "NONE"
        else:
            ph_switch = "N/A"
        
        results[com] = {
            "auc_ph6": round(auc6_mean, 4),
            "auc_ph7": round(auc7_mean, 4),
            "auc_ph6_std": round(auc6_std, 4),
            "auc_ph7_std": round(auc7_std, 4),
            "ph_ratio": round(ph_ratio, 4) if not math.isnan(ph_ratio) else "N/A",
            "od_max_ph7": round(od_max_ph7, 4),
            "od_max_ph6": round(od_max_ph6, 4),
            "binding_class": binding_class,
            "ph_switch": ph_switch,
        }
    
    return results


# ─── 4. 解析计算预测 ──────────────────────────────────
def parse_scores():
    """返回 {mut_id: {composite_score, dddG_elec, ddG_pH7.4, mutations}}"""
    scores = {}
    with open(SCORE_FILE, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mut_id = row["mut_id"].strip()
            if not mut_id:
                continue
            scores[mut_id] = {
                "composite_score": float(row["综合评分"]),
                "dddG_elec": float(row["dddG_elec"]),
                "ddG_pH7.4": float(row["ddG_pH7.4"]),
                "mutations": row["具体突变"],
            }
    return scores


# ─── 5. Spearman 相关系数（无需 scipy）────────────────
def spearman_rank_correlation(x, y):
    """计算 Spearman ρ（无需 scipy）"""
    n = len(x)
    if n < 3:
        return float('nan'), float('nan')
    
    def rankdata(vals):
        indexed = sorted(enumerate(vals), key=lambda t: t[1])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j < n - 1 and indexed[j+1][1] == indexed[j][1]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1  # 1-based average rank
            for k in range(i, j + 1):
                ranks[indexed[k][0]] = avg_rank
            i = j + 1
        return ranks
    
    rx = rankdata(x)
    ry = rankdata(y)
    
    mean_rx = sum(rx) / n
    mean_ry = sum(ry) / n
    
    num = sum((rx[i] - mean_rx) * (ry[i] - mean_ry) for i in range(n))
    den_x = sum((rx[i] - mean_rx)**2 for i in range(n)) ** 0.5
    den_y = sum((ry[i] - mean_ry)**2 for i in range(n)) ** 0.5
    
    if den_x == 0 or den_y == 0:
        return float('nan'), float('nan')
    
    rho = num / (den_x * den_y)
    
    # t-test for significance
    if abs(rho) < 1.0:
        t = rho * math.sqrt((n - 2) / (1 - rho**2))
        # Approximate p-value from t-distribution (n-2 df)
        # Using a rough approximation for small n
        p_val = 2 * (1 - t_cdf_approx(abs(t), n - 2))
    else:
        p_val = 0.0
    
    return rho, p_val


def t_cdf_approx(t, df):
    """Approximate t-distribution CDF using normal approximation for df > 3"""
    if df <= 0:
        return 0.5
    # Approximation from Abramowitz & Stegun
    x = t * (1 - 1/(4*df)) / math.sqrt(1 + t**2/(2*df))
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


# ─── 6. 主分析 ────────────────────────────────────────
def main():
    print("=" * 60)
    print("1E62_R2 Wet-Lab ELISA vs 计算预测 关联分析")
    print("=" * 60)
    
    # 1. 加载数据
    com_to_mut = parse_design_mapping()
    mut_to_com = {v: k for k, v in com_to_mut.items()}
    scores = parse_scores()
    
    print(f"\n已加载: {len(com_to_mut)} 个 variant 映射")
    print(f"已加载: {len(scores)} 个计算预测分数")
    
    # 2. 解析所有抗原的 ELISA 数据
    all_antigen_data = {}
    all_antigen_metrics = {}
    for ag_name, filepath in ELISA_FILES.items():
        print(f"\n解析 ELISA: {ag_name} ...")
        data = parse_elisa_csv(filepath)
        metrics = compute_variant_metrics(data)
        all_antigen_data[ag_name] = data
        all_antigen_metrics[ag_name] = metrics
        
        # 统计
        strong = sum(1 for v in metrics.values() if v["binding_class"] == "strong")
        weak = sum(1 for v in metrics.values() if v["binding_class"] == "weak")
        none_ = sum(1 for v in metrics.values() if v["binding_class"] == "none")
        print(f"  结合分类: strong={strong}, weak={weak}, none={none_}")
    
    # 3. 输出汇总 CSV
    summary_path = OUT_DIR / "elisa_summary.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        header = ["com_id", "mut_id", "mutations", "composite_score", "dddG_elec", "ddG_pH7.4"]
        for ag in ELISA_FILES:
            header.extend([
                f"{ag}_auc_ph6", f"{ag}_auc_ph7", f"{ag}_ph_ratio",
                f"{ag}_od_max_ph7", f"{ag}_od_max_ph6",
                f"{ag}_binding", f"{ag}_ph_switch"
            ])
        writer.writerow(header)
        
        for com_num in range(1, 21):
            com_id = f"com{com_num}"
            mut_id = com_to_mut.get(com_id, "?")
            sc = scores.get(mut_id, {})
            
            row = [
                com_id, mut_id,
                sc.get("mutations", ""),
                sc.get("composite_score", ""),
                sc.get("dddG_elec", ""),
                sc.get("ddG_pH7.4", ""),
            ]
            
            for ag in ELISA_FILES:
                m = all_antigen_metrics[ag].get(com_id, {})
                row.extend([
                    m.get("auc_ph6", ""),
                    m.get("auc_ph7", ""),
                    m.get("ph_ratio", ""),
                    m.get("od_max_ph7", ""),
                    m.get("od_max_ph6", ""),
                    m.get("binding", m.get("binding_class", "")),
                    m.get("ph_switch", ""),
                ])
            
            writer.writerow(row)
    
    print(f"\n汇总已写入: {summary_path}")
    
    # 4. 相关性分析
    print("\n" + "=" * 60)
    print("相关性分析 (Spearman ρ)")
    print("=" * 60)
    
    report_lines = []
    report_lines.append("# 1E62_R2 实验 vs 计算预测 关联分析报告\n")
    report_lines.append(f"> 分析时间: 2026-03-16\n")
    report_lines.append(f"> 样本数: 20 个 AI 设计变体, 3 种 HBsAg 基因型 (Ae, B, D1)\n")
    report_lines.append(f"> ELISA 条件: pH 6.0 / pH 7.4, 6 浓度梯度, 三重复\n\n")
    
    report_lines.append("---\n\n")

    # ─── 4a. 各抗原的结合概览 ─────────────
    report_lines.append("## 1. 实验结合概览\n\n")
    
    for ag in ELISA_FILES:
        metrics = all_antigen_metrics[ag]
        report_lines.append(f"### HBsAg-{ag}\n\n")
        report_lines.append("| com | mut_id | OD_max(pH7.4) | OD_max(pH6.0) | AUC(pH7.4) | AUC(pH6.0) | pH_ratio | 结合 | pH开关 |\n")
        report_lines.append("|-----|--------|--------------|--------------|-----------|-----------|---------|------|-----|\n")
        
        for com_num in range(1, 21):
            com_id = f"com{com_num}"
            mut_id = com_to_mut.get(com_id, "?")
            m = metrics.get(com_id, {})
            if not m:
                continue
            
            report_lines.append(
                f"| {com_id} | {mut_id} | {m.get('od_max_ph7', 'N/A')} | "
                f"{m.get('od_max_ph6', 'N/A')} | {m.get('auc_ph7', 'N/A')} | "
                f"{m.get('auc_ph6', 'N/A')} | {m.get('ph_ratio', 'N/A')} | "
                f"{m.get('binding_class', '')} | {m.get('ph_switch', '')} |\n"
            )
        report_lines.append("\n")
    
    # ─── 4b. 相关性矩阵 ─────────────────
    report_lines.append("---\n\n")
    report_lines.append("## 2. 计算预测 vs 实验的 Spearman 相关性\n\n")
    
    comp_metrics = ["composite_score", "dddG_elec", "ddG_pH7.4"]
    exp_metrics_labels = ["auc_ph7", "auc_ph6", "ph_ratio", "od_max_ph7"]
    
    for ag in ELISA_FILES:
        metrics = all_antigen_metrics[ag]
        report_lines.append(f"### HBsAg-{ag}\n\n")
        
        # Header
        report_lines.append("| 计算指标 | vs AUC(pH7.4) | vs AUC(pH6.0) | vs pH_ratio | vs OD_max(pH7.4) |\n")
        report_lines.append("|---------|--------------|--------------|------------|------------------|\n")
        
        for comp_m in comp_metrics:
            row_vals = []
            for exp_m in exp_metrics_labels:
                # Collect paired data
                comp_vals = []
                exp_vals = []
                for com_num in range(1, 21):
                    com_id = f"com{com_num}"
                    mut_id = com_to_mut.get(com_id)
                    if not mut_id or mut_id not in scores:
                        continue
                    m = metrics.get(com_id, {})
                    if not m:
                        continue
                    
                    comp_v = scores[mut_id].get(comp_m)
                    exp_v = m.get(exp_m)
                    
                    if comp_v is not None and exp_v is not None and exp_v != "N/A":
                        comp_vals.append(float(comp_v))
                        exp_vals.append(float(exp_v))
                
                if len(comp_vals) >= 5:
                    rho, pval = spearman_rank_correlation(comp_vals, exp_vals)
                    sig = "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                    row_vals.append(f"ρ={rho:.3f} (p={pval:.3f}){sig}")
                else:
                    row_vals.append("N/A")
            
            report_lines.append(f"| {comp_m} | {' | '.join(row_vals)} |\n")
        
        report_lines.append("\n")
        print(f"\n--- HBsAg-{ag} 相关性 ---")
        # Also print to console
        for comp_m in comp_metrics:
            for exp_m in exp_metrics_labels:
                comp_vals = []
                exp_vals = []
                for com_num in range(1, 21):
                    com_id = f"com{com_num}"
                    mut_id = com_to_mut.get(com_id)
                    if not mut_id or mut_id not in scores:
                        continue
                    m = metrics.get(com_id, {})
                    if not m:
                        continue
                    comp_v = scores[mut_id].get(comp_m)
                    exp_v = m.get(exp_m)
                    if comp_v is not None and exp_v is not None and exp_v != "N/A":
                        comp_vals.append(float(comp_v))
                        exp_vals.append(float(exp_v))
                
                if len(comp_vals) >= 5:
                    rho, pval = spearman_rank_correlation(comp_vals, exp_vals)
                    print(f"  {comp_m} vs {exp_m}: ρ={rho:.3f}, p={pval:.3f} (n={len(comp_vals)})")
    
    # ─── 4c. 关键发现 ────────────────────
    report_lines.append("---\n\n")
    report_lines.append("## 3. 跨抗原一致性分析\n\n")
    
    # 看哪些 variant 在所有 3 种抗原上都有 strong binding
    report_lines.append("### 跨抗原结合一致性\n\n")
    report_lines.append("| com | mut_id | composite_score | Ae结合 | B结合 | D1结合 | 全阳性? |\n")
    report_lines.append("|-----|--------|----------------|-------|------|-------|--------|\n")
    
    for com_num in range(1, 21):
        com_id = f"com{com_num}"
        mut_id = com_to_mut.get(com_id, "?")
        sc = scores.get(mut_id, {})
        
        bindings = {}
        for ag in ["Ae", "B", "D1"]:
            m = all_antigen_metrics[ag].get(com_id, {})
            bindings[ag] = m.get("binding_class", "?")
        
        all_strong = all(b == "strong" for b in bindings.values())
        
        report_lines.append(
            f"| {com_id} | {mut_id} | {sc.get('composite_score', 'N/A')} | "
            f"{bindings['Ae']} | {bindings['B']} | {bindings['D1']} | "
            f"{'✅' if all_strong else '❌'} |\n"
        )
    
    report_lines.append("\n")
    
    # ─── 4d. pH-switch 质量分析 ──────────
    report_lines.append("### pH-switch 质量 (仅 strong binder)\n\n")
    report_lines.append("| com | mut_id | Ae_ratio | B_ratio | D1_ratio | 综合pH开关 |\n")
    report_lines.append("|-----|--------|---------|--------|---------|----------|\n")
    
    for com_num in range(1, 21):
        com_id = f"com{com_num}"
        mut_id = com_to_mut.get(com_id, "?")
        
        ratios = {}
        for ag in ["Ae", "B", "D1"]:
            m = all_antigen_metrics[ag].get(com_id, {})
            if m.get("binding_class") == "strong":
                ratios[ag] = m.get("ph_ratio", "N/A")
            else:
                ratios[ag] = "-"
        
        # 判断综合 pH 开关
        valid_ratios = [r for r in ratios.values() if isinstance(r, (int, float))]
        if valid_ratios:
            avg_ratio = sum(valid_ratios) / len(valid_ratios)
            if avg_ratio < 0.7:
                switch_qual = "🟢 强"
            elif avg_ratio < 0.9:
                switch_qual = "🟡 中"
            else:
                switch_qual = "🔴 无"
        else:
            switch_qual = "-"
        
        report_lines.append(
            f"| {com_id} | {mut_id} | {ratios['Ae']} | {ratios['B']} | "
            f"{ratios['D1']} | {switch_qual} |\n"
        )
    
    report_lines.append("\n")
    
    # ─── 4e. Top 发现 & 建议 ─────────────
    report_lines.append("---\n\n")
    report_lines.append("## 4. 关键发现\n\n")
    
    # 找出最好的 pH-switcher
    best_switchers = []
    for com_num in range(1, 21):
        com_id = f"com{com_num}"
        mut_id = com_to_mut.get(com_id, "?")
        
        # 检查至少在一种抗原上是 strong binder
        has_strong = False
        ratios = []
        for ag in ["Ae", "B", "D1"]:
            m = all_antigen_metrics[ag].get(com_id, {})
            if m.get("binding_class") == "strong":
                has_strong = True
                r = m.get("ph_ratio")
                if isinstance(r, (int, float)):
                    ratios.append(r)
        
        if has_strong and ratios:
            avg_r = sum(ratios) / len(ratios)
            best_switchers.append((com_id, mut_id, avg_r, len(ratios)))
    
    best_switchers.sort(key=lambda x: x[2])
    
    if best_switchers:
        report_lines.append("### 最佳 pH-switch 候选 (按 pH_ratio 排序，越低越好)\n\n")
        for rank, (com_id, mut_id, avg_r, n_ag) in enumerate(best_switchers[:5], 1):
            sc = scores.get(mut_id, {})
            report_lines.append(
                f"{rank}. **{mut_id}** ({com_id}): pH_ratio={avg_r:.3f} "
                f"(在 {n_ag} 种抗原上), "
                f"composite={sc.get('composite_score', 'N/A')}\n"
            )
        report_lines.append("\n")
    
    # Non-binders analysis
    non_binders = []
    for com_num in range(1, 21):
        com_id = f"com{com_num}"
        mut_id = com_to_mut.get(com_id, "?")
        
        all_none = True
        for ag in ["Ae", "B", "D1"]:
            m = all_antigen_metrics[ag].get(com_id, {})
            if m.get("binding_class") != "none":
                all_none = False
                break
        
        if all_none:
            sc = scores.get(mut_id, {})
            non_binders.append((com_id, mut_id, sc.get("composite_score", 0)))
    
    if non_binders:
        report_lines.append("### 完全不结合的变体 (需关注 — 计算 false positive)\n\n")
        for com_id, mut_id, cscore in non_binders:
            sc = scores.get(mut_id, {})
            report_lines.append(
                f"- **{mut_id}** ({com_id}): composite={cscore}, "
                f"dddG_elec={sc.get('dddG_elec', 'N/A')}, "
                f"ddG_pH7.4={sc.get('ddG_pH7.4', 'N/A')}\n"
            )
        report_lines.append("\n")
        report_lines.append(
            f"> ⚠️ **{len(non_binders)}/{len(com_to_mut)}** 个变体完全不结合。"
            f"这些是计算 pipeline 的 false positives，需深入分析原因。\n\n"
        )
    
    # 5. 写报告
    report_path = OUT_DIR / "correlation_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.writelines(report_lines)
    
    print(f"\n报告已写入: {report_path}")
    print(f"汇总 CSV: {summary_path}")
    print("\n" + "=" * 60)
    print("分析完成！")


if __name__ == "__main__":
    main()
