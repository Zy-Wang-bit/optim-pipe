#!/usr/bin/env python3
"""
Step 12 (Tier 2): 门槛筛选

合并 PyRosetta 线 (pKa, dddG_elec, pH-score) 和 SimpleFold 3x 线 (CDR RMSD) 的结果。
按 pKa 相对于 WT/基底抗体排序，取 top-N；按 CDR RMSD 硬门槛过滤。
输出 tier2_candidates.csv 供 Tier 3 使用。
"""
import os, sys, yaml
import pandas as pd
import numpy as np
from pathlib import Path


def _load_tier1(cfg):
    """加载 Tier 1 候选列表。"""
    t1_out = cfg["tier1"].get("output", "results/tier1_candidates.csv")
    df = pd.read_csv(t1_out)
    # 添加 variant_id
    if "mpdb" in df.columns:
        df["variant_id"] = df["mpdb"].apply(
            lambda x: Path(str(x)).stem if pd.notna(x) else None)
    return df


def _attach_pka(df, t2_dir):
    """左连接 pKa 结果。"""
    pka_path = os.path.join(t2_dir, "pka", "pka_summary.csv")
    if not os.path.exists(pka_path):
        print(f"[Tier 2] pKa 文件不存在: {pka_path}")
        return df
    pka = pd.read_csv(pka_path)
    if "variant_name" in pka.columns and "variant_id" not in pka.columns:
        pka = pka.rename(columns={"variant_name": "variant_id"})
    pka_cols = [c for c in ["variant_id", "avg_shift_propka", "avg_shift_pkai",
                            "overall_consensus", "pKa_propka", "pKa_pkai"] if c in pka.columns]
    if "variant_id" in pka_cols:
        df = df.merge(pka[pka_cols], on="variant_id", how="left")
        print(f"[Tier 2] pKa: {pka_path} ({len(pka)} 行)")
    return df


def _attach_rosetta(df, t2_dir):
    """左连接 Rosetta 结果 (dddG_elec + pH-score)。"""
    for filename, label in [("dddg_elec.csv", "dddG_elec"), ("ph_scores.csv", "pH-score")]:
        path = os.path.join(t2_dir, "rosetta", filename)
        if not os.path.exists(path):
            continue
        sub = pd.read_csv(path)
        if "variant_name" in sub.columns and "variant_id" not in sub.columns:
            sub = sub.rename(columns={"variant_name": "variant_id"})
        cols = [c for c in sub.columns if c == "variant_id" or c in
                ["dddG_elec", "ph_score", "ddG_elec_pH7", "ddG_elec_pH5"]]
        df = df.merge(sub[cols], on="variant_id", how="left")
        print(f"[Tier 2] {label}: {path} ({len(sub)} 行)")
    return df


def _attach_rmsd(df, t2_dir):
    """左连接 RMSD 结果（优先 SimpleFold 3x）。"""
    rmsd_path = os.path.join(t2_dir, "rmsd", "rmsd_summary.csv")
    if not os.path.exists(rmsd_path):
        print(f"[Tier 2] RMSD 文件不存在: {rmsd_path}")
        return df
    rmsd = pd.read_csv(rmsd_path)
    if "variant_name" in rmsd.columns and "variant_id" not in rmsd.columns:
        rmsd = rmsd.rename(columns={"variant_name": "variant_id"})
    # 优先 simplefold_3x 来源
    if "source" in rmsd.columns:
        sf3x = rmsd[rmsd["source"] == "simplefold_3x"]
        if len(sf3x) > 0:
            rmsd = sf3x
        else:
            rmsd = rmsd.drop_duplicates(subset=["variant_id"], keep="first")
    rmsd_cols = ["variant_id"] + [c for c in rmsd.columns if c.endswith("_rmsd")]
    df = df.merge(rmsd[rmsd_cols], on="variant_id", how="left")
    print(f"[Tier 2] RMSD: {rmsd_path} ({len(rmsd)} 行)")
    return df


def _filter_pka_relative(df, t2_dir, filter_cfg):
    """按 pKa 相对于 WT 排序，取 top-N。"""
    # 尝试读取 WT pKa 基线
    wt_pka_path = os.path.join(t2_dir, "pka", "wt_pka.csv")
    wt_pka_avg = None
    if os.path.exists(wt_pka_path):
        wt_pka = pd.read_csv(wt_pka_path)
        if "pKa_propka" in wt_pka.columns:
            wt_pka_avg = wt_pka["pKa_propka"].mean()

    # 排序指标优先级: pKa_propka 相对提升 > avg_shift_propka > 跳过
    sort_col = None
    if "pKa_propka" in df.columns and df["pKa_propka"].notna().any():
        if wt_pka_avg is not None:
            df["pka_uplift"] = df["pKa_propka"] - wt_pka_avg
            sort_col = "pka_uplift"
        else:
            sort_col = "pKa_propka"
    elif "avg_shift_propka" in df.columns and df["avg_shift_propka"].notna().any():
        sort_col = "avg_shift_propka"

    if sort_col is None:
        print("[Tier 2] 无可用 pKa 数据，跳过 pKa 排序")
        return df

    top_n = filter_cfg.get("pka_top_n", 200)

    # 排除比 WT 差的
    if "pka_uplift" in df.columns:
        before = len(df)
        df = df[df["pka_uplift"] > 0]
        print(f"[Tier 2] pKa: 排除 uplift <= 0 ({before} → {len(df)})")

    df = df.sort_values(sort_col, ascending=False).head(top_n)
    print(f"[Tier 2] pKa top-{top_n}: {len(df)} 候选 (排序列: {sort_col})")
    return df


def _filter_rmsd(df, filter_cfg):
    """CDR RMSD 硬门槛。"""
    metric = filter_cfg.get("cdr_rmsd_metric", "h1_rmsd")
    threshold = float(filter_cfg.get("cdr_rmsd_max", 0.5))

    if metric not in df.columns or df[metric].isna().all():
        print(f"[Tier 2] RMSD 列 {metric} 不存在或全为空，跳过 RMSD 筛选")
        return df

    before = len(df)
    df = df[df[metric] <= threshold]
    print(f"[Tier 2] RMSD ({metric} <= {threshold}Å): {before} → {len(df)}")
    return df


def main(cfg_path):
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t2 = cfg["tier2"]
    t2_dir = t2["paths"].get("tier2_dir", t2["paths"].get("phase_c_dir", "tier2"))
    filter_cfg = t2["filter"]

    print("[Tier 2] 加载 Tier 1 候选...")
    df = _load_tier1(cfg)
    print(f"  Tier 1 候选: {len(df)}")

    # 合并评估结果
    df = _attach_pka(df, t2_dir)
    df = _attach_rosetta(df, t2_dir)
    df = _attach_rmsd(df, t2_dir)

    # 筛选
    df = _filter_pka_relative(df, t2_dir, filter_cfg)
    df = _filter_rmsd(df, filter_cfg)

    # 输出
    out_path = filter_cfg.get("output", "results/tier2_candidates.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[Tier 2] 输出: {out_path} (N={len(df)})")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)
