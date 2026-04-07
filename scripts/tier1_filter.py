#!/usr/bin/env python3
"""
Step 8 (Tier 1): 自适应门槛筛选

从 FoldX + ESM 结果中按 dG_pH7.4 和 delta 做自适应门槛筛选。
输出 tier1_candidates.csv 供 Tier 2 使用。
"""
import os, sys, glob, yaml
import pandas as pd
import numpy as np


def _load_foldx_all(foldx_dir):
    """加载所有 FoldX 批次的 summary CSV。"""
    rows = []
    for csvp in glob.glob(os.path.join(foldx_dir, "batches", "*", "batch_*", "foldx_summary.csv")):
        pid = csvp.split(os.sep)[-3]
        bid = csvp.split(os.sep)[-2]
        df = pd.read_csv(csvp)
        if "mpdb" not in df.columns:
            print(f"[WARN] {csvp} 缺少 mpdb 列，跳过")
            continue
        num_cols = ["dG_pH7_4", "dG_pH6_0", "WT_dG_pH7_4", "WT_dG_pH6_0",
                    "ddG_pH7_4", "ddG_pH6_0", "delta", "delta_wt", "delta_delta"]
        for c in num_cols:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        df["pdb_id"] = pid
        df["batch"] = bid
        rows.append(df)
    if not rows:
        raise RuntimeError("未找到 foldx_summary.csv")
    fx = pd.concat(rows, ignore_index=True)
    # 补算缺失的 delta
    if "delta" not in fx or fx["delta"].isna().all():
        if {"dG_pH6_0", "dG_pH7_4"}.issubset(fx.columns):
            fx["delta"] = fx["dG_pH6_0"] - fx["dG_pH7_4"]
    return fx


def _load_batch_meta(foldx_dir):
    """加载批次元数据（序列、突变、ESM 分数等）。"""
    rows = []
    for bcsv in glob.glob(os.path.join(foldx_dir, "batches", "*", "batch_*", "batch_seqs.csv")):
        pid = bcsv.split(os.sep)[-3]
        bid = bcsv.split(os.sep)[-2]
        df = pd.read_csv(bcsv)
        if "mpdb" not in df.columns:
            if "mpdb_base" in df.columns:
                df["mpdb"] = df["mpdb_base"].astype(str) + ".pdb"
            else:
                continue
        df["pdb_id"] = pid
        df["batch"] = bid
        rows.append(df)
    if not rows:
        return pd.DataFrame(columns=["pdb_id", "batch", "mpdb"])
    return pd.concat(rows, ignore_index=True)


def _merge_foldx_meta(fx, meta):
    """合并 FoldX 结果和元数据。"""
    keep_cols = [c for c in ["sequence", "mutations", "source", "esm_avg_logprob", "hits_hotspots"]
                 if c in meta.columns]
    merged = fx.merge(meta[["pdb_id", "batch", "mpdb"] + keep_cols],
                      on=["pdb_id", "batch", "mpdb"], how="left")
    return merged


def _flag_esm(df, percentile):
    """标记 ESM 异常值（低于指定百分位）。"""
    if "esm_avg_logprob" not in df.columns or df["esm_avg_logprob"].isna().all():
        df["esm_flag"] = False
        return df
    threshold = np.nanpercentile(df["esm_avg_logprob"].astype(float), percentile)
    df["esm_flag"] = df["esm_avg_logprob"].astype(float) < threshold
    return df


def _adaptive_filter(df, t1_cfg):
    """自适应门槛筛选：调整 dG_pH7.4 和 delta 阈值直到候选量落入目标范围。"""
    target_lo, target_hi = t1_cfg["adaptive"]["target_range"]
    start = t1_cfg["adaptive"]["start"]
    step = t1_cfg["adaptive"]["relax_step"]
    floor = t1_cfg["adaptive"]["floors"]

    dG = float(start["dG74_max"])
    dm = float(start["delta_min"])
    dG_step = float(step["dG74_max"])
    dm_step = float(step["delta_min"])
    dG_floor = float(floor["dG74_max"])
    dm_floor = float(floor["delta_min"])

    best = pd.DataFrame()
    while True:
        q = (df["dG_pH7_4"] < dG) & (df["delta"] >= dm)
        pool = df[q]

        # 去重
        dedup_key = "mutations" if "mutations" in pool.columns else "sequence" if "sequence" in pool.columns else None
        if dedup_key:
            pool = pool.drop_duplicates(subset=[dedup_key], keep="first")

        n = len(pool)
        print(f"  dG74 < {dG:.1f}, delta >= {dm:.2f} → {n} 候选")

        if target_lo <= n <= target_hi:
            best = pool
            break
        if n > target_hi:
            # 通过数太多，收紧
            if len(pool) > len(best):
                best = pool.copy()
            dG -= abs(dG_step)
            dm += abs(dm_step)
            # 防止无限收紧
            if dG < float(start["dG74_max"]) - 5 * abs(dG_step):
                best = pool.head(target_hi)
                break
        else:
            # 通过数太少，放宽
            if len(pool) > len(best):
                best = pool.copy()
            stop_dG = dG >= dG_floor - 1e-9
            stop_dm = dm <= dm_floor + 1e-9
            if stop_dG and stop_dm:
                break
            dG = min(dG_floor, dG + abs(dG_step))
            dm = max(dm_floor, dm + dm_step)

    return best


def main(cfg_path):
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t1_cfg = cfg["tier1"]
    paths = cfg.get("paths", {})
    foldx_dir = paths.get("foldx_dir", cfg.get("foldx_dir", "foldx"))

    print("[Tier 1] 加载 FoldX 结果...")
    fx = _load_foldx_all(foldx_dir)
    print(f"  FoldX 总行数: {len(fx)}")

    meta = _load_batch_meta(foldx_dir)
    merged = _merge_foldx_meta(fx, meta)

    # ESM 异常标记
    esm_pct = t1_cfg.get("esm", {}).get("flag_below_percentile", 5)
    merged = _flag_esm(merged, esm_pct)

    print("[Tier 1] 自适应筛选...")
    filtered = _adaptive_filter(merged, t1_cfg)

    # 按 delta 降序排序
    filtered = filtered.sort_values("delta", ascending=False).reset_index(drop=True)

    out_path = t1_cfg.get("output", "results/tier1_candidates.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    out_cols = [c for c in [
        "sequence", "mutations", "pdb_id", "batch", "mpdb", "source",
        "esm_avg_logprob", "esm_flag", "hits_hotspots",
        "dG_pH7_4", "dG_pH6_0", "delta",
        "WT_dG_pH7_4", "WT_dG_pH6_0",
        "ddG_pH7_4", "ddG_pH6_0",
    ] if c in filtered.columns]
    filtered[out_cols].to_csv(out_path, index=False)
    print(f"[Tier 1] 输出: {out_path} (N={len(filtered)})")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)
