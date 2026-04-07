#!/usr/bin/env python3
import os, sys, glob, json, yaml
import pandas as pd
import numpy as np
from math import ceil
from pathlib import Path

def _paths(cfg):
    p = cfg.get("paths", {})
    return {
        "foldx_dir": p.get("foldx_dir", cfg.get("foldx_dir", "foldx")),
        "results_dir": p.get("results_dir", cfg.get("results_dir", "results")),
    }

def _load_foldx_all(P):
    rows = []
    for csvp in glob.glob(os.path.join(P["foldx_dir"], "batches", "*", "batch_*", "foldx_summary.csv")):
        pid = csvp.split(os.sep)[-3]; bid = csvp.split(os.sep)[-2]
        df = pd.read_csv(csvp)
        if "mpdb" not in df.columns:
            print(f"[WARN] {csvp} 缺少 mpdb 列，跳过"); continue
        # 转数值
        num_cols = ["dG_pH7_4","dG_pH6_0","WT_dG_pH7_4","WT_dG_pH6_0",
                    "ddG_pH7_4","ddG_pH6_0","delta","delta_wt","delta_delta"]
        for c in num_cols:
            if c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
        df["pdb_id"]=pid; df["batch"]=bid
        rows.append(df)
    if not rows:
        raise RuntimeError("未找到 foldx_summary.csv")
    fx = pd.concat(rows, ignore_index=True)
    # 保险：若没有 delta/delta_wt（老文件），补算
    if "delta" not in fx or fx["delta"].isna().all():
        if {"dG_pH6_0","dG_pH7_4"}.issubset(fx.columns):
            fx["delta"] = fx["dG_pH6_0"] - fx["dG_pH7_4"]
    if "delta_wt" not in fx and {"WT_dG_pH6_0","WT_dG_pH7_4"}.issubset(fx.columns):
        fx["delta_wt"] = fx["WT_dG_pH6_0"] - fx["WT_dG_pH7_4"]
    if "delta_delta" not in fx and {"delta","delta_wt"}.issubset(fx.columns):
        fx["delta_delta"] = fx["delta"] - fx["delta_wt"]
    return fx

def _load_batch_meta(P):
    rows = []
    for bcsv in glob.glob(os.path.join(P["foldx_dir"], "batches", "*", "batch_*", "batch_seqs.csv")):
        pid = bcsv.split(os.sep)[-3]; bid = bcsv.split(os.sep)[-2]
        df = pd.read_csv(bcsv)
        if "mpdb" not in df.columns:
            if "mpdb_base" in df.columns:
                df["mpdb"] = df["mpdb_base"].astype(str) + ".pdb"
            else:
                print(f"[WARN] {bcsv} 缺少 mpdb/mpdb_base 列，跳过"); continue
        df["pdb_id"]=pid; df["batch"]=bid
        rows.append(df)
    if not rows:
        print("[WARN] 未找到 batch_seqs.csv；合并将仅包含 FoldX 指标")
        return pd.DataFrame(columns=["pdb_id","batch","mpdb"])
    return pd.concat(rows, ignore_index=True)

def _attach_meta(fx, meta, screen_dir):
    keep_cols = [c for c in ["sequence","mutations","source","esm_avg_logprob","hits_hotspots"] if c in meta.columns]
    merged = fx.merge(meta[["pdb_id","batch","mpdb"] + keep_cols], on=["pdb_id","batch","mpdb"], how="left")
    # ESM 归一化
    if "esm_avg_logprob" in merged.columns and merged["esm_avg_logprob"].notna().any():
        x = merged["esm_avg_logprob"].astype(float)
        x = x.fillna(x.mean())
        merged["esm_norm"] = (x - x.min()) / (x.max() - x.min() + 1e-8)
    else:
        merged["esm_norm"] = 0.0
    merged.to_csv(os.path.join(screen_dir,"foldx_merged.csv"), index=False)
    return merged

def _attach_phase_c(merged, cfg):
    """左连接 Phase C 结果 (pKa, Rosetta, RMSD)。未启用或无结果时返回原 df。"""
    pc = cfg.get("phase_c", {})
    if not pc.get("enabled", False):
        return merged

    pc_dir = pc.get("paths", {}).get("phase_c_dir", "phase_c")

    map_csv = os.path.join(pc_dir, "structures", "mutant_map.csv")
    if not os.path.exists(map_csv):
        print("[Phase C] mutant_map.csv 未找到，跳过 Phase C 合并")
        return merged

    # 为 merged 添加 variant_id 列 (从 mpdb 列取 stem)
    if "mpdb" in merged.columns:
        merged["variant_id"] = merged["mpdb"].apply(
            lambda x: Path(str(x)).stem if pd.notna(x) else None
        )
    else:
        print("[Phase C] merged 缺少 mpdb 列，无法建立映射")
        return merged

    def _rename_vid(df):
        if "variant_name" in df.columns and "variant_id" not in df.columns:
            df = df.rename(columns={"variant_name": "variant_id"})
        return df

    # 统一的 Phase C 数据源定义: (子目录, 文件名, 需要的列, 标签)
    joins = [
        ("pka",     "pka_summary.csv",   ["variant_id", "avg_shift_propka", "avg_shift_pkai", "overall_consensus"], "pKa"),
        ("rosetta", "ph_scores.csv",     ["variant_id", "ph_score"],        "pH-score"),
        ("rosetta", "dddg_elec.csv",     ["variant_id", "dddG_elec"],       "dddG_elec"),
    ]

    for subdir, filename, want_cols, label in joins:
        csv_path = os.path.join(pc_dir, subdir, filename)
        if not os.path.exists(csv_path):
            continue
        df = _rename_vid(pd.read_csv(csv_path))
        cols = [c for c in want_cols if c in df.columns]
        merged = merged.merge(df[cols], on="variant_id", how="left")
        print(f"[Phase C] {label}: {csv_path} ({len(df)} 行)")

    # RMSD 需要额外的 primary 过滤和动态列选择
    rmsd_path = os.path.join(pc_dir, "rmsd", "rmsd_summary.csv")
    if os.path.exists(rmsd_path):
        rmsd = _rename_vid(pd.read_csv(rmsd_path))
        primary = pc.get("structure_generation", {}).get("primary", "rosetta")
        rmsd = rmsd[rmsd["source"] == primary]
        rmsd_cols = ["variant_id", "global_rmsd"] + [
            c for c in rmsd.columns if c.endswith("_rmsd") and c != "global_rmsd"
        ]
        merged = merged.merge(rmsd[rmsd_cols], on="variant_id", how="left")
        print(f"[Phase C] RMSD: {rmsd_path} ({len(rmsd)} 行)")

    return merged

def _score(df, sel):
    # 评分：delta + α*(-dG7.4) + β*ESM + γ_his + ω*delta_delta(可选，默认0)
    alpha = float(sel.get("alpha", 0.3))
    beta  = float(sel.get("beta", 0.1))
    dd_w  = float(sel.get("delta_delta_weight", 0.0))  # 如需把 delta_delta 纳入评分，设 >0
    bonus = float(sel.get("his_bonus", 0.0))
    src = df.get("source")
    his_from_src = src.str.contains("his", case=False, na=False) if src is not None else False
    his_from_hot = df.get("hits_hotspots", 0).fillna(0).astype(bool) if "hits_hotspots" in df.columns else False
    his_bonus = np.where(his_from_src | his_from_hot, bonus, 0.0)
    dd_term = df["delta_delta"] if "delta_delta" in df.columns else 0.0
    sc = df["delta"] + alpha*(-df["dG_pH7_4"]) + beta*df["esm_norm"] + his_bonus + dd_w*dd_term
    out = df.copy()
    out["score"] = sc
    out["__is_his__"] = (his_from_src | his_from_hot)
    return out

def _fixed_filter_rank(df, sel):
    dG74_max = float(sel["dG74_max"])
    delta_min = float(sel["delta_min"])
    top_n = int(sel.get("top_n", 10000))
    # 若配置给了 delta_delta_min，则也做一下门槛（可选）
    ddd_min = sel.get("delta_delta_min", None)
    q = (df["dG_pH7_4"] < dG74_max) & (df["delta"] >= delta_min)
    if ddd_min is not None and "delta_delta" in df.columns:
        q = q & (df["delta_delta"] >= float(ddd_min))
    df = _score(df[q].copy(), sel).sort_values("score", ascending=False)
    dedup_key = "mutations" if "mutations" in df.columns else "sequence" if "sequence" in df.columns else None
    if dedup_key:
        df = df.drop_duplicates(subset=[dedup_key], keep="first")
    return df.head(top_n)

def _adaptive_pick(pool, sel):
    N = int(sel.get("target_n", sel.get("top_n", 10000)))
    per_pdb_frac = float(sel.get("per_pdb_max_frac", 1.0))
    his_min_frac = float(sel.get("his_min_frac", 0.0))
    cap = max(1, ceil(per_pdb_frac * N))
    need_his = int(np.floor(his_min_frac * N))

    pool = pool.sort_values("score", ascending=False)
    chosen = []
    per_pdb_cnt = {}
    def can_take(pid): return per_pdb_cnt.get(pid,0) < cap

    # 先保底 his
    for _, r in pool[pool["__is_his__"]].iterrows():
        if len(chosen) >= need_his: break
        pid = r["pdb_id"]
        if not can_take(pid): continue
        chosen.append(r); per_pdb_cnt[pid]=per_pdb_cnt.get(pid,0)+1

    # 再补满
    for _, r in pool.iterrows():
        if len(chosen) >= N: break
        pid = r["pdb_id"]
        if not can_take(pid): continue
        chosen.append(r); per_pdb_cnt[pid]=per_pdb_cnt.get(pid,0)+1

    out = pd.DataFrame(chosen)
    # 去重（如有序列/突变）
    dedup_key = "mutations" if "mutations" in out.columns else "sequence" if "sequence" in out.columns else None
    if dedup_key and len(out):
        out = out.drop_duplicates(subset=[dedup_key], keep="first")
    return out.head(N)

def _adaptive_filter_rank(df, sel):
    # 起点/步长/底线（与之前一致，若未配置采用默认）
    start = sel.get("start", {"dG74_max": -3.0, "delta_min": 1.0})
    step  = sel.get("relax_step", {"dG74_max": 0.5, "delta_min": -0.25})
    floor = sel.get("floors", {"dG74_max": 1.0, "delta_min": 0.0})
    # 可选：对 delta_delta 也做底线/步长
    ddd_min0 = sel.get("delta_delta_start", None)
    ddd_step = sel.get("delta_delta_step", None)
    ddd_floor= sel.get("delta_delta_floor", None)

    dG = float(start["dG74_max"]); dG_step = float(step["dG74_max"]); dG_floor = float(floor["dG74_max"])
    dm = float(start["delta_min"]); dm_step = float(step["delta_min"]); dm_floor = float(floor["delta_min"])
    ddd = float(ddd_min0) if ddd_min0 is not None else None
    ddd_s= float(ddd_step) if ddd_step is not None else None
    ddd_f= float(ddd_floor) if ddd_floor is not None else None

    # 预计算评分（不依赖过滤阈值），避免每次循环重复计算
    scored = _score(df.copy(), sel)

    best = pd.DataFrame()
    while True:
        q = (scored["dG_pH7_4"] < dG) & (scored["delta"] >= dm)
        if ddd is not None and "delta_delta" in scored.columns:
            q = q & (scored["delta_delta"] >= ddd)
        pool = scored[q]
        picked = _adaptive_pick(pool, sel)
        if len(picked) >= int(sel.get("target_n", sel.get("top_n", 10000))):
            best = picked; break
        if len(picked) > len(best): best = picked.copy()

        # 是否触底
        stop_dG = dG >= dG_floor - 1e-9
        stop_dm = dm <= dm_floor + 1e-9
        stop_dd = (ddd is None) or (ddd <= ddd_f + 1e-9)
        if stop_dG and stop_dm and stop_dd:
            break
        # 放宽
        dG  = min(dG_floor, dG + dG_step)
        dm  = max(dm_floor, dm + dm_step)
        if ddd is not None:
            ddd = max(ddd_f, ddd + (ddd_s or 0.0))

    return best

def main(cfg_path):
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    P = _paths(cfg)
    res_dir = P["results_dir"]; screen_dir = os.path.join(res_dir, "screening")
    os.makedirs(screen_dir, exist_ok=True)

    fx = _load_foldx_all(P)
    fx.to_csv(os.path.join(screen_dir,"foldx_all.csv"), index=False)
    print(f"[OK] foldx_all.csv n={len(fx)}")

    meta = _load_batch_meta(P)
    merged = _attach_meta(fx, meta, screen_dir)

    # Phase C: 左连接 pKa / Rosetta / RMSD 结果
    merged = _attach_phase_c(merged, cfg)
    if cfg.get("phase_c", {}).get("enabled"):
        merged.to_csv(os.path.join(res_dir, "merged_all.csv"), index=False)
        print(f"[OK] merged_all.csv n={len(merged)}")

    sel = cfg["selection"]; mode = str(sel.get("mode","fixed")).lower()
    if mode == "adaptive":
        final = _adaptive_filter_rank(merged, sel)
    else:
        final = _fixed_filter_rank(merged, sel)

    out_csv = os.path.join(res_dir, "final_top10k.csv")
    cols = [c for c in [
        "sequence","mutations","pdb_id","batch","mpdb","source",
        "esm_avg_logprob","esm_norm","hits_hotspots",
        "dG_pH7_4","dG_pH6_0","delta",
        "WT_dG_pH7_4","WT_dG_pH6_0","delta_wt",
        "ddG_pH7_4","ddG_pH6_0","delta_delta",
        # Phase C 列（存在时输出，缺失时忽略）
        "ph_score","dddG_elec",
        "avg_shift_propka","avg_shift_pkai","overall_consensus",
        "global_rmsd",
        "score"
    ] if c in final.columns]
    final[cols].to_csv(out_csv, index=False)
    print(f"[OK] final -> {out_csv}  (N={len(final)})")

    audit = {
        "mode": mode,
        "rows_after_foldx": int(len(merged)),
        "rows_final": int(len(final)),
        "selection": sel
    }
    os.makedirs(os.path.join(res_dir,"audit"), exist_ok=True)
    with open(os.path.join(res_dir,"audit","manifest.json"), "w") as f:
        json.dump(audit, f, indent=2)

if __name__=="__main__":
    cfg = sys.argv[1] if len(sys.argv)>1 else "configs/config.yaml"
    main(cfg)