#!/usr/bin/env python3
"""
Step 13 (Tier 3): 精排 + 最终输出

从 tier2_candidates.csv 读取候选，按 dddG_elec 排序，添加软标记。
输出 final_candidates.csv。
"""
import os, sys, json, yaml
import pandas as pd
import numpy as np


def _add_soft_flags(df, t3_cfg):
    """添加软标记列。"""
    flags = t3_cfg.get("soft_flags", {})

    # ESM 异常标记（从 Tier 1 继承）
    if flags.get("esm_flag") and "esm_flag" not in df.columns:
        df["esm_flag"] = False

    # pH-score 异常标记（z-score > 2）
    if flags.get("phscore_flag") and "ph_score" in df.columns:
        ph = pd.to_numeric(df["ph_score"], errors="coerce")
        mean, std = ph.mean(), ph.std()
        if std > 0:
            df["phscore_flag"] = (ph - mean).abs() > 2 * std
        else:
            df["phscore_flag"] = False
    else:
        df["phscore_flag"] = False

    # pKa consensus 标记
    if flags.get("consensus_flag") and "overall_consensus" in df.columns:
        df["consensus_flag"] = df["overall_consensus"].apply(
            lambda x: x not in ("agree", "neutral") if pd.notna(x) else True)
    else:
        df["consensus_flag"] = False

    return df


def main(cfg_path):
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    t2_filter = cfg["tier2"]["filter"]
    t3_cfg = cfg["tier3"]

    input_csv = t2_filter.get("output", "results/tier2_candidates.csv")
    df = pd.read_csv(input_csv)
    print(f"[Tier 3] 输入: {input_csv} ({len(df)} 候选)")

    # 排序
    rank_col = t3_cfg.get("rank_by", "dddG_elec")
    rank_asc = t3_cfg.get("rank_ascending", False)

    if rank_col in df.columns and df[rank_col].notna().any():
        df = df.sort_values(rank_col, ascending=rank_asc).reset_index(drop=True)
        print(f"[Tier 3] 按 {rank_col} {'升序' if rank_asc else '降序'}排序")
    else:
        print(f"[Tier 3] 排序列 {rank_col} 不存在或全为空，保持原顺序")

    # 软标记
    df = _add_soft_flags(df, t3_cfg)

    # 输出列
    out_cols = [c for c in [
        "variant_id", "sequence", "mutations", "pdb_id", "mpdb",
        # Tier 1 指标
        "dG_pH7_4", "ddG_pH7_4", "delta", "esm_avg_logprob",
        # Tier 2 指标
        "pKa_propka", "pKa_pkai", "avg_shift_propka", "avg_shift_pkai",
        "overall_consensus", "pka_uplift",
        "dddG_elec", "ddG_elec_pH7", "ddG_elec_pH5", "ph_score",
        "global_rmsd", "h1_rmsd", "h2_rmsd", "h3_rmsd",
        "l1_rmsd", "l2_rmsd", "l3_rmsd",
        # 标记
        "esm_flag", "phscore_flag", "consensus_flag",
    ] if c in df.columns]

    out_path = t3_cfg.get("output", "results/final_candidates.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df[out_cols].to_csv(out_path, index=False)

    # 审计日志
    audit = {
        "tier3_input": len(pd.read_csv(input_csv)),
        "tier3_output": len(df),
        "rank_by": rank_col,
        "flags": {
            "esm_flagged": int(df["esm_flag"].sum()) if "esm_flag" in df.columns else 0,
            "phscore_flagged": int(df["phscore_flag"].sum()) if "phscore_flag" in df.columns else 0,
            "consensus_flagged": int(df["consensus_flag"].sum()) if "consensus_flag" in df.columns else 0,
        }
    }
    audit_dir = os.path.join(os.path.dirname(out_path), "audit")
    os.makedirs(audit_dir, exist_ok=True)
    with open(os.path.join(audit_dir, "tier3_manifest.json"), "w") as f:
        json.dump(audit, f, indent=2)

    print(f"[Tier 3] 输出: {out_path} (N={len(df)})")
    n_flagged = sum(1 for _, r in df.iterrows()
                    if r.get("esm_flag") or r.get("phscore_flag") or r.get("consensus_flag"))
    print(f"[Tier 3] 软标记: {n_flagged} 候选至少有 1 个 flag")


if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)
