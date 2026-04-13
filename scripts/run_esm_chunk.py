#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, csv, yaml, json, glob
import pandas as pd
import torch, esm

def load_model(device):
    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    model.eval().to(device)
    bc = alphabet.get_batch_converter()
    return model, bc

def score_block(model, bc, device, seqs):
    data = [(f"id{i}", s) for i, s in enumerate(seqs)]
    _, _, toks = bc(data)
    toks = toks.to(device)
    with torch.no_grad():
        out = model(toks, repr_layers=[], return_contacts=False)
        logp = out["logits"].log_softmax(-1).cpu()
    scores = []
    for i, s in enumerate(seqs):
        L = len(s)
        idx = torch.arange(1, L + 1)
        aa = toks[i, 1:L + 1].cpu()
        scores.append(float(logp[i, idx, aa].mean()))
    return scores

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)
    return p

def gather_mpnn_all(cfg, results_screen_dir):
    """优先使用 results/screening/mpnn_all.csv；若不存在则从 mpnn_outputs/**/seqs*.csv 汇总生成。"""
    mpnn_all = os.path.join(results_screen_dir, "mpnn_all.csv")
    if os.path.exists(mpnn_all):
        return mpnn_all

    roots = [
        os.path.join(cfg["paths"]["mpnn_out_dir"], "**", "seqs*.csv"),
        os.path.join(cfg["paths"]["mpnn_out_dir"], "**", "seqs.csv"),
    ]
    hits = []
    for pat in roots:
        hits += glob.glob(pat, recursive=True)
    hits = sorted(set(hits))
    if not hits:
        return None

    rows = []
    for p in hits:
        try:
            df = pd.read_csv(p)
            if "sequence" not in df.columns:
                continue
            # 尽量保留来源信息（若已有 pdb_id/source 等）
            keep_cols = [c for c in df.columns if c != "sequence"]
            df2 = df[["sequence"] + keep_cols].copy()
            # 标注来源为 mpnn（不覆盖已有 source）
            if "source" not in df2.columns:
                df2["source"] = "mpnn"
            rows.append(df2)
        except Exception:
            continue
    if not rows:
        return None

    mpnn_df = pd.concat(rows, ignore_index=True)
    # 轻度去重：以 sequence+pdb_id（若有）为键
    keys = ["sequence"] + ([ "pdb_id"] if "pdb_id" in mpnn_df.columns else [])
    mpnn_df = mpnn_df.drop_duplicates(subset=keys, keep="first")
    ensure_dir(results_screen_dir)
    mpnn_df.to_csv(mpnn_all, index=False)
    print(f"[OK] 构建 mpnn_all.csv -> {mpnn_all}  n={len(mpnn_df)}")
    return mpnn_all

def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    device = cfg["esm"]["device"]
    bs     = int(cfg["esm"]["batch_size"])
    chunk  = int(cfg["resources"]["esm_chunk_size"])
    esm_outdir = ensure_dir(cfg["paths"]["esm_out_dir"])

    # 统一结果目录（向后兼容）
    results_dir = cfg["paths"].get("results_dir", cfg.get("results_dir", "results"))
    results_screen_dir = ensure_dir(os.path.join(results_dir, "screening"))

    # 收集输入（MPNN 汇总 + His 种子）
    inputs = []

    mpnn_all = gather_mpnn_all(cfg, results_screen_dir)
    if mpnn_all:
        inputs.append(mpnn_all)

    his_files = sorted(glob.glob(os.path.join(cfg["paths"]["his_seed_dir"], "his_seed_*.csv")))
    inputs += his_files

    if not inputs:
        raise RuntimeError("未找到可评分的输入CSV（既无 mpnn_all 也无 his_seed_*）")

    # 评分去重控制（防止重复跑）
    done_path = os.path.join(esm_outdir, "esm_done.json")
    done = set(json.load(open(done_path))) if os.path.exists(done_path) else set()

    # 载入模型
    model, bc = load_model(device)

    # 逐文件评分，并在每个 _esm.csv 中保留原始列（除了重复的 sequence）
    out_paths = []
    for path in inputs:
        base = os.path.basename(path)
        tag = os.path.splitext(base)[0]  # mpnn_all / his_seed_S1_single / his_seed_S2_pairs
        if base in done:
            print(f"[skip] {base}")
            out_paths.append(os.path.join(esm_outdir, f"{tag}_esm.csv"))
            continue

        df_in = pd.read_csv(path)
        if "sequence" not in df_in.columns:
            print(f"[WARN] {path} 缺少 sequence 列，已跳过")
            continue

        # -------- 多链格式处理：提取抗体链用于 ESM 评分，原始序列保留 --------
        has_multi = df_in["sequence"].astype(str).str.contains("/").any()
        if has_multi:
            # 保留原始多链序列
            df_in["_original_sequence"] = df_in["sequence"].copy()
            # 提取抗体链（去掉抗原）用于 ESM 评分
            design_chain = cfg["design"]["chain"]
            if isinstance(design_chain, list):
                # 多链设计：拼接所有设计链
                chain_indices = list(range(len(design_chain)))
            else:
                idx_map = {"A": 0, "B": 1, "C": 2}
                chain_indices = [idx_map.get(str(design_chain).strip().upper(), 0)]

            def extract_ab_chains(seq: str) -> str:
                parts = [p for p in str(seq).split("/") if p != ""]
                return "".join(parts[i] for i in chain_indices if i < len(parts))

            df_in["sequence"] = df_in["sequence"].astype(str).map(extract_ab_chains)
            print(f"[INFO] {os.path.basename(path)} 多链格式 → 提取链 {design_chain} 用于 ESM 评分")
        # -----------------------------------------------------------------------

        seqs = df_in["sequence"].astype(str).tolist()

        # 小块送 GPU（每块 chunk，再切 batch）
        scores = []
        for i in range(0, len(seqs), chunk):
            block = seqs[i:i + chunk]
            # 再分成批
            blk_scores = []
            for j in range(0, len(block), bs):
                part = block[j:j + bs]
                blk_scores += score_block(model, bc, device, part)
            scores += blk_scores

        # 输出 _esm.csv：保留原始列 + esm_avg_logprob + source（若已存在 source 列则不覆盖）
        out_csv = os.path.join(esm_outdir, f"{tag}_esm.csv")
        df_out = df_in.copy()
        df_out["esm_avg_logprob"] = scores
        # 恢复原始多链序列
        if "_original_sequence" in df_out.columns:
            df_out["sequence"] = df_out["_original_sequence"]
            df_out = df_out.drop(columns=["_original_sequence"])
        if "source" not in df_out.columns:
            # tag 以 his_seed* 开头的标为 his_seed，其他标为 mpnn
            src = "his_seed" if os.path.basename(path).startswith("his_seed") else "mpnn"
            df_out["source"] = src
        df_out.to_csv(out_csv, index=False)
        print(f"[OK] {base} -> {out_csv}  n={len(df_out)}")

        done.add(base)
        json.dump(sorted(done), open(done_path, "w"))
        out_paths.append(out_csv)

    # 合并为一个总表：results/screening/esm_all.csv
    rows = []
    for p in out_paths:
        if not os.path.exists(p): 
            continue
        try:
            rows.append(pd.read_csv(p))
        except Exception:
            pass
    if not rows:
        raise RuntimeError("未能生成任何 *_esm.csv，无法合并 esm_all.csv")

    esm_all = pd.concat(rows, ignore_index=True)

    # 轻度去重：尽量保留更多上下文列（sequence + pdb_id 若存在）
    keys = ["sequence"] + ([ "pdb_id"] if "pdb_id" in esm_all.columns else [])
    esm_all = esm_all.drop_duplicates(subset=keys + ["source"], keep="first")

    out_esm_all = os.path.join(results_screen_dir, "esm_all.csv")
    esm_all.to_csv(out_esm_all, index=False)
    print(f"[OK] 合并 -> {out_esm_all}  n={len(esm_all)}")

if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv)>1 else "configs/config.yaml"
    main(cfg)
