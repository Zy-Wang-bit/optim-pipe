#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, yaml
import pandas as pd
import re

AA20 = set(list("ACDEFGHIKLMNPQRSTVWY"))

def clean_seq(s: str) -> str:
    """去空白、转大写、仅保留20种标准AA。"""
    if not isinstance(s, str):
        return ""
    s = "".join([c for c in s.strip().upper() if c.isalpha()])
    s = "".join([c for c in s if c in AA20])
    return s

def resolve_wt_fasta(cfg, chain_id: str) -> str:
    """
    仅按 PDB 链 ID（A/B/C/…）解析 WT FASTA。
    - 若 paths.wt_files.<CHAIN> 是绝对路径：直接使用；
    - 若是相对路径：拼接到 paths.wt_dir 下（可含多级子目录，如 "IGH/heavy.fasta"）。
    未配置则报错，避免误判。
    """
    chain_id = str(chain_id).strip().upper()
    wt_dir = cfg["paths"]["wt_dir"]
    mapping = (cfg["paths"].get("wt_files", {}) or {})
    if chain_id in mapping:
        p = mapping[chain_id]
        # 绝对路径直接用；相对路径拼到 wt_dir
        return p if os.path.isabs(p) else os.path.join(wt_dir, p)
    raise FileNotFoundError(
        f"缺少 WT FASTA 映射：paths.wt_files.{chain_id}。"
        f"请在 configs/config.yaml 为链 {chain_id} 显式指定文件（允许相对 paths.wt_dir 的多级子目录）。"
    )

def load_wt(cfg, chain_id: str) -> str:
    """读取并校验 WT 序列（按链 ID）。"""
    path = resolve_wt_fasta(cfg, chain_id)
    if not os.path.exists(path):
        raise FileNotFoundError(f"未找到 {chain_id} 链 WT 序列文件：{path}")
    seq = ""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq += "".join([c for c in line.upper() if c.isalpha()])
    if not seq:
        raise ValueError(f"{chain_id} 链 WT 序列为空：{path}")
    bad = set(seq) - AA20
    if bad:
        raise ValueError(f"{chain_id} 链 WT 序列含非标准AA {bad} 于 {path}")
    print(f"[WT] chain={chain_id} file={path} length={len(seq)}")
    return seq

def seq_to_mutcodes(wt: str, var: str, chain_id: str, region):
    """仅在 region 内比较差异，生成 FoldX 突变码列表（例：FA102H）。"""
    s, e = region
    muts = []
    for i in range(s, e + 1):
        if wt[i - 1] != var[i - 1]:
            muts.append(f"{wt[i - 1]}{chain_id}{i}{var[i - 1]}")
    return muts

def _parse_design_chains(cfg):
    """解析 design 配置，支持单链和多链。返回 {chain: (start, end)}"""
    chain = cfg["design"]["chain"]
    if isinstance(chain, str):
        region = cfg["design"]["region"]
        return {chain: (region[0], region[1])}
    else:
        regions = cfg["design"]["regions"]
        return {ch: (regions[ch][0], regions[ch][1]) for ch in chain}


def multi_seq_to_mutcodes(wt_seqs, var_seq_str, design_chains):
    """多链版本的突变码生成。
    wt_seqs: {chain: wt_seq}
    var_seq_str: "H_SEQ/L_SEQ/..." 或单链序列
    design_chains: {chain: (start, end)}
    返回: 突变码列表 ['FA102H', 'SB42H', ...]
    """
    chain_order = sorted(wt_seqs.keys())
    parts = var_seq_str.split("/")

    muts = []
    for idx, ch in enumerate(chain_order):
        if ch not in design_chains:
            continue
        if idx >= len(parts):
            continue
        var = parts[idx]
        wt = wt_seqs[ch]
        s, e = design_chains[ch]
        if len(var) != len(wt):
            continue
        for i in range(s, e + 1):
            if i <= len(wt) and i <= len(var) and wt[i-1] != var[i-1]:
                muts.append(f"{wt[i-1]}{ch}{i}{var[i-1]}")
    return muts


def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))

    design_chains = _parse_design_chains(cfg)
    chain_ids = sorted(design_chains.keys())
    is_multi = len(chain_ids) > 1

    chunk   = int(cfg["resources"]["foldx_chunk_size"])
    root    = os.path.join(cfg["paths"]["foldx_dir"], "batches")
    os.makedirs(root, exist_ok=True)

    # 输入候选（ESM/His seeds 合并后的 for_foldx.csv）
    for_foldx_csv = os.path.join(cfg["paths"]["results_dir"], "screening", "for_foldx.csv")
    if not os.path.exists(for_foldx_csv):
        alt = os.path.join(cfg.get("results_dir", "results"), "screening", "for_foldx.csv")
        for_foldx_csv = alt if os.path.exists(alt) else for_foldx_csv
    if not os.path.exists(for_foldx_csv):
        raise FileNotFoundError(f"未找到候选列表：{for_foldx_csv}")

    # 加载 WT 序列（所有设计链）
    wt_seqs = {}
    for ch in chain_ids:
        wt_seqs[ch] = load_wt(cfg, ch)

    # 兼容单链：用第一条链的长度和 region 做基本校验
    if not is_multi:
        chain_id = chain_ids[0]
        wt = wt_seqs[chain_id]
        Lw = len(wt)
        s, e = design_chains[chain_id]
        if not (1 <= s <= e <= Lw):
            raise ValueError(f"design.region={[s,e]} 超出 WT 长度（WT={Lw} aa）。")
    else:
        # 多链：总长度用于序列校验
        Lw = sum(len(wt_seqs[ch]) for ch in chain_ids)

    df = pd.read_csv(for_foldx_csv)
    if "pdb_id" not in df.columns or "sequence" not in df.columns:
        raise ValueError("for_foldx.csv 需至少包含列：pdb_id, sequence")

    # 基本清洗（多链序列含 /，不能 clean_seq 整体处理）
    if is_multi:
        def _clean_multi(s):
            if not isinstance(s, str): return ""
            return "/".join(clean_seq(p) for p in s.split("/"))
        df["sequence"] = df["sequence"].apply(_clean_multi)
    else:
        df["sequence"] = df["sequence"].apply(clean_seq)
    df = df[df["sequence"].str.len() > 0].copy()

    # 加载所有链的 WT（包括非设计链，用于序列长度校验）
    all_wt_files = cfg["paths"].get("wt_files", {})
    all_chain_ids = sorted(all_wt_files.keys())  # 所有链 ID（如 A, B, C）

    # 过滤并收集可用行
    bad_rows, keep_rows = [], []
    for _, row in df.iterrows():
        seq = row["sequence"]
        if is_multi or "/" in seq:
            parts = seq.split("/")
            # 多链序列可能包含抗原链（如 A/B/C），检查设计链的长度匹配
            ok = True
            for idx, ch in enumerate(chain_ids):
                # 找到该设计链在完整序列中的位置
                if ch in all_chain_ids:
                    ch_pos = all_chain_ids.index(ch)
                else:
                    ch_pos = idx
                if ch_pos >= len(parts):
                    bad_rows.append({"reason": f"missing_chain:{ch}", **row.to_dict()})
                    ok = False; break
                if len(parts[ch_pos]) != len(wt_seqs[ch]):
                    bad_rows.append({"reason": f"len_mismatch:{ch}:{len(parts[ch_pos])}!=WT{len(wt_seqs[ch])}", **row.to_dict()})
                    ok = False; break
            if not ok: continue
        else:
            # 单链序列
            bad = set(seq) - AA20
            if bad:
                bad_rows.append({"reason": f"illegal_aa:{''.join(sorted(bad))}", **row.to_dict()}); continue
            if len(seq) != Lw:
                bad_rows.append({"reason": f"len_mismatch:{len(seq)}!=WT{Lw}", **row.to_dict()}); continue
        keep_rows.append(row)

    if not keep_rows:
        skip_root = os.path.join(root, "_skipped"); os.makedirs(skip_root, exist_ok=True)
        pd.DataFrame(bad_rows).to_csv(os.path.join(skip_root, "skipped_all.csv"), index=False)
        raise RuntimeError(f"所有序列被过滤。请检查 for_foldx.csv、WT长度与 design.region。跳过明细见 {skip_root}/skipped_all.csv")

    clean_df = pd.DataFrame(keep_rows)

    total_batches = 0
    for pid, sub in clean_df.groupby("pdb_id"):
        out_p = os.path.join(root, pid); os.makedirs(out_p, exist_ok=True)

        parts = [sub.iloc[i:i + chunk].copy() for i in range(0, len(sub), chunk)]
        for bi, part in enumerate(parts):
            bdir = os.path.join(out_p, f"batch_{bi:05d}"); os.makedirs(bdir, exist_ok=True)

            kept_records, skipped = [], []
            for _, r in part.iterrows():
                seq = r["sequence"]
                try:
                    if is_multi:
                        muts = multi_seq_to_mutcodes(wt_seqs, seq, design_chains)
                    else:
                        muts = seq_to_mutcodes(wt_seqs[chain_ids[0]], seq, chain_ids[0], design_chains[chain_ids[0]])
                except Exception as ex:
                    skipped.append({"reason": f"exception:{ex}", **r.to_dict()}); continue
                if not muts:
                    skipped.append({"reason": "no_mut_in_region", **r.to_dict()}); continue
                rr = r.to_dict()
                rr["mutations"] = ",".join(muts)
                kept_records.append(rr)

            # individual_list.txt（仅 kept；行号即 *_Repair_<i>.pdb 的 i）
            ind_path = os.path.join(bdir, "individual_list.txt")
            with open(ind_path, "w") as g:
                for rr in kept_records:
                    line = rr["mutations"]
                    if not line.endswith(";"):
                        line += ";"
                    g.write(line + "\n")

            # batch_seqs.csv（只包含 kept；写出 mpdb 对应关系）
            if kept_records:
                df_keep = pd.DataFrame(kept_records).reset_index(drop=True)
                df_keep["row_index"]  = range(len(df_keep))  # 0..N-1
                df_keep["mpdb_base"]  = df_keep["row_index"].apply(lambda i: f"{pid}_Repair_{i+1}")
                df_keep["mpdb"]       = df_keep["mpdb_base"] + ".pdb"
                df_keep["pdb_id"]     = pid
                df_keep["batch"]      = f"batch_{bi:05d}"
                cols_first = ["pdb_id","batch","row_index","mpdb_base","mpdb","sequence","mutations"]
                cols_more  = [c for c in ["source","esm_avg_logprob","hits_hotspots"] if c in df_keep.columns]
                out_df = df_keep[cols_first + cols_more]
                out_df.to_csv(os.path.join(bdir, "batch_seqs.csv"), index=False)
                written = len(df_keep)
            else:
                pd.DataFrame(columns=["pdb_id","batch","row_index","mpdb_base","mpdb","sequence","mutations"]).to_csv(
                    os.path.join(bdir, "batch_seqs.csv"), index=False
                )
                written = 0

            if skipped:
                pd.DataFrame(skipped).to_csv(os.path.join(bdir, "skipped.csv"), index=False)

            print(f"[{pid}][batch {bi:05d}] 写入 {written} 条；跳过 {len(skipped)} 条 -> {bdir}")
            total_batches += 1

    print(f"[OK] make_mutlist_chunk -> {root}  批次数={total_batches}")

if __name__ == "__main__":
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    main(cfg)