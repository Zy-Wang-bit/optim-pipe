#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, re, yaml, csv
import itertools
import pandas as pd

AA20 = set("ACDEFGHIKLMNPQRSTVWY")

def resolve_wt_fasta(cfg, chain):
    wt_dir = cfg["paths"]["wt_dir"]
    mapping = cfg["paths"].get("wt_files", {})
    # 1) 显式配置优先
    if mapping and chain in mapping:
        p = mapping[chain]
        return p if os.path.isabs(p) else os.path.join(wt_dir, p)
    # 2) 关键词匹配兜底
    patts = {
        "H": [r"heavy", r"\bVH\b", r"[_\-]H(\.|$)", r"\bIGH\b"],
        "L": [r"light", r"\bVL\b", r"[_\-]L(\.|$)", r"\bIG[LK]\b"],
        "A": [r"antigen", r"\bAg\b", r"[_\-]A(\.|$)"],
    }
    cand=[]
    if os.path.isdir(wt_dir):
        for fn in os.listdir(wt_dir):
            if not fn.lower().endswith((".fa",".fasta",".faa")): continue
            for pat in patts.get(chain, []):
                if re.search(pat, fn, re.IGNORECASE):
                    cand.append(os.path.join(wt_dir, fn)); break
    if cand:
        cand.sort(key=lambda x: len(os.path.basename(x)))
        return cand[0]
    # 3) 固定名最终兜底
    fallback = {"H":"heavy.fasta","L":"light.fasta","A":"antigen.fasta"}
    return os.path.join(wt_dir, fallback.get(chain, "heavy.fasta"))

def load_wt(cfg, chain):
    path = resolve_wt_fasta(cfg, chain)
    if not os.path.exists(path):
        raise FileNotFoundError(f"找不到 {chain} 链 WT fasta：{path}")
    seq=""
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith(">"): continue
            seq += "".join([c for c in line.upper() if c.isalpha()])
    if not seq: raise ValueError(f"{chain} 链 WT 序列为空：{path}")
    bad=set(seq)-AA20
    if bad: raise ValueError(f"{chain} 链 WT 含非标准AA {bad}：{path}")
    print(f"[WT] chain={chain} file={path} length={len(seq)}")
    return seq

def read_hotspots(cfg):
    """读取界面热点位点：优先 his_bias.target_positions，否则用派生文件，再否则用 his_hotspots.csv TopK"""
    ch = cfg["design"]["chain"]
    region = cfg["design"]["region"]
    prefer_k = int(cfg["his_bias"].get("prefer_sites_topk", 20))

    # 1) 显式指定
    tp = cfg["his_bias"].get("target_positions", []) or []
    if tp:
        pos=[]
        for x in tp:
            # 兼容 "H:32" 或 "32"
            if ":" in x:
                c, n = x.split(":"); 
                if c.strip()!=ch: continue
                pos.append(int(n))
            else:
                pos.append(int(x))
        # 只保留在设计窗口内
        a,b = region
        pos = [p for p in pos if a<=p<=b]
        return sorted(set(pos))

    # 2) 派生位点 txt
    dfile = os.path.join(cfg["paths"]["results_dir"], "screening", "derived_his_positions.txt")
    if os.path.exists(dfile):
        pos=[]
        with open(dfile) as f:
            for line in f:
                line=line.strip()
                if not line: continue
                if ":" in line:
                    c,n = line.split(":")
                    if c.strip()!=ch: continue
                    pos.append(int(n))
        a,b = region
        pos = [p for p in pos if a<=p<=b]
        if pos: return sorted(set(pos))[:prefer_k]

    # 3) his_hotspots.csv TopK
    hcsv = os.path.join(cfg["paths"]["results_dir"], "screening", "his_hotspots.csv")
    if not os.path.exists(hcsv):
        raise FileNotFoundError("未找到 his_hotspots.csv；请先运行 scan_interface.py")
    df = pd.read_csv(hcsv)
    df = df[(df["chain"]==ch)]
    df = df.sort_values("hotness", ascending=False)
    a,b = region
    df = df[(df["resno"]>=a) & (df["resno"]<=b)]
    return df["resno"].head(prefer_k).tolist()

def mutate_to_h(seq, positions):
    """把给定 1-based positions 窗口内的位点改为 H"""
    s = list(seq)
    for p in positions:
        if s[p-1] != "H":
            s[p-1] = "H"
    return "".join(s)

def mk_mutcodes(wt, seq, chain, positions):
    muts=[]
    for p in positions:
        if wt[p-1] != seq[p-1]:
            muts.append(f"{wt[p-1]}{chain}{p}{seq[p-1]}")
    return ",".join(muts)

def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    out_dir = cfg["paths"]["his_seed_dir"]
    os.makedirs(out_dir, exist_ok=True)

    chain  = cfg["design"]["chain"]   # 你现在就是 "H"
    region = cfg["design"]["region"]
    a,b = region

    wt = load_wt(cfg, chain)
    if b > len(wt):
        raise ValueError(f"design.region 末端 {b} 超出 {chain} 链 WT 长度 {len(wt)}")

    # 取热点位点（已确保只在 H 链 & 设计窗口内）
    sites = read_hotspots(cfg)
    if not sites:
        raise RuntimeError("没有可用的热点位点（请检查 his_hotspots 或 target_positions 与 design.region 是否匹配）")
    print(f"[His sites] {sites}")

    # S1：单点变 H
    s1_rows=[]
    for p in sites:
        if wt[p-1] == "H": 
            continue
        seq = mutate_to_h(wt, [p])
        muts = mk_mutcodes(wt, seq, chain, [p])
        s1_rows.append({"sequence": seq, "mutations": muts, "source": "his_seed_S1"})

    # S2：成对变 H（受 pair_scan_limit 控制）
    pair_limit = int(cfg["his_bias"].get("pair_scan_limit", 1000))
    s2_rows=[]
    for i, (p1,p2) in enumerate(itertools.combinations(sites, 2)):
        if i >= pair_limit: break
        pos = [p1,p2]
        seq = mutate_to_h(wt, pos)
        muts = mk_mutcodes(wt, seq, chain, pos)
        # 若两位点本来就是 H，muts 为空，则跳过
        if not muts: 
            continue
        s2_rows.append({"sequence": seq, "mutations": muts, "source": "his_seed_S2"})

    s1 = pd.DataFrame(s1_rows)
    s2 = pd.DataFrame(s2_rows)
    f1 = os.path.join(out_dir, "his_seed_S1_single.csv")
    f2 = os.path.join(out_dir, "his_seed_S2_pairs.csv")
    s1.to_csv(f1, index=False)
    s2.to_csv(f2, index=False)
    print(f"[OK] seeds -> {f1} (n={len(s1)}), {f2} (n={len(s2)})")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv)>1 else "configs/config.yaml")