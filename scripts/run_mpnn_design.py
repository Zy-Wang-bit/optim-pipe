#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, re, glob, json, yaml, argparse, subprocess
import pandas as pd
from collections import OrderedDict

def sh(cmd, cwd=None, env=None):
    env2 = dict(os.environ, **(env or {}))
    print("[RUN]", " ".join(cmd))
    p = subprocess.run(cmd, cwd=cwd, env=env2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        print(p.stdout)
        print(p.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return p

def extract_chain_res_order(pdb_path, chain_id):
    """按出现顺序收集该链的残基编号（PDB resSeq+iCode），用于把 region=76..115 映射到1..N索引。"""
    chain_id = str(chain_id).strip()
    seen = OrderedDict()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 27:
                continue
            ch = line[21].strip()
            if ch != chain_id:
                continue
            resseq = line[22:26].strip()
            icode  = line[26].strip()
            key = (resseq, icode)
            if key not in seen:
                seen[key] = None
    res_list = list(seen.keys())
    return res_list  # 长度 = 该链的残基数，1-based索引就是枚举顺序

def _parse_design_chains(cfg):
    """解析 design 配置，支持单链和多链。

    单链格式: design.chain="A", design.region=[1,115]
    多链格式: design.chain=["A","B"], design.regions={"A":[1,115],"B":[1,113]}

    返回: {chain_id: (start, end), ...}
    """
    chain = cfg["design"]["chain"]
    if isinstance(chain, str):
        # 单链
        region = cfg["design"]["region"]
        return {chain: (region[0], region[1])}
    else:
        # 多链
        regions = cfg["design"]["regions"]
        return {ch: (regions[ch][0], regions[ch][1]) for ch in chain}


def build_jsonl_for_mpnn(cfg, pdb_dir, design_chains, out_dir):
    """
    生成 MPNN 所需的 JSONL 文件，支持多链同时设计。

    design_chains: {chain_id: (start, end), ...}
    """
    import os, glob, json
    from collections import OrderedDict

    def get_chain_len_map(pdb_path):
        order = {}
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                ch = line[21].strip()
                if not ch:
                    continue
                key = (ch, line[22:26].strip(), line[26].strip())
                order.setdefault(ch, OrderedDict()).setdefault(key, None)
        return {ch: len(od) for ch, od in order.items()}

    os.makedirs(out_dir, exist_ok=True)
    chain_jsonl = os.path.join(out_dir, "chain_id.jsonl")
    fixed_jsonl = os.path.join(out_dir, "fixed_positions.jsonl")

    design_chain_ids = list(design_chains.keys())

    # 1) chain_id.jsonl
    CJ = []
    for pdb in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
        base = os.path.splitext(os.path.basename(pdb))[0]
        abs_p = os.path.abspath(pdb)
        chain_tag = "+".join(design_chain_ids)
        aliases = [
            base, f"{base}.pdb", abs_p,
            f"{base}_{chain_tag}", f"{base}.pdb_{chain_tag}", f"{abs_p}_{chain_tag}",
        ]
        for name in aliases:
            CJ.append({"name": name, "chain_ids": design_chain_ids})
    with open(chain_jsonl, "w") as f:
        for x in CJ:
            f.write(json.dumps(x) + "\n")

    # 2) fixed_positions.jsonl
    FP = {}
    for pdb in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
        base = os.path.splitext(os.path.basename(pdb))[0]
        abs_p = os.path.abspath(pdb)
        chain_tag = "+".join(design_chain_ids)
        aliases = [
            base, f"{base}.pdb", abs_p,
            f"{base}_{chain_tag}", f"{base}.pdb_{chain_tag}", f"{abs_p}_{chain_tag}",
        ]
        len_map = get_chain_len_map(pdb)

        for ch_id in design_chain_ids:
            if ch_id not in len_map:
                raise ValueError(f"[{base}] PDB 中找不到链 {ch_id}")

        fixed_map = {}
        for ch, L in len_map.items():
            if ch in design_chains:
                s, e = design_chains[ch]
                if e > L:
                    raise ValueError(f"[{base}] region={design_chains[ch]} 超出链 {ch} 长度 N={L}")
                fixed_map[ch] = list(range(1, s)) + list(range(e+1, L+1))
            else:
                fixed_map[ch] = list(range(1, L+1))

        for name in aliases:
            FP[name] = fixed_map

    with open(fixed_jsonl, "w") as f:
        f.write(json.dumps(FP))

    chain_jsonl = os.path.abspath(chain_jsonl)
    fixed_jsonl = os.path.abspath(fixed_jsonl)
    print(f"[OK] JSONL -> {chain_jsonl}  /  {fixed_jsonl}")
    print(f"     设计链: {design_chains}")
    return chain_jsonl, fixed_jsonl

def find_fasta_files(folder):
    fas = glob.glob(os.path.join(folder, "*.fa")) + \
          glob.glob(os.path.join(folder, "*.fasta")) + \
          glob.glob(os.path.join(folder, "seqs", "*.fa")) + \
          glob.glob(os.path.join(folder, "seqs", "*.fasta"))
    return fas

def load_fasta(fp):
    seqs=[]; s=""
    with open(fp) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if s: seqs.append(s); s=""
            else:
                s+=line
    if s: seqs.append(s)
    return seqs

def mpnn_to_csv_per_pid(pid_dir, pid, out_csv):
    rows=[]
    for shard_dir in sorted(glob.glob(os.path.join(pid_dir, "shard_*"))):
        shard_name = os.path.basename(shard_dir)
        # 解析温度与seed（写在 marker.json）
        meta = {}
        mpath = os.path.join(shard_dir,"marker.json")
        if os.path.exists(mpath):
            try:
                meta = json.load(open(mpath))
            except Exception:
                meta = {}
        fas = find_fasta_files(shard_dir)
        for fa in fas:
            for s in load_fasta(fa):
                rows.append({
                    "pdb_id": pid,
                    "sequence": s,
                    "source": "mpnn",
                    "shard": shard_name,
                    "sampling_temp": meta.get("sampling_temp", None),
                    "seed": meta.get("seed", None),
                })
    if rows:
        df = pd.DataFrame(rows).drop_duplicates(subset=["pdb_id","sequence"])
        df.to_csv(out_csv, index=False)
        print(f"[MPNN->CSV] {pid}: {len(df)} -> {out_csv}")
    else:
        print(f"[WARN] {pid_dir} 未找到FASTA，跳过写出")

def main():
    ap = argparse.ArgumentParser(description="Run ProteinMPNN by window and standardize outputs.")
    ap.add_argument("config", help="configs/config.yaml")
    ap.add_argument("--runner", default=None, help="ProteinMPNN 入口脚本（默认从 config.paths.mpnn_runner 或 PATH 中的 protein_mpnn_run.py）")
    ap.add_argument("--num-per-pdb", type=int, default=200000, help="每个PDB生成的序列数")
    ap.add_argument("--shards", type=int, default=10, help="分成多少个 shard 运行")
    ap.add_argument("--temps", default="0.10,0.15,0.20,0.25,0.30", help="采样温度列表，用逗号分隔")
    ap.add_argument("--seed", type=int, default=1234, help="随机种子基数（每个 shard 会加上 i）")
    ap.add_argument("--build-his-seeds", action="store_true", help="顺带执行热点扫描与 His 种子生成")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    pdb_dir   = cfg["paths"]["pdb_dir"]
    out_root  = cfg["paths"]["mpnn_out_dir"]
    res_dir   = cfg["paths"]["results_dir"]
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(os.path.join(res_dir,"screening"), exist_ok=True)

    design_chains = _parse_design_chains(cfg)        # {chain: (start, end), ...}
    chain_id  = cfg["design"]["chain"]                # 兼容旧代码（可能是 str 或 list）

    runner = args.runner or cfg["paths"].get("mpnn_runner") or "protein_mpnn_run.py"
    temps  = [float(x) for x in args.temps.split(",") if x.strip()]
    shards = int(args.shards)
    n_per  = int(args["num_per_pdb"]) if isinstance(args, dict) and "num_per_pdb" in args else int(args.num_per_pdb)
    per_shard = max(1, n_per // shards)

    # max_len 检查（多链时检查总长度）
    max_len = cfg["design"].get("max_len", 40)
    total_win = sum(e - s + 1 for s, e in design_chains.values())
    if total_win > max_len:
        raise ValueError(f"设计窗口总长度 {total_win} 超出 max_len {max_len}，"
                         f"请调大 design.max_len 或缩小 region")

    # 1) JSONL（准确链长）— 使用 mpnn_out_dir 下的子目录避免并行冲突
    jsonl_dir = os.path.join(out_root, "_mpnn_config")
    chain_jsonl, fixed_jsonl = build_jsonl_for_mpnn(cfg, pdb_dir, design_chains, jsonl_dir)

    # 2) 逐 PDB 逐 shard 跑 MPNN
    for pdb in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
        pid = os.path.splitext(os.path.basename(pdb))[0]
        pdb_abs = os.path.abspath(pdb)
        pid_out = os.path.join(out_root, pid)
        os.makedirs(pid_out, exist_ok=True)

        for i in range(shards):
            temp = temps[i % len(temps)]
            seed = args.seed + i
            shard_dir = os.path.join(pid_out, f"shard_{i:03d}")
            os.makedirs(shard_dir, exist_ok=True)

            # 记一个标记，便于之后汇总写入元信息
            json.dump({"sampling_temp":temp, "seed":seed, "per_shard":per_shard},
                      open(os.path.join(shard_dir,"marker.json"),"w"))

            cmd = [
                sys.executable, runner,
                "--pdb_path", pdb_abs,
                "--out_folder", shard_dir,
                "--chain_id_jsonl", chain_jsonl,
                "--fixed_positions_jsonl", fixed_jsonl,
                "--num_seq_per_target", str(per_shard),
                "--sampling_temp", str(temp),
                "--seed", str(seed),
            ]
            sh(cmd)

        # 3) 汇总该 PDB 的 MPNN 输出为 CSV
        out_csv = os.path.join(pid_out, f"mpnn_{pid}.csv")
        mpnn_to_csv_per_pid(pid_out, pid, out_csv)

    # 4) 汇总所有 PDB 的 CSV 为一个总表（方便后续 ESM/pick 脚本读取）
    rows=[]
    for csvp in glob.glob(os.path.join(out_root, "*", "mpnn_*.csv")):
        rows.append(pd.read_csv(csvp))
    if rows:
        raw_df = pd.concat(rows, ignore_index=True)
        n_raw = len(raw_df)
        # 全局按 sequence 去重（保留首次出现的 pdb_id）
        all_df = raw_df.drop_duplicates(subset=["sequence"], keep="first")
        n_dedup = len(all_df)
        out_all = os.path.join(res_dir, "screening", "mpnn_all.csv")
        all_df.to_csv(out_all, index=False)
        dup_pct = (1 - n_dedup / n_raw) * 100 if n_raw > 0 else 0
        print(f"[OK] mpnn_all.csv -> {out_all}  原始={n_raw}  去重后={n_dedup}  重复率={dup_pct:.1f}%")
    else:
        print("[WARN] 未合并到任何 MPNN 序列")

    # 5) （可选）顺带执行热点扫描与 His 种子生成
    if args.build_his_seeds:
        sh([sys.executable, "scripts/scan_interface.py", args.config])
        sh([sys.executable, "scripts/build_his_seeds.py", args.config])
        print("[OK] 已生成 His 种子（his_seeds/*）")

if __name__ == "__main__":
    main()