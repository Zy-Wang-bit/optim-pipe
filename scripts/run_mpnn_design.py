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

def build_jsonl_for_mpnn(cfg, pdb_dir, chain_id, region, out_dir):
    """
    生成:
      - chain_id.jsonl: 每行 {"name": <别名>, "chain_ids": [<目标链>]}
      - fixed_positions.jsonl: 单个 JSON 映射:
            { "<alias>": {"A":[...], "B":[...], "C":[...]}, ... }
        其中:
          * 目标链 (chain_id)：窗口外索引固定
          * 其余存在的链：整条链 1..N 固定（确保不被设计）
    """
    import os, glob, json
    from collections import OrderedDict

    def get_chain_len_map(pdb_path):
        """返回 {chain_letter: length}，按 PDB 中实际存在的链统计"""
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

    # 1) chain_id.jsonl：逐行写（你的分支能正确加载这类格式）
    CJ = []
    s, e = region
    for pdb in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
        base = os.path.splitext(os.path.basename(pdb))[0]
        abs_p = os.path.abspath(pdb)
        # 多种别名，最大兼容
        aliases = [
            base, f"{base}.pdb", abs_p,
            f"{base}_{chain_id}", f"{base}.pdb_{chain_id}", f"{abs_p}_{chain_id}",
        ]
        for name in aliases:
            CJ.append({"name": name, "chain_ids": [chain_id]})
    with open(chain_jsonl, "w") as f:
        for x in CJ:
            f.write(json.dumps(x) + "\n")

    # 2) fixed_positions.jsonl：写成单个 JSON 映射，覆盖 A/B/C
    FP = {}
    for pdb in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
        base = os.path.splitext(os.path.basename(pdb))[0]
        abs_p = os.path.abspath(pdb)
        aliases = [
            base, f"{base}.pdb", abs_p,
            f"{base}_{chain_id}", f"{base}.pdb_{chain_id}", f"{abs_p}_{chain_id}",
        ]
        len_map = get_chain_len_map(pdb)  # 例如 {"A":122, "B":110, "C":241}
        if chain_id not in len_map:
            raise ValueError(f"[{base}] 在 PDB 中找不到链 {chain_id}")

        # 为所有存在的链准备固定位点
        fixed_map = {}
        for ch, L in len_map.items():
            if ch == chain_id:
                if e > L:
                    raise ValueError(f"[{base}] region={region} 超出链 {chain_id} 长度 N={L}")
                # 目标链: 仅窗口外固定 -> 允许窗口内可变
                fixed_map[ch] = list(range(1, s)) + list(range(e+1, L+1))
            else:
                # 其它链: 全部固定 (1..L)，即不允许任何位置被设计
                fixed_map[ch] = list(range(1, L+1))

        for name in aliases:
            FP[name] = fixed_map

    with open(fixed_jsonl, "w") as f:
        f.write(json.dumps(FP))

    chain_jsonl = os.path.abspath(chain_jsonl)
    fixed_jsonl = os.path.abspath(fixed_jsonl)
    print("[OK] JSONL ->", chain_jsonl, " / ", fixed_jsonl)
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

    chain_id  = cfg["design"]["chain"]              # A/B/C
    region    = cfg["design"]["region"]             # [76,115]

    runner = args.runner or cfg["paths"].get("mpnn_runner") or "protein_mpnn_run.py"
    temps  = [float(x) for x in args.temps.split(",") if x.strip()]
    shards = int(args.shards)
    n_per  = int(args["num_per_pdb"]) if isinstance(args, dict) and "num_per_pdb" in args else int(args.num_per_pdb)
    per_shard = max(1, n_per // shards)

    # 1) JSONL（准确链长）
    jsonl_dir = "configs/mpnn"
    chain_jsonl, fixed_jsonl = build_jsonl_for_mpnn(cfg, pdb_dir, chain_id, region, jsonl_dir)

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
        all_df = pd.concat(rows, ignore_index=True).drop_duplicates(subset=["pdb_id","sequence"])
        out_all = os.path.join(res_dir, "screening", "mpnn_all.csv")
        all_df.to_csv(out_all, index=False)
        print("[OK] mpnn_all.csv ->", out_all, "N=", len(all_df))
    else:
        print("[WARN] 未合并到任何 MPNN 序列")

    # 5) （可选）顺带执行热点扫描与 His 种子生成
    if args.build_his_seeds:
        sh([sys.executable, "scripts/scan_interface.py", args.config])
        sh([sys.executable, "scripts/build_his_seeds.py", args.config])
        print("[OK] 已生成 His 种子（his_seeds/*）")

if __name__ == "__main__":
    main()