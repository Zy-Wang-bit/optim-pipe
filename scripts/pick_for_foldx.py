#!/usr/bin/env python3
import os, sys, glob, yaml
import pandas as pd

def list_pdb_ids(pdb_dir):
    ids=[]
    for fn in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
        ids.append(os.path.splitext(os.path.basename(fn))[0])
    return ids

def add_pdb_round_robin(df, pdb_ids):
    # 将没有 pdb_id 的行按轮询方式分配
    mask = ~df.columns.isin(["pdb_id"])
    if "pdb_id" in df.columns:
        miss = df["pdb_id"].isna() | (df["pdb_id"]=="")
    else:
        df["pdb_id"] = ""
        miss = df["pdb_id"]==""
    idx = df[miss].index.tolist()
    for i,k in enumerate(idx):
        df.at[k,"pdb_id"] = pdb_ids[i % len(pdb_ids)]
    return df

def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    esm_dir = cfg["paths"]["esm_out_dir"]
    out_csv = os.path.join(cfg["paths"]["results_dir"], "screening", "for_foldx.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    # 合并所有 ESM 结果
    rows=[]
    for fn in sorted(glob.glob(os.path.join(esm_dir, "*_esm.csv"))):
        df = pd.read_csv(fn)
        rows.append(df)
    big = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["sequence","esm_avg_logprob","source","pdb_id","mutations"])

    # 补齐/分配 pdb_id（针对 his_seeds 或旧格式）
    pdb_ids = list_pdb_ids(cfg["paths"]["pdb_dir"])
    big = add_pdb_round_robin(big, pdb_ids)

    # 每个 pdb_id 内取前 keep_top_pct_per_pdb
    keep_pct = float(cfg["screening"]["keep_top_pct_per_pdb"])
    kept=[]
    for pid, sub in big.groupby("pdb_id"):
        n = max(1, int(len(sub)*keep_pct))
        kept.append(sub.sort_values("esm_avg_logprob", ascending=False).head(n))
    kept = pd.concat(kept, ignore_index=True)
    kept.to_csv(out_csv, index=False)
    print(f"[OK] {out_csv} (rows={len(kept)})")

if __name__=="__main__":
    cfg = sys.argv[1] if len(sys.argv)>1 else "configs/config.yaml"
    main(cfg)