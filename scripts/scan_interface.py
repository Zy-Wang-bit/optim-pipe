#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
扫描 PDB 复合物的抗体-抗原接口位点，统计“热度”(hotness)。
- 从 configs/config.yaml 读取：
  - paths.pdb_dir：PDB 目录（扫描 *.pdb）
  - paths.results_dir：输出目录（results/screening）
  - interface.ab_chains：抗体链集合（默认 ["H","L"] 或根据你的 A/B）
  - interface.ag_chains：抗原链集合（默认 ["A"] 或根据你的 C）
  - interface.cutoff：距离阈值 Å（默认 5.0）
  - design.chain：后续要设计的链（如 "A"）
  - design.region：窗口 [s, e]（仅用于派生位点）
  - his_bias.target_positions：若非空，则不生成 derived_his_positions.txt
  - his_bias.prefer_sites_topk：派生位点 Top-K（默认 20）
- 输出：
  1) results/screening/his_hotspots.csv    (chain,resno,hotness)
  2) results/screening/derived_his_positions.txt  (仅当 target_positions 为空时；每行 "A:32" 这类)
"""

import os
import sys
import glob
import yaml
import pandas as pd
from collections import defaultdict
from Bio.PDB import PDBParser

def load_cfg(cfg_path: str):
    if not os.path.exists(cfg_path):
        raise FileNotFoundError(f"找不到配置文件：{cfg_path}")
    with open(cfg_path, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg

def get_cfg_val(cfg, keys, default=None):
    """安全取值：keys 为 'a.b.c' 字符串或列表"""
    if isinstance(keys, str):
        keys = keys.split(".")
    cur = cfg
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur

def is_std_residue(res):
    """仅保留标准氨基酸（排除水/配体/插入码标识 HETATM）"""
    hetflag, resseq, icode = res.get_id()
    if hetflag != " ":  # 非标准残基（含 HOH、配体等）
        return False
    # 可按需过滤非蛋白残基名；这里允许 FoldX 常见的 20AA 三字母
    return True

def residue_min_distance(res_a, res_b, cutoff=5.0):
    """返回两个残基最小原子-原子距离（若全部>cutoff，可提前返回>cutoff的值）"""
    min_d = float("inf")
    for at_a in res_a.get_atoms():
        coord_a = at_a.get_coord()
        for at_b in res_b.get_atoms():
            # 小的优化：先做一个快速 L1/L2 近似裁剪可加速；这里直接算欧氏距离
            diff = coord_a - at_b.get_coord()
            d2 = (diff * diff).sum()  # numpy array 点乘
            if d2 < min_d * min_d:
                # 只有可能比当前更小才开方
                d = float(d2) ** 0.5
                if d < min_d:
                    min_d = d
                    if min_d <= cutoff:
                        return min_d
    return min_d

def scan_interface_for_pdb(pdb_path, ab_chains, ag_chains, cutoff):
    """返回本 PDB 中满足距离阈值的抗体侧接口残基集合：{(chain_id, resno_int)}"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_path), pdb_path)

    # 只看第一个 model（AF/复合物通常只有一个）
    model = next(structure.get_models())

    # 收集链
    chain_map = {c.id: c for c in model.get_chains()}
    ab_list = [chain_map[cid] for cid in ab_chains if cid in chain_map]
    ag_list = [chain_map[cid] for cid in ag_chains if cid in chain_map]

    if not ab_list or not ag_list:
        # 链名不匹配时给出提示，但不中断
        missing_ab = [cid for cid in ab_chains if cid not in chain_map]
        missing_ag = [cid for cid in ag_chains if cid not in chain_map]
        msg = []
        if missing_ab: msg.append(f"缺少抗体链: {missing_ab}")
        if missing_ag: msg.append(f"缺少抗原链: {missing_ag}")
        sys.stderr.write(f"[WARN] {os.path.basename(pdb_path)}: {', '.join(msg)}\n")

    # 预取所有残基
    ag_residues = []
    for ch in ag_list:
        for res in ch.get_residues():
            if is_std_residue(res):
                ag_residues.append((ch.id, res))

    interface_res = set()
    for ch in ab_list:
        for res in ch.get_residues():
            if not is_std_residue(res):
                continue
            # 与所有抗原残基求最小距离（提前剪枝：若某对 <= cutoff 就判定为接口）
            for _, ag_res in ag_residues:
                dmin = residue_min_distance(res, ag_res, cutoff=cutoff)
                if dmin <= cutoff:
                    # 记录为（抗体侧）接口位点
                    # 注意：res id 为 (' ', resseq, icode)
                    resno = res.get_id()[1]  # 插入码将被忽略（例如 100A->100）
                    interface_res.add((ch.id, int(resno)))
                    break  # 该残基已被判定为接口，无需再比其它抗原残基
    return interface_res

def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "configs/config.yaml"
    cfg = load_cfg(cfg_path)

    pdb_dir    = get_cfg_val(cfg, "paths.pdb_dir")
    results_dir = get_cfg_val(cfg, "paths.results_dir", "results")
    screen_dir = os.path.join(results_dir, "screening")
    os.makedirs(screen_dir, exist_ok=True)

    iface_ab   = get_cfg_val(cfg, "interface.ab_chains", ["H","L"])
    iface_ag   = get_cfg_val(cfg, "interface.ag_chains", ["A"])
    cutoff     = float(get_cfg_val(cfg, "interface.cutoff", 5.0))

    # 设计链与窗口（用于派生位点）
    design_chain = get_cfg_val(cfg, "design.chain", "H")
    design_region = get_cfg_val(cfg, "design.region", [1, 9999])
    region_s, region_e = int(design_region[0]), int(design_region[1])

    # his_bias
    his_bias = get_cfg_val(cfg, "his_bias", {}) or {}
    target_positions = his_bias.get("target_positions", []) or []
    prefer_topk = int(his_bias.get("prefer_sites_topk", 20))

    # 扫描所有 PDB
    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")))
    if not pdb_files:
        raise FileNotFoundError(f"在 {pdb_dir} 下没有找到 .pdb 文件")

    print(f"[scan] PDB 数量: {len(pdb_files)} | 抗体链: {iface_ab} | 抗原链: {iface_ag} | cutoff={cutoff}Å")
    hotspots = defaultdict(int)  # (chain, resno) -> hotness(出现的PDB计数)

    for pdb in pdb_files:
        inter = scan_interface_for_pdb(pdb, iface_ab, iface_ag, cutoff)
        for key in inter:
            hotspots[key] += 1
        print(f"[scan] {os.path.basename(pdb)} -> 接口位点(抗体侧) {len(inter)} 个")

    # 汇总为 DataFrame
    if not hotspots:
        sys.stderr.write("[WARN] 所有 PDB 未检测到接口位点，请检查链名/距离阈值。\n")

    rows = [{"chain": ch, "resno": resno, "hotness": cnt}
            for (ch, resno), cnt in hotspots.items()]
    df = pd.DataFrame(rows, columns=["chain", "resno", "hotness"]).sort_values(
        ["hotness", "chain", "resno"], ascending=[False, True, True]
    )

    out_csv = os.path.join(screen_dir, "his_hotspots.csv")
    df.to_csv(out_csv, index=False)
    print(f"[OK] 写出热点表 -> {out_csv} (rows={len(df)})")

    # 若未显式指定 target_positions，则派生一份当前设计链+窗口的 Top-K 位点清单
    if not target_positions:
        sub = df[(df["chain"] == design_chain) &
                 (df["resno"] >= region_s) &
                 (df["resno"] <= region_e)].copy()
        sub = sub.sort_values(["hotness", "resno"], ascending=[False, True])
        top_sites = sub["resno"].head(prefer_topk).tolist()

        out_txt = os.path.join(screen_dir, "derived_his_positions.txt")
        with open(out_txt, "w") as f:
            for p in top_sites:
                f.write(f"{design_chain}:{int(p)}\n")
        print(f"[OK] 写出派生位点 -> {out_txt} (Top-{prefer_topk} within {design_chain}[{region_s},{region_e}])")
    else:
        print("[skip] his_bias.target_positions 已指定，跳过派生位点生成。")

if __name__ == "__main__":
    main()