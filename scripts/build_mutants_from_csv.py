#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import csv
import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd


def run_foldx(foldx_bin, args, cwd, env=None, log_prefix=None):
    env2 = dict(os.environ, **(env or {}))
    env2.setdefault("OMP_NUM_THREADS", "1")
    cmd = [foldx_bin] + args
    result = subprocess.run(
        cmd,
        cwd=cwd,
        env=env2,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if log_prefix:
        Path(log_prefix + ".out").write_text(result.stdout or "", encoding="utf-8")
        Path(log_prefix + ".err").write_text(result.stderr or "", encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"FoldX 失败 (code={result.returncode}): {' '.join(cmd)}\n{result.stderr}"
        )
    return result


def parse_mutation_token(token):
    token = token.strip()
    if not token:
        return None
    m = re.match(r"^([A-Za-z])_([A-Za-z])(\d+)([A-Za-z])$", token)
    if not m:
        raise ValueError(f"无法解析突变标记: {token}")
    chain, orig, pos, new = m.groups()
    return f"{orig}{chain}{pos}{new}"


def mutations_to_foldx(mutation_str):
    if mutation_str is None or (isinstance(mutation_str, float) and pd.isna(mutation_str)):
        return None
    items = [parse_mutation_token(t) for t in str(mutation_str).split(",")]
    items = [x for x in items if x]
    if not items:
        return None
    return ",".join(items) + ";"


def sort_mutant_files(files):
    def key_fn(p: Path):
        m = re.search(r"_(\d+)\.pdb$", p.name)
        return int(m.group(1)) if m else 10**9
    return sorted(files, key=key_fn)


def main():
    parser = argparse.ArgumentParser(description="用 FoldX 由 CSV 批量生成突变体 PDB")
    parser.add_argument(
        "--csv",
        default="/public/home/ziyang/code/optim-pipe/experiments/1E62/R2/1E62_R2_热点区域推荐_突变序列.csv",
        help="输入CSV路径",
    )
    parser.add_argument(
        "--wt",
        default="/public/home/ziyang/code/optim-pipe/experiments/1E62/data/AF3-abag-n1.pdb",
        help="野生型PDB路径",
    )
    parser.add_argument(
        "--foldx",
        default="/public/home/ziyang/code/optim-pipe/foldX/foldx",
        help="FoldX可执行文件路径",
    )
    parser.add_argument(
        "--out-dir",
        default="/public/home/ziyang/code/optim-pipe/experiments/1E62/R2/foldx_mutant_pdbs",
        help="输出目录",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    wt_path = Path(args.wt)
    foldx_bin = args.foldx
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV不存在: {csv_path}")
    if not wt_path.exists():
        raise FileNotFoundError(f"WT PDB不存在: {wt_path}")
    if not Path(foldx_bin).exists():
        raise FileNotFoundError(f"FoldX不存在: {foldx_bin}")

    work_dir = out_dir / "foldx_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    # 复制WT到工作目录
    wt_name = wt_path.name
    wt_in_work = work_dir / wt_name
    if not wt_in_work.exists():
        shutil.copy2(wt_path, wt_in_work)

    # 读取CSV
    df = pd.read_csv(csv_path)
    if "突变" not in df.columns:
        raise ValueError("CSV缺少列: 突变")

    # 生成FoldX突变格式
    df["foldx_mutations"] = df["突变"].apply(mutations_to_foldx)
    df_valid = df[df["foldx_mutations"].notna()].copy()

    if df_valid.empty:
        raise RuntimeError("没有可用的突变记录")

    # RepairPDB
    wt_stem = wt_path.stem
    repair_log = str(work_dir / "repair")
    run_foldx(
        foldx_bin,
        [
            "--command=RepairPDB",
            f"--pdb={wt_name}",
            f"--output-file={wt_stem}",
        ],
        cwd=str(work_dir),
        log_prefix=repair_log,
    )

    repaired_name = f"{wt_stem}_Repair.pdb"
    repaired_path = work_dir / repaired_name
    if not repaired_path.exists():
        raise FileNotFoundError(f"RepairPDB输出未找到: {repaired_path}")

    # 写入 individual_list.txt
    mut_list_path = work_dir / "individual_list.txt"
    with open(mut_list_path, "w", encoding="utf-8") as f:
        for mut in df_valid["foldx_mutations"].tolist():
            f.write(mut + "\n")

    # BuildModel
    build_log = str(work_dir / "build")
    run_foldx(
        foldx_bin,
        [
            "--command=BuildModel",
            f"--pdb={repaired_name}",
            "--mutant-file=individual_list.txt",
            "--pH=7.4",
            "--numberOfRuns=1",
        ],
        cwd=str(work_dir),
        log_prefix=build_log,
    )

    # 收集突变体PDB并按顺序重命名
    mutant_files = [
        p for p in work_dir.glob(f"{wt_stem}_Repair_*.pdb")
        if not p.name.startswith("WT_") and not p.name.endswith("_Repair.pdb")
    ]
    mutant_files = sort_mutant_files(mutant_files)

    if not mutant_files:
        raise RuntimeError("未生成任何突变体PDB，请检查FoldX输出")

    # 输出映射
    map_csv = out_dir / "mutant_map.csv"
    with open(map_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["序号", "突变", "foldx_mutations", "mutant_pdb"])

        for idx, (_, row) in enumerate(df_valid.iterrows()):
            if idx >= len(mutant_files):
                break
            src = mutant_files[idx]
            mut_id = row.get("序号", idx + 1)
            try:
                mut_id_int = int(mut_id)
            except Exception:
                mut_id_int = idx + 1
            dst_name = f"mut_{mut_id_int:03d}.pdb"
            dst = out_dir / dst_name
            shutil.copy2(src, dst)
            writer.writerow([row.get("序号", idx + 1), row.get("突变", ""), row["foldx_mutations"], dst_name])

    print(f"[OK] 生成突变体PDB: {out_dir}")
    print(f"[OK] 映射表: {map_csv}")
    print(f"[OK] 工作目录: {work_dir}")


if __name__ == "__main__":
    main()
