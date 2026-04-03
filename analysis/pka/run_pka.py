#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pKa 预测统一 wrapper — PROPKA3 + pKAI+ 双工具交叉验证

可复用模块，支持单文件和批量目录模式。
输出统一格式 CSV，包含双工具预测值、delta_pKa 和共识标注。

用法:
    # 单文件
    python run_pka.py --pdb structure.pdb

    # 批量目录
    python run_pka.py --pdb-dir ./pdbs/ -o results.csv

    # 带 WT 参考（计算 shift）
    python run_pka.py --pdb-dir ./pdbs/ --wt-pdb wt.pdb -o results.csv

    # 仅关注指定 His 位点
    python run_pka.py --pdb-dir ./pdbs/ --wt-pdb wt.pdb \
        --his-filter B:42,B:43,B:44 -o results.csv
"""

import argparse
import glob
import os
import sys
import warnings
from pathlib import Path

import pandas as pd

# His model pKa (solution value)
HIS_MODEL_PKA = 6.5

# ── PROPKA3 backend ─────────────────────────────────────────────────────────


def _run_propka(pdb_path: str) -> list[dict]:
    """Run PROPKA3 on a PDB file, return list of His pKa results."""
    try:
        from propka.run import single
    except ImportError:
        raise ImportError("propka 未安装，请运行: pip install propka")

    results = []
    # PROPKA 会在 cwd 生成 .pka 文件，切换到临时目录避免污染项目
    import tempfile
    original_cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                mol = single(pdb_path, optargs=["--quiet"])
            except Exception as e:
                os.chdir(original_cwd)
                print(f"  [PROPKA] 失败 {os.path.basename(pdb_path)}: {e}", file=sys.stderr)
                return results
        os.chdir(original_cwd)

    for conformation in mol.conformations:
        for group in mol.conformations[conformation].groups:
            if group.residue_type == "HIS":
                # label format: "HIS  61 A" → parse chain and resid
                parts = group.label.split()
                chain = parts[2] if len(parts) >= 3 else ""
                resid = int(parts[1]) if len(parts) >= 2 else 0
                results.append(
                    {
                        "chain": chain,
                        "resid": resid,
                        "resname": "HIS",
                        "pKa_propka": group.pka_value,
                    }
                )
        break  # 只取第一个构象
    return results


# ── pKAI+ backend ───────────────────────────────────────────────────────────


def _run_pkai(pdb_path: str) -> list[dict]:
    """Run pKAI+ on a PDB file, return list of His pKa results."""
    try:
        from pkai.pKAI import pKAI
    except ImportError:
        raise ImportError("pKAI 未安装，请运行: pip install pKAI")

    results = []
    try:
        # pKAI prints all residues to stdout; suppress
        import io
        import contextlib
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pks = pKAI(pdb_path)
    except Exception as e:
        print(f"  [pKAI+] 失败 {os.path.basename(pdb_path)}: {e}", file=sys.stderr)
        return results

    for chain, resnum, resname, pka_val in pks:
        if resname == "HIS":
            results.append(
                {
                    "chain": chain,
                    "resid": int(resnum),
                    "resname": "HIS",
                    "pKa_pkai": pka_val,
                }
            )
    return results


# ── 核心预测函数 ─────────────────────────────────────────────────────────────


def predict_pka(
    pdb_path: str,
    his_filter: list[tuple[str, int]] | None = None,
) -> pd.DataFrame:
    """
    对单个 PDB 文件运行 PROPKA3 + pKAI+ 双工具预测。

    Args:
        pdb_path: PDB 文件路径
        his_filter: 可选，仅保留指定 His 位点 [(chain, resid), ...]

    Returns:
        DataFrame with columns:
            pdb_name, chain, resid, resname,
            pKa_propka, pKa_pkai, pKa_model,
            delta_propka, delta_pkai
    """
    pdb_name = os.path.basename(pdb_path)

    # Run both tools
    propka_results = _run_propka(pdb_path)
    pkai_results = _run_pkai(pdb_path)

    # Merge on (chain, resid)
    propka_df = pd.DataFrame(propka_results) if propka_results else pd.DataFrame(
        columns=["chain", "resid", "resname", "pKa_propka"]
    )
    pkai_df = pd.DataFrame(pkai_results) if pkai_results else pd.DataFrame(
        columns=["chain", "resid", "resname", "pKa_pkai"]
    )

    if propka_df.empty and pkai_df.empty:
        return pd.DataFrame()

    # Outer merge to keep all His from both tools
    if not propka_df.empty and not pkai_df.empty:
        merged = pd.merge(
            propka_df, pkai_df, on=["chain", "resid", "resname"], how="outer"
        )
    elif not propka_df.empty:
        merged = propka_df
        merged["pKa_pkai"] = None
    else:
        merged = pkai_df
        merged["pKa_propka"] = None

    merged["pdb_name"] = pdb_name
    merged["pKa_model"] = HIS_MODEL_PKA
    merged["delta_propka"] = merged["pKa_propka"] - HIS_MODEL_PKA
    merged["delta_pkai"] = merged["pKa_pkai"] - HIS_MODEL_PKA

    # Apply filter
    if his_filter:
        mask = merged.apply(
            lambda r: (r["chain"], int(r["resid"])) in his_filter, axis=1
        )
        merged = merged[mask]

    cols = [
        "pdb_name", "chain", "resid", "resname",
        "pKa_propka", "pKa_pkai", "pKa_model",
        "delta_propka", "delta_pkai",
    ]
    return merged[cols].reset_index(drop=True)


def predict_pka_with_shift(
    pdb_path: str,
    wt_pka: pd.DataFrame,
    his_filter: list[tuple[str, int]] | None = None,
) -> pd.DataFrame:
    """
    预测单个突变体的 pKa，并与 WT 基线比较计算 shift。

    Args:
        pdb_path: 突变体 PDB 路径
        wt_pka: WT 的 predict_pka() 输出
        his_filter: 可选 His 位点过滤

    Returns:
        DataFrame with additional columns:
            shift_propka, shift_pkai, consensus
    """
    mut_pka = predict_pka(pdb_path, his_filter=his_filter)
    if mut_pka.empty:
        return mut_pka

    # Merge with WT on (chain, resid)
    wt_cols = wt_pka[["chain", "resid", "pKa_propka", "pKa_pkai"]].rename(
        columns={"pKa_propka": "wt_propka", "pKa_pkai": "wt_pkai"}
    )
    merged = pd.merge(mut_pka, wt_cols, on=["chain", "resid"], how="left")

    merged["shift_propka"] = merged["pKa_propka"] - merged["wt_propka"]
    merged["shift_pkai"] = merged["pKa_pkai"] - merged["wt_pkai"]

    # Consensus: both shifts in same direction
    def _consensus(row):
        sp = row["shift_propka"]
        sk = row["shift_pkai"]
        if pd.isna(sp) or pd.isna(sk):
            return "incomplete"
        if abs(sp) < 0.05 and abs(sk) < 0.05:
            return "neutral"
        if (sp > 0 and sk > 0) or (sp < 0 and sk < 0):
            return "agree"
        return "disagree"

    merged["consensus"] = merged.apply(_consensus, axis=1)

    # Drop intermediate WT columns
    merged.drop(columns=["wt_propka", "wt_pkai"], inplace=True)

    return merged


# ── 批量处理 ─────────────────────────────────────────────────────────────────


def batch_predict(
    pdb_paths: list[str],
    wt_pdb: str | None = None,
    his_filter: list[tuple[str, int]] | None = None,
) -> pd.DataFrame:
    """
    批量预测多个 PDB 文件的 His pKa。

    Args:
        pdb_paths: PDB 文件路径列表
        wt_pdb: 可选 WT PDB 路径（用于计算 shift）
        his_filter: 可选 His 位点过滤

    Returns:
        合并后的 DataFrame
    """
    # WT baseline
    wt_pka = None
    if wt_pdb:
        print(f"[WT] 计算基线: {os.path.basename(wt_pdb)}")
        wt_pka = predict_pka(wt_pdb, his_filter=his_filter)
        if wt_pka.empty:
            print("  警告: WT 未检测到 His 位点", file=sys.stderr)

    all_results = []
    for i, pdb_path in enumerate(sorted(pdb_paths), 1):
        pdb_name = os.path.basename(pdb_path)
        print(f"[{i}/{len(pdb_paths)}] {pdb_name}...", end=" ")

        if wt_pka is not None and not wt_pka.empty:
            result = predict_pka_with_shift(pdb_path, wt_pka, his_filter=his_filter)
        else:
            result = predict_pka(pdb_path, his_filter=his_filter)

        if result.empty:
            print("无 His")
        else:
            n_his = len(result)
            print(f"{n_his} His")

        all_results.append(result)

    if not all_results:
        return pd.DataFrame()

    return pd.concat(all_results, ignore_index=True)


# ── CLI ──────────────────────────────────────────────────────────────────────


def parse_his_filter(s: str) -> list[tuple[str, int]]:
    """Parse 'B:42,B:43,B:44' → [('B', 42), ('B', 43), ('B', 44)]"""
    result = []
    for item in s.split(","):
        chain, resid = item.strip().split(":")
        result.append((chain.strip(), int(resid.strip())))
    return result


def main():
    parser = argparse.ArgumentParser(
        description="pKa 预测 (PROPKA3 + pKAI+ 双工具交叉验证)"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb", help="单个 PDB 文件路径")
    group.add_argument("--pdb-dir", help="PDB 文件目录")

    parser.add_argument("--wt-pdb", help="WT PDB 路径（用于计算 shift）")
    parser.add_argument(
        "--his-filter",
        help="仅关注指定 His 位点，格式: B:42,B:43,B:44",
    )
    parser.add_argument("-o", "--output", help="输出 CSV 路径")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()

    # Collect PDB files
    if args.pdb:
        pdb_paths = [args.pdb]
    else:
        pdb_paths = sorted(glob.glob(os.path.join(args.pdb_dir, "*.pdb")))
        if not pdb_paths:
            print(f"错误: {args.pdb_dir} 下未找到 PDB 文件")
            sys.exit(1)

    his_filter = parse_his_filter(args.his_filter) if args.his_filter else None

    print(f"pKa 预测: {len(pdb_paths)} 个 PDB 文件")
    if his_filter:
        print(f"His 过滤: {his_filter}")
    if args.wt_pdb:
        print(f"WT 参考: {args.wt_pdb}")
    print()

    result = batch_predict(pdb_paths, wt_pdb=args.wt_pdb, his_filter=his_filter)

    if result.empty:
        print("未检测到任何 His 位点")
        sys.exit(0)

    # Output
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, index=False)
        print(f"\n结果已保存: {args.output}")
    else:
        print()
        print(result.to_string(index=False))


if __name__ == "__main__":
    main()
