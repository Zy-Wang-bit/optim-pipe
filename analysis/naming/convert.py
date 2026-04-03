#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
突变位点命名格式转换工具

统一格式: {Chain}{OrigAA}{Position}{NewAA}
  - Chain: H (重链/Heavy) 或 L (轻链/Light)
  - 例: LS42H = 轻链42位 Ser→His, HE1H = 重链1位 Glu→His

位点引用格式: {Chain}{Position}
  - 例: H42 = 重链42位, L42 = 轻链42位

PDB链映射 (1E62 体系):
  PDB chain A → H (重链)
  PDB chain B → L (轻链)
  PDB chain C → 靶标 (不变)
"""

import re

# 默认 PDB 链号 → H/L 映射 (1E62 体系)
DEFAULT_PDB_TO_HL = {"A": "H", "B": "L"}
DEFAULT_HL_TO_PDB = {"H": "A", "L": "B"}


def pdb_chain_to_hl(chain, mapping=None):
    """PDB 链号转为 H/L 标识。 'A' → 'H', 'B' → 'L'"""
    m = mapping or DEFAULT_PDB_TO_HL
    return m.get(chain, chain)


def hl_to_pdb_chain(chain, mapping=None):
    """H/L 标识转为 PDB 链号。 'H' → 'A', 'L' → 'B'"""
    m = mapping or DEFAULT_HL_TO_PDB
    return m.get(chain, chain)


def position_label(chain, resid, mapping=None):
    """生成统一位点标签。 ('B', 42) → 'L42', ('A', 48) → 'H48'"""
    return f"{pdb_chain_to_hl(chain, mapping)}{resid}"


# ── R1 格式转换 ─────────────────────────────────────────────────────────────

def r1_to_unified(s):
    """R1 旧格式 → 统一格式。 'A_E1H' → 'HE1H', 'B_L52H' → 'LL52H'"""
    m = re.match(r"^([AB])_([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析 R1 格式: {s}")
    chain, orig, pos, new = m.groups()
    return f"{pdb_chain_to_hl(chain)}{orig}{pos}{new}"


def unified_to_r1(s):
    """统一格式 → R1 旧格式。 'HE1H' → 'A_E1H'"""
    m = re.match(r"^([HL])([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析统一格式: {s}")
    chain, orig, pos, new = m.groups()
    return f"{hl_to_pdb_chain(chain)}_{orig}{pos}{new}"


# ── R3 箭头格式转换 ──────────────────────────────────────────────────────────

def r3_arrow_to_unified(s):
    """R3 箭头格式 → 统一格式。 'H18L>I' → 'HL18I', 'L42S>H' → 'LS42H'"""
    m = re.match(r"^([HL])(\d+)([A-Z])>([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析 R3 箭头格式: {s}")
    chain, pos, orig, new = m.groups()
    return f"{chain}{orig}{pos}{new}"


def unified_to_r3_arrow(s):
    """统一格式 → R3 箭头格式。 'HL18I' → 'H18L>I'"""
    m = re.match(r"^([HL])([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析统一格式: {s}")
    chain, orig, pos, new = m.groups()
    return f"{chain}{pos}{orig}>{new}"


# ── FoldX 格式转换 ───────────────────────────────────────────────────────────

def foldx_to_unified(s, mapping=None):
    """FoldX 格式 → 统一格式。 'SA40A' → 'HS40A' (A链=重链)"""
    m = re.match(r"^([A-Z])([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析 FoldX 格式: {s}")
    orig, chain, pos, new = m.groups()
    return f"{pdb_chain_to_hl(chain, mapping)}{orig}{pos}{new}"


def unified_to_foldx(s, mapping=None):
    """统一格式 → FoldX 格式。 'HS40A' → 'SA40A'"""
    m = re.match(r"^([HL])([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析统一格式: {s}")
    chain, orig, pos, new = m.groups()
    pdb_chain = hl_to_pdb_chain(chain, mapping)
    return f"{orig}{pdb_chain}{pos}{new}"


# ── R2 下划线格式转换 ────────────────────────────────────────────────────────

def r2_underscore_to_unified(s):
    """R2 下划线格式 → 统一格式。 'H_S40A' → 'HS40A', 'L_S42H' → 'LS42H'"""
    m = re.match(r"^([HL])_([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析 R2 下划线格式: {s}")
    chain, orig, pos, new = m.groups()
    return f"{chain}{orig}{pos}{new}"


def unified_to_r2_underscore(s):
    """统一格式 → R2 下划线格式。 'HS40A' → 'H_S40A'"""
    m = re.match(r"^([HL])([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析统一格式: {s}")
    chain, orig, pos, new = m.groups()
    return f"{chain}_{orig}{pos}{new}"


# ── 批量转换工具 ─────────────────────────────────────────────────────────────

def convert_mutation_string(s, converter, sep=";"):
    """对分号/逗号分隔的突变串逐个转换。支持 'none(WT)' 透传。"""
    parts = s.split(sep)
    result = []
    for p in parts:
        p = p.strip()
        if not p or p.startswith("none"):
            result.append(p)
            continue
        result.append(converter(p))
    return sep.join(result)


def parse_unified(s):
    """解析统一格式。 'HS40A' → ('H', 'S', 40, 'A')"""
    m = re.match(r"^([HL])([A-Z])(\d+)([A-Z])$", s)
    if not m:
        raise ValueError(f"无法解析统一格式: {s}")
    chain, orig, pos, new = m.groups()
    return chain, orig, int(pos), new


# ── pKa 列名转换 ─────────────────────────────────────────────────────────────

def rename_pka_columns(columns, mapping=None):
    """将 pKa CSV 列名中的 PDB 链号替换为 H/L。
    'B42_propka' → 'L42_propka', 'A48_shift_pkai' → 'H48_shift_pkai'
    """
    pat = re.compile(r"^([AB])(\d+)(_.+)$")
    new_cols = []
    for c in columns:
        m = pat.match(c)
        if m:
            chain, resid, suffix = m.groups()
            new_cols.append(f"{pdb_chain_to_hl(chain, mapping)}{resid}{suffix}")
        else:
            new_cols.append(c)
    return new_cols


# ── CLI 入口 ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="突变命名格式转换工具")
    sub = parser.add_subparsers(dest="cmd")

    p_conv = sub.add_parser("convert", help="转换单个突变标识")
    p_conv.add_argument("mutation", help="突变标识字符串")
    p_conv.add_argument("--from-format", choices=["r1", "r3_arrow", "foldx", "r2_underscore", "unified"],
                        required=True, help="源格式")
    p_conv.add_argument("--to-format", choices=["unified", "foldx", "r1", "r3_arrow", "r2_underscore"],
                        default="unified", help="目标格式")

    args = parser.parse_args()

    if args.cmd == "convert":
        converters_to_unified = {
            "r1": r1_to_unified,
            "r3_arrow": r3_arrow_to_unified,
            "foldx": foldx_to_unified,
            "r2_underscore": r2_underscore_to_unified,
            "unified": lambda x: x,
        }
        converters_from_unified = {
            "unified": lambda x: x,
            "foldx": unified_to_foldx,
            "r1": unified_to_r1,
            "r3_arrow": unified_to_r3_arrow,
            "r2_underscore": unified_to_r2_underscore,
        }
        mid = converters_to_unified[args.from_format](args.mutation)
        result = converters_from_unified[args.to_format](mid)
        print(result)
