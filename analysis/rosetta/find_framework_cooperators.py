#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Framework Cooperation Site Analysis for pH-Switch Design

在 sdab nanobody 的框架区寻找可与 CDR His 突变协同的位点，
通过分析 His 微环境中的框架残基距离、电荷和 pH 响应特性，
输出候选列表用于优化组合设计。

输出：
  - framework_cooperators.csv   排序候选列表
  - cooperators_overview.pml    PyMOL 可视化
  - cooperation_report.md       分析报告

需要 pyrosetta conda 环境。
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

BACKBONE_ATOMS = {"N", "CA", "C", "O", "H", "HA", "HA2", "HA3", "1H", "2H", "3H", "OXT"}

# His 质子化位点（电荷中心）
HIS_CHARGE_ATOMS = ["ND1", "NE2"]

# 带电残基的功能基团原子（电荷中心）
CHARGED_FUNC_ATOMS = {
    "ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"],
    "ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"],
}

# 不可突变的结构关键残基
EXCLUDED_AA = {"CYS", "PRO"}


def load_cdr_regions(config_path):
    """从 config YAML 读取 CDR 定义，返回 {region_name: (start, end)}。"""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    cdr = {}
    for name, v in cfg["tier2"]["cdr_regions"].items():
        cdr[name] = (v["start"], v["end"])
    return cdr


def load_his_targets(config_path):
    """从 config 读取 His 靶位和 ELISA ratio，返回 [(resid, ratio)]。"""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    targets = []
    for entry in cfg["his_bias"]["target_positions"]:
        # format: "A:27"  # R27H CDR1, ratio 1.45
        resid = int(entry.split(":")[1].split()[0])
        # 从注释提取 ratio
        if "ratio" in entry:
            ratio = float(entry.split("ratio")[1].strip())
        else:
            ratio = 1.0
        targets.append((resid, ratio))
    return targets


def classify_residue(resno, cdr_regions):
    """将残基编号分类为 FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4。"""
    for name, (start, end) in cdr_regions.items():
        if start <= resno <= end:
            return name
    # 框架区
    h1_start = cdr_regions["H1"][0]
    h1_end = cdr_regions["H1"][1]
    h2_start = cdr_regions["H2"][0]
    h2_end = cdr_regions["H2"][1]
    h3_start = cdr_regions["H3"][0]
    h3_end = cdr_regions["H3"][1]

    if resno < h1_start:
        return "FR1"
    elif h1_end < resno < h2_start:
        return "FR2"
    elif h2_end < resno < h3_start:
        return "FR3"
    else:
        return "FR4"


class FrameworkCooperationAnalyzer:

    def __init__(self, wt_pdb, binder_chains=("A",), target_chains=("B",),
                 interface_dist=8.0, cooperation_radius=8.0, verbose=False):
        self.wt_pdb = wt_pdb
        self.binder_chains = set(binder_chains)
        self.target_chains = set(target_chains)
        self.interface_dist = interface_dist
        self.coop_radius = cooperation_radius
        self.verbose = verbose

        if not hasattr(pyrosetta, "_initialized"):
            pyrosetta.init("-mute all")
            pyrosetta._initialized = True

        self.pose = pyrosetta.pose_from_pdb(wt_pdb)
        self.pdb_info = self.pose.pdb_info()
        self.chain_map = {}
        self.resname_map = {}
        self.resno_to_idx = {}
        for i in range(1, self.pose.total_residue() + 1):
            c = self.pdb_info.chain(i)
            r = self.pdb_info.number(i)
            self.chain_map[i] = c
            self.resname_map[i] = self.pose.residue(i).name3().strip()
            if c in self.binder_chains:
                self.resno_to_idx[r] = i

    def _res_label(self, idx):
        return f"{self.chain_map[idx]}:{self.resname_map[idx]}{self.pdb_info.number(idx)}"

    # ------------------------------------------------------------------ #
    #  Phase 1: WT Baseline
    # ------------------------------------------------------------------ #
    def get_interface_residues(self):
        binder_idxs = [i for i, c in self.chain_map.items() if c in self.binder_chains]
        target_idxs = [i for i, c in self.chain_map.items() if c in self.target_chains]
        binder_if = set()
        for bi in binder_idxs:
            bxyz = self.pose.residue(bi).nbr_atom_xyz()
            for ti in target_idxs:
                if self.pose.residue(ti).nbr_atom_xyz().distance(bxyz) <= self.interface_dist:
                    binder_if.add(bi)
                    break
        return binder_if

    @staticmethod
    def _get_framework_atoms(res):
        """返回框架残基用于距离测量的原子列表和类型标签。

        - 带电残基（Arg/Lys/Asp/Glu）：返回功能基团原子（电荷中心）
        - 中性残基：返回 CB（Lys NZ 从 CB 延伸约 5Å，CB 指示侧链方向）
        - Gly：返回 CA
        """
        resname = res.name3().strip()
        if resname in CHARGED_FUNC_ATOMS:
            return CHARGED_FUNC_ATOMS[resname], "functional"
        if resname == "GLY":
            return ["CA"], "CA_proxy"
        return ["CB"], "CB_proxy"

    @staticmethod
    def _atom_distance(res1, atoms1, res2, atoms2):
        """计算两组原子间的最短距离。"""
        min_d = 999.0
        best = ("", "")
        for a1 in atoms1:
            try:
                xyz1 = res1.atom(res1.atom_index(a1)).xyz()
            except Exception:
                continue
            for a2 in atoms2:
                try:
                    xyz2 = res2.atom(res2.atom_index(a2)).xyz()
                except Exception:
                    continue
                d = xyz1.distance(xyz2)
                if d < min_d:
                    min_d = d
                    best = (a1, a2)
        return min_d, best

    def classify_framework(self, cdr_regions):
        """Phase 1: 分类框架残基，识别界面，返回基础信息。

        不测量 His 距离（WT 中 His 靶位不是 His，距离无意义）。
        His 距离在 Phase 2 中从变体 PDB 获取。
        """
        interface_set = self.get_interface_residues()
        interface_resnos = {self.pdb_info.number(i) for i in interface_set}

        rows = []
        for resno, idx in sorted(self.resno_to_idx.items()):
            region = classify_residue(resno, cdr_regions)
            if not region.startswith("FR"):
                continue
            resname = self.resname_map[idx]
            charge = "+"  if resname in ("ARG", "LYS") else \
                     "-"  if resname in ("ASP", "GLU") else \
                     "0"

            rows.append({
                "position": resno,
                "wt_aa": resname,
                "region": region,
                "at_interface": resno in interface_resnos,
                "charge": charge,
            })

        if self.verbose:
            n_if = sum(1 for r in rows if r["at_interface"])
            print(f"框架残基: {len(rows)} 个 ({n_if} 个在界面)")
        return rows

    # ------------------------------------------------------------------ #
    #  Phase 2: Per-Variant Analysis
    # ------------------------------------------------------------------ #
    def _simulate_ph5(self, pose):
        from pyrosetta.rosetta.core.pose import replace_pose_residue_copying_existing_coordinates
        pose_copy = Pose()
        pose_copy.assign(pose)
        rsd_set = pose_copy.residue_type_set_for_pose()
        target_name = None
        for name in ["HIS_D", "HIP"]:
            if rsd_set.has_name(name):
                target_name = name
                break
        if not target_name:
            return pose_copy

        target_type = rsd_set.name_map(target_name)
        for i in range(1, pose_copy.total_residue() + 1):
            if pose_copy.residue(i).name3().strip().startswith("HIS"):
                try:
                    replace_pose_residue_copying_existing_coordinates(pose_copy, i, target_type)
                except Exception:
                    pass
        try:
            pose_copy.conformation().detect_disulfides()
            pose_copy.update_residue_neighbors()
        except Exception:
            pass
        return pose_copy

    def _repack_around(self, pose, center_indices, radius=8.0):
        pose_copy = Pose()
        pose_copy.assign(pose)
        # 找邻域
        focus = set()
        for ci in center_indices:
            cxyz = pose.residue(ci).nbr_atom_xyz()
            for i in range(1, pose.total_residue() + 1):
                if pose.residue(i).nbr_atom_xyz().distance(cxyz) <= radius:
                    focus.add(i)
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        task = tf.create_task_and_apply_taskoperations(pose_copy)
        for i in range(1, pose_copy.total_residue() + 1):
            if i not in focus:
                task.nonconst_residue_task(i).prevent_repacking()
        sf = pyrosetta.get_score_function()
        packer = PackRotamersMover(sf)
        packer.task_factory(tf)
        packer.apply(pose_copy)
        return pose_copy

    def analyze_variant(self, variant_pdb, his_resno, cdr_regions):
        """对单个 His 变体分析框架残基与 His 的功能基团距离。

        距离定义（电荷中心到电荷中心）：
        - His 侧: ND1/NE2（质子化位点）
        - 框架侧:
          - Arg/Lys/Asp/Glu: 功能基团原子（NH1/NH2/NE, NZ, OD1/OD2, OE1/OE2）
          - 中性残基: CB（代理：若突变为 Lys，NZ 从 CB 延伸约 5Å）
        """
        pose = pyrosetta.pose_from_pdb(variant_pdb)
        pdb_info = pose.pdb_info()

        # 找到突变 His 的 Rosetta index
        his_idx = None
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.chain(i) in self.binder_chains and pdb_info.number(i) == his_resno:
                his_idx = i
                break
        if his_idx is None:
            return {}

        # pH 7.4: repack His 周围
        pose_ph7 = self._repack_around(pose, [his_idx])

        # pH 5.0: His→HIS_D + repack
        pose_ph5 = self._simulate_ph5(pose_ph7)
        pose_ph5 = self._repack_around(pose_ph5, [his_idx])

        # 框架残基与 His 的功能基团距离
        results = {}
        for resno, idx_wt in self.resno_to_idx.items():
            region = classify_residue(resno, cdr_regions)
            if not region.startswith("FR"):
                continue

            # 在变体 pose 中找对应残基
            idx_var = None
            for i in range(1, pose.total_residue() + 1):
                if pdb_info.chain(i) in self.binder_chains and pdb_info.number(i) == resno:
                    idx_var = i
                    break
            if idx_var is None:
                continue

            fr_atoms, fr_type = self._get_framework_atoms(pose_ph7.residue(idx_var))

            d7, pair7 = self._atom_distance(
                pose_ph7.residue(idx_var), fr_atoms,
                pose_ph7.residue(his_idx), HIS_CHARGE_ATOMS)
            d5, pair5 = self._atom_distance(
                pose_ph5.residue(idx_var), fr_atoms,
                pose_ph5.residue(his_idx), HIS_CHARGE_ATOMS)

            results[resno] = {
                "dist_pH7": round(d7, 2),
                "dist_pH5": round(d5, 2),
                "atoms_pH7": pair7,
                "atoms_pH5": pair5,
                "fr_atom_type": fr_type,
            }

        return results

    # ------------------------------------------------------------------ #
    #  Phase 3: Scoring
    # ------------------------------------------------------------------ #
    def score_cooperators(self, baseline_rows, variant_results, his_targets):
        """对每个框架残基计算协同评分。

        评分完全基于 Phase 2 变体 PDB 的功能基团距离。
        不使用 Phase 1 WT 距离（WT 中 His 靶位不是 His，距离无意义）。
        """
        max_ratio = max(t[1] for t in his_targets)
        his_ratio = {t[0]: t[1] / max_ratio for t in his_targets}

        scored = []
        for row in baseline_rows:
            pos = row["position"]
            if row["at_interface"]:
                continue
            if row["wt_aa"] in EXCLUDED_AA:
                continue

            total_score = 0.0
            cooperating_his = []

            for h_resno, norm_ratio in his_ratio.items():
                vdata = variant_results.get(h_resno, {}).get(pos)
                if not vdata:
                    continue
                dist = vdata["dist_pH7"]
                fr_type = vdata.get("fr_atom_type", "unknown")

                # 距离阈值根据测量类型调整：
                # functional: 电荷中心距离，8Å 是合理的静电作用范围
                # CB_proxy: CB 到 His ND1/NE2 的距离；若突变为 Lys，
                #           NZ 从 CB 延伸约 5Å，所以有效范围 ≈ dist - 5Å
                #           用 12Å 作为 CB 阈值（≈ 7Å functional + 5Å 延伸）
                if fr_type == "functional":
                    cutoff = self.coop_radius
                else:
                    cutoff = self.coop_radius + 4.0  # CB proxy 放宽

                if dist > cutoff:
                    continue

                # 对 CB_proxy 估算功能基团距离 ≈ dist - 3Å (保守估计)
                effective_dist = dist if fr_type == "functional" else max(dist - 3.0, 2.0)

                dist_score = 1.0 if effective_dist <= 6.0 else 0.5
                charge = row["charge"]
                if charge == "+":
                    charge_factor = 1.0
                elif charge == "-":
                    charge_factor = -0.5
                else:
                    charge_factor = 0.5  # 中性：可突变

                contrib = dist_score * norm_ratio * charge_factor
                total_score += contrib
                cooperating_his.append((h_resno, round(dist, 1), fr_type))

            if total_score <= 0 and not cooperating_his:
                continue

            # 最近 His 和距离（从变体数据）
            if cooperating_his:
                nearest = min(cooperating_his, key=lambda x: x[1])
                nearest_his = nearest[0]
                min_dist = nearest[1]
            else:
                nearest_his = None
                min_dist = 999

            # 建议突变
            if row["charge"] == "+":
                proposed = "keep (already +)"
            elif row["charge"] == "-":
                proposed = "→ Lys/Arg (flip charge)"
            elif row["wt_aa"] == "GLY":
                proposed = "→ Lys (caution: Gly flexibility loss)"
            else:
                proposed = "→ Lys (add + charge)"

            scored.append({
                "position": pos,
                "wt_aa": row["wt_aa"],
                "region": row["region"],
                "charge": row["charge"],
                "nearest_his": nearest_his,
                "min_distance": min_dist,
                "n_his_within_cutoff": len(cooperating_his),
                "cooperation_score": round(total_score, 3),
                "proposed_mutation": proposed,
                "cooperating_his": cooperating_his,
                "at_interface": False,
            })

        scored.sort(key=lambda x: x["cooperation_score"], reverse=True)
        return scored

    # ------------------------------------------------------------------ #
    #  Phase 4: Output
    # ------------------------------------------------------------------ #
    def generate_csv(self, scored, output_dir):
        rows = []
        for s in scored:
            his_str = "; ".join(f"{h}({d}Å,{t})" for h, d, t in s["cooperating_his"])
            rows.append({
                "position": s["position"],
                "wt_aa": s["wt_aa"],
                "region": s["region"],
                "current_charge": s["charge"],
                "nearest_his": s["nearest_his"],
                "min_distance": s["min_distance"],
                "n_his_within_cutoff": s["n_his_within_cutoff"],
                "cooperation_score": s["cooperation_score"],
                "proposed_mutation": s["proposed_mutation"],
                "cooperating_his_detail": his_str,
            })
        df = pd.DataFrame(rows)
        path = os.path.join(output_dir, "framework_cooperators.csv")
        df.to_csv(path, index=False)
        if self.verbose:
            print(f"  CSV: {path} ({len(rows)} candidates)")

    def generate_pymol(self, scored, his_targets, interface_set, output_dir):
        lines = []
        pdb_basename = os.path.basename(self.wt_pdb)
        lines.append(f"# Framework Cooperation Analysis: {pdb_basename}")
        lines.append(f"load {pdb_basename}, sdab")
        lines.append("")
        lines.append("hide all")
        lines.append("show cartoon, sdab")
        lines.append("color lightblue, chain A")
        lines.append("color gray80, chain B")
        lines.append("")

        # His targets
        his_sels = [f"(chain A and resi {t[0]})" for t in his_targets]
        lines.append(f"select his_targets, {' or '.join(his_sels)}")
        lines.append("show sticks, his_targets")
        lines.append("color yellow, his_targets and elem C")
        lines.append("")

        # Interface surface
        if_sels = [f"(chain A and resi {self.pdb_info.number(i)})" for i in sorted(interface_set)]
        if if_sels:
            lines.append(f"select interface, {' or '.join(if_sels)}")
            lines.append("show surface, interface")
            lines.append("set transparency, 0.7, interface")
            lines.append("")

        # Top candidates
        top = scored[:15]
        cand_sels = []
        dash_count = 0
        for s in top:
            pos = s["position"]
            cand_sels.append(f"(chain A and resi {pos})")
            # dashes to cooperating His
            for h_resno, dist, _ftype in s["cooperating_his"][:3]:
                dash_count += 1
                lines.append(f"distance coop_{dash_count}, "
                             f"/sdab//A/{pos}/CA, /sdab//A/{h_resno}/CA")

        if cand_sels:
            lines.append(f"select candidates, {' or '.join(cand_sels)}")
            lines.append("show sticks, candidates")
            lines.append("color red, candidates and elem C")
            lines.append("")

        if dash_count:
            lines.append("color magenta, coop_*")
            lines.append("hide labels, coop_*")
            lines.append("set dash_gap, 0.4, coop_*")
            lines.append("group cooperation_dashes, coop_*")
            lines.append("")

        lines.append("set dash_radius, 0.06")
        lines.append("bg_color white")
        lines.append("zoom chain A")
        lines.append("")
        lines.append("# Red sticks = top framework candidates")
        lines.append("# Yellow sticks = His target positions")
        lines.append("# Magenta dashes = candidate-His proximity")

        path = os.path.join(output_dir, "cooperators_overview.pml")
        with open(path, "w") as f:
            f.write("\n".join(lines))
        if self.verbose:
            print(f"  PML: {path}")

    def generate_report(self, scored, baseline_rows, his_targets, variant_results,
                        interface_set, cdr_regions, output_dir):
        lines = []
        pdb_name = os.path.basename(self.wt_pdb).replace(".pdb", "")
        lines.append(f"# Framework Cooperation Analysis: {pdb_name}")
        lines.append("")
        lines.append("## 方法说明")
        lines.append("")
        lines.append("对 sdab WT 结构和 14 个单 His 突变体进行 PyRosetta 分析，"
                     "寻找框架区中与 CDR His 突变空间邻近（≤8Å）的残基，"
                     "评估其通过静电协同增强 pH-switch 效果的潜力。")
        lines.append("")
        lines.append("**注意：** 协同评分基于空间距离和电荷类型的启发式评估，"
                     "不能替代 MD 模拟或实验验证。评分高的候选表示**值得进一步测试**，"
                     "而非已证明有效。")
        lines.append("")

        # Structure overview
        n_binder = sum(1 for c in self.chain_map.values() if c in self.binder_chains)
        n_fr = len(baseline_rows)
        n_if = sum(1 for r in baseline_rows if r["at_interface"])
        lines.append("## 1. 结构概览")
        lines.append("")
        lines.append(f"- sdab (chain A): {n_binder} 残基")
        lines.append(f"- 框架残基: {n_fr} 个（{n_if} 个在界面，已排除）")
        lines.append(f"- His 靶位: {len(his_targets)} 个")
        lines.append(f"- 协同半径: {self.coop_radius} Å")
        lines.append("")

        # CDR regions
        lines.append("**区域划分：**")
        lines.append("")
        lines.append("| Region | Range | His targets |")
        lines.append("|--------|-------|-------------|")
        his_resnos = {t[0] for t in his_targets}
        for region_name in ["FR1", "H1", "FR2", "H2", "FR3", "H3", "FR4"]:
            if region_name in cdr_regions:
                s, e = cdr_regions[region_name]
                display = f"CDR-{region_name}"
            else:
                # 计算 FR 范围
                if region_name == "FR1":
                    s, e = 1, cdr_regions["H1"][0] - 1
                elif region_name == "FR2":
                    s, e = cdr_regions["H1"][1] + 1, cdr_regions["H2"][0] - 1
                elif region_name == "FR3":
                    s, e = cdr_regions["H2"][1] + 1, cdr_regions["H3"][0] - 1
                else:
                    s, e = cdr_regions["H3"][1] + 1, 122
                display = region_name
            his_in = [str(h) for h in sorted(his_resnos) if s <= h <= e]
            lines.append(f"| {display} | {s}-{e} | {', '.join(his_in) or '-'} |")
        lines.append("")

        # Top candidates
        lines.append("## 2. Top 候选位点")
        lines.append("")
        if not scored:
            lines.append("*未找到候选位点。*")
        else:
            lines.append("| Rank | Position | WT AA | Region | Charge | Score | Nearest His | Min Dist | Proposed | Cooperating His |")
            lines.append("|:----:|----------|-------|--------|:------:|:-----:|-------------|:--------:|----------|----------------|")
            for i, s in enumerate(scored[:20], 1):
                his_str = ", ".join(f"{h}({d}Å)" for h, d, _t in s["cooperating_his"][:3])
                lines.append(f"| {i} | **{s['position']}** | {s['wt_aa']} | {s['region']} | "
                             f"{s['charge']} | {s['cooperation_score']:.2f} | "
                             f"{s['nearest_his']} | {s['min_distance']}Å | "
                             f"{s['proposed_mutation']} | {his_str} |")
            lines.append("")

        # Per-candidate analysis for top 10
        lines.append("## 3. Top 10 详细分析")
        lines.append("")
        for s in scored[:10]:
            pos = s["position"]
            lines.append(f"### A:{s['wt_aa']}{pos} ({s['region']})")
            lines.append("")
            lines.append(f"- 当前电荷: {s['charge']}")
            lines.append(f"- 协同评分: {s['cooperation_score']:.3f}")
            lines.append(f"- 建议: {s['proposed_mutation']}")
            lines.append("")

            if s["cooperating_his"]:
                lines.append("| His 靶位 | Func Dist pH7 | Func Dist pH5 | Atom type | pH effect |")
                lines.append("|---------|:---:|:---:|-----------|-----------|")
                for h_resno, dist_val, fr_type in s["cooperating_his"]:
                    vdata = variant_results.get(h_resno, {}).get(pos, {})
                    d7 = vdata.get("dist_pH7", "-")
                    d5 = vdata.get("dist_pH5", "-")
                    if s["charge"] == "+":
                        effect = "His(+) vs (+) → repulsion at pH5"
                    elif s["charge"] == "-":
                        effect = "His(+) vs (-) → attraction at pH5"
                    else:
                        effect = "neutral → mutate to + for repulsion"
                    lines.append(f"| {h_resno} | {d7}Å | {d5}Å | {fr_type} | {effect} |")
                lines.append("")

        # Design recommendations
        lines.append("## 4. 设计建议")
        lines.append("")
        pos_charged = [s for s in scored[:10] if s["charge"] == "+"]
        neutral = [s for s in scored[:10] if s["charge"] == "0"]
        neg_charged = [s for s in scored[:10] if s["charge"] == "-"]

        if pos_charged:
            lines.append("### 已有正电荷（保持或验证）")
            lines.append("")
            for s in pos_charged:
                lines.append(f"- **A:{s['wt_aa']}{s['position']}** ({s['region']}): "
                             f"已带正电，与 {len(s['cooperating_his'])} 个 His 靶位邻近。"
                             f"若在组合设计中此位置被 MPNN 替换为中性残基，可能削弱 pH-switch。"
                             f"**建议在 MPNN 设计中固定此位置。**")
            lines.append("")

        if neutral:
            lines.append("### 中性残基（候选突变为 Lys）")
            lines.append("")
            for s in neutral:
                lines.append(f"- **A:{s['wt_aa']}{s['position']}** ({s['region']}): "
                             f"当前中性，距最近 His 靶位 {s['min_distance']}Å。"
                             f"突变为 Lys 可在 pH5 下与质子化 His 产生正-正排斥。"
                             f"**候选框架突变。**")
            lines.append("")

        if neg_charged:
            lines.append("### 负电荷（需关注）")
            lines.append("")
            for s in neg_charged:
                lines.append(f"- **A:{s['wt_aa']}{s['position']}** ({s['region']}): "
                             f"带负电，在 pH5 下会与 His(+) 吸引，可能**稳定** His 而非促进释放。"
                             f"考虑是否突变为正电或中性。")
            lines.append("")

        path = os.path.join(output_dir, "cooperation_report.md")
        with open(path, "w") as f:
            f.write("\n".join(lines))
        if self.verbose:
            print(f"  Report: {path}")

    # ------------------------------------------------------------------ #
    #  Main
    # ------------------------------------------------------------------ #
    def run(self, his_targets, cdr_regions, variant_pdbs, output_dir):
        os.makedirs(output_dir, exist_ok=True)

        # Phase 1
        if self.verbose:
            print("=== Phase 1: WT Baseline ===")
        baseline = self.classify_framework(cdr_regions)
        interface_set = self.get_interface_residues()

        # Phase 2
        if self.verbose:
            print(f"=== Phase 2: Analyzing {len(variant_pdbs)} variants ===")
        variant_results = {}
        for h_resno, pdb_path in sorted(variant_pdbs.items()):
            if self.verbose:
                print(f"  Variant: His{h_resno} ({os.path.basename(pdb_path)})")
            variant_results[h_resno] = self.analyze_variant(pdb_path, h_resno, cdr_regions)

        # Phase 3
        if self.verbose:
            print("=== Phase 3: Scoring ===")
        scored = self.score_cooperators(baseline, variant_results, his_targets)
        if self.verbose:
            print(f"  {len(scored)} candidates scored")

        # Phase 4
        if self.verbose:
            print("=== Phase 4: Output ===")
        self.generate_csv(scored, output_dir)
        self.generate_pymol(scored, his_targets, interface_set, output_dir)
        self.generate_report(scored, baseline, his_targets, variant_results,
                             interface_set, cdr_regions, output_dir)

        print(f"输出目录: {output_dir}")
        return scored


# ====================================================================== #
#  CLI
# ====================================================================== #
def main():
    parser = argparse.ArgumentParser(
        description="Framework Cooperation Site Analysis for pH-Switch Design")
    parser.add_argument("--wt-pdb", required=True, help="WT PDB file")
    parser.add_argument("--variant-dir", required=True, help="Directory with His variant PDBs")
    parser.add_argument("--variants", nargs="+", required=True,
                        help="Variant names (e.g. A_N52H A_V105H ...)")
    parser.add_argument("--config", required=True, help="Pipeline config YAML")
    parser.add_argument("-o", "--output-dir", default=None)
    parser.add_argument("--binder", nargs="+", default=["A"])
    parser.add_argument("--target", nargs="+", default=["B"])
    parser.add_argument("--interface-dist", type=float, default=8.0)
    parser.add_argument("--coop-radius", type=float, default=8.0)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(os.path.dirname(args.wt_pdb), "framework_cooperation")

    # Load config
    cdr_regions = load_cdr_regions(args.config)
    his_targets = load_his_targets(args.config)

    # Map variant names to PDB files and His residue numbers
    variant_pdbs = {}
    for vname in args.variants:
        pdb_path = os.path.join(args.variant_dir, f"{vname}.pdb")
        if not os.path.exists(pdb_path):
            print(f"WARNING: {pdb_path} not found, skipping")
            continue
        # Parse resno from name: A_N52H → 52
        import re
        m = re.search(r'[A-Z]_[A-Z]?(\d+)H', vname)
        if m:
            resno = int(m.group(1))
            variant_pdbs[resno] = pdb_path

    analyzer = FrameworkCooperationAnalyzer(
        wt_pdb=args.wt_pdb,
        binder_chains=args.binder,
        target_chains=args.target,
        interface_dist=args.interface_dist,
        cooperation_radius=args.coop_radius,
        verbose=args.verbose,
    )
    analyzer.run(his_targets, cdr_regions, variant_pdbs, args.output_dir)


if __name__ == "__main__":
    main()
