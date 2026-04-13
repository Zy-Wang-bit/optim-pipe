#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pH-Dependent Interface Interaction Analysis

对抗体-抗原复合物做双 pH 状态（pH 7.4 vs pH 5.0）的静态分析，
检测氢键和盐桥在 His 质子化前后的变化，生成：
  - interactions_all.csv    全部相互作用清单
  - visualize_interactions.pml  PyMOL 可视化脚本
  - interaction_report.md   结构化分析报告

需要 pyrosetta conda 环境。
"""

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

# 盐桥相关原子定义
POSITIVE_ATOMS = {
    "ARG": ["NH1", "NH2", "NE"],
    "LYS": ["NZ"],
}
POSITIVE_ATOMS_HIP = {"ND1", "NE2"}  # HIP (protonated His) 的带正电原子

NEGATIVE_ATOMS = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}

BACKBONE_ATOMS = {"N", "CA", "C", "O", "H", "HA", "HA2", "HA3", "1H", "2H", "3H", "OXT"}


class InterfaceInteractionAnalyzer:

    def __init__(self, pdb_path, binder_chains=("A", "B"), target_chains=("C"),
                 interface_dist=8.0, salt_bridge_dist=4.0, verbose=False):
        self.pdb_path = pdb_path
        self.pdb_name = os.path.basename(pdb_path).replace(".pdb", "")
        self.binder_chains = set(binder_chains)
        self.target_chains = set(target_chains)
        self.interface_dist = interface_dist
        self.salt_bridge_dist = salt_bridge_dist
        self.verbose = verbose

        if not hasattr(pyrosetta, "_initialized"):
            pyrosetta.init("-mute all")
            pyrosetta._initialized = True

        self.pose = pyrosetta.pose_from_pdb(pdb_path)
        self.pdb_info = self.pose.pdb_info()
        self.chain_map = self._build_chain_map()
        self.resname_map = self._build_resname_map()

    # ------------------------------------------------------------------ #
    #  基础设施
    # ------------------------------------------------------------------ #
    def _build_chain_map(self):
        chain_map = {}
        for i in range(1, self.pose.total_residue() + 1):
            chain_map[i] = self.pdb_info.chain(i)
        return chain_map

    def _build_resname_map(self):
        resname_map = {}
        for i in range(1, self.pose.total_residue() + 1):
            resname_map[i] = self.pose.residue(i).name3().strip()
        return resname_map

    def _resno(self, idx):
        """Rosetta index → PDB residue number."""
        return self.pdb_info.number(idx)

    def _res_label(self, idx):
        """Human-readable label: 'A:Asp33'."""
        chain = self.chain_map[idx]
        resname = self.resname_map[idx]
        resno = self._resno(idx)
        return f"{chain}:{resname}{resno}"

    def _is_backbone_atom(self, pose, rsd_idx, atm_idx):
        try:
            atom_name = pose.residue(rsd_idx).atom_name(atm_idx).strip()
            return atom_name in BACKBONE_ATOMS
        except Exception:
            return False

    def _get_atom_name(self, pose, rsd_idx, atm_idx):
        try:
            return pose.residue(rsd_idx).atom_name(atm_idx).strip()
        except Exception:
            return "UNK"

    def _get_neighbors(self, pose, center_indices, radius_A=8.0):
        neighbor_set = set()
        for cidx in center_indices:
            center_xyz = pose.residue(cidx).nbr_atom_xyz()
            for i in range(1, pose.total_residue() + 1):
                if pose.residue(i).nbr_atom_xyz().distance(center_xyz) <= radius_A:
                    neighbor_set.add(i)
        return sorted(neighbor_set)

    def _get_histidine_residues(self, pose):
        his = []
        for i in range(1, pose.total_residue() + 1):
            if pose.residue(i).name3().strip().startswith("HIS"):
                his.append(i)
        return his

    # ------------------------------------------------------------------ #
    #  界面残基识别
    # ------------------------------------------------------------------ #
    def _get_interface_residues(self, pose):
        """返回 (binder_interface_set, target_interface_set)."""
        binder_idxs = [i for i in range(1, pose.total_residue() + 1)
                       if self.chain_map[i] in self.binder_chains]
        target_idxs = [i for i in range(1, pose.total_residue() + 1)
                       if self.chain_map[i] in self.target_chains]

        binder_if, target_if = set(), set()
        for bi in binder_idxs:
            bxyz = pose.residue(bi).nbr_atom_xyz()
            for ti in target_idxs:
                if pose.residue(ti).nbr_atom_xyz().distance(bxyz) <= self.interface_dist:
                    binder_if.add(bi)
                    target_if.add(ti)
        return binder_if, target_if

    # ------------------------------------------------------------------ #
    #  Repack & pH 模拟
    # ------------------------------------------------------------------ #
    def _repack_focused(self, pose, center_indices):
        focus = self._get_neighbors(pose, center_indices, radius_A=8.0)
        pose_copy = Pose()
        pose_copy.assign(pose)

        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        task = tf.create_task_and_apply_taskoperations(pose_copy)
        focus_set = set(focus)
        for i in range(1, pose_copy.total_residue() + 1):
            if i not in focus_set:
                task.nonconst_residue_task(i).prevent_repacking()

        sf = pyrosetta.get_score_function()
        packer = PackRotamersMover(sf)
        packer.task_factory(tf)
        packer.apply(pose_copy)
        return pose_copy

    def _simulate_ph5_all_his(self, pose):
        from pyrosetta.rosetta.core.pose import replace_pose_residue_copying_existing_coordinates

        pose_copy = Pose()
        pose_copy.assign(pose)
        his_residues = self._get_histidine_residues(pose_copy)

        if not his_residues:
            return pose_copy

        # 获取 HIP/HIS_D 残基类型（pose 级别的 residue type set）
        rsd_set = pose_copy.residue_type_set_for_pose()
        target_name = None
        if rsd_set.has_name("HIP"):
            target_name = "HIP"
        elif rsd_set.has_name("HIS_D"):
            target_name = "HIS_D"
        else:
            print("  WARNING: 数据库中未找到 HIP 或 HIS_D，跳过质子化")
            return pose_copy

        target_type = rsd_set.name_map(target_name)
        changed = 0
        for his_idx in his_residues:
            try:
                replace_pose_residue_copying_existing_coordinates(pose_copy, his_idx, target_type)
                changed += 1
            except Exception as e:
                if self.verbose:
                    print(f"  WARNING: His idx={his_idx} 替换失败: {e}")

        try:
            pose_copy.conformation().detect_disulfides()
            pose_copy.update_residue_neighbors()
        except Exception:
            pass

        if self.verbose:
            print(f"  pH5 模拟: {changed}/{len(his_residues)} His 替换为 {target_name}")
        return pose_copy

    # ------------------------------------------------------------------ #
    #  氢键检测
    # ------------------------------------------------------------------ #
    def _analyze_hbonds(self, pose):
        """返回氢键列表 [{donor_idx, acceptor_idx, donor_atom, acceptor_atom, energy, ...}]."""
        pose.update_residue_neighbors()
        sf = pyrosetta.get_score_function()
        sf(pose)
        hbond_set = pose.get_hbonds(False, False, False, False, False)

        hbonds = []
        for i in range(1, hbond_set.nhbonds() + 1):
            hb = hbond_set.hbond(i)
            don_rsd = hb.don_res()
            acc_rsd = hb.acc_res()
            don_chain = self.chain_map.get(don_rsd)
            acc_chain = self.chain_map.get(acc_rsd)
            if don_chain is None or acc_chain is None:
                continue

            don_hatm = hb.don_hatm()
            acc_atm = hb.acc_atm()
            don_atom_name = self._get_atom_name(pose, don_rsd, don_hatm)
            acc_atom_name = self._get_atom_name(pose, acc_rsd, acc_atm)

            # 计算距离 (donor heavy atom → acceptor atom)
            # donor的H原子对应的heavy atom大约是 don_hatm - 1
            # 但更可靠的是直接用 don_res 的 N/O → acc atom 距离
            try:
                don_xyz = pose.residue(don_rsd).atom(don_hatm).xyz()
                acc_xyz = pose.residue(acc_rsd).atom(acc_atm).xyz()
                dist = don_xyz.distance(acc_xyz)
            except Exception:
                dist = 0.0

            don_bb = self._is_backbone_atom(pose, don_rsd, don_hatm)
            acc_bb = self._is_backbone_atom(pose, acc_rsd, acc_atm)
            if don_bb and acc_bb:
                subtype = "backbone-backbone"
            elif don_bb or acc_bb:
                subtype = "backbone-sidechain"
            else:
                subtype = "sidechain-sidechain"

            don_resname = self.resname_map.get(don_rsd, "")
            acc_resname = self.resname_map.get(acc_rsd, "")
            involves_his = don_resname.startswith("HIS") or acc_resname.startswith("HIS")

            hbonds.append({
                "donor_idx": don_rsd,
                "acceptor_idx": acc_rsd,
                "donor_atom": don_atom_name,
                "acceptor_atom": acc_atom_name,
                "distance": round(dist, 2),
                "energy": round(hb.energy(), 3),
                "subtype": subtype,
                "involves_his": involves_his,
                "type": "hbond",
            })
        return hbonds

    # ------------------------------------------------------------------ #
    #  盐桥检测
    # ------------------------------------------------------------------ #
    def _detect_salt_bridges(self, pose, include_hip=False):
        """检测盐桥。include_hip=True 时把 HIS_D/HIP 视为正电残基。"""
        pos_residues = []  # (idx, atom_names)
        neg_residues = []

        for i in range(1, pose.total_residue() + 1):
            resname = pose.residue(i).name3().strip()
            if resname in POSITIVE_ATOMS:
                pos_residues.append((i, POSITIVE_ATOMS[resname]))
            if include_hip and resname.startswith("HIS"):
                pos_residues.append((i, list(POSITIVE_ATOMS_HIP)))
            if resname in NEGATIVE_ATOMS:
                neg_residues.append((i, NEGATIVE_ATOMS[resname]))

        bridges = []
        for pi, p_atoms in pos_residues:
            res_p = pose.residue(pi)
            for ni, n_atoms in neg_residues:
                if pi == ni:
                    continue
                res_n = pose.residue(ni)
                min_dist = 999.0
                best_pair = (None, None)
                for pa in p_atoms:
                    try:
                        pa_idx = res_p.atom_index(pa)
                    except Exception:
                        continue
                    pa_xyz = res_p.atom(pa_idx).xyz()
                    for na in n_atoms:
                        try:
                            na_idx = res_n.atom_index(na)
                        except Exception:
                            continue
                        d = pa_xyz.distance(res_n.atom(na_idx).xyz())
                        if d < min_dist:
                            min_dist = d
                            best_pair = (pa, na)
                if min_dist <= self.salt_bridge_dist:
                    p_resname = self.resname_map.get(pi, "")
                    n_resname = self.resname_map.get(ni, "")
                    involves_his = p_resname.startswith("HIS") or n_resname.startswith("HIS")
                    bridges.append({
                        "donor_idx": pi,
                        "acceptor_idx": ni,
                        "donor_atom": best_pair[0],
                        "acceptor_atom": best_pair[1],
                        "distance": round(min_dist, 2),
                        "energy": None,
                        "subtype": "sidechain-sidechain",
                        "involves_his": involves_his,
                        "type": "salt_bridge",
                    })
        return bridges

    # ------------------------------------------------------------------ #
    #  His 微环境分析
    # ------------------------------------------------------------------ #
    def _analyze_his_environment(self, pose, his_residues, binder_if, radius=6.0):
        """分析每个 His 侧链 6Å 内的带电残基，返回 {his_idx: [contacts]}。

        距离定义：His 电荷中心（ND1/NE2）到 partner 功能基团原子的最短距离。
        """
        his_sc_atoms = ["ND1", "NE2"]
        charged_types = {
            "ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"],
            "ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"],
        }

        env = {}
        for hidx in his_residues:
            contacts = []
            res_h = pose.residue(hidx)
            h_chain = self.chain_map[hidx]
            h_resno = self._resno(hidx)

            for j in range(1, pose.total_residue() + 1):
                resname_j = self.resname_map.get(j, "")
                if resname_j not in charged_types:
                    continue
                res_j = pose.residue(j)
                j_chain = self.chain_map[j]
                j_resno = self._resno(j)

                min_d = 999.0
                best_pair = ("", "")
                for ha in his_sc_atoms:
                    try:
                        ha_idx = res_h.atom_index(ha)
                        ha_xyz = res_h.atom(ha_idx).xyz()
                    except Exception:
                        continue
                    for ta in charged_types[resname_j]:
                        try:
                            ta_idx = res_j.atom_index(ta)
                            d = ha_xyz.distance(res_j.atom(ta_idx).xyz())
                            if d < min_d:
                                min_d = d
                                best_pair = (ha, ta)
                        except Exception:
                            continue

                if min_d <= radius:
                    charge = "+" if resname_j in ("ARG", "LYS") else "-"
                    cross_interface = ((h_chain in self.binder_chains and j_chain in self.target_chains) or
                                       (h_chain in self.target_chains and j_chain in self.binder_chains))
                    contacts.append({
                        "partner": f"{j_chain}:{resname_j}{j_resno}",
                        "partner_charge": charge,
                        "his_atom": best_pair[0],
                        "partner_atom": best_pair[1],
                        "distance": round(min_d, 2),
                        "cross_interface": cross_interface,
                        "at_interface": hidx in binder_if,
                    })
            env[hidx] = sorted(contacts, key=lambda c: c["distance"])
        return env

    # ------------------------------------------------------------------ #
    #  交叉匹配
    # ------------------------------------------------------------------ #
    def _interaction_key(self, ix):
        """生成匹配 key: (chain1:resno1:atom, chain2:resno2:atom, type)."""
        d_chain = self.chain_map[ix["donor_idx"]]
        a_chain = self.chain_map[ix["acceptor_idx"]]
        d_resno = self._resno(ix["donor_idx"])
        a_resno = self._resno(ix["acceptor_idx"])
        return (d_chain, d_resno, ix["donor_atom"],
                a_chain, a_resno, ix["acceptor_atom"],
                ix["type"])

    def _classify_interaction(self, ix, binder_if, target_if):
        """判断 interface/internal_antibody。"""
        d_chain = self.chain_map[ix["donor_idx"]]
        a_chain = self.chain_map[ix["acceptor_idx"]]
        d_in_binder = d_chain in self.binder_chains
        a_in_binder = a_chain in self.binder_chains
        d_in_target = d_chain in self.target_chains
        a_in_target = a_chain in self.target_chains

        if (d_in_binder and a_in_target) or (d_in_target and a_in_binder):
            return "interface"
        if d_in_binder and a_in_binder:
            return "internal_antibody"
        return "other"

    def _merge_interactions(self, ph7_list, ph5_list, binder_if, target_if):
        """合并两个 pH 状态的相互作用列表。"""
        ph7_map = {}
        for ix in ph7_list:
            key = self._interaction_key(ix)
            ph7_map[key] = ix

        ph5_map = {}
        for ix in ph5_list:
            key = self._interaction_key(ix)
            ph5_map[key] = ix

        all_keys = set(ph7_map.keys()) | set(ph5_map.keys())
        merged = []
        for key in sorted(all_keys):
            in7 = key in ph7_map
            in5 = key in ph5_map
            ref = ph7_map.get(key) or ph5_map[key]

            if in7 and in5:
                status = "stable"
            elif in7 and not in5:
                status = "lost"
            else:
                status = "gained"

            category = self._classify_interaction(ref, binder_if, target_if)

            merged.append({
                "donor": self._res_label(ref["donor_idx"]),
                "donor_chain": self.chain_map[ref["donor_idx"]],
                "donor_resno": self._resno(ref["donor_idx"]),
                "donor_atom": ref["donor_atom"],
                "acceptor": self._res_label(ref["acceptor_idx"]),
                "acceptor_chain": self.chain_map[ref["acceptor_idx"]],
                "acceptor_resno": self._resno(ref["acceptor_idx"]),
                "acceptor_atom": ref["acceptor_atom"],
                "type": ref["type"],
                "subtype": ref["subtype"],
                "distance_pH7": ph7_map[key]["distance"] if in7 else None,
                "distance_pH5": ph5_map[key]["distance"] if in5 else None,
                "energy_pH7": ph7_map[key]["energy"] if in7 else None,
                "energy_pH5": ph5_map[key]["energy"] if in5 else None,
                "present_pH7": in7,
                "present_pH5": in5,
                "involves_his": ref["involves_his"],
                "category": category,
                "status": status,
            })
        return merged

    # ------------------------------------------------------------------ #
    #  主流程
    # ------------------------------------------------------------------ #
    def run(self, output_dir):
        os.makedirs(output_dir, exist_ok=True)

        # 1. 界面残基
        binder_if, target_if = self._get_interface_residues(self.pose)
        if self.verbose:
            print(f"界面残基: {len(binder_if)} binder, {len(target_if)} target")

        # 2. His 残基
        his_residues = self._get_histidine_residues(self.pose)
        if self.verbose:
            labels = [self._res_label(h) for h in his_residues]
            print(f"His 残基 ({len(his_residues)}): {', '.join(labels)}")

        # 3. pH 7.4 分析
        if self.verbose:
            print("--- pH 7.4 分析 ---")
        repack_shell = self._get_neighbors(self.pose, his_residues, 8.0) if his_residues else None
        pose_ph7 = self._repack_focused(self.pose, his_residues) if his_residues else Pose(self.pose)
        hbonds_ph7 = self._analyze_hbonds(pose_ph7)
        salt_ph7 = self._detect_salt_bridges(pose_ph7, include_hip=False)
        all_ph7 = hbonds_ph7 + salt_ph7
        if self.verbose:
            print(f"  pH7.4: {len(hbonds_ph7)} hbonds, {len(salt_ph7)} salt bridges")

        # 4. pH 5.0 分析
        if self.verbose:
            print("--- pH 5.0 分析 ---")
        pose_ph5 = self._simulate_ph5_all_his(pose_ph7)
        pose_ph5 = self._repack_focused(pose_ph5, his_residues) if his_residues else pose_ph5
        hbonds_ph5 = self._analyze_hbonds(pose_ph5)
        salt_ph5 = self._detect_salt_bridges(pose_ph5, include_hip=True)
        all_ph5 = hbonds_ph5 + salt_ph5
        if self.verbose:
            print(f"  pH5.0: {len(hbonds_ph5)} hbonds, {len(salt_ph5)} salt bridges")

        # 5. 合并
        merged = self._merge_interactions(all_ph7, all_ph5, binder_if, target_if)
        if self.verbose:
            n_stable = sum(1 for m in merged if m["status"] == "stable")
            n_lost = sum(1 for m in merged if m["status"] == "lost")
            n_gained = sum(1 for m in merged if m["status"] == "gained")
            print(f"合并: {len(merged)} 条 (stable={n_stable}, lost={n_lost}, gained={n_gained})")

        # 5b. His 微环境（基于原始结构，不受 repack rotamer 偏移影响）
        his_env_ph7 = self._analyze_his_environment(self.pose, his_residues, binder_if)
        his_env_ph5 = self._analyze_his_environment(pose_ph5, his_residues, binder_if)

        # 6. 输出
        self._generate_csv(merged, output_dir)
        self._generate_pymol_script(merged, binder_if, target_if, his_residues,
                                    his_env_ph7, output_dir)
        self._generate_report(merged, binder_if, target_if, his_residues,
                              his_env_ph7, his_env_ph5, output_dir)

        print(f"输出目录: {output_dir}")
        return merged

    # ------------------------------------------------------------------ #
    #  CSV 输出
    # ------------------------------------------------------------------ #
    def _generate_csv(self, merged, output_dir):
        df = pd.DataFrame(merged)
        path = os.path.join(output_dir, "interactions_all.csv")
        df.to_csv(path, index=False)
        if self.verbose:
            print(f"  CSV: {path}")

    # ------------------------------------------------------------------ #
    #  PyMOL 脚本
    # ------------------------------------------------------------------ #
    def _generate_pymol_script(self, merged, binder_if, target_if, his_residues,
                               his_env_ph7, output_dir):
        lines = []
        lines.append(f"# pH-Dependent Interaction Visualization: {self.pdb_name}")
        lines.append(f"# Generated by analyze_interface_interactions.py")
        lines.append("")

        # Load
        pdb_basename = os.path.basename(self.pdb_path)
        lines.append(f"load {pdb_basename}, complex")
        lines.append("")

        # Display
        lines.append("# Basic display")
        lines.append("hide all")
        lines.append("show cartoon, complex")
        lines.append("color lightblue, chain A")
        lines.append("color lightorange, chain B")
        lines.append("color gray80, chain C")
        lines.append("")

        # Interface residues
        binder_sels = []
        for idx in sorted(binder_if):
            c = self.chain_map[idx]
            r = self._resno(idx)
            binder_sels.append(f"(chain {c} and resi {r})")
        target_sels = []
        for idx in sorted(target_if):
            c = self.chain_map[idx]
            r = self._resno(idx)
            target_sels.append(f"(chain {c} and resi {r})")

        if binder_sels:
            lines.append(f"select interface_ab, {' or '.join(binder_sels)}")
            lines.append("show sticks, interface_ab")
        if target_sels:
            lines.append(f"select interface_ag, {' or '.join(target_sels)}")
            lines.append("show sticks, interface_ag")
        lines.append("")

        # His residues
        his_sels = []
        for idx in his_residues:
            c = self.chain_map[idx]
            r = self._resno(idx)
            his_sels.append(f"(chain {c} and resi {r})")
        if his_sels:
            lines.append(f"select his_residues, {' or '.join(his_sels)}")
            lines.append("color yellow, his_residues and elem C")
            lines.append("show sticks, his_residues")
            lines.append("")

        # Interaction dashes
        counters = {"hbond_stable": 0, "salt_stable": 0, "lost": 0, "gained": 0}

        for ix in merged:
            # 显示: 界面相互作用、His 相关、或 pH 变化的相互作用
            if ix["status"] == "stable":
                if ix["category"] == "other":
                    continue
                if ix["category"] == "internal_antibody" and not ix["involves_his"]:
                    continue

            dc = ix["donor_chain"]
            dr = ix["donor_resno"]
            da = ix["donor_atom"]
            ac = ix["acceptor_chain"]
            ar = ix["acceptor_resno"]
            aa = ix["acceptor_atom"]

            sel1 = f"/complex//{dc}/{dr}/{da}"
            sel2 = f"/complex//{ac}/{ar}/{aa}"

            if ix["status"] == "stable":
                if ix["type"] == "hbond":
                    counters["hbond_stable"] += 1
                    name = f"hbond_stable_{counters['hbond_stable']}"
                else:
                    counters["salt_stable"] += 1
                    name = f"salt_stable_{counters['salt_stable']}"
            elif ix["status"] == "lost":
                counters["lost"] += 1
                name = f"lost_{counters['lost']}"
            else:
                counters["gained"] += 1
                name = f"gained_{counters['gained']}"

            lines.append(f"distance {name}, {sel1}, {sel2}")

        lines.append("")

        # Color dashes
        if counters["hbond_stable"]:
            lines.append("color cyan, hbond_stable_*")
            lines.append("hide labels, hbond_stable_*")
        if counters["salt_stable"]:
            lines.append("color blue, salt_stable_*")
            lines.append("hide labels, salt_stable_*")
        if counters["lost"]:
            lines.append("color red, lost_*")
            lines.append("hide labels, lost_*")
        if counters["gained"]:
            lines.append("color orange, gained_*")
            lines.append("hide labels, gained_*")
        lines.append("")

        # Groups
        if counters["hbond_stable"]:
            lines.append("group hbonds_stable, hbond_stable_*")
        if counters["salt_stable"]:
            lines.append("group salt_bridges_stable, salt_stable_*")
        if counters["lost"]:
            lines.append("group interactions_lost, lost_*")
        if counters["gained"]:
            lines.append("group interactions_gained, gained_*")
        lines.append("")

        # His microenvironment: charged residue proximity dashes
        env_counter = 0
        for hidx in his_residues:
            contacts = his_env_ph7.get(hidx, [])
            h_chain = self.chain_map[hidx]
            h_resno = self._resno(hidx)
            for c in contacts:
                env_counter += 1
                p_parts = c["partner"].split(":")
                p_chain = p_parts[0]
                p_resno = ''.join(filter(str.isdigit, p_parts[1]))
                sel1 = f"/complex//{h_chain}/{h_resno}/{c['his_atom']}"
                sel2 = f"/complex//{p_chain}/{p_resno}/{c['partner_atom']}"
                name = f"his_env_{env_counter}"
                lines.append(f"distance {name}, {sel1}, {sel2}")

        if env_counter:
            lines.append("color magenta, his_env_*")
            lines.append("set dash_width, 1.5, his_env_*")
            lines.append("set dash_gap, 0.4, his_env_*")
            lines.append("hide labels, his_env_*")
            lines.append("group his_environment, his_env_*")
            lines.append("")

        # Display settings
        lines.append("# Display settings")
        lines.append("set dash_gap, 0.2")
        lines.append("set dash_radius, 0.08")
        lines.append("bg_color white")
        lines.append("set ray_opaque_background, 1")
        if binder_sels or target_sels:
            lines.append("zoom interface_ab or interface_ag")
        lines.append("")
        lines.append("# Tip: toggle groups in PyMOL object panel to show/hide interaction types")

        pml_path = os.path.join(output_dir, "visualize_interactions.pml")
        with open(pml_path, "w") as f:
            f.write("\n".join(lines))
        if self.verbose:
            print(f"  PML: {pml_path}")

    # ------------------------------------------------------------------ #
    #  Markdown 报告
    # ------------------------------------------------------------------ #
    def _generate_report(self, merged, binder_if, target_if, his_residues,
                         his_env_ph7, his_env_ph5, output_dir):
        lines = []
        lines.append(f"# pH-Dependent Interaction Analysis: {self.pdb_name}")
        lines.append("")

        # 1. Structure overview
        n_total = self.pose.total_residue()
        chain_counts = defaultdict(int)
        for i in range(1, n_total + 1):
            chain_counts[self.chain_map[i]] += 1

        lines.append("## 1. Structure Overview")
        lines.append("")
        for c in sorted(chain_counts):
            role = "Heavy chain" if c == "A" else "Light chain" if c == "B" else "Antigen"
            lines.append(f"- Chain {c} ({role}): {chain_counts[c]} residues")
        lines.append(f"- His residues ({len(his_residues)}): "
                     + ", ".join(self._res_label(h) for h in his_residues))
        lines.append("")

        # 2. Interface residues
        lines.append("## 2. Interface Residues")
        lines.append("")
        lines.append(f"Cutoff: {self.interface_dist} Å")
        lines.append("")
        lines.append(f"- Antibody interface: {len(binder_if)} residues")
        binder_labels = sorted([self._res_label(i) for i in binder_if])
        lines.append(f"  - {', '.join(binder_labels)}")
        lines.append(f"- Antigen interface: {len(target_if)} residues")
        target_labels = sorted([self._res_label(i) for i in target_if])
        lines.append(f"  - {', '.join(target_labels)}")
        lines.append("")

        # Helper: format table
        def _table(rows, category_filter=None, type_filter=None, his_only=False):
            filtered = [r for r in merged
                        if (category_filter is None or r["category"] == category_filter)
                        and (type_filter is None or r["type"] == type_filter)
                        and (not his_only or r["involves_his"])]
            if not filtered:
                lines.append("*None detected.*")
                lines.append("")
                return
            lines.append("| Donor | Acceptor | Type | Dist pH7 | Dist pH5 | Energy pH7 | Status |")
            lines.append("|-------|----------|------|----------|----------|------------|--------|")
            for r in filtered:
                d7 = f"{r['distance_pH7']}" if r["distance_pH7"] is not None else "-"
                d5 = f"{r['distance_pH5']}" if r["distance_pH5"] is not None else "-"
                e7 = f"{r['energy_pH7']}" if r["energy_pH7"] is not None else "-"
                st = r["status"]
                if st == "lost":
                    st = "**LOST**"
                elif st == "gained":
                    st = "**GAINED**"
                lines.append(f"| {r['donor']}:{r['donor_atom']} | {r['acceptor']}:{r['acceptor_atom']} "
                             f"| {r['subtype']} | {d7} | {d5} | {e7} | {st} |")
            lines.append("")

        # 3. Hydrogen bonds
        lines.append("## 3. Hydrogen Bonds")
        lines.append("")
        lines.append("### 3.1 Interface H-bonds")
        lines.append("")
        _table(merged, category_filter="interface", type_filter="hbond")

        lines.append("### 3.2 Internal Antibody H-bonds (His-involved)")
        lines.append("")
        _table(merged, category_filter="internal_antibody", type_filter="hbond", his_only=True)

        # 4. Salt bridges
        lines.append("## 4. Salt Bridges")
        lines.append("")
        lines.append("### 4.1 Interface Salt Bridges")
        lines.append("")
        _table(merged, category_filter="interface", type_filter="salt_bridge")

        lines.append("### 4.2 Internal Salt Bridges (His-involved)")
        lines.append("")
        _table(merged, category_filter="internal_antibody", type_filter="salt_bridge", his_only=True)

        # 5. pH-dependent key findings
        lines.append("## 5. pH-Dependent Key Findings")
        lines.append("")

        lost = [r for r in merged if r["status"] == "lost"]
        gained = [r for r in merged if r["status"] == "gained"]

        lines.append(f"### 5.1 Interactions Lost at pH 5.0 ({len(lost)})")
        lines.append("")
        if lost:
            for r in lost:
                his_tag = " [His-involved]" if r["involves_his"] else ""
                lines.append(f"- **{r['donor']}:{r['donor_atom']} → {r['acceptor']}:{r['acceptor_atom']}** "
                             f"({r['type']}, {r['category']}{his_tag})")
            lines.append("")
            lines.append("His 质子化（pH 5.0）后，N-delta 或 N-epsilon 由氢键受体变为供体，"
                         "导致原有氢键断裂。若 His 附近有同号电荷残基（Arg/Lys），"
                         "质子化后的电荷排斥也会破坏原有接触。")
        else:
            lines.append("*No interactions lost.*")
        lines.append("")

        lines.append(f"### 5.2 Interactions Gained at pH 5.0 ({len(gained)})")
        lines.append("")
        if gained:
            for r in gained:
                his_tag = " [His-involved]" if r["involves_his"] else ""
                lines.append(f"- **{r['donor']}:{r['donor_atom']} → {r['acceptor']}:{r['acceptor_atom']}** "
                             f"({r['type']}, {r['category']}{his_tag})")
            lines.append("")
            lines.append("His 质子化后带 +1 正电，可与 Asp/Glu 的负电基团形成新的盐桥或氢键。"
                         "这些新形成的相互作用可能改变抗体构象或影响界面稳定性。")
        else:
            lines.append("*No interactions gained.*")
        lines.append("")

        # 5.3 Per-His analysis with environment
        lines.append("### 5.3 Per-His Residue Microenvironment")
        lines.append("")
        lines.append("His 侧链 6Å 内的带电残基（潜在 pH 敏感接触）：")
        lines.append("")
        for hidx in his_residues:
            label = self._res_label(hidx)
            chain = self.chain_map[hidx]
            location = "interface" if hidx in binder_if else "internal"
            lines.append(f"#### {label} ({location})")
            lines.append("")

            # H-bonds / salt bridges involving this His
            his_interactions = [r for r in merged
                                if r["involves_his"]
                                and (self._resno(hidx) in (r["donor_resno"], r["acceptor_resno"]))
                                and (chain in (r["donor_chain"], r["acceptor_chain"]))]
            if his_interactions:
                lines.append("**Direct bonds:**")
                for r in his_interactions:
                    partner = r["acceptor"] if r["donor"].startswith(f"{chain}:HIS") else r["donor"]
                    lines.append(f"- {r['type']} with {partner} "
                                 f"— status: **{r['status']}**"
                                 f" (pH7: {r['distance_pH7'] or '-'} Å, pH5: {r['distance_pH5'] or '-'} Å)")
                lines.append("")

            # His environment contacts
            env7 = his_env_ph7.get(hidx, [])
            env5 = his_env_ph5.get(hidx, [])
            if env7 or env5:
                lines.append("**Charged residues within 6Å of sidechain:**")
                lines.append("")
                lines.append("| Partner | Charge | His atom | Partner atom | Dist pH7 | Dist pH5 | Cross-interface | pH effect |")
                lines.append("|---------|--------|----------|-------------|----------|----------|----------------|-----------|")

                # Merge pH7 and pH5 environments by partner
                all_partners = {}
                for c in env7:
                    all_partners[c["partner"]] = {"ph7": c, "ph5": None}
                for c in env5:
                    if c["partner"] in all_partners:
                        all_partners[c["partner"]]["ph5"] = c
                    else:
                        all_partners[c["partner"]] = {"ph7": None, "ph5": c}

                for partner, data in sorted(all_partners.items(), key=lambda x: (x[1]["ph7"] or x[1]["ph5"])["distance"]):
                    c7 = data["ph7"]
                    c5 = data["ph5"]
                    ref = c7 or c5
                    charge = ref["partner_charge"]
                    his_atom = ref["his_atom"]
                    p_atom = ref["partner_atom"]
                    d7 = f"{c7['distance']}" if c7 else "-"
                    d5 = f"{c5['distance']}" if c5 else "-"
                    cross = "YES" if ref["cross_interface"] else "no"

                    # pH effect annotation
                    if charge == "+":
                        effect = "His(+) + (+) → repulsion at pH5"
                    else:
                        effect = "His(+) + (-) → attraction at pH5"

                    lines.append(f"| {partner} | {charge} | {his_atom} | {p_atom} | {d7} | {d5} | {cross} | {effect} |")
                lines.append("")
            else:
                lines.append("*No charged residues within 6Å of sidechain.*")
                lines.append("")

        # 6. Summary
        lines.append("## 6. Summary")
        lines.append("")
        n_if_ph7 = sum(1 for r in merged if r["category"] == "interface" and r["present_pH7"])
        n_if_ph5 = sum(1 for r in merged if r["category"] == "interface" and r["present_pH5"])
        n_int_ph7 = sum(1 for r in merged if r["category"] == "internal_antibody" and r["present_pH7"])
        n_int_ph5 = sum(1 for r in merged if r["category"] == "internal_antibody" and r["present_pH5"])
        lines.append(f"| Metric | pH 7.4 | pH 5.0 | Change |")
        lines.append(f"|--------|--------|--------|--------|")
        lines.append(f"| Interface interactions | {n_if_ph7} | {n_if_ph5} | {n_if_ph5 - n_if_ph7:+d} |")
        lines.append(f"| Internal Ab interactions | {n_int_ph7} | {n_int_ph5} | {n_int_ph5 - n_int_ph7:+d} |")
        lines.append(f"| Interactions lost | - | - | {len(lost)} |")
        lines.append(f"| Interactions gained | - | - | {len(gained)} |")
        lines.append("")

        report_path = os.path.join(output_dir, "interaction_report.md")
        with open(report_path, "w") as f:
            f.write("\n".join(lines))
        if self.verbose:
            print(f"  Report: {report_path}")


# ====================================================================== #
#  CLI
# ====================================================================== #
def main():
    parser = argparse.ArgumentParser(
        description="pH-Dependent Interface Interaction Analysis (PyRosetta)")
    parser.add_argument("pdb_path", help="PDB file path")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Output directory (default: same dir as PDB)")
    parser.add_argument("--binder", nargs="+", default=["A", "B"],
                        help="Binder chain IDs (default: A B)")
    parser.add_argument("--target", nargs="+", default=["C"],
                        help="Target chain IDs (default: C)")
    parser.add_argument("--interface-dist", type=float, default=8.0,
                        help="Interface distance cutoff in Å (default: 8.0)")
    parser.add_argument("--salt-bridge-dist", type=float, default=4.0,
                        help="Salt bridge distance cutoff in Å (default: 4.0)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(os.path.dirname(args.pdb_path), "interaction_analysis")

    analyzer = InterfaceInteractionAnalyzer(
        pdb_path=args.pdb_path,
        binder_chains=args.binder,
        target_chains=args.target,
        interface_dist=args.interface_dist,
        salt_bridge_dist=args.salt_bridge_dist,
        verbose=args.verbose,
    )
    analyzer.run(args.output_dir)


if __name__ == "__main__":
    main()
