#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量计算 dddG_elec (Debug Version)
修复：只对突变引入的His进行pH5质子化，而非所有天然His

原理：
dddG_elec = ddG_elec(pH5) - ddG_elec(pH7)
只对突变位点的His进行HIP替换，评估该位点的静电效应。
"""

import sys
import os
import re
import glob
import pandas as pd
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring import ScoreType, ScoreFunction
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pose import replace_pose_residue_copying_existing_coordinates

# ================= 用户配置区域 =================
DEFAULT_BINDER_CHAINS = ["A", "B"]
DEFAULT_TARGET_CHAINS = ["C"]
# ===============================================


class dddGElecCalculator:
    def __init__(self, pdb_path, binder_chains, target_chains, verbose=False):
        self.pdb_path = pdb_path
        self.binder_chains = set(binder_chains)
        self.target_chains = set(target_chains)
        self.verbose = verbose

        if not hasattr(pyrosetta, '_initialized'):
            pyrosetta.init("-mute all")
            pyrosetta._initialized = True

        try:
            self.pose = pyrosetta.pose_from_pdb(pdb_path)
        except Exception as e:
            raise RuntimeError(f"无法加载 PDB: {pdb_path} | 错误: {e}")

        self.chain_map = self._build_chain_map()
        self._validate_chains()

        self.sf_elec = ScoreFunction()
        self.sf_elec.set_weight(ScoreType.fa_elec, 1.0)
        self.sf_full = pyrosetta.get_score_function()

    def _build_chain_map(self):
        chain_map = {}
        pdb_info = self.pose.pdb_info()
        for i in range(1, self.pose.total_residue() + 1):
            chain_map[i] = pdb_info.chain(i)
        return chain_map

    def _validate_chains(self):
        existing = set(self.chain_map.values())
        if not self.binder_chains.issubset(existing):
            raise ValueError(f"Binder 链缺失. 需: {self.binder_chains}, 存: {existing}")
        if not self.target_chains.issubset(existing):
            raise ValueError(f"Target 链缺失. 需: {self.target_chains}, 存: {existing}")

    def _parse_mutation_site(self, pdb_name):
        """
        从文件名解析突变位点
        支持格式: A_K50H.pdb, B_Y31H.pdb 等
        返回: (chain, resno) 或 (None, None)
        """
        name = pdb_name.replace('.pdb', '')
        # 匹配: 链_野生型残基+位置+H
        m = re.match(r'^([A-Z])_([A-Z])(\d+)H$', name)
        if m:
            chain, wt_aa, resno = m.groups()
            return chain, int(resno)
        return None, None

    def _find_rosetta_index(self, chain, resno):
        """将PDB链+残基号转换为Rosetta内部索引"""
        pdb_info = self.pose.pdb_info()
        for i in range(1, self.pose.total_residue() + 1):
            if pdb_info.chain(i) == chain and pdb_info.number(i) == resno:
                return i
        return None

    def _get_all_histidine_residues(self, pose):
        """获取所有His残基（用于统计）"""
        his_indices = []
        for i in range(1, pose.total_residue() + 1):
            n = pose.residue(i).name3()
            if n.startswith("HIS") or n in ("HIP", "HID", "HIE"):
                his_indices.append(i)
        return his_indices

    def _extract_chains(self, pose, chain_ids):
        from pyrosetta.rosetta.core.pose import append_subpose_to_pose
        new_pose = Pose()
        pdb_info = pose.pdb_info()
        residues = [i for i in range(1, pose.total_residue() + 1) if pdb_info.chain(i) in chain_ids]

        if not residues:
            return new_pose

        segments = []
        start = residues[0]
        end = start
        for r in residues[1:]:
            if r == end + 1:
                end = r
            else:
                segments.append((start, end))
                start = r
                end = r
        segments.append((start, end))

        if segments:
            s0, e0 = segments[0]
            new_pose = Pose(pose, s0, e0)
            for s, e in segments[1:]:
                sub = Pose(pose, s, e)
                append_subpose_to_pose(new_pose, sub, 1, sub.total_residue())
        return new_pose

    def _get_neighbors(self, pose, center_idx, radius_A=8.0):
        """获取 center_idx 周围 radius_A 内的所有残基索引（含自身）"""
        center_res = pose.residue(center_idx)
        center_xyz = center_res.nbr_atom_xyz()
        neighbors = []
        for i in range(1, pose.total_residue() + 1):
            res_xyz = pose.residue(i).nbr_atom_xyz()
            if center_xyz.distance(res_xyz) <= radius_A:
                neighbors.append(i)
        return neighbors

    def _repack_focused(self, pose, focus_residues=None):
        pose_copy = Pose()
        pose_copy.assign(pose)

        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        packer_task = tf.create_task_and_apply_taskoperations(pose_copy)

        if focus_residues is not None:
            focus_set = set(focus_residues)
            for i in range(1, pose_copy.total_residue() + 1):
                if i not in focus_set:
                    packer_task.nonconst_residue_task(i).prevent_repacking()

        packer = PackRotamersMover(self.sf_full)
        packer.task_factory(tf)
        packer.apply(pose_copy)
        return pose_copy

    def _simulate_ph5_his_charge(self, pose, target_residues=None):
        """
        只对指定残基进行HIP替换（模拟pH5质子化）
        
        Args:
            pose: 输入pose
            target_residues: 要质子化的残基索引列表。如果None则不处理任何残基。
        """
        pose_copy = Pose()
        pose_copy.assign(pose)

        if target_residues is None or len(target_residues) == 0:
            if self.verbose:
                print("  > 无目标残基需要质子化")
            return pose_copy, None

        # 验证目标残基确实是His
        valid_targets = []
        for idx in target_residues:
            if idx is None:
                continue
            n = pose_copy.residue(idx).name3()
            if n.startswith("HIS") or n in ("HIP", "HID", "HIE"):
                valid_targets.append(idx)
            elif self.verbose:
                print(f"  > Warning: 位置 {idx} 不是His ({n})，跳过")

        if not valid_targets:
            if self.verbose:
                print("  > 无有效His残基需要质子化")
            return pose_copy, None

        # 获取HIP残基类型
        rsd_set = pose_copy.residue_type_set_for_pose()
        target_name = None
        if rsd_set.has_name("HIP"):
            target_name = "HIP"
        elif rsd_set.has_name("HIS_D"):
            target_name = "HIS_D"
        else:
            raise RuntimeError("数据库中未找到 HIP 或 HIS_D")

        target_type = rsd_set.name_map(target_name)

        # 记录替换前后的状态用于诊断
        diag_info = {}
        
        for idx in valid_targets:
            # 替换前状态
            res_before = pose_copy.residue(idx)
            name_before = res_before.name3()
            natoms_before = res_before.natoms()
            
            try:
                replace_pose_residue_copying_existing_coordinates(pose_copy, idx, target_type)
                
                # 替换后状态
                res_after = pose_copy.residue(idx)
                name_after = res_after.name3()
                natoms_after = res_after.natoms()
                
                diag_info[idx] = {
                    "before": name_before,
                    "after": name_after,
                    "atoms_before": natoms_before,
                    "atoms_after": natoms_after,
                    "success": name_after != name_before
                }
                
            except Exception as e:
                diag_info[idx] = {"error": str(e)}

        # 更新Pose状态
        try:
            pose_copy.conformation().detect_disulfides()
            pose_copy.update_residue_neighbors()
        except:
            pass

        return pose_copy, diag_info

    def _calc_components(self, pose):
        E_c = self.sf_elec(pose)
        binder = self._extract_chains(pose, self.binder_chains)
        target = self._extract_chains(pose, self.target_chains)

        if binder.total_residue() == 0 or target.total_residue() == 0:
            return None, None, None

        return E_c, self.sf_elec(binder), self.sf_elec(target)

    def run(self):
        pdb_name = os.path.basename(self.pdb_path)
        if self.verbose:
            print(f"Processing: {pdb_name}")

        # 解析突变位点
        mut_chain, mut_resno = self._parse_mutation_site(pdb_name)
        mutation_idx = None
        
        if mut_chain and mut_resno:
            mutation_idx = self._find_rosetta_index(mut_chain, mut_resno)
            if self.verbose:
                print(f"  > 突变位点: {mut_chain}_{mut_resno} -> Rosetta idx: {mutation_idx}")
        else:
            if self.verbose:
                print(f"  > Warning: 无法从文件名解析突变位点")

        # 统计总His数
        all_his = self._get_all_histidine_residues(self.pose)
        n_his_total = len(all_his)

        # 提前校验突变位点（避免无用 repack）
        if mutation_idx is None:
            return {
                "pdb_name": pdb_name,
                "status": "error",
                "error": "无法定位突变位点",
                "n_his_total": n_his_total
            }

        res_at_site = self.pose.residue(mutation_idx).name3()
        is_his = res_at_site.startswith("HIS") or res_at_site in ("HIP", "HID", "HIE")

        if not is_his:
            return {
                "pdb_name": pdb_name,
                "status": "error",
                "error": f"突变位点{mutation_idx}不是His，而是{res_at_site}（PDB可能是野生型）",
                "res_at_site": res_at_site,
                "rosetta_idx": mutation_idx
            }

        # 1. 预处理 (focused repack — 只 repack 突变位点周围 8Å)
        repack_shell = self._get_neighbors(self.pose, mutation_idx, radius_A=8.0)
        if self.verbose:
            print(f"  > Repack shell: {len(repack_shell)} residues within 8Å of idx {mutation_idx}")
        pose_ph7 = self._repack_focused(self.pose, focus_residues=repack_shell)

        # 2. pH 7 计算
        ec7, eb7, et7 = self._calc_components(pose_ph7)
        if ec7 is None:
            return {"pdb_name": pdb_name, "status": "error", "error": "pH7 extraction failed"}
        ddG_ph7 = ec7 - eb7 - et7
        
        try:
            pose_ph5_raw, diag_info = self._simulate_ph5_his_charge(pose_ph7, target_residues=[mutation_idx])
        except RuntimeError as e:
            return {"pdb_name": pdb_name, "status": "error", "error": str(e)}

        # 提取诊断信息
        replace_success = False
        res_before = res_at_site
        res_after = res_at_site
        if diag_info and mutation_idx in diag_info:
            d = diag_info[mutation_idx]
            if "error" not in d:
                res_before = d.get("before", "?")
                res_after = d.get("after", "?")
                replace_success = d.get("success", False)

        # 4. pH 5 Repack (same shell as pH7 for consistency)
        try:
            pose_ph5 = self._repack_focused(pose_ph5_raw, focus_residues=repack_shell)
        except Exception as e:
            return {"pdb_name": pdb_name, "status": "error", "error": f"Repack failed: {e}"}

        # 5. pH 5 计算
        ec5, eb5, et5 = self._calc_components(pose_ph5)
        if ec5 is None:
            return {"pdb_name": pdb_name, "status": "error", "error": "pH5 extraction failed"}
        ddG_ph5 = ec5 - eb5 - et5

        return {
            "pdb_name": pdb_name,
            "status": "success",
            "mutation": f"{mut_chain}_{mut_resno}H",
            "rosetta_idx": mutation_idx,
            "res_before": res_before,
            "res_after": res_after,
            "replace_ok": replace_success,
            "dddG_elec": ddG_ph5 - ddG_ph7,
            "ddG_elec_pH7": ddG_ph7,
            "ddG_elec_pH5": ddG_ph5,
            "delta_E_complex": ec5 - ec7,
            "n_his_total": n_his_total,
            "error": None
        }


def batch_process(pdb_dir, output_csv):
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    if not pdb_files:
        print(f"错误: {pdb_dir} 为空")
        return

    print(f"开始处理 {len(pdb_files)} 个文件 | Binder:{DEFAULT_BINDER_CHAINS} Target:{DEFAULT_TARGET_CHAINS}")
    print("=" * 70)
    print("【Debug模式】只对突变位点进行pH5质子化，输出诊断信息")
    print("=" * 70)

    results = []
    for i, f in enumerate(sorted(pdb_files), 1):
        name = os.path.basename(f)
        print(f"[{i}/{len(pdb_files)}] {name}...", end="", flush=True)

        try:
            calc = dddGElecCalculator(f, DEFAULT_BINDER_CHAINS, DEFAULT_TARGET_CHAINS)
            res = calc.run()

            if res.get("error"):
                print(f" ❌ {res['error']}")
            else:
                # 显示诊断信息
                replace_ok = "✓" if res.get('replace_ok') else "✗"
                print(f" {res['res_before']}->{res['res_after']}[{replace_ok}] dddG={res['dddG_elec']:+.4f}")
            results.append(res)

        except Exception as e:
            print(f" ❌ Error: {e}")
            results.append({"pdb_name": name, "status": "error", "error": str(e)})

    df = pd.DataFrame(results)
    cols = ["pdb_name", "status", "mutation", "rosetta_idx", 
            "res_before", "res_after", "replace_ok",
            "dddG_elec", "ddG_elec_pH7", "ddG_elec_pH5", 
            "delta_E_complex", "n_his_total", "error"]
    df = df[[c for c in cols if c in df.columns]]
    df.to_csv(output_csv, index=False)
    print(f"\n完成。结果已保存: {output_csv}")
    
    # 诊断统计
    success = df[df['status'] == 'success']
    if len(success) > 0:
        print(f"\n===== 诊断统计 (n={len(success)}) =====")
        
        # 替换成功率
        replace_ok_count = success['replace_ok'].sum() if 'replace_ok' in success.columns else 0
        print(f"残基替换成功: {replace_ok_count}/{len(success)}")
        
        # res_before/after 分布
        if 'res_before' in success.columns:
            print(f"替换前残基类型: {success['res_before'].value_counts().to_dict()}")
        if 'res_after' in success.columns:
            print(f"替换后残基类型: {success['res_after'].value_counts().to_dict()}")
        
        # dddG统计
        print(f"\ndddG_elec 范围: [{success['dddG_elec'].min():.4f}, {success['dddG_elec'].max():.4f}]")
        print(f"dddG_elec 唯一值: {success['dddG_elec'].round(4).nunique()}")
        near_zero = (success['dddG_elec'].abs() < 1e-6).sum()
        print(f"dddG_elec ≈ 0: {near_zero}/{len(success)}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python calc_dddg_elec_debug.py <pdb_dir> [output.csv]")
        sys.exit(1)
    batch_process(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "dddg_elec_debug_results.csv")