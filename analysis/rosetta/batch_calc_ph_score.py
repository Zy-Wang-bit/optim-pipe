#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量计算pH score
对目录下所有PDB文件计算pH score并输出CSV
"""

import sys
import os
import glob
from collections import defaultdict
import pandas as pd
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring.sasa import SasaCalc


class pHScoreCalculator:
    """计算pH score"""
    
    def __init__(self, pdb_path, binder_chains=["A", "B"], target_chains=["C"], 
                 buried_cutoff=0.20, verbose=False):
        """
        Args:
            pdb_path: PDB文件路径
            binder_chains: binder链ID列表(抗体重链和轻链)
            target_chains: target链ID列表(抗原)
            buried_cutoff: buried判断阈值(SASA相对值)
            verbose: 是否打印详细信息
        """
        self.pdb_path = pdb_path
        self.binder_chains = set(binder_chains)
        self.target_chains = set(target_chains)
        self.buried_cutoff = buried_cutoff
        self.verbose = verbose
        
        # 初始化PyRosetta（静默模式）
        if not hasattr(pyrosetta, '_initialized'):
            pyrosetta.init("-mute all")
            pyrosetta._initialized = True
        
        # 加载pose
        self.pose = pyrosetta.pose_from_pdb(pdb_path)
        
        # 构建链ID映射(Rosetta内部编号 -> PDB链ID)
        self.chain_map = self._build_chain_map()
        
    def _build_chain_map(self):
        """构建残基序号到PDB链ID的映射"""
        chain_map = {}
        pdb_info = self.pose.pdb_info()
        for i in range(1, self.pose.total_residue() + 1):
            chain_id = pdb_info.chain(i)
            chain_map[i] = chain_id
        return chain_map
    
    def _get_histidines(self):
        """获取所有His残基"""
        histidines = []
        for i in range(1, self.pose.total_residue() + 1):
            if self.pose.residue(i).name3().startswith("HIS"):
                chain = self.chain_map[i]
                histidines.append((i, chain))
        return histidines
    
    def _set_his_protonation(self, pose, ph_value):
        """
        实现逻辑：
        1. 状态切换：
           - pH < 6.0: 强制突变为 HIS_D（双质子化态，带正电）。这等效于论文中 pH=0 的物理环境。
           - pH >= 6.0: 强制突变为 HIS（中性态）。
        2. 几何优化 (Repack)：
           - 仅对 His 残基进行侧链重排。
           - 锁定周围所有其他残基（"Control Variable"），确保能量/氢键的变化完全来源于 His 的质子化改变。
           - 使用 ref2015 全原子打分函数，确保生成的构象在几何上没有碰撞（Clash），从而能准确统计氢键。
        """
        # 在函数内部引入必要组件，保持模块独立性
        from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
        from pyrosetta.rosetta.core.pack.task import TaskFactory
        from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
        
        # 1. 深度复制 Pose，确保不污染原始结构
        new_pose = pose.clone()
        
        # 2. 确定目标残基类型
        # HIS_D: Rosetta 内部定义的双质子化组氨酸 (delta H + epsilon H)，净电荷 +1
        # HIS:   标准中性组氨酸，Rosetta 会在后续 Repack 中自动采样 tautomers (HID/HIE)
        target_his_type = "HIS_D" if ph_value < 6.0 else "HIS"
        
        # 3. 收集所有 His 的残基索引
        # 注意：Rosetta 索引从 1 开始
        his_indices = []
        for i in range(1, new_pose.total_residue() + 1):
            name3 = new_pose.residue(i).name3()
            if name3.startswith("HIS"):
                his_indices.append(i)

        # 如果没有 His，直接返回副本
        if not his_indices:
            return new_pose

        # 4. 执行突变 (状态切换)
        for i in his_indices:
            current_name = new_pose.residue(i).name3()
            # 只有当类型不一致时才突变 (避免不必要的计算)
            if current_name != target_his_type:
                # MutateResidue 会将残基替换为目标类型
                # set_preserve_atom_coords(False) 允许侧链坐标重置，因为质子化改变了化学拓扑
                mutator = MutateResidue(i, target_his_type)
                mutator.set_preserve_atom_coords(False) 
                mutator.apply(new_pose)
        
        # 5. 配置侧链重排任务 (Repacking)
        # 目标：对应论文 "repacking all histidines"，严格限制优化范围
        tf = TaskFactory()
        task = tf.create_packer_task(new_pose)
        
        # 5.1 初始化：先冻结整个蛋白 (Prevent Repacking)
        task.restrict_to_repacking() # 全局禁止设计(改变氨基酸种类)
        for i in range(1, new_pose.total_residue() + 1):
            task.nonconst_residue_task(i).prevent_repacking() # 禁止重排侧链
            
        # 5.2 解冻：仅允许 His 残基重排侧链
        for i in his_indices:
            task.nonconst_residue_task(i).restrict_to_repacking()
            # 此时，非 His 残基完全静止，His 残基可以旋转 Chi 角寻找最佳氢键位置

        # 6. 执行优化
        # 使用 ref2015 标准打分函数，包含范德华力(fa_rep)和静电(fa_elec)
        # 这比论文中只用 e_pH 更科学，能避免质子化后的侧链与其他原子发生碰撞
        sfxn = pyrosetta.create_score_function("ref2015")
        packer = PackRotamersMover(sfxn, task)
        packer.apply(new_pose)
        
        return new_pose
        
    def _calculate_sasa(self, pose):
        """计算每个残基的SASA"""
        sasa_calc = SasaCalc()
        sasa_calc.calculate(pose)
        
        sasa_values = {}
        for i in range(1, pose.total_residue() + 1):
            rsd_sasa = sasa_calc.get_residue_sasa()[i]
            sasa_values[i] = rsd_sasa
        return sasa_values
    
    def _is_buried(self, rsd_idx, sasa_values):
        """判断残基是否buried"""
        return sasa_values.get(rsd_idx, 100) < 40
    
    def _is_backbone_atom(self, pose, rsd_idx, atm_idx):
        """判断原子是否为主链原子（兼容不同PyRosetta版本）"""
        BACKBONE_ATOMS = {"N", "CA", "C", "O", "H", "HA"}
        try:
            atm_name = pose.residue(rsd_idx).atom_name(atm_idx).strip()
            return atm_name in BACKBONE_ATOMS
        except Exception:
            return False

    def _analyze_hbonds(self, pose):
        """分析氢键网络"""
        pose.update_residue_neighbors()
        hbond_set = pose.get_hbonds(False, False, False, False, False)

        total_hbonds = hbond_set.nhbonds()
        his_any = 0
        his_cross = 0

        his_hbonds = defaultdict(lambda: {"acceptor": [], "donor": []})

        for i in range(1, total_hbonds + 1):
            hbond = hbond_set.hbond(i)

            don_rsd = hbond.don_res()
            acc_rsd = hbond.acc_res()

            don_chain = self.chain_map.get(don_rsd)
            acc_chain = self.chain_map.get(acc_rsd)

            don_res_name = pose.residue(don_rsd).name3()
            acc_res_name = pose.residue(acc_rsd).name3()

            involves_his = (don_res_name.startswith("HIS") or
                            acc_res_name.startswith("HIS"))
            if involves_his:
                his_any += 1

            is_cross_chain = (
                (don_chain in self.binder_chains and acc_chain in self.target_chains) or
                (don_chain in self.target_chains and acc_chain in self.binder_chains)
            )

            if involves_his and is_cross_chain:
                his_cross += 1

            if not is_cross_chain:
                continue

            # His 作为 donor
            if don_res_name.startswith("HIS") and don_chain in self.binder_chains:
                is_backbone = self._is_backbone_atom(pose, don_rsd, hbond.don_hatm())
                his_hbonds[don_rsd]["donor"].append({
                    "partner": acc_rsd,
                    "backbone": is_backbone,
                    "energy": hbond.energy()
                })

            # His 作为 acceptor
            if acc_res_name.startswith("HIS") and acc_chain in self.binder_chains:
                is_backbone = self._is_backbone_atom(pose, acc_rsd, hbond.acc_atm())
                his_hbonds[acc_rsd]["acceptor"].append({
                    "partner": don_rsd,
                    "backbone": is_backbone,
                    "energy": hbond.energy()
                })

        if self.verbose:
            print(f"  总氢键: {total_hbonds}, His氢键: {his_any}, 跨链His氢键: {his_cross}")

        return his_hbonds    
    
    def _classify_his_pattern(self, his_idx, hbonds):
        """分类His的氢键模式"""
        n_acceptor = len(hbonds["acceptor"])
        n_donor = len(hbonds["donor"])
        
        has_bb_acceptor = any(h["backbone"] for h in hbonds["acceptor"])
        has_bb_donor = any(h["backbone"] for h in hbonds["donor"])
        
        patterns = []
        
        if n_acceptor > 0 and n_donor == 0:
            patterns.append("A")
            if has_bb_acceptor:
                patterns.append("A_bb")
        elif n_donor > 0 and n_acceptor == 0:
            patterns.append("D")
            if has_bb_donor:
                patterns.append("D_bb")
        elif n_acceptor > 0 and n_donor > 0:
            patterns.append("AD")
            if has_bb_acceptor or has_bb_donor:
                patterns.append("AD_bb")
        
        if n_donor >= 2:
            patterns.append("DD")
            if has_bb_donor:
                patterns.append("DD_bb")
        
        return patterns
    
    def calculate_ph_score(self):
        """计算pH score"""
        if self.verbose:
            print(f"处理: {os.path.basename(self.pdb_path)}")
        
        # 生成两个pH条件的pose
        pose_ph7 = self._set_his_protonation(self.pose, 7.0)
        pose_ph5 = self._set_his_protonation(self.pose, 5.0)
        
        # 计算SASA
        sasa_ph7 = self._calculate_sasa(pose_ph7)
        sasa_ph5 = self._calculate_sasa(pose_ph5)
        
        # 分析氢键
        hbonds_ph7 = self._analyze_hbonds(pose_ph7)
        hbonds_ph5 = self._analyze_hbonds(pose_ph5)
        
        # 统计各项term
        terms = defaultdict(int)
        
        # 获取所有His
        histidines = self._get_histidines()
        
        # 只关注binder上的His
        binder_his = [(idx, chain) for idx, chain in histidines 
                      if chain in self.binder_chains]
        
        # pH 5条件
        for his_idx, chain in binder_his:
            if his_idx not in hbonds_ph5:
                continue
            
            is_buried = self._is_buried(his_idx, sasa_ph5)
            location = "core" if is_buried else "surf"
            
            patterns = self._classify_his_pattern(his_idx, hbonds_ph5[his_idx])
            
            for pattern in patterns:
                term_name = f"low_pH_{location}_{pattern}"
                terms[term_name] += 1
        
        # pH 7条件
        for his_idx, chain in binder_his:
            if his_idx not in hbonds_ph7:
                continue
            
            is_buried = self._is_buried(his_idx, sasa_ph7)
            location = "core" if is_buried else "surf"
            
            patterns = self._classify_his_pattern(his_idx, hbonds_ph7[his_idx])
            
            for pattern in patterns:
                term_name = f"high_pH_{location}_{pattern}"
                terms[term_name] += 1
        
        # 计算pH score(根据论文的公式)
        weights = {
            "low_pH_surf_D": 1.0,
            "low_pH_surf_DD": 3.0,
            "low_pH_surf_D_bb": 2.0,
            "low_pH_surf_DD_bb": 6.0,
            "low_pH_core_D": 3.0,
            "low_pH_core_DD": 9.0,
            "low_pH_core_D_bb": 6.0,
            "low_pH_core_DD_bb": 18.0,
            "high_pH_surf_A": 1.0,
            "high_pH_surf_D": 1.0,
            "high_pH_surf_AD": 3.0,
            "high_pH_surf_A_bb": 2.0,
            "high_pH_surf_D_bb": 2.0,
            "high_pH_surf_AD_bb": 6.0,
            "high_pH_core_A": 3.0,
            "high_pH_core_D": 3.0,
            "high_pH_core_AD": 9.0,
            "high_pH_core_A_bb": 6.0,
            "high_pH_core_D_bb": 6.0,
            "high_pH_core_AD_bb": 18.0,
        }
        
        ph_score = 0.0
        for term, weight in weights.items():
            count = terms[term]
            contribution = weight * count
            ph_score += contribution
        
        return {
            "ph_score": ph_score,
            "terms": dict(terms),
            "n_histidines_total": len(histidines),
            "n_histidines_binder": len(binder_his),
        }


def batch_calculate_ph_score(pdb_dir, output_csv=None, verbose=False):
    """
    批量计算目录下所有PDB的pH score
    
    Args:
        pdb_dir: PDB文件目录
        output_csv: 输出CSV路径（默认为代码所在目录下的ph_scores.csv）
        verbose: 是否打印详细信息
    """
    # 查找所有PDB文件
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    if not pdb_files:
        print(f"错误: 在 {pdb_dir} 下没有找到PDB文件")
        return
    
    print(f"找到 {len(pdb_files)} 个PDB文件")
    
    # 准备结果列表
    results = []
    
    # 所有可能的term
    all_terms = [
        "low_pH_surf_D", "low_pH_surf_DD", "low_pH_surf_D_bb", "low_pH_surf_DD_bb",
        "low_pH_core_D", "low_pH_core_DD", "low_pH_core_D_bb", "low_pH_core_DD_bb",
        "high_pH_surf_A", "high_pH_surf_D", "high_pH_surf_AD", 
        "high_pH_surf_A_bb", "high_pH_surf_D_bb", "high_pH_surf_AD_bb",
        "high_pH_core_A", "high_pH_core_D", "high_pH_core_AD",
        "high_pH_core_A_bb", "high_pH_core_D_bb", "high_pH_core_AD_bb",
    ]
    
    # 逐个处理
    for i, pdb_path in enumerate(sorted(pdb_files), 1):
        pdb_name = os.path.basename(pdb_path)
        print(f"[{i}/{len(pdb_files)}] 处理 {pdb_name}...", end=" ")
        
        try:
            calc = pHScoreCalculator(
                pdb_path=pdb_path,
                # binder_chains=["A", "B"],
                # target_chains=["C"],
                binder_chains=["A"],
                target_chains=["B"],
                buried_cutoff=0.20,
                verbose=verbose
            )
            
            result = calc.calculate_ph_score()
            
            # 构建一行数据
            row = {
                "pdb_name": pdb_name,
                "ph_score": result["ph_score"],
                "n_histidines_total": result["n_histidines_total"],
                "n_histidines_binder": result["n_histidines_binder"],
                "status": "success"
            }
            
            # 添加所有term的count
            for term in all_terms:
                row[term] = result["terms"].get(term, 0)
            
            results.append(row)
            print(f"✓ pH_score={result['ph_score']:.1f}")
            
        except Exception as e:
            print(f"✗ 失败: {str(e)}")
            row = {
                "pdb_name": pdb_name,
                "ph_score": None,
                "n_histidines_total": None,
                "n_histidines_binder": None,
                "status": f"error: {str(e)}"
            }
            for term in all_terms:
                row[term] = None
            results.append(row)
    
    # 转为DataFrame
    df = pd.DataFrame(results)
    
    # 调整列顺序
    cols = ["pdb_name", "status", "ph_score", "n_histidines_total", "n_histidines_binder"] + all_terms
    df = df[cols]
    
    # 输出CSV
    if output_csv is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_csv = os.path.join(script_dir, "ph_scores.csv")
    
    df.to_csv(output_csv, index=False)
    
    print(f"\n{'='*60}")
    print(f"完成！共处理 {len(pdb_files)} 个PDB文件")
    print(f"成功: {len([r for r in results if r['status'] == 'success'])} 个")
    print(f"失败: {len([r for r in results if r['status'] != 'success'])} 个")
    print(f"结果已保存到: {output_csv}")
    print(f"{'='*60}")
    
    return df


def main():
    if len(sys.argv) < 2:
        print("用法: python batch_calc_ph_score.py <pdb_directory> [output_csv]")
        print("示例: python batch_calc_ph_score.py ./pdbs/")
        print("      python batch_calc_ph_score.py ./pdbs/ results.csv")
        sys.exit(1)
    
    pdb_dir = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.isdir(pdb_dir):
        print(f"错误: {pdb_dir} 不是有效的目录")
        sys.exit(1)
    
    batch_calculate_ph_score(pdb_dir, output_csv, verbose=False)


if __name__ == "__main__":
    main()