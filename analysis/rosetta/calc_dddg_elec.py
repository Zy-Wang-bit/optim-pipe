#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量计算dddG_elec（delta delta electrostatic ddG）
基于论文方法：计算pH 5和pH 7条件下的静电ddG差异

dddG_elec = ddG_elec(pH5) - ddG_elec(pH7)

其中ddG_elec = E_elec(complex) - E_elec(binder) - E_elec(target)
"""

import sys
import os
import glob
from collections import defaultdict
import pandas as pd
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring import ScoreType, ScoreFunction
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover


class dddGElecCalculator:
    """计算dddG_elec"""
    
    def __init__(self, pdb_path, binder_chains=["A", "B"], target_chains=["C"], verbose=False):
        """
        Args:
            pdb_path: PDB文件路径
            binder_chains: binder链ID列表(抗体重链和轻链)
            target_chains: target链ID列表(抗原)
            verbose: 是否打印详细信息
        """
        self.pdb_path = pdb_path
        self.binder_chains = set(binder_chains)
        self.target_chains = set(target_chains)
        self.verbose = verbose
        
        # 初始化PyRosetta（静默模式）
        if not hasattr(pyrosetta, '_initialized'):
            pyrosetta.init("-mute all")
            pyrosetta._initialized = True
        
        # 加载pose
        self.pose = pyrosetta.pose_from_pdb(pdb_path)
        
        # 构建链ID映射
        self.chain_map = self._build_chain_map()
        
        # 创建只含fa_elec的ScoreFunction
        self.sf_elec = ScoreFunction()
        self.sf_elec.set_weight(ScoreType.fa_elec, 1.0)
        
    def _build_chain_map(self):
        """构建残基序号到PDB链ID的映射"""
        chain_map = {}
        pdb_info = self.pose.pdb_info()
        for i in range(1, self.pose.total_residue() + 1):
            chain_id = pdb_info.chain(i)
            chain_map[i] = chain_id
        return chain_map
    
    def _get_chain_residues(self, pose, chain_ids):
        """获取指定链的残基索引列表"""
        residues = []
        pdb_info = pose.pdb_info()
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.chain(i) in chain_ids:
                residues.append(i)
        return residues
    
    def _get_histidine_residues(self, pose):
        """获取所有His残基的索引"""
        his_residues = []
        for i in range(1, pose.total_residue() + 1):
            if pose.residue(i).name3().startswith("HIS"):
                his_residues.append(i)
        return his_residues
    
    def _extract_chains(self, pose, chain_ids):
        """提取指定链创建新pose"""
        from pyrosetta.rosetta.core.pose import Pose as RosettaPose
        from pyrosetta.rosetta.core.pose import append_subpose_to_pose
        
        new_pose = Pose()
        pdb_info = pose.pdb_info()
        
        # 找到要保留的残基
        residues_to_keep = []
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.chain(i) in chain_ids:
                residues_to_keep.append(i)
        
        if not residues_to_keep:
            return new_pose
        
        # 使用slice方式提取
        # 找连续区段
        segments = []
        start = residues_to_keep[0]
        end = start
        for r in residues_to_keep[1:]:
            if r == end + 1:
                end = r
            else:
                segments.append((start, end))
                start = r
                end = r
        segments.append((start, end))
        
        # 提取第一个segment
        first_start, first_end = segments[0]
        new_pose = Pose(pose, first_start, first_end)
        
        # 追加其他segments
        for seg_start, seg_end in segments[1:]:
            temp_pose = Pose(pose, seg_start, seg_end)
            append_subpose_to_pose(new_pose, temp_pose, 1, temp_pose.total_residue())
        
        return new_pose
    
    def _repack_pose(self, pose):
        """对pose进行repack"""
        pose_copy = Pose()
        pose_copy.assign(pose)
        
        # 创建TaskFactory，只允许repack不允许design
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        
        packer_task = tf.create_task_and_apply_taskoperations(pose_copy)
        
        # 使用标准ScoreFunction进行repack
        sf_full = pyrosetta.get_score_function()
        
        packer = PackRotamersMover(sf_full)
        packer.task_factory(tf)
        packer.apply(pose_copy)
        
        return pose_copy
    
    def _repack_his_only(self, pose, ph_value):
        """
        只对His残基进行repack
        pH 5时His按双质子化处理（通过设置pH模式）
        """
        pose_copy = Pose()
        pose_copy.assign(pose)
        
        # 获取His残基
        his_residues = self._get_histidine_residues(pose_copy)
        
        if not his_residues:
            return pose_copy
        
        # 创建TaskFactory
        from pyrosetta.rosetta.core.pack.task.operation import PreventRepacking
        
        tf = TaskFactory()
        tf.push_back(RestrictToRepacking())
        
        packer_task = tf.create_task_and_apply_taskoperations(pose_copy)
        
        # 只允许His残基repack
        for i in range(1, pose_copy.total_residue() + 1):
            if i not in his_residues:
                packer_task.nonconst_residue_task(i).prevent_repacking()
        
        # 使用标准ScoreFunction进行repack
        sf_full = pyrosetta.get_score_function()
        
        packer = PackRotamersMover(sf_full)
        packer.apply(pose_copy)
        
        return pose_copy
    
    def _calc_elec_energy(self, pose):
        """计算fa_elec能量"""
        return self.sf_elec(pose)
    
    def _calc_ddg_elec(self, pose):
        """
        计算静电ddG
        ddG_elec = E_elec(complex) - E_elec(binder) - E_elec(target)
        """
        # 复合物能量
        E_complex = self._calc_elec_energy(pose)
        
        # 提取binder和target
        binder_pose = self._extract_chains(pose, self.binder_chains)
        target_pose = self._extract_chains(pose, self.target_chains)
        
        if binder_pose.total_residue() == 0 or target_pose.total_residue() == 0:
            return None, None, None
        
        # 分别计算能量
        E_binder = self._calc_elec_energy(binder_pose)
        E_target = self._calc_elec_energy(target_pose)
        
        ddG_elec = E_complex - E_binder - E_target
        
        return ddG_elec, E_complex, E_binder, E_target
    
    def _simulate_ph5_his_charge(self, pose):
        """
        模拟pH 5条件下His的质子化
        通过尝试将HIS替换为HIS_D（双质子化）
        如果失败则返回原pose
        """
        pose_copy = Pose()
        pose_copy.assign(pose)
        
        his_residues = self._get_histidine_residues(pose_copy)
        
        for his_idx in his_residues:
            try:
                # 尝试替换为双质子化His
                res_type_set = pose_copy.residue(his_idx).residue_type_set()
                
                # 尝试不同的双质子化His名称
                for his_type in ["HIS_D", "HIP", "HSP", "HIS:protonated"]:
                    try:
                        new_type = res_type_set.name_map(his_type)
                        pyrosetta.rosetta.core.pose.replace_pose_residue_copying_existing_coordinates(
                            pose_copy, his_idx, new_type
                        )
                        break
                    except:
                        continue
            except:
                # 如果替换失败，继续使用原类型
                pass
        
        return pose_copy
    
    def calculate_dddg_elec(self):
        """
        计算dddG_elec
        
        Returns:
            dict: 包含各项计算结果
        """
        if self.verbose:
            print(f"处理: {os.path.basename(self.pdb_path)}")
        
        # 先对整体结构进行repack
        pose_repacked = self._repack_pose(self.pose)
        
        # pH 7.4条件：使用默认His状态
        pose_ph7 = Pose()
        pose_ph7.assign(pose_repacked)
        
        # 对His进行repack（pH 7.4条件）
        pose_ph7 = self._repack_his_only(pose_ph7, 7.4)
        
        # 计算pH 7.4的ddG_elec
        result_ph7 = self._calc_ddg_elec(pose_ph7)
        if result_ph7[0] is None:
            return {
                "dddG_elec": None,
                "ddG_elec_pH7": None,
                "ddG_elec_pH5": None,
                "error": "Failed to extract chains"
            }
        ddG_elec_pH7, E_complex_pH7, E_binder_pH7, E_target_pH7 = result_ph7
        
        # pH 5.0条件：尝试模拟His双质子化
        pose_ph5 = self._simulate_ph5_his_charge(pose_repacked)
        
        # 对His进行repack（pH 5.0条件）
        pose_ph5 = self._repack_his_only(pose_ph5, 5.0)
        
        # 计算pH 5.0的ddG_elec
        result_ph5 = self._calc_ddg_elec(pose_ph5)
        if result_ph5[0] is None:
            return {
                "dddG_elec": None,
                "ddG_elec_pH7": ddG_elec_pH7,
                "ddG_elec_pH5": None,
                "error": "Failed to extract chains at pH5"
            }
        ddG_elec_pH5, E_complex_pH5, E_binder_pH5, E_target_pH5 = result_ph5
        
        # 计算dddG_elec
        dddG_elec = ddG_elec_pH5 - ddG_elec_pH7
        
        # 统计His数量
        his_residues = self._get_histidine_residues(self.pose)
        binder_his = [h for h in his_residues if self.chain_map.get(h) in self.binder_chains]
        
        return {
            "dddG_elec": dddG_elec,
            "ddG_elec_pH7": ddG_elec_pH7,
            "ddG_elec_pH5": ddG_elec_pH5,
            "E_complex_pH7": E_complex_pH7,
            "E_binder_pH7": E_binder_pH7,
            "E_target_pH7": E_target_pH7,
            "E_complex_pH5": E_complex_pH5,
            "E_binder_pH5": E_binder_pH5,
            "E_target_pH5": E_target_pH5,
            "n_histidines_total": len(his_residues),
            "n_histidines_binder": len(binder_his),
            "error": None
        }


def batch_calculate_dddg_elec(pdb_dir, output_csv=None, verbose=False):
    """
    批量计算目录下所有PDB的dddG_elec
    
    Args:
        pdb_dir: PDB文件目录
        output_csv: 输出CSV路径（默认为代码所在目录下的dddg_elec_scores.csv）
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
    
    # 逐个处理
    for i, pdb_path in enumerate(sorted(pdb_files), 1):
        pdb_name = os.path.basename(pdb_path)
        print(f"[{i}/{len(pdb_files)}] 处理 {pdb_name}...", end=" ")
        
        try:
            calc = dddGElecCalculator(
                pdb_path=pdb_path,
                binder_chains=["A", "B"],
                target_chains=["C"],
                # binder_chains=["A"],
                # target_chains=["B"],
                verbose=verbose
            )
            
            result = calc.calculate_dddg_elec()
            
            # 构建一行数据
            row = {
                "pdb_name": pdb_name,
                "status": "success" if result["error"] is None else f"error: {result['error']}",
                "dddG_elec": result["dddG_elec"],
                "ddG_elec_pH7": result["ddG_elec_pH7"],
                "ddG_elec_pH5": result["ddG_elec_pH5"],
                "E_complex_pH7": result.get("E_complex_pH7"),
                "E_binder_pH7": result.get("E_binder_pH7"),
                "E_target_pH7": result.get("E_target_pH7"),
                "E_complex_pH5": result.get("E_complex_pH5"),
                "E_binder_pH5": result.get("E_binder_pH5"),
                "E_target_pH5": result.get("E_target_pH5"),
                "n_histidines_total": result.get("n_histidines_total"),
                "n_histidines_binder": result.get("n_histidines_binder"),
            }
            
            results.append(row)
            
            if result["dddG_elec"] is not None:
                print(f"✓ dddG_elec={result['dddG_elec']:.3f}")
            else:
                print(f"✗ {result['error']}")
            
        except Exception as e:
            print(f"✗ 失败: {str(e)}")
            row = {
                "pdb_name": pdb_name,
                "status": f"error: {str(e)}",
                "dddG_elec": None,
                "ddG_elec_pH7": None,
                "ddG_elec_pH5": None,
                "E_complex_pH7": None,
                "E_binder_pH7": None,
                "E_target_pH7": None,
                "E_complex_pH5": None,
                "E_binder_pH5": None,
                "E_target_pH5": None,
                "n_histidines_total": None,
                "n_histidines_binder": None,
            }
            results.append(row)
    
    # 转为DataFrame
    df = pd.DataFrame(results)
    
    # 调整列顺序
    cols = [
        "pdb_name", "status", 
        "dddG_elec", "ddG_elec_pH7", "ddG_elec_pH5",
        "E_complex_pH7", "E_binder_pH7", "E_target_pH7",
        "E_complex_pH5", "E_binder_pH5", "E_target_pH5",
        "n_histidines_total", "n_histidines_binder"
    ]
    df = df[cols]
    
    # 输出CSV
    if output_csv is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_csv = os.path.join(script_dir, "dddg_elec_scores.csv")
    
    df.to_csv(output_csv, index=False)
    
    print(f"\n{'='*60}")
    print(f"完成！共处理 {len(pdb_files)} 个PDB文件")
    print(f"成功: {len([r for r in results if 'success' in r['status']])} 个")
    print(f"失败: {len([r for r in results if 'success' not in r['status']])} 个")
    print(f"结果已保存到: {output_csv}")
    print(f"{'='*60}")
    
    return df


def main():
    if len(sys.argv) < 2:
        print("用法: python batch_calc_dddg_elec.py <pdb_directory> [output_csv]")
        print("示例: python batch_calc_dddg_elec.py ./pdbs/")
        print("      python batch_calc_dddg_elec.py ./pdbs/ results.csv")
        sys.exit(1)
    
    pdb_dir = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.isdir(pdb_dir):
        print(f"错误: {pdb_dir} 不是有效的目录")
        sys.exit(1)
    
    batch_calculate_dddg_elec(pdb_dir, output_csv, verbose=False)


if __name__ == "__main__":
    main()