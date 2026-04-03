#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
计算抗体内部稳定性的 pH score (v2)

与原版的区别：
1. 检测范围：binder内部氢键（链内 + VH-VL界面），排除binder-target界面
2. Partner类型：区分 charged (His/Arg/Lys) vs non-charged
3. Buried判断：使用相对SASA (< 20%)
4. 权重体系：专注于内部稳定性，surf权重大幅降低
"""

import sys
import os
import glob
from collections import defaultdict
import pandas as pd
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring.sasa import SasaCalc


# ============================================================================
# 残基最大SASA参考值 (Gly-X-Gly三肽, 单位 Å²)
# 来源: Miller et al. / Sander & Rost
# ============================================================================
MAX_SASA = {
    "ALA": 113, "ARG": 241, "ASN": 158, "ASP": 151, "CYS": 140,
    "GLN": 189, "GLU": 183, "GLY": 85,  "HIS": 194, "ILE": 182,
    "LEU": 180, "LYS": 211, "MET": 204, "PHE": 218, "PRO": 143,
    "SER": 122, "THR": 146, "TRP": 259, "TYR": 229, "VAL": 160,
}

# Charged residues (论文核心关注的partner类型)
CHARGED_RESIDUES = {"HIS", "ARG", "LYS"}


class pHScoreCalculatorInternal:
    """计算抗体内部稳定性的pH score"""
    
    def __init__(self, pdb_path, binder_chains=["A", "B"], target_chains=["C"], 
                 buried_cutoff=0.20, verbose=False):
        """
        Args:
            pdb_path: PDB文件路径
            binder_chains: binder链ID列表(抗体重链和轻链)
            target_chains: target链ID列表(抗原) - 用于排除界面氢键
            buried_cutoff: buried判断阈值(相对SASA)
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
        
        # 构建映射
        self.chain_map = self._build_chain_map()
        self.resname_map = self._build_resname_map()
        
    def _build_chain_map(self):
        """构建残基序号到PDB链ID的映射"""
        chain_map = {}
        pdb_info = self.pose.pdb_info()
        for i in range(1, self.pose.total_residue() + 1):
            chain_id = pdb_info.chain(i)
            chain_map[i] = chain_id
        return chain_map
    
    def _build_resname_map(self):
        """构建残基序号到三字母名称的映射"""
        resname_map = {}
        for i in range(1, self.pose.total_residue() + 1):
            resname_map[i] = self.pose.residue(i).name3().strip()
        return resname_map
    
    def _get_histidines(self):
        """获取所有His残基"""
        histidines = []
        for i in range(1, self.pose.total_residue() + 1):
            if self.resname_map[i].startswith("HIS"):
                chain = self.chain_map[i]
                histidines.append((i, chain))
        return histidines
    
    def _calculate_relative_sasa(self, pose):
        """计算每个残基的相对SASA"""
        sasa_calc = SasaCalc()
        sasa_calc.calculate(pose)
        
        relative_sasa = {}
        for i in range(1, pose.total_residue() + 1):
            abs_sasa = sasa_calc.get_residue_sasa()[i]
            resname = self.resname_map[i]
            # 获取最大SASA，未知残基用200作为默认值
            max_sasa = MAX_SASA.get(resname, 200)
            relative_sasa[i] = abs_sasa / max_sasa if max_sasa > 0 else 1.0
        return relative_sasa
    
    def _is_buried(self, rsd_idx, relative_sasa_values):
        """判断残基是否buried (相对SASA < cutoff)"""
        return relative_sasa_values.get(rsd_idx, 1.0) < self.buried_cutoff
    
    def _is_charged_residue(self, rsd_idx):
        """判断残基是否是charged residue (His/Arg/Lys)"""
        resname = self.resname_map.get(rsd_idx, "")
        # 处理HIS的各种变体 (HIS, HIS_D, HIP等)
        if resname.startswith("HIS"):
            return True
        return resname in CHARGED_RESIDUES
    
    def _is_backbone_atom(self, pose, rsd_idx, atm_idx):
        """
        判断原子是否是backbone原子
        Backbone原子名称: N, CA, C, O, H, HA (以及一些变体)
        """
        try:
            residue = pose.residue(rsd_idx)
            atom_name = residue.atom_name(atm_idx).strip()
            # Backbone原子名称集合
            backbone_atoms = {"N", "CA", "C", "O", "H", "HA", "HA2", "HA3", 
                            "1H", "2H", "3H", "OXT"}
            return atom_name in backbone_atoms
        except Exception:
            return False
    
    def _analyze_hbonds_internal(self, pose):
        """
        分析抗体内部的氢键网络
        
        返回: his_hbonds dict
            key: his残基索引
            value: {
                "acceptor": [{partner, partner_charged, backbone, energy}, ...],
                "donor": [{partner, partner_charged, backbone, energy}, ...]
            }
        """
        pose.update_residue_neighbors()
        hbond_set = pose.get_hbonds(False, False, False, False, False)

        total_hbonds = hbond_set.nhbonds()
        
        # 统计信息
        stats = {
            "total": total_hbonds,
            "his_any": 0,
            "his_internal": 0,
            "his_internal_charged": 0,
        }

        his_hbonds = defaultdict(lambda: {"acceptor": [], "donor": []})

        for i in range(1, total_hbonds + 1):
            hbond = hbond_set.hbond(i)

            don_rsd = hbond.don_res()
            acc_rsd = hbond.acc_res()

            don_chain = self.chain_map.get(don_rsd)
            acc_chain = self.chain_map.get(acc_rsd)

            don_res_name = self.resname_map.get(don_rsd, "")
            acc_res_name = self.resname_map.get(acc_rsd, "")

            involves_his = (don_res_name.startswith("HIS") or
                            acc_res_name.startswith("HIS"))
            
            if involves_his:
                stats["his_any"] += 1

            # ============================================================
            # 核心修改：定义"内部氢键"
            # ============================================================
            # 内部氢键条件：
            # 1. 双方都在binder链上 (A-A, B-B, 或 A-B)
            # 2. 不能涉及target链
            
            don_in_binder = don_chain in self.binder_chains
            acc_in_binder = acc_chain in self.binder_chains
            don_in_target = don_chain in self.target_chains
            acc_in_target = acc_chain in self.target_chains
            
            # 必须双方都在binder，且不涉及target
            is_internal = (don_in_binder and acc_in_binder and 
                          not don_in_target and not acc_in_target)
            
            if not is_internal:
                continue
                
            if not involves_his:
                continue
                
            stats["his_internal"] += 1
            
            # 获取原子索引用于判断backbone
            don_hatm = hbond.don_hatm()  # donor的H原子
            acc_atm = hbond.acc_atm()    # acceptor的重原子
            
            # 检查partner是否是charged residue
            # His作为donor时，partner是acceptor残基
            # His作为acceptor时，partner是donor残基
            
            # His 作为 donor
            if don_res_name.startswith("HIS"):
                partner_charged = self._is_charged_residue(acc_rsd)
                # 检查donor的H是否是backbone H
                is_backbone = self._is_backbone_atom(pose, don_rsd, don_hatm)
                his_hbonds[don_rsd]["donor"].append({
                    "partner": acc_rsd,
                    "partner_name": acc_res_name,
                    "partner_charged": partner_charged,
                    "backbone": is_backbone,
                    "energy": hbond.energy()
                })
                if partner_charged:
                    stats["his_internal_charged"] += 1

            # His 作为 acceptor
            if acc_res_name.startswith("HIS"):
                partner_charged = self._is_charged_residue(don_rsd)
                # 检查acceptor原子是否是backbone原子
                is_backbone = self._is_backbone_atom(pose, acc_rsd, acc_atm)
                his_hbonds[acc_rsd]["acceptor"].append({
                    "partner": don_rsd,
                    "partner_name": don_res_name,
                    "partner_charged": partner_charged,
                    "backbone": is_backbone,
                    "energy": hbond.energy()
                })
                if partner_charged:
                    stats["his_internal_charged"] += 1

        if self.verbose:
            print(f"  总氢键: {stats['total']}, "
                  f"His相关: {stats['his_any']}, "
                  f"His内部: {stats['his_internal']}, "
                  f"His-charged: {stats['his_internal_charged']}")

        return his_hbonds, stats
    
    def _classify_his_pattern_extended(self, his_idx, hbonds):
        """
        分类His的氢键模式（扩展版，区分charged partner）
        
        返回: list of (pattern_name, has_charged_partner)
        """
        acceptors = hbonds["acceptor"]
        donors = hbonds["donor"]
        
        n_acceptor = len(acceptors)
        n_donor = len(donors)
        
        # 检查是否有charged partner
        has_charged_acceptor = any(h["partner_charged"] for h in acceptors)
        has_charged_donor = any(h["partner_charged"] for h in donors)
        has_any_charged = has_charged_acceptor or has_charged_donor
        
        # 检查是否涉及backbone
        has_bb_acceptor = any(h["backbone"] for h in acceptors)
        has_bb_donor = any(h["backbone"] for h in donors)
        
        patterns = []
        
        # 基础模式分类
        if n_acceptor > 0 and n_donor == 0:
            patterns.append(("A", has_charged_acceptor))
            if has_bb_acceptor:
                patterns.append(("A_bb", has_charged_acceptor))
                
        elif n_donor > 0 and n_acceptor == 0:
            patterns.append(("D", has_charged_donor))
            if has_bb_donor:
                patterns.append(("D_bb", has_charged_donor))
                
        elif n_acceptor > 0 and n_donor > 0:
            patterns.append(("AD", has_any_charged))
            if has_bb_acceptor or has_bb_donor:
                patterns.append(("AD_bb", has_any_charged))
        
        # DD模式（双donor，论文中特别关注）
        if n_donor >= 2:
            patterns.append(("DD", has_charged_donor))
            if has_bb_donor:
                patterns.append(("DD_bb", has_charged_donor))
        
        return patterns
    
    def calculate_ph_score(self):
        """计算pH score"""
        if self.verbose:
            print(f"处理: {os.path.basename(self.pdb_path)}")
        
        # 计算相对SASA
        relative_sasa = self._calculate_relative_sasa(self.pose)
        
        # 分析内部氢键
        his_hbonds, hbond_stats = self._analyze_hbonds_internal(self.pose)
        
        # 统计各项term
        terms = defaultdict(int)
        
        # 获取所有His
        histidines = self._get_histidines()
        
        # 只关注binder上的His
        binder_his = [(idx, chain) for idx, chain in histidines 
                      if chain in self.binder_chains]
        
        # 详细记录每个His的情况（用于调试）
        his_details = []
        
        for his_idx, chain in binder_his:
            if his_idx not in his_hbonds:
                continue
            
            hb = his_hbonds[his_idx]
            is_buried = self._is_buried(his_idx, relative_sasa)
            location = "core" if is_buried else "surf"
            
            patterns = self._classify_his_pattern_extended(his_idx, hb)
            
            # 记录详情
            pdb_info = self.pose.pdb_info()
            resno = pdb_info.number(his_idx)
            detail = {
                "his_idx": his_idx,
                "chain": chain,
                "resno": resno,
                "location": location,
                "rel_sasa": relative_sasa.get(his_idx, -1),
                "n_acceptor": len(hb["acceptor"]),
                "n_donor": len(hb["donor"]),
                "patterns": [p[0] for p in patterns],
                "has_charged": any(p[1] for p in patterns),
                "partners": [h["partner_name"] for h in hb["acceptor"] + hb["donor"]]
            }
            his_details.append(detail)
            
            for pattern, has_charged in patterns:
                # 构建term名称
                # 格式: {pH}_{location}_{pattern}[_chg]
                # 由于我们关注的是内部稳定性，pH标签反映的是：
                # - low_pH: His质子化时的行为（donor模式更强）
                # - high_pH: His未质子化时的行为（可以acceptor）
                
                # 简化处理：由于内部稳定性主要关注His-charged network
                # 我们直接记录pattern和charged状态
                
                base_term = f"{location}_{pattern}"
                terms[base_term] += 1
                
                if has_charged:
                    charged_term = f"{location}_{pattern}_chg"
                    terms[charged_term] += 1
        
        # ================================================================
        # 新的权重体系（专注于内部稳定性）
        # ================================================================
        # 核心原则：
        # 1. Core >> Surf（埋藏的氢键对稳定性影响更大）
        # 2. Charged partner >> Non-charged（His-His/Arg/Lys网络是关键）
        # 3. DD模式是强信号（双donor = 质子化状态）
        
        weights = {
            # Core + Charged: 最高权重
            "core_D_chg": 6.0,
            "core_DD_chg": 18.0,
            "core_D_bb_chg": 9.0,
            "core_DD_bb_chg": 27.0,
            "core_A_chg": 6.0,
            "core_AD_chg": 12.0,
            "core_A_bb_chg": 9.0,
            "core_AD_bb_chg": 18.0,
            
            # Core + Non-charged: 中等权重
            "core_D": 3.0,
            "core_DD": 9.0,
            "core_D_bb": 6.0,
            "core_DD_bb": 18.0,
            "core_A": 3.0,
            "core_AD": 6.0,
            "core_A_bb": 6.0,
            "core_AD_bb": 12.0,
            
            # Surf + Charged: 低权重（表面charged仍有一定意义）
            "surf_D_chg": 1.0,
            "surf_DD_chg": 3.0,
            "surf_D_bb_chg": 2.0,
            "surf_DD_bb_chg": 6.0,
            "surf_A_chg": 1.0,
            "surf_AD_chg": 2.0,
            "surf_A_bb_chg": 2.0,
            "surf_AD_bb_chg": 4.0,
            
            # Surf + Non-charged: 极低权重（几乎忽略）
            "surf_D": 0.1,
            "surf_DD": 0.3,
            "surf_D_bb": 0.2,
            "surf_DD_bb": 0.6,
            "surf_A": 0.1,
            "surf_AD": 0.2,
            "surf_A_bb": 0.2,
            "surf_AD_bb": 0.4,
        }
        
        # 计算总分
        ph_score = 0.0
        score_breakdown = {}
        for term, count in terms.items():
            weight = weights.get(term, 0.0)
            contribution = weight * count
            ph_score += contribution
            if contribution > 0:
                score_breakdown[term] = {"count": count, "weight": weight, "contrib": contribution}
        
        return {
            "ph_score": ph_score,
            "terms": dict(terms),
            "score_breakdown": score_breakdown,
            "n_histidines_total": len(histidines),
            "n_histidines_binder": len(binder_his),
            "n_his_with_internal_hbonds": len(his_hbonds),
            "hbond_stats": hbond_stats,
            "his_details": his_details,
        }


def batch_calculate_ph_score(pdb_dir, output_csv=None, binder_chains=["A", "B"], 
                              target_chains=["C"], verbose=False):
    """
    批量计算目录下所有PDB的pH score
    """
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    if not pdb_files:
        print(f"错误: 在 {pdb_dir} 下没有找到PDB文件")
        return
    
    print(f"找到 {len(pdb_files)} 个PDB文件")
    print(f"Binder链: {binder_chains}, Target链: {target_chains}")
    
    results = []
    
    # 所有可能的term（新版）
    locations = ["core", "surf"]
    patterns = ["D", "DD", "D_bb", "DD_bb", "A", "AD", "A_bb", "AD_bb"]
    all_terms = []
    for loc in locations:
        for pat in patterns:
            all_terms.append(f"{loc}_{pat}")
            all_terms.append(f"{loc}_{pat}_chg")
    
    for i, pdb_path in enumerate(sorted(pdb_files), 1):
        pdb_name = os.path.basename(pdb_path)
        print(f"[{i}/{len(pdb_files)}] 处理 {pdb_name}...", end=" ")
        
        try:
            calc = pHScoreCalculatorInternal(
                pdb_path=pdb_path,
                binder_chains=binder_chains,
                target_chains=target_chains,
                buried_cutoff=0.20,
                verbose=verbose
            )
            
            result = calc.calculate_ph_score()
            
            # 构建一行数据
            row = {
                "pdb_name": pdb_name,
                "ph_score": result["ph_score"],
                "n_his_total": result["n_histidines_total"],
                "n_his_binder": result["n_histidines_binder"],
                "n_his_with_hbonds": result["n_his_with_internal_hbonds"],
                "status": "success"
            }
            
            # 添加所有term的count
            for term in all_terms:
                row[term] = result["terms"].get(term, 0)
            
            results.append(row)
            
            # 打印简要信息
            n_core = sum(v for k, v in result["terms"].items() if k.startswith("core"))
            n_chg = sum(v for k, v in result["terms"].items() if k.endswith("_chg"))
            print(f"✓ score={result['ph_score']:.1f} "
                  f"(core={n_core}, charged={n_chg})")
            
            # 详细输出
            if verbose and result["his_details"]:
                for d in result["his_details"]:
                    print(f"    His {d['chain']}{d['resno']}: "
                          f"{d['location']}, patterns={d['patterns']}, "
                          f"charged={d['has_charged']}, "
                          f"partners={d['partners']}")
            
        except Exception as e:
            print(f"✗ 失败: {str(e)}")
            import traceback
            if verbose:
                traceback.print_exc()
            row = {
                "pdb_name": pdb_name,
                "ph_score": None,
                "n_his_total": None,
                "n_his_binder": None,
                "n_his_with_hbonds": None,
                "status": f"error: {str(e)}"
            }
            for term in all_terms:
                row[term] = None
            results.append(row)
    
    df = pd.DataFrame(results)
    
    # 调整列顺序
    info_cols = ["pdb_name", "status", "ph_score", "n_his_total", 
                 "n_his_binder", "n_his_with_hbonds"]
    # 只保留有非零值的term列
    term_cols = [t for t in all_terms if df[t].sum() > 0] if len(results) > 0 else all_terms
    df = df[info_cols + term_cols]
    
    # 输出CSV
    if output_csv is None:
        script_dir = os.path.dirname(os.path.abspath(__file__)) or "."
        output_csv = os.path.join(script_dir, "ph_scores_internal.csv")
    
    df.to_csv(output_csv, index=False)
    
    print(f"\n{'='*60}")
    print(f"完成！共处理 {len(pdb_files)} 个PDB文件")
    print(f"成功: {len([r for r in results if r['status'] == 'success'])} 个")
    print(f"失败: {len([r for r in results if r['status'] != 'success'])} 个")
    print(f"结果已保存到: {output_csv}")
    print(f"{'='*60}")
    
    # 统计摘要
    success_df = df[df["status"] == "success"]
    if len(success_df) > 0:
        print(f"\n统计摘要:")
        print(f"  pH score 范围: {success_df['ph_score'].min():.1f} - {success_df['ph_score'].max():.1f}")
        print(f"  pH score 平均: {success_df['ph_score'].mean():.1f}")
        print(f"  非零 score 数量: {(success_df['ph_score'] > 0).sum()}")
    
    return df


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="计算抗体内部稳定性的pH score (v2)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python calc_ph_score_internal.py ./pdbs/
  python calc_ph_score_internal.py ./pdbs/ -o results.csv -v
  python calc_ph_score_internal.py ./pdbs/ --binder A B --target C
        """
    )
    parser.add_argument("pdb_dir", help="PDB文件目录")
    parser.add_argument("-o", "--output", help="输出CSV路径")
    parser.add_argument("--binder", nargs="+", default=["A"],
                        help="Binder链ID (默认: A B)")
    parser.add_argument("--target", nargs="+", default=["B"],
                        help="Target链ID (默认: C)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="打印详细信息")
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.pdb_dir):
        print(f"错误: {args.pdb_dir} 不是有效的目录")
        sys.exit(1)
    
    batch_calculate_ph_score(
        args.pdb_dir, 
        output_csv=args.output,
        binder_chains=args.binder,
        target_chains=args.target,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()