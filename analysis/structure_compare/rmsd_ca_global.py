#!/usr/bin/env python3
"""
rmsd_ca_global.py - 计算野生型和突变体抗体的全局Cα-RMSD
使用: python rmsd_ca_global.py <wt_pdb_path> <mutant_pdb_path>
"""

import sys
import numpy as np
from Bio.PDB import PDBParser, Superimposer

def extract_ca_atoms(structure):
    """提取所有链的Cα原子，返回{(chain_id, resseq, icode): atom}"""
    ca_dict = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ':  # 跳过HETATM
                    continue
                if 'CA' in residue:
                    # residue.id = (' ', resseq, icode)
                    key = (chain.id, residue.id[1], residue.id[2])
                    ca_dict[key] = residue['CA']
        break  # 只处理第一个model
    return ca_dict

def calculate_rmsd(wt_pdb, mut_pdb):
    """计算两个结构的全局Cα-RMSD"""
    parser = PDBParser(QUIET=True)
    
    # 解析PDB
    wt_structure = parser.get_structure('wt', wt_pdb)
    mut_structure = parser.get_structure('mut', mut_pdb)
    
    # 提取Cα原子
    wt_ca = extract_ca_atoms(wt_structure)
    mut_ca = extract_ca_atoms(mut_structure)
    
    # 找共同残基
    common_keys = sorted(set(wt_ca.keys()) & set(mut_ca.keys()))
    
    if len(common_keys) < 3:
        return float('nan')
    
    # 提取坐标
    wt_coords = np.array([wt_ca[key].coord for key in common_keys])
    mut_coords = np.array([mut_ca[key].coord for key in common_keys])
    
    # Kabsch对齐并计算RMSD
    super_imposer = Superimposer()
    super_imposer.set_atoms(
        [wt_ca[key] for key in common_keys],
        [mut_ca[key] for key in common_keys]
    )
    
    return super_imposer.rms

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <wt_pdb_path> <mutant_pdb_path>")
        sys.exit(1)
    
    # 修改这里，直接指定你的PDB文件绝对路径
    wt_pdb = sys.argv[1]     # 例如: "/path/to/your/wt.pdb"
    mut_pdb = sys.argv[2]    # 例如: "/path/to/your/mutant.pdb"
    
    rmsd = calculate_rmsd(wt_pdb, mut_pdb)
    print(f"{rmsd:.3f}")

if __name__ == "__main__":
    main()