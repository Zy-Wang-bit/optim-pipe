# third_party/molecular_dynamics/lib/protonation.py
"""His 质子化态管理，用于 CpHMD λ-dynamics 配置。

职责：
1. 从 PDB 文件中检测所有 His 残基位置
2. 根据 config 过滤要参与 CpHMD 的残基列表
3. 生成 CpHMD 所需的 λ-dynamics 组参数
"""

import logging
from pathlib import Path

from Bio.PDB import PDBParser

logger = logging.getLogger(__name__)


def detect_his_residues(pdb_path: Path) -> list[dict]:
    """从 PDB 文件中检测所有 His 残基。

    Returns
    -------
    list[dict]  每个元素 {"chain": str, "resid": int, "resname": str}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", str(pdb_path))
    his_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue  # skip HETATM
                if residue.resname in ("HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP"):
                    his_list.append({
                        "chain": chain.id,
                        "resid": residue.id[1],
                        "resname": residue.resname,
                    })
        break  # only first model
    logger.info("Detected %d His residues in %s", len(his_list), pdb_path.name)
    return his_list


def filter_titratable(
    all_his: list[dict],
    config_residues: list[dict],
) -> list[dict]:
    """过滤参与 CpHMD 的 His 残基。

    Parameters
    ----------
    all_his : detect_his_residues 返回的列表
    config_residues : md_config.yaml 中 cphmd.titratable_residues
                     空列表 = 使用全部检测到的 His

    Returns
    -------
    list[dict]  过滤后的 His 列表
    """
    if not config_residues:
        return all_his
    keep = {(r["chain"], r["resid"]) for r in config_residues}
    return [h for h in all_his if (h["chain"], h["resid"]) in keep]


def build_cphmd_groups(
    titratable_his: list[dict],
    barrier_height: float = 5.0,
) -> list[dict]:
    """为每个可滴定 His 生成 CpHMD λ-dynamics 组参数。

    Parameters
    ----------
    titratable_his : filter_titratable 返回的列表
    barrier_height : λ 势垒高度 (kJ/mol)

    Returns
    -------
    list[dict]  CpHMD 组参数，可直接传给 production.mdp.j2 的 cphmd_groups
    """
    groups = []
    for his in titratable_his:
        name = f"HIS_{his['chain']}{his['resid']}"
        groups.append({
            "type": "histidine",
            "name": name,
            "index_group": name,
            "barrier": barrier_height,
            "init_lambda": 0.5,  # 从中间态开始
        })
    return groups
