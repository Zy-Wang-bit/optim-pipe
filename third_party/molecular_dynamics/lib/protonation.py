# third_party/molecular_dynamics/lib/protonation.py
"""His 质子化态管理。

职责：
1. 从 PDB 文件中检测所有 His 残基位置
2. 根据 config 过滤要参与 CpHMD 的残基列表
3. 生成 CpHMD 所需的 λ-dynamics 组参数
4. 固定质子态模式：生成 HIS→HIE/HIP 重命名后的 PDB
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


def prepare_fixed_protonation(
    pdb_path: Path,
    ph: float,
    output_path: Path,
) -> Path:
    """生成固定 His 质子态的 PDB 文件。

    pH >= 7.0 时 HIS → HIE（中性，Nε 质子化）；
    pH <  7.0 时 HIS → HIP（双质子化，+1 电荷）。

    GROMACS pdb2gmx 会根据残基名自动选择正确的质子态，
    无需交互式输入。

    Parameters
    ----------
    pdb_path : 输入 PDB 文件
    ph : 目标 pH 值
    output_path : 输出 PDB 文件

    Returns
    -------
    Path  输出文件路径
    """
    target = "HIP" if ph < 7.0 else "HIE"
    his_variants = {"HIS", "HID", "HIE", "HIP", "HSE", "HSD", "HSP"}
    n_renamed = 0

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(pdb_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM  ", "HETATM")):
                resname = line[17:20].strip()
                if resname in his_variants and resname != target:
                    line = line[:17] + f"{target:>3}" + line[20:]
                    n_renamed += 1
            fout.write(line)

    logger.info(
        "Fixed protonation (pH=%.1f): %s → %s, %d atom lines renamed in %s",
        ph, "HIS/*", target, n_renamed, output_path.name,
    )
    return output_path


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
