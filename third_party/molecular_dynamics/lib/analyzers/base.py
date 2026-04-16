# third_party/molecular_dynamics/lib/analyzers/base.py
"""分析器基类：统一接口 + 收敛判断工具。"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def resolve_chain_segids(u, chains: list[str]) -> list[str]:
    """将用户配置的链名 (如 "A", "B") 解析为实际 segid。

    GROMACS 输出的 segid 通常是 "seg_0_Protein_chain_A" 而非简单的 "A"。
    本函数依次尝试：精确匹配 → GROMACS 格式匹配 → 返回原始值。
    """
    resolved = []
    all_segids = [seg.segid for seg in u.segments]
    for c in chains:
        # 精确匹配
        if c in all_segids:
            resolved.append(c)
            continue
        # GROMACS 格式: seg_N_Protein_chain_X
        matches = [s for s in all_segids if s.endswith(f"_chain_{c}")]
        if matches:
            resolved.extend(matches)
            continue
        # 回退：保留原始值
        resolved.append(c)
    return resolved


def select_chains(u, chains: list[str], extra_sel: str = "") -> "mda.AtomGroup":
    """按链名选择原子，自动处理 GROMACS segid 格式。

    Parameters
    ----------
    u : MDAnalysis Universe
    chains : 链名列表，如 ["A", "B"]
    extra_sel : 额外的选择字符串，如 "not name H*" 或 "name CA"
    """
    import MDAnalysis as mda

    # 先尝试 segid (解析后)
    segids = resolve_chain_segids(u, chains)
    parts = " or ".join(f"segid {s}" for s in segids)
    sel = f"({parts})"
    if extra_sel:
        sel = f"{extra_sel} and {sel}"
    grp = u.select_atoms(sel)
    if len(grp) > 0:
        return grp

    # 回退到 chainID
    parts = " or ".join(f"chainID {c}" for c in chains)
    sel = f"({parts})"
    if extra_sel:
        sel = f"{extra_sel} and {sel}"
    grp = u.select_atoms(sel)
    return grp


class BaseAnalyzer(ABC):
    """轨迹分析器基类。

    Parameters
    ----------
    traj_dir : Path
        轨迹目录（包含 production.xtc 和 production.tpr）
    cfg : dict
        完整配置
    """

    name: str = "base"

    def __init__(self, traj_dir: Path, cfg: dict):
        self.traj_dir = Path(traj_dir)
        self.cfg = cfg
        self.analysis_cfg = cfg["analysis"]
        self.tpr = self.traj_dir / "production.tpr"
        # 优先使用 PBC unwrap 后的轨迹
        unwrap_xtc = self.traj_dir / "production_unwrap.xtc"
        self.xtc = unwrap_xtc if unwrap_xtc.exists() else self.traj_dir / "production.xtc"
        self.output_dir = self.traj_dir / "analysis"
        self.output_dir.mkdir(exist_ok=True)

    def check_inputs(self) -> bool:
        """检查输入文件是否存在。"""
        if not self.tpr.exists():
            logger.error("TPR file not found: %s", self.tpr)
            return False
        if not self.xtc.exists():
            logger.error("XTC file not found: %s", self.xtc)
            return False
        return True

    def load_universe(self):
        """加载 MDAnalysis Universe 并处理 PBC unwrap。

        GROMACS xtc 轨迹中分子可能跨越周期性边界，导致坐标跳变。
        统一做 unwrap + center 处理，所有分析器应使用此方法加载轨迹。
        """
        import MDAnalysis as mda
        from MDAnalysis.transformations import unwrap, center_in_box

        u = mda.Universe(str(self.tpr), str(self.xtc))
        protein = u.select_atoms("protein")
        if len(protein) > 0:
            transforms = [
                unwrap(protein),
                center_in_box(protein, center="mass"),
            ]
            u.trajectory.add_transformations(*transforms)
            logger.debug("Applied PBC unwrap + centering (%d protein atoms)", len(protein))
        return u

    @abstractmethod
    def run(self) -> dict[str, Any]:
        """执行分析，返回摘要指标。"""

    @abstractmethod
    def save(self, results: dict[str, Any]) -> Path:
        """保存详细结果到 CSV，返回文件路径。"""


def detect_convergence(
    time_ns: np.ndarray,
    values: np.ndarray,
    window_ns: float = 5.0,
    variance_threshold: float = 0.05,
) -> float:
    """检测 RMSD 时间序列的收敛起始时间。

    使用 sliding window 方差判断：当连续窗口的方差 < 阈值时认为收敛。

    Parameters
    ----------
    time_ns : 时间点（ns）
    values : 对应的值（如 RMSD in Å）
    window_ns : 滑动窗口大小（ns）
    variance_threshold : 方差阈值（Å²）

    Returns
    -------
    float  收敛起始时间（ns），未收敛返回 -1.0
    """
    if len(time_ns) < 2:
        return -1.0

    dt = time_ns[1] - time_ns[0]
    window_frames = max(1, int(window_ns / dt))

    for i in range(len(values) - window_frames):
        window = values[i : i + window_frames]
        if np.var(window) < variance_threshold:
            return float(time_ns[i])

    return -1.0
