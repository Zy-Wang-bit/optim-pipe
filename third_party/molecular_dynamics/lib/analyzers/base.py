# third_party/molecular_dynamics/lib/analyzers/base.py
"""分析器基类：统一接口 + 收敛判断工具。"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


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
        self.xtc = self.traj_dir / "production.xtc"
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
