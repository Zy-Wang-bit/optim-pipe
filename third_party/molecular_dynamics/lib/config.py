"""MD 模块配置加载与校验。"""

from pathlib import Path
from typing import Any

import yaml

_MODULE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = _MODULE_ROOT / "configs" / "md_config.yaml"


def load_config(path: Path | str | None = None) -> dict[str, Any]:
    """加载并校验 md_config.yaml。

    Parameters
    ----------
    path : 配置文件路径，默认 configs/md_config.yaml

    Returns
    -------
    dict  完整配置字典
    """
    path = Path(path) if path else DEFAULT_CONFIG
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with open(path) as f:
        cfg = yaml.safe_load(f)
    _validate(cfg)
    return cfg


def _validate(cfg: dict) -> None:
    """基本校验：必需字段存在且类型正确。"""
    md = cfg.get("md")
    if md is None:
        raise ValueError("Config missing top-level 'md' section")

    required_md = ["force_field", "water_model", "temperature", "stages"]
    for key in required_md:
        if key not in md:
            raise ValueError(f"Config md.{key} is required")

    stages = md["stages"]
    for stage in ("minimization", "nvt", "npt", "production"):
        if stage not in stages:
            raise ValueError(f"Config md.stages.{stage} is required")

    analysis = cfg.get("analysis")
    if analysis is None:
        raise ValueError("Config missing top-level 'analysis' section")

    for key in ("antibody_chains", "antigen_chains", "cdr_regions"):
        if key not in analysis:
            raise ValueError(f"Config analysis.{key} is required")


def get_production_nsteps(cfg: dict) -> int:
    """从配置计算 production 阶段的总步数。"""
    prod = cfg["md"]["stages"]["production"]
    duration_ps = prod["duration_ns"] * 1000
    dt_ps = prod["dt_ps"]
    return int(duration_ps / dt_ps)


def get_save_nsteps(cfg: dict) -> int:
    """从配置计算轨迹保存间隔的步数。"""
    prod = cfg["md"]["stages"]["production"]
    return int(prod["save_interval_ps"] / prod["dt_ps"])
