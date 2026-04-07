#!/usr/bin/env python3
# third_party/molecular_dynamics/run_md.py
"""MD 模拟执行入口。

用法：
  # 单 PDB，标准 MD
  python run_md.py --pdb input.pdb

  # 单 PDB，CpHMD 双 pH
  python run_md.py --pdb input.pdb --ph 7.4 6.0

  # 批量执行
  python run_md.py --pdb-dir variants/ --ph 7.4 6.0

  # 指定输出目录
  python run_md.py --pdb input.pdb --ph 6.0 --output-dir experiments/1E62/R3/md/

  # 指定配置文件
  python run_md.py --pdb input.pdb --config configs/my_config.yaml
"""

import argparse
import copy
import logging
import sys
from pathlib import Path

# 将模块根目录加入 sys.path
_MODULE_ROOT = Path(__file__).resolve().parent
if str(_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(_MODULE_ROOT))

from lib.config import load_config
from lib.gromacs_wrapper import GromacsWrapper
from lib.protonation import detect_his_residues, filter_titratable, build_cphmd_groups

logger = logging.getLogger(__name__)


def run_single(
    pdb_path: Path,
    ph_values: list[float],
    cfg: dict,
    output_dir: Path,
) -> None:
    """对单个 PDB 运行 MD 模拟。

    如果 ph_values 为空，运行一次标准 MD。
    如果提供 ph_values，为每个 pH 运行一次。
    """
    variant_name = pdb_path.stem
    pdb_path = pdb_path.resolve()

    if not ph_values:
        # 标准 MD，无 pH 指定
        work_dir = output_dir / variant_name
        logger.info("Running standard MD for %s", variant_name)
        wrapper = GromacsWrapper(cfg, work_dir)
        wrapper.run_full_pipeline(pdb_path)
        return

    for ph in ph_values:
        work_dir = output_dir / variant_name / f"pH_{ph}"
        logger.info("Running MD for %s at pH %.1f", variant_name, ph)

        # CpHMD 设置
        cphmd_cfg = cfg["md"].get("cphmd", {})
        if cphmd_cfg.get("enabled", False):
            his_list = detect_his_residues(pdb_path)
            titratable = filter_titratable(
                his_list, cphmd_cfg.get("titratable_residues", [])
            )
            groups = build_cphmd_groups(
                titratable,
                barrier_height=cphmd_cfg["lambda_dynamics"]["barrier_height"],
            )
            logger.info(
                "CpHMD: %d titratable His at pH %.1f", len(titratable), ph
            )
            # 将 groups 注入 config 的副本
            cfg_copy = copy.deepcopy(cfg)
            cfg_copy["md"]["cphmd"]["_runtime_groups"] = groups
            wrapper = GromacsWrapper(cfg_copy, work_dir)
        else:
            wrapper = GromacsWrapper(cfg, work_dir)

        wrapper.run_full_pipeline(pdb_path, ph=ph)


def main():
    parser = argparse.ArgumentParser(description="MD 模拟执行")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb", type=Path, help="单个 PDB 文件路径")
    group.add_argument("--pdb-dir", type=Path, help="PDB 文件目录（批量执行）")
    parser.add_argument(
        "--ph", type=float, nargs="+", default=[],
        help="模拟的 pH 值列表（如 7.4 6.0）"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="输出目录（默认从 config 读取）"
    )
    parser.add_argument(
        "--config", type=Path, default=None,
        help="配置文件路径（默认 configs/md_config.yaml）"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg = load_config(args.config)
    output_dir = args.output_dir or Path(cfg["md"]["output_dir"])
    output_dir = output_dir.resolve()

    # 收集 PDB 文件列表
    if args.pdb:
        pdb_files = [args.pdb]
    else:
        pdb_files = sorted(args.pdb_dir.glob("*.pdb"))
        if not pdb_files:
            logger.error("No PDB files found in %s", args.pdb_dir)
            sys.exit(1)
        logger.info("Found %d PDB files in %s", len(pdb_files), args.pdb_dir)

    for pdb in pdb_files:
        run_single(pdb, args.ph, cfg, output_dir)

    logger.info("All MD simulations complete. Output: %s", output_dir)


if __name__ == "__main__":
    main()
