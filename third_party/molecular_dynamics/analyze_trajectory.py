#!/usr/bin/env python3
# third_party/molecular_dynamics/analyze_trajectory.py
"""轨迹分析入口。

用法：
  python analyze_trajectory.py --traj experiments/1E62/R3/md/HE1H/pH_7.4/
  python analyze_trajectory.py --traj path/ --analyses rmsd rmsf hbond
  python analyze_trajectory.py --traj path/ --config configs/my_config.yaml
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

_MODULE_ROOT = Path(__file__).resolve().parent
if str(_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(_MODULE_ROOT))

from lib.config import load_config
from lib.analyzers import ANALYZERS

logger = logging.getLogger(__name__)


def run_analyses(
    traj_dir: Path,
    cfg: dict,
    analyses: list[str] | None = None,
) -> dict[str, dict[str, Any]]:
    """运行指定的分析项。"""
    if analyses is None:
        analyses = list(ANALYZERS.keys())

    # mmpbsa 默认不运行除非显式指定
    if "mmpbsa" in analyses and not cfg["analysis"].get("mmpbsa", {}).get("enabled"):
        if analyses == list(ANALYZERS.keys()):
            analyses.remove("mmpbsa")

    all_results = {}
    for name in analyses:
        if name not in ANALYZERS:
            logger.warning("Unknown analyzer: %s (available: %s)", name, list(ANALYZERS.keys()))
            continue
        logger.info("--- Running %s analysis ---", name)
        analyzer = ANALYZERS[name](traj_dir, cfg)
        try:
            results = analyzer.run()
            all_results[name] = results
        except Exception:
            logger.exception("Analyzer %s failed", name)
            all_results[name] = {"summary": {"error": "failed"}}

    # 保存汇总摘要
    summary = {}
    for name, results in all_results.items():
        if "summary" in results:
            for k, v in results["summary"].items():
                summary[f"{name}_{k}"] = v

    summary_path = traj_dir / "analysis" / "summary.json"
    summary_path.parent.mkdir(exist_ok=True)
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    logger.info("Summary saved: %s", summary_path)

    return all_results


def main():
    parser = argparse.ArgumentParser(description="MD 轨迹分析")
    parser.add_argument("--traj", type=Path, required=True, help="轨迹目录")
    parser.add_argument(
        "--analyses", nargs="+", default=None,
        help=f"分析项（可选：{', '.join(ANALYZERS.keys())}）"
    )
    parser.add_argument("--config", type=Path, default=None, help="配置文件路径")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg = load_config(args.config)
    run_analyses(args.traj, cfg, args.analyses)


if __name__ == "__main__":
    main()
