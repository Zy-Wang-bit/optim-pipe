# third_party/molecular_dynamics/lib/analyzers/pairwise_dist.py
"""残基对距离时间序列分析。

跟踪指定残基对之间的最短重原子距离随时间变化，
用于监控盐桥/接触的动态行为。

配置示例（md_config.yaml）：
  analysis:
    pairwise_dist:
      pairs:
        - { chain1: "B", resid1: 42, chain2: "C", resid2: 122, label: "HIS42-LYS122" }
        - { chain1: "C", resid1: 122, chain2: "A", resid2: 103, label: "LYS122-ASP103" }
"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import numpy as np
import pandas as pd

from .base import BaseAnalyzer

logger = logging.getLogger(__name__)


class PairwiseDistAnalyzer(BaseAnalyzer):
    """指定残基对的最短重原子距离时间序列。"""

    name = "pairwise_dist"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        pair_cfg = self.analysis_cfg.get("pairwise_dist", {})
        pairs = pair_cfg.get("pairs", [])
        if not pairs:
            logger.warning("No residue pairs configured in analysis.pairwise_dist.pairs")
            return {"summary": {}, "timeseries_df": pd.DataFrame()}

        u = mda.Universe(str(self.tpr), str(self.xtc))

        # 为每个 pair 准备 atom groups
        pair_selections = []
        for p in pairs:
            label = p.get("label", f"{p['chain1']}:{p['resid1']}-{p['chain2']}:{p['resid2']}")
            sel1 = self._select_residue_heavy(u, p["chain1"], p["resid1"])
            sel2 = self._select_residue_heavy(u, p["chain2"], p["resid2"])
            if len(sel1) == 0 or len(sel2) == 0:
                logger.warning("Empty selection for pair %s: sel1=%d, sel2=%d",
                               label, len(sel1), len(sel2))
                continue
            pair_selections.append((label, sel1, sel2))

        if not pair_selections:
            logger.error("No valid pairs found")
            return {"summary": {}, "timeseries_df": pd.DataFrame()}

        # 遍历轨迹
        rows = []
        for ts in u.trajectory:
            time_ns = ts.time / 1000.0
            row = {"time_ns": time_ns}
            for label, sel1, sel2 in pair_selections:
                dists = mda.analysis.distances.distance_array(
                    sel1.positions, sel2.positions)
                row[label] = float(np.min(dists))
            rows.append(row)

        df = pd.DataFrame(rows)

        # Summary statistics
        summary = {}
        for label, _, _ in pair_selections:
            vals = df[label].values
            summary[f"{label}_mean"] = float(np.mean(vals))
            summary[f"{label}_std"] = float(np.std(vals))
            summary[f"{label}_min"] = float(np.min(vals))
            summary[f"{label}_max"] = float(np.max(vals))
            # 盐桥/接触占据率 (< 4.5Å)
            summary[f"{label}_occupancy_4.5A"] = float(np.mean(vals < 4.5))

        results = {
            "timeseries_df": df,
            "summary": summary,
        }
        self.save(results)
        for label, _, _ in pair_selections:
            logger.info("  %s: mean=%.1fÅ, occ(<4.5Å)=%.1f%%",
                        label, summary[f"{label}_mean"],
                        summary[f"{label}_occupancy_4.5A"] * 100)
        return results

    def save(self, results: dict[str, Any]) -> Path:
        path = self.output_dir / "pairwise_dist_timeseries.csv"
        results["timeseries_df"].to_csv(path, index=False, float_format="%.3f")
        logger.info("Saved: %s", path)
        return path

    @staticmethod
    def _select_residue_heavy(u: mda.Universe, chain: str, resid: int) -> mda.AtomGroup:
        """选择指定链和残基号的重原子。先尝试 segid，再尝试 chainID。"""
        for kw in ("segid", "chainID"):
            sel = f"not name H* and {kw} {chain} and resid {resid}"
            grp = u.select_atoms(sel)
            if len(grp) > 0:
                return grp
        return mda.AtomGroup([], u)
