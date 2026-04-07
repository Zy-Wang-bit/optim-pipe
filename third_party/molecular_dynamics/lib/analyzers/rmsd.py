# third_party/molecular_dynamics/lib/analyzers/rmsd.py
"""RMSD 时间序列分析 + 自动收敛检测。"""

import logging
from pathlib import Path
from typing import Any

import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD as MDA_RMSD
import numpy as np
import pandas as pd

from .base import BaseAnalyzer, detect_convergence

logger = logging.getLogger(__name__)


class RMSDAnalyzer(BaseAnalyzer):
    """计算全局 + CDR 区域 RMSD 时间序列。"""

    name = "rmsd"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        ca = u.select_atoms("name CA")

        # 全局 RMSD
        rmsd_analysis = MDA_RMSD(ca, ca, ref_frame=0)
        rmsd_analysis.run()
        # results.rmsd: shape (n_frames, 3) — [frame, time(ps), rmsd(Å)]
        time_ns = rmsd_analysis.results.rmsd[:, 1] / 1000.0
        global_rmsd = rmsd_analysis.results.rmsd[:, 2]

        # CDR RMSD
        cdr_rmsds = {}
        cdr_regions = self.analysis_cfg["cdr_regions"]
        for cdr_name, cdr_def in cdr_regions.items():
            sel = (
                f"name CA and segid {cdr_def['chain']} "
                f"and resid {cdr_def['start']}:{cdr_def['end']}"
            )
            cdr_atoms = u.select_atoms(sel)
            if len(cdr_atoms) == 0:
                # 回退到 chainID 选择
                sel = (
                    f"name CA and chainID {cdr_def['chain']} "
                    f"and resid {cdr_def['start']}:{cdr_def['end']}"
                )
                cdr_atoms = u.select_atoms(sel)

            if len(cdr_atoms) == 0:
                logger.warning("CDR %s: no atoms found, skipping", cdr_name)
                cdr_rmsds[cdr_name] = np.full(len(time_ns), np.nan)
                continue

            cdr_rmsd = MDA_RMSD(cdr_atoms, cdr_atoms, ref_frame=0)
            cdr_rmsd.run()
            cdr_rmsds[cdr_name] = cdr_rmsd.results.rmsd[:, 2]

        # 收敛检测
        conv_cfg = self.analysis_cfg["convergence"]
        converge_ns = detect_convergence(
            time_ns, global_rmsd,
            window_ns=conv_cfg["window_ns"],
            variance_threshold=conv_cfg["variance_threshold"],
        )

        # 构建结果
        results = {
            "time_ns": time_ns,
            "global_rmsd": global_rmsd,
            "cdr_rmsds": cdr_rmsds,
            "converge_ns": converge_ns,
            "summary": {
                "rmsd_mean": float(np.mean(global_rmsd)),
                "rmsd_std": float(np.std(global_rmsd)),
                "converge_ns": converge_ns,
            },
        }
        # CDR 摘要
        for name, vals in cdr_rmsds.items():
            key = f"{name.lower()}_rmsd_mean"
            results["summary"][key] = float(np.nanmean(vals))

        self.save(results)
        logger.info(
            "RMSD analysis: mean=%.3f Å, converge=%.1f ns",
            results["summary"]["rmsd_mean"],
            converge_ns,
        )
        return results

    def save(self, results: dict[str, Any]) -> Path:
        """保存 RMSD 时间序列到 CSV。"""
        data = {"time_ns": results["time_ns"], "global": results["global_rmsd"]}
        for name, vals in results["cdr_rmsds"].items():
            data[name] = vals
        df = pd.DataFrame(data)
        out = self.output_dir / "rmsd_timeseries.csv"
        df.to_csv(out, index=False, float_format="%.4f")
        logger.info("Saved: %s", out)
        return out
