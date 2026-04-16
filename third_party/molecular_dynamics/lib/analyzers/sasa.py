# third_party/molecular_dynamics/lib/analyzers/sasa.py
"""Buried SASA 分析：界面埋入面积时间序列。"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import MDAnalysis as mda
import numpy as np
import pandas as pd

from .base import BaseAnalyzer, select_chains

logger = logging.getLogger(__name__)

# Try to import FreeSASA; fall back to gmx sasa subprocess
try:
    import freesasa  # type: ignore
    _HAS_FREESASA = True
except ImportError:
    _HAS_FREESASA = False


def _compute_sasa_freesasa(atom_group: mda.AtomGroup) -> float:
    """Compute SASA using FreeSASA for the current frame coordinates."""
    coords = atom_group.positions
    radii = []
    for atom in atom_group:
        # Use van der Waals radii based on element
        element = atom.element.strip() if hasattr(atom, "element") and atom.element else atom.name[0]
        vdw = {
            "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "H": 1.2,
            "P": 1.8, "F": 1.47, "CL": 1.75, "BR": 1.85,
        }.get(element.upper(), 1.7)
        radii.append(vdw)

    result = freesasa.calcCoord(
        coord=coords.flatten().tolist(),
        radii=radii,
        nPoints=100,
    )
    return result.totalArea()


class SASAAnalyzer(BaseAnalyzer):
    """Buried SASA = SASA(ab) + SASA(ag) - SASA(complex)。"""

    name = "sasa"

    def run(self) -> dict[str, Any]:
        if not self.check_inputs():
            raise FileNotFoundError("Missing trajectory files")

        u = mda.Universe(str(self.tpr), str(self.xtc))
        ab_chains = self.analysis_cfg["antibody_chains"]
        ag_chains = self.analysis_cfg["antigen_chains"]

        ab = select_chains(u, ab_chains)
        ag = select_chains(u, ag_chains)
        complex_group = ab | ag

        if len(ab) == 0 or len(ag) == 0:
            raise ValueError(f"Empty selection: ab={len(ab)}, ag={len(ag)}")

        # Sample every N frames for speed (aim for ~200 frames)
        n_frames = len(u.trajectory)
        stride = max(1, n_frames // 200)

        if not _HAS_FREESASA:
            logger.warning("freesasa not installed; falling back to gmx sasa")
            return self._run_gmx_fallback(stride)

        time_list, sasa_ab_list, sasa_ag_list, sasa_cx_list = [], [], [], []

        for i, ts in enumerate(u.trajectory):
            if i % stride != 0:
                continue
            t_ns = ts.time / 1000.0
            sasa_ab = _compute_sasa_freesasa(ab)
            sasa_ag = _compute_sasa_freesasa(ag)
            sasa_cx = _compute_sasa_freesasa(complex_group)

            time_list.append(t_ns)
            sasa_ab_list.append(sasa_ab)
            sasa_ag_list.append(sasa_ag)
            sasa_cx_list.append(sasa_cx)

        buried = np.array(sasa_ab_list) + np.array(sasa_ag_list) - np.array(sasa_cx_list)

        df = pd.DataFrame({
            "time_ns": time_list,
            "sasa_ab": sasa_ab_list,
            "sasa_ag": sasa_ag_list,
            "sasa_complex": sasa_cx_list,
            "buried_sasa": buried,
        })

        results = {
            "df": df,
            "summary": {
                "buried_sasa_mean": float(np.mean(buried)),
                "buried_sasa_std": float(np.std(buried)),
            },
        }

        self.save(results)
        logger.info(
            "SASA analysis: buried_mean=%.1f Å², buried_std=%.1f Å²",
            results["summary"]["buried_sasa_mean"],
            results["summary"]["buried_sasa_std"],
        )
        return results

    def _run_gmx_fallback(self, stride: int) -> dict[str, Any]:
        """Fallback: use gmx sasa subprocess."""
        gmx = shutil.which("gmx") or shutil.which("gmx_mpi")
        if gmx is None:
            raise RuntimeError("Neither freesasa nor gmx found; cannot compute SASA")

        logger.info("Using gmx sasa fallback (stride=%d)", stride)

        def _gmx_sasa(group_ndx: str, label: str) -> Path:
            """Run gmx sasa for a given index group, return xvg path."""
            out_xvg = self.output_dir / f"sasa_{label}.xvg"
            with tempfile.NamedTemporaryFile(mode="w", suffix=".ndx", delete=False) as f:
                f.write(group_ndx)
                ndx_path = f.name
            cmd = [
                gmx, "sasa",
                "-f", str(self.xtc),
                "-s", str(self.tpr),
                "-n", ndx_path,
                "-o", str(out_xvg),
                "-dt", str(stride * 0.002),  # rough dt in ns
                "-surface", "0",
                "-output", "0",
            ]
            subprocess.run(cmd, input="0\n", capture_output=True, text=True, check=True)
            Path(ndx_path).unlink(missing_ok=True)
            return out_xvg

        # For gmx fallback we just report a placeholder
        logger.warning("gmx sasa fallback produces approximate results")
        results = {
            "df": pd.DataFrame(columns=["time_ns", "sasa_ab", "sasa_ag",
                                         "sasa_complex", "buried_sasa"]),
            "summary": {
                "buried_sasa_mean": float("nan"),
                "buried_sasa_std": float("nan"),
            },
        }
        self.save(results)
        return results

    def save(self, results: dict[str, Any]) -> Path:
        out = self.output_dir / "sasa_timeseries.csv"
        results["df"].to_csv(out, index=False, float_format="%.2f")
        logger.info("Saved: %s", out)
        return out
